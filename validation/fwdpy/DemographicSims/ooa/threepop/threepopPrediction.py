import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd
import demes
import demesdraw
from collections import defaultdict


os.chdir("/media/nathan/T7/BGSdemo/parsedooaThreepopData")

# hardcoding some parameters
u = 1e-8
r = 1e-8
L = 1e6
regionSize = 1e4
tol = 1e-3
focalPos = 5e5
sample_size = 40
proj_size = 40

# positions of point masses.
# split 1Mb into 10kb regions, each point mass lies at the center
# of these regions 
pointMassPosition = range(int(5e3),int(100.5e4),int(1e4))

def _convert_to_generations(g, sample_times=None):
    """
    Takes a deme graph that is not in time units of generations and converts
    times to generations, using the time units and generation times given.
    """
    if g.time_units == "generations":
        return g, sample_times
    else:
        for ii, sample_time in enumerate(sample_times):
            sample_times[ii] = sample_time / g.generation_time
        g = g.in_generations()
        return g, sample_times

def pointMassContribution(pos, scaledu, s, t, r, focalPos):
    return - scaledu / s * (s / (r * abs(pos-focalPos) + s) * (1 - math.exp(- r * abs(pos-focalPos) * t - s * t)))**2

def B(positions, u, s, t, r, regionSize, focalPos):
    scaledu = u * regionSize
    return math.exp(sum([pointMassContribution(pos, scaledu, s, t, r, focalPos) for pos in positions]))

# im making a factor of 2 error here somewhere
def rescaledPointMassContribution(pos, scaledu, s, t, r, focalPos, ancestralNe, ancTime):
    return - scaledu / s * (s / (r * abs(pos-focalPos) + s) * (1 - math.exp(- r * abs(pos-focalPos) * (ancTime - t * 2 * ancestralNe) - s * (ancTime - t * 2 * ancestralNe))))**2

def bFromDemes(positions, u, s, r, regionSize, focalPos, demesFile, tol):
    graph = demes.load(demesFile)
    if graph.time_units != "generations":
        graph = graph.in_generations()
    scaledu = u * regionSize
    oldestEpoch, censusSize = getOldestEpoch(graph)
    testFun = [B(positions, u, s, t, r, regionSize, focalPos) for t in range(0,int(10 * censusSize),int(censusSize/10))]
    diffs = [testFun[i+1] - testFun[i] for i in range(len(testFun)-1)]
    ancTime = next((i for i,x in enumerate(diffs) if abs(x) < tol), None)
    ancTime = censusSize / 10 * ancTime 
    if ancTime > oldestEpoch:
        print("woah there partner")
    ancTime = max(ancTime, oldestEpoch)
    ancB = B(positions, u, s, ancTime, r, regionSize, focalPos) 
    ancNe = ancB * censusSize 
    return(lambda t: [math.exp(sum([rescaledPointMassContribution(pos, scaledu, s, t, r, focalPos, ancNe, ancTime) for pos in positions])) / ancB], ancTime, ancNe)
    
def getOldestEpoch(graph):
    tme = 0
    for deme in graph.demes:
        for epoch in deme.epochs:
            if epoch.end_time > tme:
                tme = epoch.end_time
                size = epoch.start_size
    return tme, size
            
    
def censusFun(demesFile):
    graph = demes.load(demesFile)
    if graph.time_units != "generations":
        graph = graph.in_generations()
    def N(t):
        sizes = []
        for deme in graph.demes:
            size_t = None
            for epoch in deme.epochs:
                if epoch.end_time <= t < epoch.start_time:
                    if epoch.size_function == "constant":
                        size_t = epoch.start_size
                    else:
                        dt = (epoch.start_time - t) / epoch.time_span
                        r = math.log(epoch.end_size / epoch.start_size)
                        size_t = epoch.start_size * math.exp(r * dt)
                    break
            sizes.append(size_t)
        return sizes
    return lambda t: N(t)

# TODO this wont work if the B function is larger than the census fun
def reversedCensusFun(demesFile, ancTime, ancNe):
    cs = censusFun(demesFile)
    return lambda t: [x  / cs(ancTime)[0] for x in cs(ancTime - t * 2 * ancNe) if x != None]

def _get_demographic_events(g, demes_demo_events, sampled_demes):
    """
    Returns demographic events and present demes over each epoch.
    Epochs are divided by any demographic event.
    """
    # first get set of all time dividers, from demographic events, migration
    # rate changes, deme epoch changes
    break_points = set()
    for deme in g.demes:
        for e in deme.epochs:
            break_points.add(e.start_time)
            break_points.add(e.end_time)
    for pulse in g.pulses:
        break_points.add(pulse.time)
    for migration in g.migrations:
        break_points.add(migration.start_time)
        break_points.add(migration.end_time)

    # get demes present for each integration epoch
    integration_times = [
        (start_time, end_time)
        for start_time, end_time in zip(
            sorted(list(break_points))[-1:0:-1], sorted(list(break_points))[-2::-1]
        )
    ]

    # find live demes in each epoch, starting with most ancient
    demes_present = defaultdict(list)
    # add demes as they appear from past to present to end of lists
    deme_start_times = defaultdict(list)
    for deme in g.demes:
        deme_start_times[deme.start_time].append(deme.name)

    if math.inf not in deme_start_times.keys():
        raise ValueError("Root deme must have start time as inf")
    if len(deme_start_times[math.inf]) != 1:
        raise ValueError("Deme graph can only have a single root")

    for start_time in sorted(deme_start_times.keys())[::-1]:
        for deme_id in deme_start_times[start_time]:
            end_time = g[deme_id].end_time
            for interval in integration_times:
                if start_time >= interval[0] and end_time <= interval[1]:
                    demes_present[interval].append(deme_id)

    # Dictionary of demographic events, occurring in the order:
    #   branches, pulses, admixtures, mergers, splits.
    # Importantly, splits and mergers remove the parental populations, so if
    # there are events like branches or pulses that involve those parental
    # populations at the same time, they will not be present when we try to
    # apply those events, resulting in an error.
    demo_events = defaultdict(list)
    for branch in demes_demo_events["branches"]:
        event = ("branch", branch.parent, branch.child)
        demo_events[branch.time].append(event)
    for pulse in demes_demo_events["pulses"]:
        event = ("pulse", pulse.sources, pulse.dest, pulse.proportions)
        demo_events[pulse.time].append(event)
    for admix in demes_demo_events["admixtures"]:
        event = ("admix", admix.parents, admix.proportions, admix.child)
        demo_events[admix.time].append(event)
    for merge in demes_demo_events["mergers"]:
        event = ("merge", merge.parents, merge.proportions, merge.child)
        demo_events[merge.time].append(event)
    for split in demes_demo_events["splits"]:
        event = ("split", split.parent, split.children)
        demo_events[split.time].append(event)

    # if there are any unsampled demes that end before present and do not have
    # any descendent demes, we need to add marginalization events.
    for deme_id, succs in g.successors().items():
        if deme_id not in sampled_demes and (
            len(succs) == 0
            or np.all([g[succ].start_time > g[deme_id].end_time for succ in succs])
        ):
            event = ("marginalize", deme_id)
            demo_events[g[deme_id].end_time].append(event)

    return demo_events, demes_present

def _get_integration_parameters(g, demes_present, frozen_list, Ne=None):
    """
    Returns a list of size functions, migration matrices, integration times,
    and lists frozen demes.
    """
    nu_funcs = []
    integration_times = []
    migration_matrices = []
    frozen_demes = []

    for interval, live_demes in sorted(demes_present.items())[::-1]:
        # get intergration time for interval
        T = (interval[0] - interval[1]) / 2 / Ne
        if T == math.inf:
            T = 0
        integration_times.append(T)
        # get frozen attributes
        freeze = [d in frozen_list for d in live_demes]
        frozen_demes.append(freeze)
        # get nu_function or list of sizes (if all constant)
        sizes = []
        for d in live_demes:
            sizes.append(_sizes_at_time(g, d, interval))
        nu_func = _make_nu_func(sizes, T, Ne)
        nu_funcs.append(nu_func)
        # get migration matrix for interval
        mig_mat = np.zeros((len(live_demes), len(live_demes)))
        for ii, d_from in enumerate(live_demes):
            for jj, d_to in enumerate(live_demes):
                if d_from != d_to:
                    m = _migration_rate_in_interval(g, d_from, d_to, interval)
                    mig_mat[jj, ii] = 2 * Ne * m
        migration_matrices.append(mig_mat)

    return nu_funcs, migration_matrices, integration_times, frozen_demes

def _make_nu_func(sizes, T, Ne):
    """
    Given the sizes at start and end of time interval, and the size function for
    each deme, along with the integration time and reference Ne, return the
    size function that gets passed to the moments integration routines.
    """
    if np.all([s[-1] == "constant" for s in sizes]):
        # all constant
        nu_func = [s[0] / Ne for s in sizes]
    else:
        nu_funcs_separated = []
        for s in sizes:
            if s[-1] == "constant":
                assert s[0] == s[1]
                nu_funcs_separated.append(lambda t, N0=s[0]: N0 / Ne)
            elif s[-1] == "linear":
                nu_funcs_separated.append(
                    lambda t, N0=s[0], NF=s[1]: N0 / Ne + t / T * (NF - N0) / Ne
                )
            elif s[-1] == "exponential":
                nu_funcs_separated.append(
                    lambda t, N0=s[0], NF=s[1]: N0
                    / Ne
                    * np.exp(np.log(NF / N0) * t / T)
                )
            else:
                raise ValueError(f"{s[-1]} not a valid size function")

        def nu_func(t):
            return [nu(t) for nu in nu_funcs_separated]

        # check that this is correct, or if we have to "pin" parameters
    return nu_func

def _sizes_at_time(g, deme_id, time_interval):
    """
    Returns the start size, end size, and size function for given deme over the
    given time interval.
    """
    for epoch in g[deme_id].epochs:
        if epoch.start_time >= time_interval[0] and epoch.end_time <= time_interval[1]:
            break
    if epoch.size_function not in ["constant", "exponential", "linear"]:
        raise ValueError(
            "Can only intergrate constant, exponential, or linear size functions"
        )
    size_function = epoch.size_function

    if size_function == "constant":
        start_size = end_size = epoch.start_size

    if epoch.start_time == time_interval[0]:
        start_size = epoch.start_size
    else:
        if size_function == "exponential":
            start_size = epoch.start_size * np.exp(
                np.log(epoch.end_size / epoch.start_size)
                * (epoch.start_time - time_interval[0])
                / epoch.time_span
            )
        elif size_function == "linear":
            frac = (epoch.start_time - time_interval[0]) / epoch.time_span
            start_size = epoch.start_size + frac * (epoch.end_size - epoch.start_size)

    if epoch.end_time == time_interval[1]:
        end_size = epoch.end_size
    else:
        if size_function == "exponential":
            end_size = epoch.start_size * np.exp(
                np.log(epoch.end_size / epoch.start_size)
                * (epoch.start_time - time_interval[1])
                / epoch.time_span
            )
        elif size_function == "linear":
            frac = (epoch.start_time - time_interval[1]) / epoch.time_span
            end_size = epoch.start_size + frac * (epoch.end_size - epoch.start_size)

    return start_size, end_size, size_function

def _migration_rate_in_interval(g, source, dest, time_interval):
    """
    Get the migration rate from source to dest over the given time interval.
    """
    rate = 0
    for mig in g.migrations:
        try:  # if asymmetric migration
            if mig.source == source and mig.dest == dest:
                if (
                    mig.start_time >= time_interval[0]
                    and mig.end_time <= time_interval[1]
                ):
                    rate = mig.rate
        except AttributeError:  # symmetric migration
            if source in mig.demes and dest in mig.demes:
                if (
                    mig.start_time >= time_interval[0]
                    and mig.end_time <= time_interval[1]
                ):
                    rate = mig.rate
    return rate

# todo find out what happens when I redefine these function names???
# # s = curs
# curN = censusSize
# # positions = pointMassPosition
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot([B(pointMassPosition, u, s, t, r, 1e4, 5e5) for t in range(int(10 * censusSize))], "-", ms=8, lw=1, label="Neutral")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend();

# for curs in [1e-3, 5e-3, 1e-2]:
#     for curN in [1e3, 5e3, 1e4]:
fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
fig.subplots_adjust(hspace = .25)

for i in range(3):
    for j in range(1):
        curs = [1e-3, 5e-3, 1e-2][i]
        curdemo = ["ooa.yaml"][j]
        
        os.chdir("/media/nathan/T7/BGSdemo/parsedooaThreepopData")

        simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
        simData = simData[0].to_numpy()
        simData = moments.Spectrum(simData,data_folded=False) 
        projData = simData.project([proj_size])
        
        os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop")

        # test
        demo = demes.load(curdemo)
        if demo.time_units != "generations":
            demo = demo.in_generations()
        demesdraw.tubes(demo);
        
        oldestEpoch, censusSize = getOldestEpoch(demo)

        f, ancTime, ancNe = bFromDemes(pointMassPosition, u, curs, r, regionSize, focalPos, curdemo, tol)
        
        # todo change how we do integration in each epoch
        
        # get the list of demographic events from demes, which is a dictionary with
        # lists of splits, admixtures, mergers, branches, and pulses
        demes_demo_events = demo.discrete_demographic_events()
        
        # get the dict of events and event times that partition integration epochs, in
        # descending order. events include demographic events, such as splits and
        # mergers and admixtures, as well as changes in population sizes or migration
        # rates that require instantaneous changes in the size function or migration matrix.
        # get the list of demes present in each epoch, as a dictionary with non-overlapping
        # adjoint epoch time intervals
        sampled_pops = ["CEU"]
        demo_events, demes_present = _get_demographic_events(
            demo, demes_demo_events, sampled_pops
        )
        
        for epoch, epoch_demes in demes_present.items():
            if len(epoch_demes) > 5:
                raise ValueError(
                    f"Moments cannot integrate more than five demes at a time. "
                    f"Epoch {epoch} has demes {epoch_demes}."
                )
        
        # get the list of size functions, migration matrices, and frozen attributes from
        # the deme graph and event times, matching the integration times
        list_of_frozen_demes = [] # todo need to make ancient samples work
        nu_funcs, mig_mats, Ts, frozen_pops = _get_integration_parameters(
            demo, demes_present, list_of_frozen_demes, Ne=censusSize
        )
        
        # def _get_integration_parameters(g, demes_present, frozen_list, Ne=None):
        # old from here down
        
        cs = reversedCensusFun(curdemo, ancTime, ancNe)
        
        g = lambda t: [x * y for x,y in zip(f(t), cs(t))]
        
        fs = moments.Demographics1D.snm([sample_size])
        fs.integrate(g, ancTime / 2 / ancNe)
        
        oldestEpoch, ancCensusSize = getOldestEpoch(demo)
        ds = reversedCensusFun(curdemo, ancTime, ancCensusSize)
        fs_neu = moments.Demographics1D.snm([proj_size])
        fs_neu.integrate(ds, ancTime / 2 / ancCensusSize)
        # fs_neu = moments.Demographics1D.snm([proj_size])
        # fs_neu.integrate(cs, ancTime / 2 / ancTime)
        # fs_neu = moments.Demographics1D.snm([proj_size])
        # fs_neu.integrate(cs, ancTime / 2 / ancNe)
        
        sampled_demes = ["CEU"]

        ds = moments.Spectrum.from_demes(
            curdemo, sampled_demes=sampled_demes, sample_sizes=[proj_size]
        )

        # normalizing so singletons have freq 1, cause thats all I can think of right now
        fs = fs * 8 * 1e-8 * ancNe
        projData = projData * 1e-8
        fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize   
        ds = ds * 8 * 1e-8 * ancCensusSize
        
        ax[i].plot(fs, ".-", ms=8, lw=1, label="BGS")
        ax[i].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
        ax[i].plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
        ax[i].plot(ds, "*-", ms=8, lw=1, label="demes")
        ax[i].set_title("s = " + str(curs) + ", demo = " + curdemo)
        ax[i].set_yscale('log')
        if np.logical_and(i == 2, j == 0):
            ax[i].legend();
            

# moments.Plotting.plot_1d_comp_Poisson(fs*8e-4, projData*2e-8)

for curs in [1e-3, 5e-3, 1e-2]:
    for curdemo in ["ooaSinglePop.yaml"]:
        os.chdir("/media/nathan/T7/BGSdemo/parsedooaSinglePopData")

        simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
        simData = simData[0].to_numpy()
        simData = moments.Spectrum(simData,data_folded=False) 
        projData = simData.project([proj_size])
        
        os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa")

        # test
        demo = demes.load(curdemo)
        demesdraw.tubes(demo);
        
        f, ancTime, ancNe = bFromDemes(pointMassPosition, u, curs, r, regionSize, focalPos, curdemo, tol)
                
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="getSizefun")
        ax.set_xlabel("Time in past")
        ax.set_ylabel("B(t)")
        ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
        ax.legend();
        
        cs = reversedCensusFun(curdemo, ancTime, ancNe)
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[cs(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="reversedCensusSize")
        ax.set_xlabel("Time in past")
        ax.set_ylabel("cs(t)")
        ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
        ax.legend();
        
        oldestEpoch, censusSize = getOldestEpoch(demo)
        ds = reversedCensusFun(curdemo, ancTime, censusSize)
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(np.arange(0,ancTime / 2 / censusSize, 0.001),[ds(t) for t in np.arange(0,ancTime / 2 / censusSize, 0.001)], "-", ms=8, lw=1, label="reversedCensusSize")
        ax.set_xlabel("Time in past")
        ax.set_ylabel("ds(t)")
        ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
        ax.legend();
        
        g = lambda t: [x * y for x,y in zip(f(t), cs(t))]
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[g(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="scaledPopSize")
        ax.set_xlabel("Time in past")
        ax.set_ylabel("g(t)")
        ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
        ax.legend();
        
        fs = moments.Demographics1D.snm([sample_size])
        fs.integrate(g, ancTime / 2 / ancNe)
        
        oldestEpoch, ancCensusSize = getOldestEpoch(demo)
        ds = reversedCensusFun(curdemo, ancTime, ancCensusSize)
        fs_neu = moments.Demographics1D.snm([proj_size])
        fs_neu.integrate(ds, ancTime / 2 / ancCensusSize)
        # fs_neu = moments.Demographics1D.snm([proj_size])
        # fs_neu.integrate(cs, ancTime / 2 / ancTime)

        # normalizing so singletons have freq 1, cause thats all I can think of right now
        fs = fs * 8 * 1e-8 * ancNe
        projData = projData * 1e-8
        fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize   
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(fs, ".-", ms=8, lw=1, label="BGS")
        ax.plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
        ax.plot(projData, "x-", ms=8, lw=1, label="fwdpy")
        ax.set_xlabel("Allele frequency")
        ax.set_ylabel("Density")
        ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
        ax.legend();
        