import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd
import demes
import demesdraw
from collections import defaultdict


# os.chdir("/media/nathan/T7/BGSdemo/parsedooaThreepopData")

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
            
def _get_root_Ne(g):
    # get root population and set Ne to root size
    for deme_id, preds in g.predecessors().items():
        if len(preds) == 0:
            root_deme = deme_id
            break
    Ne = g[root_deme].epochs[0].start_size
    return Ne
    
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

g = demo
sampled_demes=['OOA']
sample_sizes=[40]
sample_times=None
samples=None
unsampled_n=4
gamma=None
s=None
h=None
theta=1
u=None
reversible=False
L=1

import copy
import warnings

# new params
bgs_Ne = ancNe

def SFS_bgs(
    g,
    sampled_demes=None,
    sample_sizes=None,
    sample_times=None,
    samples=None,
    unsampled_n=4,
    gamma=None,
    s=None,
    h=None,
    theta=1,
    u=None,
    reversible=False,
    L=1,
):
    """
    Compute the SFS from a ``demes``-specified demographic model.
    ``demes`` is a package for specifying demographic models in a
    user-friendly, human-readable YAML format. This function
    automatically parses the demographic model and returns the SFS
    for the specified populations, sample sizes, and (optionally)
    sampling times.

    Selection and dominance can be specified as a single value for
    all populations, or on a per-deme basis using a dictionary
    mapping deme name to the coefficient (defaults can also be set
    if multiple demes have the same selection or dominance
    coefficient). The mutation rate can be given as either a
    per-base rate (possibly multiplied by the sequence length),
    or as a population size-scaled rate. If mutation rates are not
    given, the SFS is scaled by ``4*N_e``, so that multiplying the
    output SFS by ``u`` results in a properly scaled SFS.

    :param g: A ``demes`` DemeGraph from which to compute the SFS.
    :type g: :class:`demes.DemeGraph`
    :param sampled_demes: A list of deme IDs to take samples from. We can repeat
        demes, as long as the sampling of repeated deme IDs occurs at distinct
        times.
    :type sampled_demes: list of strings
    :param sample_sizes: A list of the same length as ``sampled_demes``,
        giving the sample sizes for each sampled deme.
    :type sample_sizes: list of ints
    :param sample_times: If None, assumes all sampling occurs at the end of the
        existence of the sampled deme. If there are
        ancient samples, ``sample_times`` must be a list of same length as
        ``sampled_demes``, giving the sampling times for each sampled
        deme. Sampling times are given in time units of the original deme graph,
        so might not necessarily be generations (e.g. if ``g.time_units`` is years)
    :type sample_times: list of scalars, optional
    :param unsampled_n: The default sample size of unsampled demes, which must be
        greater than or equal to 4.
    :type unsampled_n: int, optional
    :param gamma: The scaled selection coefficient(s), ``2*Ne*s``. Defaults to None,
        which implies neutrality. Can be given as a scalar value, in which case
        all populations have the same selection coefficient. Alternatively, can
        be given as a dictionary, with keys given as population names in the
        input Demes model. Any population missing from this dictionary will be
        assigned a selection coefficient of zero. A non-zero default selection
        coefficient can be provided, using the key ``_default``. See the Demes
        exension documentation for more details and examples.
    :type gamma: scalar or dict
    :param h: The dominance coefficient(s). Defaults to additivity (or genic
        selection). Can be given as a scalar value, in which case all populations
        have the same dominance coefficient. Alternatively, can be given as a
        dictionary, with keys given as population names in the input Demes model.
        Any population missing from this dictionary will be assigned a dominance
        coefficient of ``1/2`` (additivity). A different default dominance
        coefficient can be provided, using the key ``_default``. See the Demes
        exension documentation for more details and examples.
    :type h: scalar or dict
    :param theta: The scaled mutation rate(s), 4*Ne*u. When simulating under the
        infinite sites model (the default mutation model), ``theta`` should be given
        as a scalar value greater than zero. If it is not provided, it is computed
        using the input value of ``u`` as ``4*Ne*u``. If ``u`` is not provided, then
        the SFS is scaled by ``4*Ne``, and the user can recover a properly scaled SFS
        by multiplying it by ``u`` or ``u*L``. When simulating under the reversible
        mutation model (with ``reversible=True``), ``theta`` may be a list of length
        2 and both the forward and backward scaled mutation rates must be less
        than 1.
    :type theta: scalar or list of length 2
    :param u: The per-base mutation rate. When simulating under the infinite sites
        model (the default mutation model), ``u`` should be a scalar. When simulating
        under the reversible mutation model (with ``reversible=True``), ``u`` may
        be a list of length 2, and mutation rate(s) must be small enough so that
        the product of ``4*Ne*u`` is less than 1.
    :type u: scalar or list of length 2
    :param L: The effective sequence length, which may be used along with ``u`` to
        set the total mutation rate. Defaults to 1, and it must be 1 when using
        the reversible mutation model.
    :type L: scalar
    :return: A ``moments`` site frequency spectrum, with dimension equal to the
        length of ``sampled_demes``, and shape equal to ``sample_sizes`` plus one
        in each dimension, indexing the allele frequency in each deme from 0
        to n[i], where i is the deme index.
    :rtype: :class:`moments.Spectrum`
    """
    # could specify samples as a dict instead of sampled_demes and sample_sizes
    if samples is None:
        if sampled_demes is None or sample_sizes is None:
            raise ValueError(
                "must specify either samples (as a dict mapping demes to sample sizes,"
                " or specify both sampled_demes and sample_times"
            )
    else:
        if type(samples) is not dict:
            raise ValueError("samples must be a dict mapping demes to sample sizes")
        if sampled_demes is not None or sample_sizes is not None:
            raise ValueError(
                "if samples is given as dict, cannot "
                "specify sampled_demes or sample_sizes"
            )
        if sample_times is not None:
            raise ValueError("if samples is given as dict, cannot specify sample times")
        sampled_demes = list(samples.keys())
        sample_sizes = list(samples.values())

    if len(sampled_demes) != len(sample_sizes):
        raise ValueError("sampled_demes and sample_sizes must be same length")
    if sample_times is not None and len(sampled_demes) != len(sample_times):
        raise ValueError("sample_times must have same length as sampled_demes")
    for deme in sampled_demes:
        if deme not in g:
            raise ValueError(f"deme {deme} is not in demography")

    # we need to copy these to new variable names
    # so they don't get updated during optimization
    sampled_pops = copy.copy(sampled_demes)
    deme_sample_times = copy.copy(sample_times)

    if unsampled_n < 4:
        raise ValueError("unsampled_n must be greater than 3")

    sampled_deme_end_times = [g[d].end_time for d in sampled_pops]
    if deme_sample_times is None:
        deme_sample_times = sampled_deme_end_times

    for t, d in zip(deme_sample_times, sampled_pops):
        if t < g[d].end_time or t >= g[d].start_time:
            raise ValueError(f"sample time {t} is outside of deme {d}'s time span")

    # for any ancient samples, we need to add frozen branches
    # with this, all "sample times" are at time 0, and ancient sampled demes are frozen
    if np.any(np.array(deme_sample_times) != 0):
        g, sampled_pops, list_of_frozen_demes = _augment_with_ancient_samples(
            g, sampled_pops, deme_sample_times
        )
        deme_sample_times = [0 for _ in deme_sample_times]
    else:
        list_of_frozen_demes = []

    if g.time_units != "generations":
        g, deme_sample_times = _convert_to_generations(g, deme_sample_times)

    # if any sample sizes are less than unsample_n, we increase and project after
    sim_sample_sizes = []
    for d, n, t in zip(sampled_pops, sample_sizes, deme_sample_times):
        sim_sample_sizes.append(max(n, unsampled_n))
        if t < g[d].end_time or t >= g[d].start_time:
            raise ValueError("sample time for {deme} must be within its time span")

    # get reference Ne from demes model
    cs_Ne = _get_root_Ne(g)

    # if (unscaled) s is provided, convert into (scaled) gamma selection coefficients
    if s is not None:
        if gamma is not None:
            raise ValueError("Cannot specify both gamma and s")
        if isinstance(s, (int, float)):
            gamma = 2 * Ne * s
        elif type(s) is dict:
            gamma = {}
            for k, v in s.items():
                gamma[k] = 2 * Ne * v
        else:
            raise ValueError("Selection coefficient must be a scalar value or dict")

    # check selection and dominance inputs
    if gamma is not None:
        if "_default" in g:
            raise ValueError(
                "Cannot use `_default` as a deme name when selection is specified"
            )
        if isinstance(gamma, (int, float)):
            if not np.isfinite(gamma):
                raise ValueError("Selection coefficient must be a finite number")
        elif type(gamma) is dict:
            for k in gamma.keys():
                if k != "_default" and k not in g:
                    raise ValueError(f"Deme {k} in gamma, but {k} not in input graph")
        else:
            raise ValueError("Selection coefficient must be a scalar value or dict")
    if h is not None:
        if type(h) is dict:
            for k in h.keys():
                if k != "_default" and k not in g:
                    raise ValueError(f"Deme {k} in h, but {k} not in input graph")

    # set up the mutation rates as needed
    if theta is None:
        if u is None:
            u = 1
        if isinstance(u, (int, float)):
            theta = 4 * Ne * u * L
        else:
            if np.ndim(u) != 1 or len(u) != 2:
                raise ValueError(
                    "Mutation rates must be a list of length 2 when using "
                    "the reversible mutation model"
                )
            theta = [4 * Ne * u[0], 4 * Ne * u[1]]
    else:
        if u is not None:
            raise ValueError("Only one of u or theta may be specified")
        if isinstance(theta, (int, float)):
            theta *= L
        else:
            if np.ndim(theta) != 1 or len(theta) != 2:
                raise ValueError(
                    "Mutation rates must be a list of length 2 when using "
                    "the reversible mutation model"
                )
            theta[0] *= L
            theta[1] *= L

    # if a scalar, must be positive; if list-like, must be length 2 and both positive
    if not reversible:
        if not isinstance(theta, (int, float)):
            raise ValueError(
                "Mutation rate must be a scalar value for the default ISM model"
            )
        if theta <= 0:
            raise ValueError("Mutation rate must be positive")
    if reversible:
        if L != 1:
            raise ValueError(
                "Sequence length L must be 1 when using the reversible mutation model"
            )
        if isinstance(theta, (int, float)):
            theta = [theta, theta]
        if theta[0] <= 0 or theta[1] <= 0:
            raise ValueError("Mutation rates must be positive")
        if theta[0] >= 1 or theta[1] >= 1:
            raise ValueError("Mutation rates too large for reversible mutation model")

    # get the list of demographic events from demes, which is a dictionary with
    # lists of splits, admixtures, mergers, branches, and pulses
    demes_demo_events = g.discrete_demographic_events()

    # get the dict of events and event times that partition integration epochs, in
    # descending order. events include demographic events, such as splits and
    # mergers and admixtures, as well as changes in population sizes or migration
    # rates that require instantaneous changes in the size function or migration matrix.
    # get the list of demes present in each epoch, as a dictionary with non-overlapping
    # adjoint epoch time intervals
    demo_events, demes_present = _get_demographic_events(
        g, demes_demo_events, sampled_pops
    )

    for epoch, epoch_demes in demes_present.items():
        if len(epoch_demes) > 5:
            raise ValueError(
                f"Moments cannot integrate more than five demes at a time. "
                f"Epoch {epoch} has demes {epoch_demes}."
            )


    # get the list of size functions, migration matrices, and frozen attributes from
    # the deme graph and event times, matching the integration times
    # nu_funcs, mig_mats, Ts, frozen_pops = _get_integration_parameters(
    #     g, demes_present, list_of_frozen_demes, Ne=cs_Ne
    # )
    
    # nu_funcs_bgs, mig_mats_bgs, Ts_bgs, frozen_pops_bgs = _get_integration_parameters(
    #     g, demes_present, list_of_frozen_demes, Ne=ancNe
    # )
    
    # nu_funcs_test, mig_mats_test, Ts_test, frozen_pops_test 
    nu_funcs, mig_mats, Ts, frozen_pops = _get_integration_parameters_bgs(
        g, demes_present, list_of_frozen_demes, bgs_Ne
    )
    # get the sample sizes within each deme, given sample sizes
    deme_sample_sizes = _get_deme_sample_sizes(
        g,
        demo_events,
        sampled_pops,
        sim_sample_sizes,
        demes_present,
        unsampled_n=unsampled_n,
    )

    # compute the SFS
    fs = _compute_sfs(
        demo_events,
        demes_present,
        deme_sample_sizes,
        nu_funcs,
        mig_mats,
        Ts,
        frozen_pops,
        theta=theta,
        gamma=gamma,
        h=h,
        reversible=reversible,
    )

    fs = _reorder_fs(fs, sampled_pops)

    # project down to desired sample sizes, if needed
    fs = fs.project(sample_sizes)
    # simplify pop id name if ancient sample at end time of that deme
    for ii, pid in enumerate(fs.pop_ids):
        if "_sampled_" in pid:
            p, t = pid.split("_sampled_")
            t = float(t.replace("_", "."))
            if t == sampled_deme_end_times[ii]:
                fs.pop_ids[ii] = p

    return fs

bgs_Ne = ancNe
frozen_list = list_of_frozen_demes

def _get_integration_parameters_bgs(g, demes_present, frozen_list, bgs_Ne, cs_Ne=None):
    """
    Returns a list of size functions, migration matrices, integration times,
    and lists frozen demes.
    """
    nu_funcs = []
    integration_times = []
    migration_matrices = []
    frozen_demes = []

    if cs_Ne is None:
        cs_Ne = _get_root_Ne(g)
    else:
        if cs_Ne != _get_root_Ne(g):
            warnings.warn(
                "Input cs_Ne is different from root population initial size, "
                "subsequent population size scaling may be incorrect"
            )
    
    T_elapsed = 0
    for interval, live_demes in sorted(demes_present.items())[::-1]:
        # get intergration time for interval
        T = (interval[0] - interval[1]) / 2 / bgs_Ne
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
        # nu_func = _make_nu_func(sizes, T, cs_Ne)
        nu_func = _make_nu_func_bgs(sizes, T, cs_Ne, T_elapsed)
        T_elapsed += T

        nu_funcs.append(nu_func)
        # get migration matrix for interval
        mig_mat = np.zeros((len(live_demes), len(live_demes)))
        for ii, d_from in enumerate(live_demes):
            for jj, d_to in enumerate(live_demes):
                if d_from != d_to:
                    m = _migration_rate_in_interval(g, d_from, d_to, interval)
                    mig_mat[jj, ii] = 2 * bgs_Ne * m
        migration_matrices.append(mig_mat)

    return nu_funcs, migration_matrices, integration_times, frozen_demes

def _make_nu_func_bgs(sizes, T, Ne, T_elapsed):
    """
    Given the sizes at start and end of time interval, and the size function for
    each deme, along with the integration time and reference Ne, return the
    size function that gets passed to the moments integration routines.
    """
    if np.all([s[-1] == "constant" for s in sizes]):
        # all constant
        if T == 0:
            nu_func = [s[0] / Ne for s in sizes]
        else:
            nu_funcs_separated = []
            for s in sizes:
                assert s[0] == s[1]
                nu_funcs_separated.append(lambda t, N0=s[0]: (N0 / Ne) * f(T_elapsed + t)[0])
                
            def nu_func(t):
                return [nu(t) for nu in nu_funcs_separated]
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

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(np.arange(0,Ts_test[1], 0.001),[nu_funcs_test[1](t) for t in np.arange(0,Ts_test[1], 0.001)], "-", ms=8, lw=1, label="nu_func1")
ax.plot(np.arange(0,Ts_test[2], 0.001),[nu_funcs_test[2](t) for t in np.arange(0,Ts_test[2], 0.001)], "-", ms=8, lw=1, label="nu_func2")
# ax.plot(np.arange(0,T, 0.001),[f(t) for t in np.arange(0,T, 0.001)], "-", ms=8, lw=1, label="f")
ax.set_xlabel("Time in past")
ax.set_ylabel("value")
ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
ax.legend();

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(np.arange(0,T, 0.001),[test(t) for t in np.arange(0,T, 0.001)], "-", ms=8, lw=1, label="nu_func")
# ax.plot(np.arange(0,T, 0.001),[f(t) for t in np.arange(0,T, 0.001)], "-", ms=8, lw=1, label="f")
ax.set_xlabel("Time in past")
ax.set_ylabel("value")
ax.set_title("s = " + str(curs) + ", demo = " + curdemo)
ax.legend();

# todo find out what happens when I redefine these function names???
s = curs
curN = censusSize
# positions = pointMassPosition
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot([B(pointMassPosition, u, s, t, r, 1e4, 5e5) for t in range(int(10 * censusSize))], "-", ms=8, lw=1, label="Neutral")
ax.set_xlabel("Time in past")
ax.set_ylabel("B(t)")
ax.legend();

curs = 1e-3
curdemo = 'ooaTwoPop.yaml'
# for curs in [1e-3, 5e-3, 1e-2]:
#     for curN in [1e3, 5e3, 1e4]:
fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
fig.subplots_adjust(hspace = .25)

for i in range(3):
    for j in range(1):
        curs = [1e-3, 5e-3, 1e-2][i]
        curdemo = ["ooaTwoPop.yaml"][j]
        
        # os.chdir("/media/nathan/T7/BGSdemo/parsedooaThreepopData")

        # simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
        # simData = simData[0].to_numpy()
        # simData = moments.Spectrum(simData,data_folded=False) 
        # projData = simData.project([proj_size])
        
        os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/twopop")

        # test
        demo = demes.load(curdemo)
        if demo.time_units != "generations":
            demo = demo.in_generations()
        demesdraw.tubes(demo);
        
        oldestEpoch, censusSize = getOldestEpoch(demo)

        f, ancTime, ancNe = bFromDemes(pointMassPosition, u, curs, r, regionSize, focalPos, curdemo, tol)
        
        SFS_bgs(
           demo,
           sampled_demes=["CEU"],
           sample_sizes=[sample_size],
           theta = 1
       )
       
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
        