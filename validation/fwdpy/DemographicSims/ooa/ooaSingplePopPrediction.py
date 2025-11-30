import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd
import demes
import demesdraw

os.chdir("/media/nathan/T7/BGSdemo/parsedooaSinglePopData")

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

def getSizeFun(positions, u, s, r, regionSize, focalPos, censusSize, tol):
    # rescale u for each region, eventually region size will be a vector, u could also change??
    scaledu = u * regionSize
    # finding B(t) at several time points
    # TODO need to do something to deal with longer times to reach equil
    testFun = [B(positions, u, s, t, r, regionSize, focalPos) for t in range(0,int(10 * censusSize),int(censusSize/10))]
    # finding when equilibrium is reached
    diffs = [testFun[i+1] - testFun[i] for i in range(len(testFun)-1)]
    ancTime = next((i for i,x in enumerate(diffs) if abs(x) < tol), None)
    # time to equil B(t) in generations
    ancTime = censusSize / 10 * ancTime 
    # need to think about what happens if ancTime exceeds censusSize (the largest time considered above)
    ancB = B(positions, u, s, ancTime, r, regionSize, focalPos) 
    ancNe = ancB * censusSize 
    return(lambda t: [math.exp(sum([rescaledPointMassContribution(pos, scaledu, s, t, r, focalPos, ancNe, ancTime) for pos in positions])) / ancB], ancTime, ancNe)

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

# s = curs
curN = censusSize
# positions = pointMassPosition
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot([B(pointMassPosition, u, s, t, r, 1e4, 5e5) for t in range(int(10 * curN))], "-", ms=8, lw=1, label="Neutral")
ax.set_xlabel("Time in past")
ax.set_ylabel("B(t)")
ax.legend();

# for curs in [1e-3, 5e-3, 1e-2]:
#     for curN in [1e3, 5e3, 1e4]:
fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
fig.subplots_adjust(hspace = .25)

for i in range(3):
    for j in range(1):
        curs = [1e-3, 5e-3, 1e-2][i]
        curdemo = ["ooaSinglePop.yaml"][j]
        
        os.chdir("/media/nathan/T7/BGSdemo/parsedooaSinglePopData")

        simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
        simData = simData[0].to_numpy()
        simData = moments.Spectrum(simData,data_folded=False) 
        projData = simData.project([proj_size])
        
        os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa")

        # test
        demo = demes.load(curdemo)
        # demesdraw.tubes(demo);
        
        f, ancTime, ancNe = bFromDemes(pointMassPosition, u, curs, r, regionSize, focalPos, curdemo, tol)
        
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
        