import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd

os.chdir("/media/nathan/T7/BGSdemo/equilAFS")

# hardcoding some parameters
u = 1e-8
r = 1e-8
L = 1e6
regionSize = 1e4
tol = 1e-3
focalPos = 5e5
sample_size = 40

# positions of point masses.
# split 1Mb into 10kb regions, each point mass lies at the center
# of these regions 
pointMassPosition = range(int(5e3),int(100.5e4),int(1e4))

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

# s = curs
# censusSize = curN
# positions = pointMassPosition
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot([B(pointMassPosition, u, s, t, r, 1e4, 5e5) for t in range(int(10 * curN))], "-", ms=8, lw=1, label="Neutral")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend();

# for curs in [1e-3, 5e-3, 1e-2]:
#     for curN in [1e3, 5e3, 1e4]:
fig, ax = plt.subplots(3, 3, figsize=(16, 8), sharex=True, sharey=False)
fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
fig.subplots_adjust(hspace = .25)

for i in range(3):
    for j in range(3):
        curs = [1e-3, 5e-3, 1e-2][i]
        curN = [1e3, 5e3, 1e4][j]
        simData = pd.read_csv(str(curs) + "_" + str(int(curN)) + ".csv", header = None)
        simData = simData[0].to_numpy()
        simData = moments.Spectrum(simData,data_folded=False)
        projData = simData.project([sample_size])
        fs = moments.Demographics1D.snm([sample_size])
        f, ancTime, ancNe = getSizeFun(pointMassPosition, u, curs, r, regionSize, focalPos, curN, tol)
        
        # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        # ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="getSizefun")
        # ax.set_xlabel("Time in past")
        # ax.set_ylabel("B(t)")
        # ax.set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        # ax.legend();
                
        fs.integrate(f, ancTime / 2 / ancNe)
        
        fs_neu = moments.Demographics1D.snm([sample_size])
        # normalizing so singletons have freq 1, cause thats all I can think of right now
        # fs = fs / fs[1]
        # projData = projData / projData[1]
        
        # normalizing based on N and mu 
        # should it be ancNe
        fs = fs * 8 * 1e-8 * ancNe
        projData = projData * 1e-8
        fs_neu = fs_neu * 8 * 1e-8 * curN
        
        # todo i think the correct thing is tp divide the previous lines by 2
        # regular theta for a single site and span normalized projData
        
        ax[i,j].plot(fs, ".-", ms=8, lw=1, label="BGS")
        ax[i,j].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
        ax[i,j].plot(fs_neu, "+-", ms=8, lw=1, label="SNM")
        ax[i,j].set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        if np.logical_and(i == 2, j == 2):
            ax[i,j].legend();

# moments.Plotting.plot_1d_comp_Poisson(fs*8e-4, projData*2e-8)
        