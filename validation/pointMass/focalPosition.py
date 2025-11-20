import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd

os.chdir("/media/nathan/T7/BGSdemo/parsedEquilSelData")

# hardcoding some parameters
u = 1e-8
r = 1e-8
L = 1e6
regionSize = 1e4
tol = 1e-3
sample_size = 500
proj_size = 40

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
fig, ax = plt.subplots(3, 3, figsize=(16, 8), sharex=False, sharey=False)
fig.text(0.5, 0.04, 'time', ha='center')
fig.text(0.04, 0.5, 'B(t)', va='center', rotation='vertical')
fig.subplots_adjust(hspace = .25)

for i in range(3):
    for j in range(3):
        curs = [1e-3, 5e-3, 1e-2][i]
        curN = [1e3, 5e3, 1e4][j]
        ax[i,j].set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        for pos in  range(int(3e5),int(5.5e5),int(5e4)):# range(int(4.5e5),int(5.1e5),int(1e4)): 
            f, ancTime, ancNe = getSizeFun(pointMassPosition, u, curs, r, regionSize, pos, curN, tol)
            
            ax[i,j].plot(np.arange(0,ancTime / 2 / ancNe, 0.001),
                    [f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], 
                    "-", ms=8, lw=1, label=str(pos))
            if np.logical_and(i == 2, j == 2):
                ax[i,j].legend();
            

        
