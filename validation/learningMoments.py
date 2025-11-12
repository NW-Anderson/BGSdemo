import moments
import numpy as np
import math
import matplotlib.pylab as plt

# hardcoding some parameters
u = 1e-8
r = 1e-8
s = 1e-3
censusSize = 1e4

# positions of point masses.
# split 1Mb into 10kb regions, each point mass lies at the center
# of these regions 
pointMassPosition = range(int(5e3),int(100.5e4),int(1e4))
print(pointMassPosition[0]) # first point mass position
print(pointMassPosition[len(pointMassPosition)- 1]) # final position
    
# this comes from the final line in the eq before eq (3) in ND13
# recombination fraction is computed r |x_pointMass - x_focal|
# in general, s and u can be altered to reflect the local
# constraint and mutation rate
def pointMassContribution(pos, scaledu, s, t, r, focalPos):
    return - scaledu / s * (s / (r * abs(pos-focalPos) + s) * (1 - math.exp(- r * abs(pos-focalPos) * t - s * t)))**2

# wrapper for the above function, sums the contribution 
# of each point mass and puts them in the exponent
# this is what is multiplied by N(t) to compute Ne(t)
def B(positions, u, s, t, r, regionSize, focalPos):
    scaledu = u * regionSize
    return math.exp(sum([pointMassContribution(pos, scaledu, s, t, r, focalPos) for pos in positions]))

# plotting B(t) for a site at the center of the 1Mb region
testFun = [B(pointMassPosition, u, s, t, r, 1e4, 5e5) for t in range(int(1e4))]
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(testFun, "-", ms=8, lw=1, label="Neutral")
ax.set_xlabel("Time in past")
ax.set_ylabel("B(t)")
ax.legend();

# need to rescale and reverse time
# B(t) needs to be in terms of the ancestral population size (i.e >1)
# first rescale so that B(t) shows growth from an ancestral scale of 1
# instead of a reduction of Ne in the past
rescaledFun = [x / testFun[len(testFun) - 1] for x in testFun] 
# reversing time
rescaledFun.reverse() 
# computing the ancestral effective population size
ancestralNe = censusSize * testFun[len(testFun) - 1]
# rescaled time generations -> 2 Ne(anc) generations
tms = [t / 2 / ancestralNe for t in range(int(1e4))]

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(tms, rescaledFun, "-", ms=8, lw=1, label="Neutral")
ax.set_xlabel("Time since equilibrium")
ax.set_ylabel("B(t)")
ax.legend();

# this is our size function, but it needs to be automated
# for moments.integrate we need a function that takes in t and only t
# and returns the corresponding value of rescaledFun

# same as above but time is now in units of 2 Ne(anc) generations
# and time is reversed relative to ancTime (when B(t) is flat) so we start at equil
def rescaledPointMassContribution(pos, scaledu, s, t, r, focalPos, ancestralNe, ancTime):
    return - scaledu / s * (s / (r * abs(pos-focalPos) + s) * (1 - math.exp(- r * abs(pos-focalPos) * (ancTime - t * 2 * ancestralNe) - s * (ancTime - t * 2 * ancestralNe))))**2

# this function automates the earlier work going from testFun to rescaledFun. only valid for times in [0,ancTime]
regionSize = 1e4
tol = 1e-3
focalPos = 5e5
def getSizeFun(positions, u, s, r, regionSize, focalPos, censusSize, tol):
    # rescale u for each region, eventually region size will be a vector, u could also change??
    scaledu = u * regionSize
    # finding B(t) at several time points
    testFun = [B(positions, u, s, t, r, regionSize, focalPos) for t in range(0,int(censusSize),int(censusSize/10))]
    # finding when equilibrium is reached
    diffs = [testFun[i+1] - testFun[i] for i in range(len(testFun)-1)]
    ancTime = next((i for i,x in enumerate(diffs) if abs(x) < tol), None)
    # time to equil B(t) in generations
    # ancTime = 2 * censusSize / 10 * ancTime
    ancTime = censusSize / 10 * ancTime
    # need to think about what happens if ancTime exceeds censusSize (the largest time considered above)
    ancB = B(positions, u, s, ancTime, r, regionSize, focalPos) 
    ancNe = ancB * censusSize 
    return(lambda t: [math.exp(sum([rescaledPointMassContribution(pos, scaledu, s, t, r, focalPos, ancNe, ancTime) for pos in positions])) / ancB], ancTime, ancNe)

f, ancTime, ancNe = getSizeFun(pointMassPosition, u, s, r, regionSize, focalPos, censusSize, tol)

# confirming we start at the equilibrium size
f(0)
# confirming we recover the census size at current
f(ancTime / 2 / ancNe) [0] * ancNe
# going too far in the future (beyond the present) results in extinction
f(ancTime / ancNe)

# plotting size function
# note the difference in where time starts using f and rescaledFun
# ancTime is computed automatically in the former
# chosen arbitrarily censusSize generations in the past in the latter
test = [f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)]
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(np.arange(ancTime / 2 / ancNe, 1e4 / 2 / ancNe, 0.001), test, "-", ms=8, lw=1, label="getSizefun")
ax.plot(tms, rescaledFun, "-", ms=8, lw=1, label="rescaledFun")
ax.set_xlabel("Time in past")
ax.set_ylabel("B(t)")
ax.legend();

# getting sfs under several conditions    
fs = moments.Demographics1D.snm([10])
fs.integrate(f, ancTime / 2 / ancNe)
fs_neu = moments.Demographics1D.snm([10])
g, gancTime, gancNe = getSizeFun(pointMassPosition, u, s, 4e-8, regionSize, focalPos, censusSize, tol)
gs = moments.Demographics1D.snm([10])
gs.integrate(g, gancTime / 2 / gancNe)


fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(moments.Spectrum(fs) * ancNe / censusSize, ".-", ms=8, lw=1, label="BGS: r=1e-8")
ax.plot(moments.Spectrum(gs)* gancNe / censusSize, "+-", ms=8, lw=1, label="BGS: r=4e-8")
ax.plot(fs_neu, "x-", ms=8, lw=1, label="Neutral")
ax.set_xlabel("Allele frequency")
ax.set_ylabel("Density")
ax.legend();

# comparison to constant census or consant ancestral Ne to size function
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(moments.Spectrum(gs) * gancNe / censusSize, "+-", ms=8, lw=1, label="BGS: r=4e-8")
ax.plot(fs_neu, "x-", ms=8, lw=1, label="Neutral, Census")
ax.plot(fs_neu * gancNe / censusSize, "+-", ms=8, lw=1, label="Neutral, Ne(anc)")
ax.set_xlabel("Allele frequency")
ax.set_ylabel("Density")
ax.legend();
