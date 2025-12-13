import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd

# hardcoding some parameters
tol = 1e-4
sample_size = 400

def pointMassContribution(pos, scaledu, s, t, r, focalPos):
    return - scaledu / s * (s / (r * abs(pos-focalPos) + s) * (1 - math.exp(- r * abs(pos-focalPos) * t - s * t)))**2

def B(positions, u, s, t, r, regionSize, focalPos):
    scaledu = u * regionSize
    return math.exp(sum([pointMassContribution(pos, scaledu, s, t, r, focalPos) for pos in positions]))

# im making a factor of 2 error here somewhere # am I? did i fix and forget to remove?
def rescaledPointMassContribution(pos, scaledu, s, t, r, focalPos, ancestralNe, ancTime):
    return - scaledu / s * (s / (r * abs(pos-focalPos) + s) * (1 - math.exp(- r * abs(pos-focalPos) * (ancTime - t * 2 * ancestralNe) - s * (ancTime - t * 2 * ancestralNe))))**2

def getSizeFun(positions, u, s, r, regionSize, focalPos, censusSize, tol):
    if any([type(u) == list, type(r) == list]):
        getSizeFun_maps()
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

def getSizeFun_maps(positions, u, s, r, focalPos, censusSize, tol):
    scaledu = [z * (y - x) for x,y,z in u]
    recDist = [getRecDist(pos, r, focalPos) for pos in positions]
    testFun = [B_maps(scaledu, s, t, recDist) for t in range(0,int(10 * censusSize),int(censusSize/10))]  
    diffs = [testFun[i+1] - testFun[i] for i in range(len(testFun)-1)]
    ancTime = next((i for i,x in enumerate(diffs) if abs(x) < tol), None)
    ancTime = censusSize / 10 * ancTime 
    ancB = B_maps(scaledu, s, ancTime, recDist)
    ancNe = ancB * censusSize 
    return(lambda t: [math.exp(sum([rescaledPointMassContribution_map(u, s, t, r, ancNe, ancTime) for u,r in zip(scaledu, recDist)])) / ancB], ancTime, ancNe)

def rescaledPointMassContribution_map(scaledu, s, t, r, ancestralNe, ancTime):
    return - scaledu / s * (s / (r + s) * (1 - math.exp(- r * (ancTime - t * 2 * ancestralNe) - s * (ancTime - t * 2 * ancestralNe))))**2

    
def B_maps(scaledu, s, t, recDist):
    return math.exp(sum([pointMassContribution_map(u, s, t, r) for u,r in zip(scaledu, recDist)]))

def pointMassContribution_map(u, s, t, r):
    return - u / s * (s / (r + s) * (1 - math.exp(- r * t - s * t)))**2

def getRecDist(pos, r, focalPos):
    left = min(pos, focalPos)
    right = max(pos, focalPos)
    leftIndex = [i for i,x in enumerate(r) if x[1] > left and x[0] < left][0] # todo could probably speed this up
    rightIndex = [i for i,x in enumerate(r) if x[1] > right and x[0] < right][0]
    rleft = r[leftIndex]
    rright = r[rightIndex]
    intermediate = r[(leftIndex + 1):(rightIndex-1)]    
    bigR = [(y - x) * z for x,y,z in intermediate]
    bigR.append((rleft[1]-left) * rleft[2]) # todo need to deal with possibility left and right are inside the same window
    bigR.append((right-rright[0]) * rright[2])
    bigR = np.sum(bigR)
    return (1 - np.exp(- 2 * bigR)) / 2

# f, ancTime, ancNe = getSizeFun_maps(midPoints, exonMutMap, curs, recMap, focalPos, curN, tol)


# s = curs
# censusSize = curN
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot([t for t in range(0,int(3*censusSize),10)],[B_maps(scaledu, s, t, recDist) for t in range(0,int(3 * censusSize), 10)], "-", ms=8, lw=1, label="B(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend();
 


# need exonMutMap, need rec map
def read_rec_map(bed_path):
    intervals = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split()
            start = int(fields[1])
            end   = int(fields[2])
            rate  = float(fields[3])
            intervals.append([start, end, rate])
    return intervals

def read_mut_rates(bed_mut):
    muts = []
    with open(bed_mut) as f:
        for line in f:
            if line.startswith("chrom") or line.strip() == "":
                continue
            fields = line.strip().split(",")
            start = fields[1]
            end = fields[2]
            mu = fields[4]
            start, end, mu = int(start), int(end), float(mu)
            muts.append([start, end, mu])
    return muts

def read_exon_map(bed_exons):
    exons = []
    with open(bed_exons) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            exons.append([start, end])
    return exons

def simplify_rate_map(data):
    merged  = []
    
    prev_start = None
    prev_end = None
    prev_rate = None
        
    for x in data:
        start,end,rate = x

        # First valid line initializes the merge
        if prev_start is None:
            prev_start = start
            prev_end = end
            prev_rate = rate
            continue
        
        # adjacent, and same rate
        if (start == prev_end
            and rate == prev_rate      # exact match
        ):
            # Extend the previous interval
            prev_end = end
            continue

        # Otherwise: flush previous, start new
        merged.append([prev_start, prev_end, prev_rate])

        prev_start = start
        prev_end = end
        prev_rate = rate
    
    # Flush last interval at EOF
    if prev_start is not None:
        merged.append([prev_start, prev_end, prev_rate])
        
    return merged

def make_exon_only_mutmap(mutMap, exonMap):
    exonMutMap = []
    for m_start, m_end, mu in mutMap:
        for e_start, e_end in exonMap:
            inter = intersect(m_start, m_end, e_start, e_end)
            if inter is not None:
                beg, end = inter
                exonMutMap.append([beg, end, mu])
    totalRate = get_total_rate(exonMutMap)
    return exonMutMap, totalRate

def intersect(a_start, a_end, b_start, b_end):
    beg = max(a_start, b_start)
    end = min(a_end, b_end)
    if beg < end:
        return beg, end
    return None

def get_total_rate(xMap):
    tmp = [rate * (end - start) for start, end, rate in xMap]
    return np.sum(tmp)

def get_max_positions(recMap, mutMap):
    return max(recMap[-1][1], mutMap[-1][1])

os.chdir("/media/nathan/T7/BGSdemo/humanData")
recName = "YRI_recombination_map_hg38_chr_22.bed"
mutName = "roulette_tbl_chr22.csv"
exonName = "exons_chr22.bed"
recMap = read_rec_map(recName)
mutMap = read_mut_rates(mutName)
exonMap = read_exon_map(exonName)

recMap = simplify_rate_map(recMap)
mutMap = simplify_rate_map(mutMap)

# todo fix demo for B ancTime > oldest demo event
# todo change test fun
# todo what if a region is too large
exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap) # todo fix slow

L = get_max_positions(recMap, mutMap)
focalPos = L/2

regionSize = [y-x for x,y,z in exonMutMap]
midPoints = [(x + y)/2 for x,y,z in exonMutMap]

# np.savetxt("exonMutMap.csv", exonMutMap, delimiter = ",")

focalPos = 50270000
# f, ancTime, ancNe = getSizeFun_maps(midPoints, exonMutMap, s, recMap, focalPos, censusSize, tol)
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="getSizefun")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.set_title("s = " + str(s) + ", N = " + str(int(censusSize)))
# ax.legend();

# for curs in [1e-3, 5e-3, 1e-2]:
#     for curN in [1e3, 5e3, 1e4]:
os.chdir("/media/nathan/T7/BGSdemo/parsedequilHuman")
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
        f, ancTime, ancNe = getSizeFun_maps(midPoints, exonMutMap, curs, recMap, focalPos, curN, tol)
        
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
        
             
        # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        # ax.plot(fs, "-", ms=8, lw=1, label="BGS")
        # ax.plot(fs_neu, "-", ms=8, lw=1, label="neutral")
        # # ax.plot(projData, "x-", ms=8, lw=1, label="fwdpy")
        # ax.set_xlabel("Allele frequency")
        # ax.set_ylabel("Density")
        # ax.set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        # ax.set_yscale("log")
        # ax.annotate("B=" + str(round(fs.pi()/fs_neu.pi(),3)),xy=(sample_size / 10, max(fs[1:-1]) * 0.8),size="x-large")
        # ax.legend();
        
        
        ax[i,j].plot(fs, "-", ms=8, lw=1, label="BGS")
        # ax[i,j].plot(projData, "-", ms=8, lw=1, label="fwdpy")
        ax[i,j].plot(fs_neu, "-", ms=8, lw=1, label="SNM")
        ax[i,j].set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        ax[i,j].set_yscale('log')
        ax[i,j].annotate("B=" + str(round(fs.pi()/fs_neu.pi(),3)),xy=(sample_size / 10, max(fs[1:-1]) * 0.8),size="x-large")
        if np.logical_and(i == 2, j == 2):
            ax[i,j].legend();

x = np.arange(10)
y = np.array([5, 3, 4, 2, 7, 5, 4, 6, 3, 2])

fig, ax = plt.subplots() # Get the figure and axes objects
ax.plot(x, y, 'o-')

for i, j in zip(x, y):
    ax.annotate(str(j), xy=(i, j), xytext=(i + 0.2, j + 0.8),
                arrowprops=dict(facecolor='black', shrink=0.05))

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_title("Plot with annotated values")
plt.show()

# moments.Plotting.plot_1d_comp_Poisson(fs*8e-4, projData*2e-8)

for curs in [1e-3, 5e-3, 1e-2]:
    for curN in [1e3, 5e3, 1e4]:
        simData = pd.read_csv(str(curs) + "_" + str(int(curN)) + ".csv", header = None)
        simData = simData[0].to_numpy()
        simData = moments.Spectrum(simData,data_folded=False) / 11 / 1000
        projData = simData.project([sample_size])
        fs = moments.Demographics1D.snm([sample_size])
        
        # fs_neu = moments.Spectrum(moments.Demographics1D.snm([sample_size]))
        # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        # ax.plot(fs, ".-", ms=8, lw=1, label="neg sel")
        # ax.plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
        # ax.set_xlabel("Allele frequency")
        # ax.set_ylabel("Density")
        # ax.set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        # ax.legend();
        
        f, ancTime, ancNe = getSizeFun(pointMassPosition, u, curs, r, regionSize, focalPos, curN, tol)
        
        # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        # ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="getSizefun")
        # ax.set_xlabel("Time in past")
        # ax.set_ylabel("B(t)")
        # ax.set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        # ax.legend();
        
        # gamma_f = lambda t: [- 2 * ancNe * curs * x for x in f(t)]
        # fs = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sample_size, gamma = - 2 * ancNe * curs))
        # fs.integrate(f, ancTime / 2 / ancNe, gamma = gamma_f, h = 1/2, adapt_dt=True)
        
        gamma_f = lambda t: [- 2 * ancNe * curs * x for x in f(t)]
        fs = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sample_size, gamma = - 2 * ancNe * curs))
        fs.integrate(f, ancTime / 2 / ancNe, gamma = - 2 * ancNe * curs, h = 1/2, adapt_dt=True)
        
        fs_indep = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sample_size, gamma = - 2 * curN * curs))

        if curwi == "wi5e4.csv":
            fs = fs *  4 * 2 * 5e4 * 1e-8 * ancNe
            fs_indep = fs_indep * 4 * 2 * 5e4 * 1e-8 * curN
            projData = projData 
        if curwi == "wi1e5.csv":
            fs = fs *  4 * 2 * 1e5 * 1e-8 * ancNe
            fs_indep = fs_indep * 4 * 2 * 1e5 * 1e-8 * curN
            projData = projData 
        
        
        # fs = fs / fs[1]
        # fs_indep = fs_indep / fs_indep[1]
        # projData = projData / projData[1]
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(fs, ".-", ms=8, lw=1, label="BGS")
        ax.plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
        ax.plot(projData, "x-", ms=8, lw=1, label="fwdpy")
        ax.set_xlabel("Allele frequency")
        ax.set_ylabel("Density")
        ax.set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        ax.annotate("B=" + str(round(fs.pi()/fs_neu.pi(),3)),xy=(max(fs[1:-1]),sample_size * 0.9))
        ax.legend();
        
        