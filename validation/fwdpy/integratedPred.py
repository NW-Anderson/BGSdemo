import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd

# need to account for maps vs constant values
# selected vs neutral
# demes vs nu_func split into two functions?

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
    # startTime = datetime.now()
    exonMutMap = []
    i = 0
    for e_start, e_end in exonMap:
        # if e_start == prob:
        #     break
        inter = None
        i -= 1 
        while inter is None:
            i += 1
            m_start, m_end, mu = mutMap[i]
            inter = intersect(m_start, m_end, e_start, e_end)
        i -= 1
        while inter is not None:
            i += 1
            m_start, m_end, mu = mutMap[i]
            inter = intersect(m_start, m_end, e_start, e_end)
            if inter is not None:
                beg, end = inter
                exonMutMap.append([beg, end, mu])
        i -= 1 # want to repeat same window of mutmap for next exon

    # endTime = datetime.now()
    # endTime - startTime
    # test = exonMutMap
    
    # startTime = datetime.now()
    # exonMutMap = []
    # for m_start, m_end, mu in mutMap:
    #     for e_start, e_end in exonMap:
    #         inter = intersect(m_start, m_end, e_start, e_end)
    #         if inter is not None:
    #             beg, end = inter
    #             exonMutMap.append([beg, end, mu])
    # bench = exonMutMap
    # endTime = datetime.now()
    # endTime - startTime            
    
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

# todo still need to split but none of them are that big
def combine_and_split_regions(exonMutMap, targetSize = 1e4):
    combined = []
    tmp = []
    tmp_start = exonMutMap[0][0]
    for start, stop, mu in exonMutMap:
        if (abs(start - tmp_start) < targetSize):
            tmp.append([start, stop, mu])
        else:
            new_start = tmp[0][0]
            new_end = tmp[-1][1]
            new_L = new_end - new_start
            old_U = np.sum([(y-x)*z for x,y,z in tmp])
            new_mu = old_U / new_L

            combined.append([new_start, new_end, new_mu])
            tmp = [[start, stop, mu]]
            tmp_start = start
         
    new_start = tmp[0][0]
    new_end = tmp[-1][1]
    new_L = new_end - new_start
    old_U = np.sum([(y-x)*z for x,y,z in tmp])
    new_mu = old_U / new_L
    combined.append([new_start, new_end, new_mu])
    
    return combined

def getRecDist(pos, r, focalPos):
    left = min(pos, focalPos)
    right = max(pos, focalPos)
    
    # startTime = datetime.now()
    i = 0
    done = True
    leftCount = 0
    rightCount = 0
    while done:
        x = r[i]
        if x[1] >= left and x[0] <= left:
            leftCount += 1
            leftIndex = i
        if x[1] >= right and x[0] <= right:
            rightCount += 1
            rightIndex = i
            done = False
        i += 1
    # endTime = datetime.now()
    # endTime - startTime
    
    # startTime = datetime.now()
    # leftIndex = [i for i,x in enumerate(r) if x[1] > left and x[0] < left][0] # todo could probably speed this up
    # rightIndex = [i for i,x in enumerate(r) if x[1] > right and x[0] < right][0] # todo missing equals
    # endTime = datetime.now()
    # endTime - startTime
    
    rleft = r[leftIndex]
    rright = r[rightIndex]
    
    if leftIndex == rightIndex:
        bigR = [(right - left) * rleft[2]]
    else:
        intermediate = r[(leftIndex + 1):(rightIndex-1)]    
        bigR = [(y - x) * z for x,y,z in intermediate]
        bigR.append((rleft[1]-left) * rleft[2]) # todo need to deal with possibility left and right are inside the same window
        bigR.append((right-rright[0]) * rright[2])
        
    bigR = np.sum(bigR)
    return (1 - np.exp(- 2 * bigR)) / 2

def pointMassContribution(u, s, t, r):
    return - u / s * (s / (r + s) * (1 - math.exp(- r * t - s * t)))**2

def B(scaledu, ss, ps, t, recDist):
    return math.exp(sum([p * pointMassContribution(u, s, t, r) for u,r in zip(scaledu, recDist) for p,s in zip(ps, ss)]))

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot([B(scaledu, ss, ps, t, recDist) for t in range(int(censusSize))], "-", ms=8, lw=1, label="B(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend();

def rescaledPointMassContribution(scaledu, s, t, r, ancestralNe, ancTime):
    return - scaledu / s * (s / (r + s) * (1 - math.exp(- r * (ancTime - t * 2 * ancestralNe) - s * (ancTime - t * 2 * ancestralNe))))**2

def get_scaling_fun(positions, u, ss, ps, r, focalPos, censusSize, totalT, tol):
    scaledu = [z * (y - x) for x,y,z in u]
    recDist = [getRecDist(pos, r, focalPos) for pos in positions]
    
    diff = 100
    start = 1
    i = 0
    while abs(diff) > tol:
        end = B(scaledu, ss, ps, (i + 1) * censusSize / 10, recDist)
        diff = end - start
        start = end
        i += 1
        
    ancTime = censusSize / 10 * i
    ancTime = max(ancTime, totalT * 2 * censusSize)
    ancB = B(scaledu, ss, ps, ancTime, recDist)
    ancNe = ancB * censusSize
    return(lambda t: [math.exp(sum([p * rescaledPointMassContribution(u, s, t, r, ancNe, ancTime) for u,r in zip(scaledu, recDist) for p,s in zip(ps, ss)])) / ancB], ancTime, ancNe)

# human maps
os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/EquilibriumSims/humanMaps")
recName = "YRI_recombination_map_hg38_chr_22.bed"
mutName = "roulette_tbl_chr22.csv"
exonName = "exons_chr22.bed"
recMap = read_rec_map(recName)
mutMap = read_mut_rates(mutName)
exonMap = read_exon_map(exonName)

recMap = simplify_rate_map(recMap)
mutMap = simplify_rate_map(mutMap)

exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap) # TODO need to impliment splitting large regions

# equil
# hardcoding some parameters
u = 1e-8
r = 1e-8
L = 1e6
# regionSize = 1e4
# tol = 1e-3
focalPos = 5e5
sample_size = [40]
ss = [1e-2]

cs = lambda t: [1e3 + 2 * 1e3 * t]
totalT = 1

g = None
ps = None
targetSize = 1e4
tol = 1e-4

# positions of point masses.
# split 1Mb into 10kb regions, each point mass lies at the center
# of these regions 
# pointMassPosition = range(int(5e3),int(100.5e4),int(1e4))

#############################
######## integrated #########
#############################
def bgs_wrapper(u,
                r,
                focalPos,
                sample_size,
                ss,
                cs = None,
                g = None,
                totalT = None,
                L = None,
                ps = None,
                targetSize = 1e4,
                tol = 1e-4,
                ):
    if type(u) is not list and type(r) is not list and L is None:
            raise ValueError("If u and r are constant than the chrom. size, L, must be specified.")
    if type(u) is list: 
        if np.shape(u)[1] != 3:
            raise ValueError("If u is list it must be of the form [[start, stop, mu per bp]]")
        diffs = []
        i = 0
        while i < len(u) - 1:
            i += 1
            diffs.append(u[i][0] - u[i-1][1])
        if any([x < 0 for x in diffs]):
            raise ValueError("List u must be in increasing order by position and not overlap: u[i][0] >= u[i-1][1]")            
    if type(r) is list: 
        if np.shape(r)[1] != 3:
            raise ValueError("If r is list it must be of the form [[start, stop, r per bp]]")
        diffs = []
        i = 0
        while i < len(r) - 1:
            i += 1
            diffs.append(r[i][0] - r[i-1][1])
        if any([x <= 0 for x in diffs]): # TODO it may work with gaps
            raise ValueError("List u must be in increasing order by position with no gaps and no overlap: r[i][0] == r[i-1][1]")   
    if L is not None and type(r) is list:
        if r[-1][1] != L:
            raise ValueError("r must span the entire chromosome")
    if type(r) is list and type(u) is list:
        if r[-1][1] < u[1][1]:
            raise ValueError("r must span the entire chromosome")
    if type(focalPos) is not float and type(focalPos) is not int:
        raise ValueError("focalPos must be float or int")
    # TODO test samplesize, l, target size, and tol, ss, ps, cs, g, totalT
    
    if L is None:
        L = r[-1][1] 
        
    if type(u) is list:
        u = combine_and_split_regions(u)
    else:
        nregions = math.ceil(L/targetSize)
        regionSize = L / nregions
        u = [[regionSize * i, regionSize * (i+1), u] for i in range(nregions)]
        
    if ps is None:
        ps = [(i + 1) / len(ss) for i in range(len(ss))]
        
    positions = [(x + y)/2 for x,y,z in u]

    if type(r) is not list:
        r = [[0,L,r]]
    
    censusSize = cs(0)[0]
    # oldestEpoch, censusSize = getOldestEpoch(graph)

        
    fs = moments.Demographics1D.snm(sample_size)
    f, ancTime, ancNe = getSizeFun(positions, u, ss, r, focalPos, censusSize, tol)
    
            

        