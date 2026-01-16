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

exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap) # todo need to impliment splitting large regions

# equil
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

#############################
######## integrated #########
#############################
def bgs_wrapper(u,
                r,
                focalPos,
                sample_size,
                L = None,
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
            
    if type(r) is list: 
        if np.shape(r)[1] != 3:
            raise ValueError("If r is list it must be of the form [[start, stop, r per bp]]")

    
            

        