import numpy as np
import fwdpy11
# import time
from dataclasses import dataclass
import sys
# from typing import List
# from collections import defaultdict
# import pickle
from datetime import datetime
import argparse
# import pandas as pd
#test
# import os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.flush()


def current_time():
    return " [" + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S") + "]"


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("flank_simulation.py", formatter_class=ADHF)
    parser.add_argument("--seed", required=True, type=int)
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "--population_size",
        "-N",
        type=int,
        default=10000,
        help="Diploid population size, defaults to 10,000.",
    )
    optional.add_argument(
        "--mean_sel_coef",
        "-means",
        type=float,
        default = -0.002,
        help="Mean of gamma dfe",
    )
    optional.add_argument(
        "--exon_map",
        "-eMap",
        type=str,
        default = "eMap.bed",
        help="name of bed file for exon positions",
    )
    optional.add_argument(
        "--mutation_map",
        "-mutMap",
        type=str,
        default = "mutmap.csv",
        help="name of csv describing mutation rate across the chromosome",
    )
    optional.add_argument(
        "--recombination_map",
        "-recMap",
        type=str,
        default = "recmap.csv",
        help="name of bed file describing recombination rate along the chromosome",
    )
    return parser

def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()] 
    n = np.array([len(x) for x in sample_sets])
    def f(x):
       return x / n
    return ts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode='site', strict=False, span_normalise=False)

@dataclass
class Recorder:
    burnin: int
    population_size: int
    
    
    def __call__(self, pop, sampler):
        if pop.generation % 1000 == 0:
            eprint(current_time(), f"at generation {pop.generation}")
        if np.logical_and(pop.generation >= self.burnin, pop.generation % population_size == 0):
            sampler.assign(range(pop.N))
            
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

def make_rec_regions(recMap):
    intervals = []
    for x in recMap:
        start, end, rate = x
        intervals.append(
            fwdpy11.PoissonInterval(start, end, rate * (end - start))
        )
    return intervals

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

def make_sel_regions(mutMap, exonMap, mean, scaling):
    exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap)
    sel_regions = []
    for start, end, mu in exonMutMap:
        sel_regions.append(
            fwdpy11.ConstantS(
                beg=start, end=end, weight=mu, s=mean * scaling, h=1
            )
        )
        
    return sel_regions,U

def get_max_positions(recMap, mutMap):
    return max(recMap[-1][1], mutMap[-1][1])

def get_total_rate(xMap):
    tmp = [rate * (end - start) for start, end, rate in xMap]
    return np.sum(tmp)
    
def runsim(args):
    """
    args: The parsed arguments
    Ne: The effective population size

    Returns a tree sequence
    """
    # Set the rng with the given seed
    rng = fwdpy11.GSLrng(args.seed)
    mean = - args.mean_sel_coef
    population_size = args.population_size
    recName = args.recombination_map
    mutName = args.mutation_map
    exonName = args.exon_map
    scaling = 1
    # L = 1e6
    # r = 1e-8
    # u = 1e-8
    
    
    # test params
    # os.chdir("/media/nathan/T7/BGSdemo/humanData")
    # rng = fwdpy11.GSLrng(1)
    # mean = - 0.01
    # population_size = 1e3
    # recName = "YRI_recombination_map_hg38_chr_22.bed"
    # mutName = "roulette_tbl_chr22.csv"
    # exonName = "exons_chr22.bed"
    
    recMap = read_rec_map(recName)
    mutMap = read_mut_rates(mutName)
    exonMap = read_exon_map(exonName)
    
    recMap = simplify_rate_map(recMap)
    mutMap = simplify_rate_map(mutMap)
    
    rec_regions = make_rec_regions(recMap)
    sel_regions,U = make_sel_regions(mutMap, exonMap, mean, scaling)

    # Initialize the population
    Ne = int(population_size / scaling)
    
    L = get_max_positions(recMap, mutMap)
    pop = fwdpy11.DiploidPopulation(Ne, L)

    burnin = 20 * Ne
    sampling = 10 * Ne  # number of sampling generations
    # sampling = 100
    simlen = burnin + sampling
    eprint(current_time(), "total simulation length:", simlen)

    pdict = {
        # Multiplicative selection model
        "gvalue": fwdpy11.Multiplicative(2.0),
        # Rates: (neutral, selection, recombination)
        "rates": (0.0, U, None),
        "nregions": [],
        # Selection within the non-recombining locus
        "sregions": sel_regions,
        # Recombination to the right of this locus
        "recregions": rec_regions,
        # Evolve a single deme of size N for 20*N generations
        "demography": fwdpy11.ForwardDemesGraph.tubes([Ne], simlen),
        "simlen": simlen,
    }
    params = fwdpy11.ModelParams(**pdict)

    recorder = Recorder(burnin=burnin, population_size=int(population_size))

    eprint(current_time(), "starting simulation")
    fwdpy11.evolvets(
        rng,
        pop,
        params,
        100,
        recorder=recorder,
        suppress_table_indexing=True,
        preserve_first_generation=False,
    )
    eprint(current_time(), "finished simulation")

    ts = pop.dump_tables_to_tskit()
    
    # todo add pseudo replicates
    return ts, L

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    seed = args.seed
    population_size = args.population_size
    
    # test params
    # seed = 1
    
    eprint(
        current_time(),
        f"starting simulation for seed {args.seed}",
    )

    ts,L = runsim(args)
    
    times = ts.nodes_time
    sampleTimes = times[np.logical_and(times % population_size == 0, times <= 10 * population_size)]
    times = np.unique(sampleTimes)
    
    # for branch lengths SFS
    
    data = np.empty((len(times), 1 + 2 * int(population_size)))
    
    for j,curTime in enumerate(times):    
        sampleIndex = [i for i, x in enumerate(sampleTimes == curTime) if x]
        
        afs = ts.allele_frequency_spectrum(sample_sets=[sampleIndex],
                                           windows=[0,L/2-1,L/2+1,L],
                                           mode="branch", 
                                           polarised=True, 
                                           span_normalise=False)
        
        midAfs = afs[1]
            
        data[j,:] = midAfs
    
    np.savetxt(str(seed) + ".csv", data, delimiter=",")
    
    # for selected SFS
    
    # closeIndex = [i for i,x in enumerate(ts.sites_position) if abs(x - 5e5)<5e4]

    # data = np.empty((len(times), len(closeIndex)))
    
    # for j,curTime in enumerate(times):
    #     sampleIndex = [i for i,x in enumerate(sampleTimes == curTime) if x]
    #     frq = allele_frequencies(ts, sample_sets=[sampleIndex])[:,0][closeIndex]
        
    #     data[j,:] = frq
        
    # np.savetxt(str(seed) + "_wi5e4.csv", data, delimiter = ",")
    
    # closeIndex = [i for i,x in enumerate(ts.sites_position) if abs(x - 5e5)<1e5]

    # data = np.empty((len(times), len(closeIndex)))
    
    # for j,curTime in enumerate(times):
    #     sampleIndex = [i for i,x in enumerate(sampleTimes == curTime) if x]
    #     frq = allele_frequencies(ts, sample_sets=[sampleIndex])[:,0][closeIndex]
        
    #     data[j,:] = frq
        
    # np.savetxt(str(seed) + "_wi1e5.csv", data, delimiter = ",")
    
    # for neutral sites SFS
    
    # neutMuts = msprime.sim_mutations(ts, 
    #                                  rate=1e-6,
    #                                  random_seed = seed,
    #                                  model="binary",
    #                                  discrete_genome=False,
    #                                  keep=False)    
    # # neutMuts.num_mutations
    # closeIndex = [i for i,x in enumerate(neutMuts.sites_position) if abs(x - 5e5)<5e3]
    # frq = allele_frequencies(neutMuts)[closeIndex]
    # # min(abs(neutMuts.sites_position - 5e5))
    
    # import csv
    
    # with open(str(seed), 'w') as f:
     
    #     # using csv.writer method from CSV package
    #     write = csv.writer(f)
    #     write.writerows(map(lambda x: x, frq))
        

