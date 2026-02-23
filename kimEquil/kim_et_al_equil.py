import numpy as np
import fwdpy11
# import time
import moments
from dataclasses import dataclass
import sys
# from typing import List
# from collections import defaultdict
# import pickle
from datetime import datetime
import argparse
import msprime


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
        "--exon_map",
        "-eMap",
        type=str,
        default="exons_chr22.bed",
        help="name of bed file for exon positions",
    )
    optional.add_argument(
        "--mutation_map",
        "-mutMap",
        type=str,
        default="roulette_tbl_chr22.csv",
        help="name of csv describing mutation rate across the chromosome",
    )
    optional.add_argument(
        "--recombination_map",
        "-recMap",
        type=str,
        default="YRI_recombination_map_hg38_chr_22.bed",
        help="name of bed file describing recombination rate along the chromosome",
    )
    optional.add_argument(
        "--sample_size",
        "-n_sample",
        type=int,
        default="200",
        help="number of samples to project SFS down to",
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

def simplify_rate_map(data):
    merged = []

    prev_start = None
    prev_end = None
    prev_rate = None

    for x in data:
        start, end, rate = x

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

def read_rec_map(bed_path):
    intervals = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split()
            start = int(fields[1])
            end = int(fields[2])
            rate = float(fields[3])
            intervals.append([start, end, rate])
    return intervals

def intersect(a_start, a_end, b_start, b_end):
    beg = max(a_start, b_start)
    end = min(a_end, b_end)
    if beg < end:
        return beg, end
    return None

def get_max_positions(recMap, mutMap):
    return max(recMap[-1][1], mutMap[-1][1])

def get_total_rate(xMap):
    tmp = [rate * (end - start) for start, end, rate in xMap]
    return np.sum(tmp)

def make_sel_regions(mean_mu, exonMap, mean, scaling, shape = None):
    exonMutMap = [[x,y,mean_mu] for x,y in exonMap]
    U = get_total_rate(exonMutMap)
    sel_regions = []
    if shape is None:
        for start, end, mu in exonMutMap:
            sel_regions.append(
                fwdpy11.ConstantS(
                    beg=start, end=end, weight=mu, s=mean * scaling, h=1
                )
            )
    else:
        for start, end, mu in exonMutMap:
            sel_regions.append(
                fwdpy11.GammaS(
                    beg = start, end = end, weight = mu, mean = mean, shape_parameter = shape
               )
            )

    return sel_regions, U

def get_exon_windows(args):
    recName = args.recombination_map
    mutName = args.mutation_map
    exonName = args.exon_map

    # test
    # os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop/weakSelection")
    # recName = "YRI_recombination_map_hg38_chr_22.bed"
    # mutName = "roulette_tbl_chr22.csv"
    # exonName = "exons_chr22.bed"

    recMap = read_rec_map(recName)
    recMap = simplify_rate_map(recMap)

    mutMap = read_mut_rates(mutName)
    mutMap = simplify_rate_map(mutMap)
    
    exonMap = read_exon_map(exonName)
    
    chrom_end = get_max_positions(recMap, mutMap)
    
    windows = [x for w in exonMap for x in w]
    windows.insert(0, 0)
    windows.append(chrom_end)

    return windows 

def get_muts_in_windows(pos, w):
    idx = [i for i,x in enumerate(pos) if x >= w[0] and x <= w[1]]
    return idx

def runsim(args):
    """
    args: The parsed arguments
    Ne: The effective population size

    Returns a tree sequence
    """
    # Set the rng with the given seed
    rng = fwdpy11.GSLrng(args.seed)
    population_size = args.population_size
    recName = args.recombination_map
    mutName = args.mutation_map
    exonName = args.exon_map
    scaling = 1

    # test params
    # rng = fwdpy11.GSLrng(1)
    # population_size = 1e3    
    # recName = "YRI_recombination_map_hg38_chr_22.bed"
    # mutName = "roulette_tbl_chr22.csv"
    # exonName = "exons_chr22.bed"
    # scaling = 1

    recMap = read_rec_map(recName)
    mutMap = read_mut_rates(mutName)
    exonMap = read_exon_map(exonName)

    recMap = simplify_rate_map(recMap)
    
    mean_mu = np.mean([z for x,y,z in mutMap]) * 2.31 / (1+2.31) 
        
    rec_regions = make_rec_regions(recMap)

    means = -0.00657416 / 2
    shape = 0.186
    sel_regions,U = make_sel_regions(mean_mu, exonMap, means, scaling, shape)
    
    L = get_max_positions(recMap, mutMap)

    # Initialize the population
    Ne = int(population_size / scaling)
    # Ne = int(5e2)
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
    return ts

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    seed = args.seed
    population_size = args.population_size
    sample_size = args.sample_size
    exonMap = read_exon_map(args.exon_map)

    
    # test params
    # seed = 1
    
    eprint(
        current_time(),
        f"starting simulation for seed {args.seed}",
    )

    ts = runsim(args)
    
    times = ts.nodes_time
    sampleTimes = times[np.logical_and(times % population_size == 0, times <= 10 * population_size)]
    times = np.unique(sampleTimes)
    
    windows = get_exon_windows(args)
    
    # for branch lengths SFS
    
    combined_data = np.zeros((len(exonMap), 1 + int(sample_size)))

    for j,curTime in enumerate(times):    
        sampleIndex = [i for i, x in enumerate(sampleTimes == curTime) if x]
        
        afs = ts.allele_frequency_spectrum(sample_sets=[sampleIndex],
                                           windows=windows,
                                           mode="branch", 
                                           polarised=True, 
                                           span_normalise=False)
        
        midAfs = afs[1::2]
        
        data = np.empty((len(exonMap), 1 + int(sample_size)))
        for k,x in enumerate(midAfs):
            y = moments.Spectrum(x, data_folded = False)
            y = y.project([sample_size])
            data[k,:] = y
         # TODO is this right?
        combined_data = np.add(data, combined_data)
    
    np.savez_compressed(str(seed)+'_branch', exons = combined_data)

    # for selected SFS
    
    combined_data = np.zeros((len(exonMap), 1 + int(sample_size)))
    for j, curTime in enumerate(times):
        sampleIndex = [i for i,x in enumerate(sampleTimes == curTime) if x]
        frq = allele_frequencies(ts, sample_sets=[sampleIndex])[:,0]
        data = np.empty((len(exonMap), 1 + int(sample_size)))
        for k,w in enumerate(exonMap):
            wIndex = get_muts_in_windows(ts.sites_position, w)
            f = frq[wIndex]
            x = np.zeros(1 + 2 * int(population_size))
            for ff in f:
                # if ff != 0:
                #     raise ValueError("woah there, partner!")
                if ff * 2 * population_size % 1 != 0:
                    raise ValueError('not a count')
                x[int(ff * 2 * population_size)] += 1
            y = moments.Spectrum(x, data_folded = False)
            y = y.project([sample_size])
            data[k,:] = y
        combined_data = np.add(data, combined_data)
    
    np.savez_compressed(str(seed) + '_selSites', exons = combined_data)
    
    # for neutral sites SFS
    
    mutMap = read_mut_rates(args.mutation_map)
    neu_mu = np.mean([z for x,y,z in mutMap]) / (1+2.31) 
    
    neutMuts = msprime.sim_mutations(ts, 
                                     rate=neu_mu,
                                     random_seed = seed,
                                     model="binary",
                                     discrete_genome=False,
                                     keep=False)    
    # neutMuts.num_mutations
    
    combined_data = np.zeros((len(exonMap), 1 + int(sample_size)))
    for j, curTime in enumerate(times):
        sampleIndex = [i for i,x in enumerate(sampleTimes == curTime) if x]
        frq = allele_frequencies(neutMuts, sample_sets=[sampleIndex])[:,0]
        data = np.empty((len(exonMap), 1 + int(sample_size)))
        for k,w in enumerate(exonMap):
            wIndex = get_muts_in_windows(neutMuts.sites_position, w)
            f = frq[wIndex]
            x = np.zeros(1 + 2 * int(population_size))
            for ff in f:
                # if ff != 0:
                #     raise ValueError("woah there, partner!")
                if ff * 2 * population_size % 1 != 0:
                    raise ValueError('not a count')
                x[int(ff * 2 * population_size)] += 1
            y = moments.Spectrum(x, data_folded = False)
            y = y.project([sample_size])
            data[k,:] = y
        combined_data = np.add(data, combined_data)
    
    np.savez_compressed(str(seed) + '_neutSites', exons = combined_data)
        

