import numpy as np
import fwdpy11
import time
from dataclasses import dataclass
import sys
from typing import List
from collections import defaultdict
import pickle
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
        "--mean_sel_coef",
        "-means",
        type=float,
        default = -0.002,
        help="Mean of gamma dfe",
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
    scaling = 1
    L = 1e6
    r = 1e-8
    u = 1e-8
    
    
    # test params
    rng = fwdpy11.GSLrng(1)
    mean = - 0.01
    population_size = 1e3

    U = u * L
    R = r * L

    rec_regions = [fwdpy11.PoissonInterval(0, L, L * r * scaling)]
    sel_regions = [
        fwdpy11.ConstantS(
            beg=0, end=L, weight=1, s=mean * scaling, h=1
        )
    ]

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
    
    # test params
    seed = 1
    
    eprint(
        current_time(),
        f"starting simulation for seed {args.seed}",
    )

    ts = runsim(args)
    
    times = ts.nodes_time
    sampleTimes = times[np.logical_and(times % population_size == 0, times <= 10 * population_size)]
    times = np.unique(sampleTimes)
    
    data = np.empty((len(times), 1 + 2 * int(population_size)))
    
    for j,curTime in enumerate(times):    
        sampleIndex = [i for i, x in enumerate(sampleTimes == curTime) if x]
        
        afs = ts.allele_frequency_spectrum(sample_sets=[sampleIndex],
                                           windows=[0,5e5-1,5e5+1,1e6],
                                           mode="branch", 
                                           polarised=True, 
                                           span_normalise=False)
        
        midAfs = afs[1]
            
        data[j,:] = midAfs
    
    np.savetxt(str(seed) + ".csv", data, delimiter=",")
        
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
        

