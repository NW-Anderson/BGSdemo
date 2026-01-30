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
# import msprime
import os
import demes
# test
# import demesdraw

# test
# os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa")


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


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.flush()


def current_time():
    return " [" + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S") + "]"


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(
        "flank_simulation.py", formatter_class=ADHF)
    parser.add_argument("--seed", required=True, type=int)
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "--demes_graph",
        "-yaml",
        type=str,
        default="demo.yaml",
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

    def __call__(self, pop, sampler):
        if pop.generation % 1000 == 0:
            eprint(current_time(), f"at generation {pop.generation}")


def getRecDist(pos, r, focalPos):
    left = min(pos, focalPos)
    right = max(pos, focalPos)

    bigR = r(right) - r(left)
    return (1 - np.exp(- 2 * bigR)) / 2


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
        i -= 1  # want to repeat same window of mutmap for next exon

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


def make_sel_regions(mutMap, exonMap, mean, scaling, shape = None):
    exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap)
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


def get_max_positions(recMap, mutMap):
    return max(recMap[-1][1], mutMap[-1][1])


def get_total_rate(xMap):
    tmp = [rate * (end - start) for start, end, rate in xMap]
    return np.sum(tmp)


def get_windows(args):
    recName = args.recombination_map
    mutName = args.mutation_map

    # test
    # os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop/weakSelection")
    # recName = "YRI_recombination_map_hg38_chr_22.bed"
    # mutName = "roulette_tbl_chr22.csv"

    recMap = read_rec_map(recName)
    recMap = simplify_rate_map(recMap)

    mutMap = read_mut_rates(mutName)
    mutMap = simplify_rate_map(mutMap)

    chrom_start = recMap[0][0]
    chrom_end = get_max_positions(recMap, mutMap)

    pos = [chrom_start + 5e5]
    done = False
    while not done:
        new_pos = pos[-1] + 1e6
        if new_pos <= chrom_end:
            pos.append(new_pos)
        else:
            done = True

    windows = [[p - 1, p + 1] for p in pos]
    windows = [x for w in windows for x in w]
    windows.insert(0, 0)
    windows.append(chrom_end)

    return windows


def runsim(args):
    """
    args: The parsed arguments
    Ne: The effective population size

    Returns a tree sequence
    """

    # Set the rng with the given seed
    rng = fwdpy11.GSLrng(args.seed)
    yaml = args.demes_graph
    recName = args.recombination_map
    mutName = args.mutation_map
    exonName = args.exon_map
    scaling = 1

    # test params
    # os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop/weakSelection")
    # rng = fwdpy11.GSLrng(1)
    # yaml = "ooa.yaml"
    # recName = "YRI_recombination_map_hg38_chr_22.bed"
    # mutName = "roulette_tbl_chr22.csv"
    # exonName = "exons_chr22.bed"
    # scaling = 1

    recMap = read_rec_map(recName)
    mutMap = read_mut_rates(mutName)
    exonMap = read_exon_map(exonName)

    recMap = simplify_rate_map(recMap)
    mutMap = simplify_rate_map(mutMap)

    graph = demes.load(yaml)

    # test
    # demesdraw.tubes(graph);

    rec_regions = make_rec_regions(recMap)

    means = -0.00657416 / 2
    shape = 0.186
    sel_regions, U = make_sel_regions(mutMap, exonMap, means, scaling, shape)

    # // m1 mutation type: gamma
    # // note: some rescaling is done since we use 1,1+2sh,1+s and SLiM uses 1,1+hs,1+s
    # initializeMutationType("m1", 0.5, "g", -0.00657402090856, 0.186); 1,1+hs,1+s

    demography = fwdpy11.ForwardDemesGraph.from_demes(
        graph, 20, burnin_is_exact=False)
    L = get_max_positions(recMap, mutMap)
    pop = fwdpy11.DiploidPopulation(demography.initial_sizes, L)

    # burnin = 20 * Ne
    # sampling = 10 * Ne  # number of sampling generations
    # sampling = 100
    # simlen = burnin + sampling
    # eprint(current_time(), "total simulation length:", simlen)
    eprint(current_time(), "total simulation length:",
           demography.final_generation)

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
        "demography": demography,
        "simlen": demography.final_generation,
    }
    params = fwdpy11.ModelParams(**pdict)

    recorder = Recorder()

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

    return ts, demography


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    seed = args.seed

    # test params
    # seed = 1

    eprint(
        current_time(),
        f"starting simulation for seed {args.seed}",
    )

    ts, demography = runsim(args)

    ceuIndex = None
    for deme in demography.demes_at_final_generation:
        # print(demography.deme_labels[deme])
        if demography.deme_labels[deme] == "CEU":
            ceuIndex = deme

    windows = get_windows(args)

    afs = ts.allele_frequency_spectrum(sample_sets=[ts.samples(population_id=ceuIndex)],
                                       windows=windows,
                                       mode="branch",
                                       polarised=True,
                                       span_normalise=False)

    midAfs = afs[1::2]

    # np.save("ceuData_" + str(seed) + ".npy", midAfs)
    # np.savetxt(str(seed) + "_ceu.csv", midAfs, delimiter=",")

    import random
    ss = [random.sample(sorted(ts.samples(population_id=d)), int(40))
          for d in demography.demes_at_final_generation]

    afs = ts.allele_frequency_spectrum(sample_sets=ss,
                                       windows=windows,
                                       mode="branch",
                                       polarised=True,
                                       span_normalise=False)

    np.savez_compressed(str(seed), ceu=midAfs, joint=afs[1::2])
