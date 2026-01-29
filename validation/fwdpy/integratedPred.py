import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd
import types
import demes
import demesdraw
import copy
from scipy import interpolate
import warnings

os.chdir('/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy')
from deme_funs import _get_demographic_events, _get_deme_sample_sizes, _get_root_Ne, _sizes_at_time, _migration_rate_in_interval, _compute_sfs, _reorder_fs

def make_cum_map(r, isCum = False):
    if isCum:
        pos = [x for x,y in r]
        cum = [y for x,y in r]
    else:
        cum = np.cumsum([(y-x)*z for x,y,z in r])
        pos = [y for x,y,z in r]
        
        cum = np.insert(cum,0,0)
        pos = np.insert(pos, 0, r[0][0])
    
    return interpolate.interp1d(pos, cum)

def rate_diff(thing):
    rates = [x[-1] for x in thing]
    return np.diff(rates)
    
def unique(thing):
    tmp = []
    for x in thing:
        if x not in tmp:
            tmp.append(x)
    return tmp

def parseJointData():
    os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/twopop")
    
    file_path = "ooaTwoPop.txt"
    data_list = []
    
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line by whitespace and append the list of elements
            # You may want to convert elements to their correct types (e.g., int, float)
            elements = line.strip().split() 
            if elements: # Avoid processing empty lines
                data_list.append(elements)
                
    os.chdir("/media/nathan/T7/BGSdemo/ooaTwoPopData/joint")
    
    pdata = np.zeros((3,101,101))
    numSum = np.zeros(3)
    for s, d, seed in data_list:
        simData = np.load(seed +'.npy')
        
        if s == '0.001':
            i = 0
        if s == '0.005':
            i = 1
        if s == '0.01':
            i = 2
            
        pdata[i] = np.add(pdata[i], simData)
        numSum[i] += 1
    
    parsedData = [p / n for p,n in zip(pdata, numSum)]
    return parsedData

def censusFun(demo):
    def N(t):
        sizes = []
        for deme in demo.demes:
            size_t = None
            for epoch in deme.epochs:
                if epoch.end_time <= t < epoch.start_time:
                    if epoch.size_function == "constant":
                        size_t = epoch.start_size
                    else:
                        dt = (epoch.start_time - t) / epoch.time_span
                        r = math.log(epoch.end_size / epoch.start_size)
                        size_t = epoch.start_size * math.exp(r * dt)
                    break
            sizes.append(size_t)
        return sizes
    return lambda t: N(t)

def reversedCensusFun(demo, ancTime, ancNe):
    tmp_fun = censusFun(demo)
    denom = tmp_fun(ancTime)[0]

    def ds(t):
        # map coalescent time t (in units of ancNe) back to demes time
        tt = ancTime - t * 2 * ancNe
        # clamp into [0, ancTime]
        if tt < 0:
            tt = 0.0
            # print("overshoot:", t, "tt:", tt)
        elif tt > ancTime:
            tt = ancTime
        return [x / denom for x in tmp_fun(tt) if x is not None]
    return ds

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

def combine_and_split_regions(exonMutMap, targetSize = 1e4):
    # first split
    split = []
    for start, stop, mu in exonMutMap:
        if stop - start > targetSize:
            i = 2
            done = True
            while done:
                new_stop = start + (stop - start) / i
                if new_stop - start < targetSize:
                    new_size = new_stop - start
                    done = False
                else:
                    i += 1
            tmp = [[start + new_size * j, start + new_size * (j + 1), mu] for j in range(i)]
            split.extend(tmp)
        else:
            split.append([start, stop, mu])
            
    combined = []
    tmp = []
    tmp_start = split[0][0]
    for start, stop, mu in split:        
        if (abs(stop - tmp_start) < targetSize):
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

def getRecDist(pos,r,focalPos):
    left = min(pos, focalPos)
    right = max(pos, focalPos)
    
    bigR = r(right) - r(left)
    return (1 - np.exp(- 2 * bigR)) / 2    
    
# getRecDist(pos, recMap, focalPos)
# def getRecDist(pos, r, focalPos):
#     left = min(pos, focalPos)
#     right = max(pos, focalPos)
    
#     # startTime = datetime.now()
#     i = 0
#     done = True
#     leftCount = 0
#     rightCount = 0
#     while done:
#         x = r[i]
#         if x[1] >= left and x[0] <= left:
#             leftCount += 1
#             leftIndex = i
#         if x[1] >= right and x[0] <= right:
#             rightCount += 1 
#             rightIndex = i
#             done = False
#         i += 1
#     # endTime = datetime.now()
#     # endTime - startTime
    
#     # startTime = datetime.now()
#     # leftIndex = [i for i,x in enumerate(r) if x[1] > left and x[0] < left][0] # todo could probably speed this up
#     # rightIndex = [i for i,x in enumerate(r) if x[1] > right and x[0] < right][0] # todo missing equals
#     # endTime = datetime.now()
#     # endTime - startTime
    
#     rleft = r[leftIndex]
#     rright = r[rightIndex]
    
#     if leftIndex == rightIndex:
#         bigR = [(right - left) * rleft[2]]
#     else:
#         intermediate = r[(leftIndex + 1):(rightIndex-1)]    
#         bigR = [(y - x) * z for x,y,z in intermediate]
#         bigR.append((rleft[1]-left) * rleft[2]) 
#         bigR.append((right-rright[0]) * rright[2])
        
#     bigR = np.sum(bigR)
#     return (1 - np.exp(- 2 * bigR)) / 2

def pointMassContribution(u, s, t, r):
    return - u / s * (s / (r + s) * (1 - math.exp(- r * t - s * t)))**2

def B(scaledu, ss, ps, t, recDist):
    return math.exp(sum([p * pointMassContribution(u, s, t, r) for u,r in zip(scaledu, recDist) for p,s in zip(ps, ss)]))

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
  
def getOldestEpoch(graph):
    tme = 0
    for deme in graph.demes:
        for epoch in deme.epochs:
            if epoch.end_time > tme:
                tme = epoch.end_time
                size = epoch.start_size
    return tme, size
     
def add_new_ancestral_root(graph,
                          split_time,
                          new_name='buf',
                          size=None) -> demes.Graph:
    """
    Add a new ancestral deme above the current root(s).

    Parameters
    ----------
    graph : demes.Graph
        Existing demes graph.
    new_name : str
        Name for the new ancestral population.
    split_time : float
        Time (in generations before present) when the new ancestral ends and the
        old root(s) begin. Must be > 0 and finite.
    size : float | None
        Constant size for the new ancestral. If None, use the first root's
        initial size (its size at its earliest time) as a default.

    Returns
    -------
    demes.Graph
        New graph with the added ancestral root.
    """
    if not (math.isfinite(split_time) and split_time > 0):
        raise ValueError("split_time must be finite and > 0.")
    if any(d.name == new_name for d in graph.demes):
        raise ValueError(f"Deme '{new_name}' already exists.")

    # Convert to a mutable representation
    d = copy.deepcopy(graph.asdict())

    # Find root demes (no ancestors)
    root_demes = [deme for deme in d["demes"] if len(deme.get("ancestors", [])) == 0]
    if len(root_demes) == 0:
        raise ValueError("No root demes found (unexpected).")

    # Choose a default size if not provided
    if size is None:
        # Use the first root's earliest (ancestral) epoch start_size
        # (In demes dict form, epochs[0]["start_size"] is typically the size at deme start_time.)
        size = root_demes[0]["epochs"][0]["start_size"]

    # Add the new ancestral deme
    new_deme = {
        "name": new_name,
        "start_time": math.inf,
        "epochs": [
            {"end_time": float(split_time), "start_size": float(size)}
        ],
    }
    # d["demes"].append(new_deme)
    d["demes"].insert(0, new_deme)

    # Re-root each old root under the new ancestral
    for deme in root_demes:
        # Ensure the old root starts at split_time (it cannot still start at inf if it has an ancestor).
        deme["start_time"] = float(split_time)
        deme["ancestors"] = [new_name]
        deme["proportions"] = [1.0]

    # Rebuild graph
    return demes.Graph.fromdict(d)

def SFS_bgs(
    g,
    scaling_fun,
    bgs_Ne,
    anc_gen,
    sampled_demes=None,
    sample_sizes=None,
    sample_times=None,
    samples=None,
    unsampled_n=4,
    gamma=None,
    s=None,
    h=None,
    theta=1,
    u=None,
    reversible=False,
    L=1,
):
    """
    Compute the SFS from a ``demes``-specified demographic model.
    ``demes`` is a package for specifying demographic models in a
    user-friendly, human-readable YAML format. This function
    automatically parses the demographic model and returns the SFS
    for the specified populations, sample sizes, and (optionally)
    sampling times.

    Selection and dominance can be specified as a single value for
    all populations, or on a per-deme basis using a dictionary
    mapping deme name to the coefficient (defaults can also be set
    if multiple demes have the same selection or dominance
    coefficient). The mutation rate can be given as either a
    per-base rate (possibly multiplied by the sequence length),
    or as a population size-scaled rate. If mutation rates are not
    given, the SFS is scaled by ``4*N_e``, so that multiplying the
    output SFS by ``u`` results in a properly scaled SFS.

    :param g: A ``demes`` DemeGraph from which to compute the SFS.
    :type g: :class:`demes.DemeGraph`
    :param sampled_demes: A list of deme IDs to take samples from. We can repeat
        demes, as long as the sampling of repeated deme IDs occurs at distinct
        times.
    :type sampled_demes: list of strings
    :param sample_sizes: A list of the same length as ``sampled_demes``,
        giving the sample sizes for each sampled deme.
    :type sample_sizes: list of ints
    :param sample_times: If None, assumes all sampling occurs at the end of the
        existence of the sampled deme. If there are
        ancient samples, ``sample_times`` must be a list of same length as
        ``sampled_demes``, giving the sampling times for each sampled
        deme. Sampling times are given in time units of the original deme graph,
        so might not necessarily be generations (e.g. if ``g.time_units`` is years)
    :type sample_times: list of scalars, optional
    :param unsampled_n: The default sample size of unsampled demes, which must be
        greater than or equal to 4.
    :type unsampled_n: int, optional
    :param gamma: The scaled selection coefficient(s), ``2*Ne*s``. Defaults to None,
        which implies neutrality. Can be given as a scalar value, in which case
        all populations have the same selection coefficient. Alternatively, can
        be given as a dictionary, with keys given as population names in the
        input Demes model. Any population missing from this dictionary will be
        assigned a selection coefficient of zero. A non-zero default selection
        coefficient can be provided, using the key ``_default``. See the Demes
        exension documentation for more details and examples.
    :type gamma: scalar or dict
    :param h: The dominance coefficient(s). Defaults to additivity (or genic
        selection). Can be given as a scalar value, in which case all populations
        have the same dominance coefficient. Alternatively, can be given as a
        dictionary, with keys given as population names in the input Demes model.
        Any population missing from this dictionary will be assigned a dominance
        coefficient of ``1/2`` (additivity). A different default dominance
        coefficient can be provided, using the key ``_default``. See the Demes
        exension documentation for more details and examples.
    :type h: scalar or dict
    :param theta: The scaled mutation rate(s), 4*Ne*u. When simulating under the
        infinite sites model (the default mutation model), ``theta`` should be given
        as a scalar value greater than zero. If it is not provided, it is computed
        using the input value of ``u`` as ``4*Ne*u``. If ``u`` is not provided, then
        the SFS is scaled by ``4*Ne``, and the user can recover a properly scaled SFS
        by multiplying it by ``u`` or ``u*L``. When simulating under the reversible
        mutation model (with ``reversible=True``), ``theta`` may be a list of length
        2 and both the forward and backward scaled mutation rates must be less
        than 1.
    :type theta: scalar or list of length 2
    :param u: The per-base mutation rate. When simulating under the infinite sites
        model (the default mutation model), ``u`` should be a scalar. When simulating
        under the reversible mutation model (with ``reversible=True``), ``u`` may
        be a list of length 2, and mutation rate(s) must be small enough so that
        the product of ``4*Ne*u`` is less than 1.
    :type u: scalar or list of length 2
    :param L: The effective sequence length, which may be used along with ``u`` to
        set the total mutation rate. Defaults to 1, and it must be 1 when using
        the reversible mutation model.
    :type L: scalar
    :return: A ``moments`` site frequency spectrum, with dimension equal to the
        length of ``sampled_demes``, and shape equal to ``sample_sizes`` plus one
        in each dimension, indexing the allele frequency in each deme from 0
        to n[i], where i is the deme index.
    :rtype: :class:`moments.Spectrum`
    """
    
    # get reference Ne from demes model
    # cs_Ne = _get_root_Ne(g)
    oldest_epoch, cs_Ne  = getOldestEpoch(g)
    if anc_gen > oldest_epoch:
        g = add_new_ancestral_root(g, anc_gen)
    # demesdraw.tubes(g);
    # demesdraw.tubes(test);
    
    # could specify samples as a dict instead of sampled_demes and sample_sizes
    if samples is None:
        if sampled_demes is None or sample_sizes is None:
            raise ValueError(
                "must specify either samples (as a dict mapping demes to sample sizes,"
                " or specify both sampled_demes and sample_times"
            )
    else:
        if type(samples) is not dict:
            raise ValueError("samples must be a dict mapping demes to sample sizes")
        if sampled_demes is not None or sample_sizes is not None:
            raise ValueError(
                "if samples is given as dict, cannot "
                "specify sampled_demes or sample_sizes"
            )
        if sample_times is not None:
            raise ValueError("if samples is given as dict, cannot specify sample times")
        sampled_demes = list(samples.keys())
        sample_sizes = list(samples.values())

    if len(sampled_demes) != len(sample_sizes):
        raise ValueError("sampled_demes and sample_sizes must be same length")
    if sample_times is not None and len(sampled_demes) != len(sample_times):
        raise ValueError("sample_times must have same length as sampled_demes")
    for deme in sampled_demes:
        if deme not in g:
            raise ValueError(f"deme {deme} is not in demography")

    # we need to copy these to new variable names
    # so they don't get updated during optimization
    sampled_pops = copy.copy(sampled_demes)
    deme_sample_times = copy.copy(sample_times)

    if unsampled_n < 4:
        raise ValueError("unsampled_n must be greater than 3")

    sampled_deme_end_times = [g[d].end_time for d in sampled_pops]
    if deme_sample_times is None:
        deme_sample_times = sampled_deme_end_times

    for t, d in zip(deme_sample_times, sampled_pops):
        if t < g[d].end_time or t >= g[d].start_time:
            raise ValueError(f"sample time {t} is outside of deme {d}'s time span")

    # for any ancient samples, we need to add frozen branches
    # with this, all "sample times" are at time 0, and ancient sampled demes are frozen
    if np.any(np.array(deme_sample_times) != 0):
        g, sampled_pops, list_of_frozen_demes = _augment_with_ancient_samples(
            g, sampled_pops, deme_sample_times
        )
        deme_sample_times = [0 for _ in deme_sample_times]
    else:
        list_of_frozen_demes = []

    if g.time_units != "generations":
        g, deme_sample_times = _convert_to_generations(g, deme_sample_times)

    # if any sample sizes are less than unsample_n, we increase and project after
    sim_sample_sizes = []
    for d, n, t in zip(sampled_pops, sample_sizes, deme_sample_times):
        sim_sample_sizes.append(max(n, unsampled_n))
        if t < g[d].end_time or t >= g[d].start_time:
            raise ValueError("sample time for {deme} must be within its time span")

    # if (unscaled) s is provided, convert into (scaled) gamma selection coefficients
    if s is not None:
        if gamma is not None:
            raise ValueError("Cannot specify both gamma and s")
        if isinstance(s, (int, float)):
            gamma = 2 * Ne * s
        elif type(s) is dict:
            gamma = {}
            for k, v in s.items():
                gamma[k] = 2 * Ne * v
        else:
            raise ValueError("Selection coefficient must be a scalar value or dict")

    # check selection and dominance inputs
    if gamma is not None:
        if "_default" in g:
            raise ValueError(
                "Cannot use `_default` as a deme name when selection is specified"
            )
        if isinstance(gamma, (int, float)):
            if not np.isfinite(gamma):
                raise ValueError("Selection coefficient must be a finite number")
        elif type(gamma) is dict:
            for k in gamma.keys():
                if k != "_default" and k not in g:
                    raise ValueError(f"Deme {k} in gamma, but {k} not in input graph")
        else:
            raise ValueError("Selection coefficient must be a scalar value or dict")
    if h is not None:
        if type(h) is dict:
            for k in h.keys():
                if k != "_default" and k not in g:
                    raise ValueError(f"Deme {k} in h, but {k} not in input graph")

    # set up the mutation rates as needed
    if theta is None:
        if u is None:
            u = 1
        if isinstance(u, (int, float)):
            theta = 4 * Ne * u * L
        else:
            if np.ndim(u) != 1 or len(u) != 2: 
                raise ValueError(
                    "Mutation rates must be a list of length 2 when using "
                    "the reversible mutation model"
                )
            theta = [4 * Ne * u[0], 4 * Ne * u[1]]
    else:
        if u is not None:
            raise ValueError("Only one of u or theta may be specified")
        if isinstance(theta, (int, float)):
            theta *= L
        else:
            if np.ndim(theta) != 1 or len(theta) != 2:
                raise ValueError(
                    "Mutation rates must be a list of length 2 when using "
                    "the reversible mutation model"
                )
            theta[0] *= L
            theta[1] *= L

    # if a scalar, must be positive; if list-like, must be length 2 and both positive
    if not reversible:
        if not isinstance(theta, (int, float)):
            raise ValueError(
                "Mutation rate must be a scalar value for the default ISM model"
            )
        if theta <= 0:
            raise ValueError("Mutation rate must be positive")
    if reversible:
        if L != 1:
            raise ValueError(
                "Sequence length L must be 1 when using the reversible mutation model"
            )
        if isinstance(theta, (int, float)):
            theta = [theta, theta]
        if theta[0] <= 0 or theta[1] <= 0:
            raise ValueError("Mutation rates must be positive")
        if theta[0] >= 1 or theta[1] >= 1:
            raise ValueError("Mutation rates too large for reversible mutation model")

    # get the list of demographic events from demes, which is a dictionary with
    # lists of splits, admixtures, mergers, branches, and pulses
    demes_demo_events = g.discrete_demographic_events()

    # get the dict of events and event times that partition integration epochs, in
    # descending order. events include demographic events, such as splits and
    # mergers and admixtures, as well as changes in population sizes or migration
    # rates that require instantaneous changes in the size function or migration matrix.
    # get the list of demes present in each epoch, as a dictionary with non-overlapping
    # adjoint epoch time intervals
    demo_events, demes_present = _get_demographic_events(
        g, demes_demo_events, sampled_pops
    )

    for epoch, epoch_demes in demes_present.items():
        if len(epoch_demes) > 5:
            raise ValueError(
                f"Moments cannot integrate more than five demes at a time. "
                f"Epoch {epoch} has demes {epoch_demes}."
            )


    # get the list of size functions, migration matrices, and frozen attributes from
    # the deme graph and event times, matching the integration times
    # nu_funcs, mig_mats, Ts, frozen_pops = _get_integration_parameters(
    #     g, demes_present, list_of_frozen_demes, Ne=cs_Ne
    # )
    
    # nu_funcs_bgs, mig_mats_bgs, Ts_bgs, frozen_pops_bgs = _get_integration_parameters(
    #     g, demes_present, list_of_frozen_demes, Ne=ancNe
    # )
    
    # nu_funcs_test, mig_mats_test, Ts_test, frozen_pops_test 
    nu_funcs, mig_mats, Ts, frozen_pops = _get_integration_parameters_bgs(
        g, demes_present, list_of_frozen_demes, bgs_Ne, scaling_fun
    )
    
    # get the sample sizes within each deme, given sample sizes
    deme_sample_sizes = _get_deme_sample_sizes(
        g,
        demo_events,
        sampled_pops,
        sim_sample_sizes,
        demes_present,
        unsampled_n=unsampled_n,
    )

    # compute the SFS
    fs = _compute_sfs(
        demo_events,
        demes_present,
        deme_sample_sizes,
        nu_funcs,
        mig_mats,
        Ts,
        frozen_pops,
        theta=theta,
        gamma=gamma,
        h=h,
        reversible=reversible,
    )

    fs = _reorder_fs(fs, sampled_pops)

    # project down to desired sample sizes, if needed
    fs = fs.project(sample_sizes)
    # simplify pop id name if ancient sample at end time of that deme
    for ii, pid in enumerate(fs.pop_ids):
        if "_sampled_" in pid:
            p, t = pid.split("_sampled_")
            t = float(t.replace("_", "."))
            if t == sampled_deme_end_times[ii]:
                fs.pop_ids[ii] = p

    return fs

def _get_integration_parameters_bgs(g, demes_present, frozen_list, bgs_Ne, scaling_fun, cs_Ne=None):
    """
    Returns a list of size functions, migration matrices, integration times,
    and lists frozen demes.
    """
    nu_funcs = []
    integration_times = []
    migration_matrices = []
    frozen_demes = []

    if cs_Ne is None:
        cs_Ne = _get_root_Ne(g)
    else:
        if cs_Ne != _get_root_Ne(g):
            warnings.warn(
                "Input cs_Ne is different from root population initial size, "
                "subsequent population size scaling may be incorrect"
            )
    
    T_elapsed = 0
    for interval, live_demes in sorted(demes_present.items())[::-1]:
        # get intergration time for interval
        T = (interval[0] - interval[1]) / 2 / bgs_Ne
        if T == math.inf:
            T = 0
        integration_times.append(T)
        # get frozen attributes
        freeze = [d in frozen_list for d in live_demes]
        frozen_demes.append(freeze)
        # get nu_function or list of sizes (if all constant)
        sizes = []
        for d in live_demes:
            sizes.append(_sizes_at_time(g, d, interval))
        # nu_func = _make_nu_func(sizes, T, cs_Ne)
        nu_func = _make_nu_func_bgs(sizes, T, cs_Ne, T_elapsed, scaling_fun)
        T_elapsed += T

        nu_funcs.append(nu_func)
        # get migration matrix for interval
        mig_mat = np.zeros((len(live_demes), len(live_demes)))
        for ii, d_from in enumerate(live_demes):
            for jj, d_to in enumerate(live_demes):
                if d_from != d_to:
                    m = _migration_rate_in_interval(g, d_from, d_to, interval)
                    mig_mat[jj, ii] = 2 * bgs_Ne * m
        migration_matrices.append(mig_mat)

    return nu_funcs, migration_matrices, integration_times, frozen_demes

def _make_nu_func_bgs(sizes, T, Ne, T_elapsed, scaling_fun):
    """
    Given the sizes at start and end of time interval, and the size function for
    each deme, along with the integration time and reference Ne, return the
    size function that gets passed to the moments integration routines.
    """
    if np.all([s[-1] == "constant" for s in sizes]):
        # all constant
        if T == 0:
            nu_func = [s[0] / Ne for s in sizes]
        else:
            nu_funcs_separated = []
            for s in sizes:
                assert s[0] == s[1]
                nu_funcs_separated.append(lambda t, N0=s[0]: (N0 / Ne) * scaling_fun(T_elapsed + t)[0])
                
            def nu_func(t):
                return [nu(t) for nu in nu_funcs_separated]
    else:
        nu_funcs_separated = []
        for s in sizes:
            if s[-1] == "constant":
                assert s[0] == s[1]
                nu_funcs_separated.append(lambda t, N0=s[0]: (N0 / Ne) * scaling_fun(T_elapsed + t)[0])
            elif s[-1] == "linear":
                nu_funcs_separated.append(
                    lambda t, N0=s[0], NF=s[1]: (N0 / Ne + t / T * (NF - N0) / Ne) * scaling_fun(T_elapsed + t)[0] 
                )
            elif s[-1] == "exponential":
                nu_funcs_separated.append(
                    lambda t, N0=s[0], NF=s[1]: (N0
                    / Ne
                    * np.exp(np.log(NF / N0) * t / T)) * scaling_fun(T_elapsed + t)[0]
                )
            else:
                raise ValueError(f"{s[-1]} not a valid size function")

        def nu_func(t):
            return [nu(t) for nu in nu_funcs_separated]

        # check that this is correct, or if we have to "pin" parameters
    return nu_func

def rescale_cs(cs, totalT, ancTime, ancNe, censusSize): 
    if totalT * 2 * censusSize == ancTime:
        return lambda t: [cs(t / censusSize * ancNe)[0] / censusSize] 
    elif totalT * 2 * censusSize < ancTime:
        # time to add before demographic event
        buffer = ancTime - 2 * censusSize * totalT
        def ds(t):
            # convert time to generations
            gen = t * 2 * ancNe
            # if input time is before the start of demography
            if gen <= buffer:
                tt = 0
            else: 
                tt = (gen - buffer) / 2 / censusSize
            return [cs(tt)[0] / censusSize]
        return ds
    else: 
        raise ValueError("ancTime < totalT provided. This should never happen") 
    
def get_map_diffs(thing):
    diffs = []
    i = 0
    while i < len(thing) - 1:
        i += 1
        diffs.append(thing[i][0] - thing[i-1][1])
    return diffs

# human maps
# os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/EquilibriumSims/humanMaps")
# recName = "YRI_recombination_map_hg38_chr_22.bed"
# mutName = "roulette_tbl_chr22.csv"
# exonName = "exons_chr22.bed"
# recMap = read_rec_map(recName)
# mutMap = read_mut_rates(mutName)
# exonMap = read_exon_map(exonName)

# recMap = simplify_rate_map(recMap)
# mutMap = simplify_rate_map(mutMap)

# exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap) 

# # hardcoding some parameters
# # u = 1e-8
# u = exonMutMap
# # r = 1e-8
# r = recMap
# # cum = np.cumsum([(y-x)*z for x,y,z in r])
# # pos = [y for x,y,z in r]

# # cum = np.insert(cum,0,0)
# # pos = np.insert(pos, 0, r[0][0])

# # r = [[x,y] for x,y in zip(pos,cum)]

# # L = 1e6
# L = None
# # focalPos = 5e5
# focalPos = r[-1][1] / 2
# sample_size = [40]
# ss = [1e-2]
# # ss = [1e-2, 5e-3]

# # cs = lambda t: [1e3 + 2 * 1e3 * t]
# # cs = [1e3]
# # totalT = 0.1
# # totalT = 1
# # totalT = 0
# cs = None
# totalT = None

# os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop")

# # g = None
# # sampled_demes = None
# g = demes.load('ooa.yaml')
# sampled_demes = ["CEU"]

# ps = None
# targetSize = 1e4
# tol = 1e-4
# minPos = None
# # focal_s = 1e-3
# focal_s = None

# todo make theta != 1 support is that just setting theta != 1???

#############################
######## integrated #########
#############################

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[g(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="g(t)")
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="f(t)")
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[rescaledcs(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="rescaledcs(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("value")
# ax.legend(); 
     
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[rescaledcs(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="rescaledcs(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("rescaledcs(t)")
# ax.legend();     

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,totalT, 0.001),[cs(t) for t in np.arange(0,totalT, 0.001)], "-", ms=8, lw=1, label="cs(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("cs(t)")
# ax.legend();   

# [cs(0),rescaledcs(0)]

# [cs(totalT), rescaledcs(ancTime / 2/ ancNe)]

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="f(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("f(t)")
# ax.legend();     

def bgs_wrapper(u,
                r,
                focalPos,
                sample_size,
                ss,
                cs = None,
                g = None,
                sampled_demes = None,
                totalT = None,
                L = None,
                ps = None,
                targetSize = 1e4,
                tol = 1e-4,
                minPos = None,
                focal_s = None,
                ):
    # check u is correct shape
    if type(u) is list:
        if np.ndim(u) != 2 or np.shape(u)[1] != 3:
            raise ValueError("If u is list it must be of the form [[start, stop, mu per bp]]")
        else: 
            if any([x < 0 for x in get_map_diffs(u)]):
                raise ValueError('List u must be in increasing order by position and not overlap: u[i][0] >= u[i-1][1]')
    elif isinstance(u, (int, float)):
        if u < 0:
            raise ValueError('u must be greater than 0.')
    else:
        raise ValueError("u must be a float or integer or list of the form [[start, stop, mu per bp]]")
        
    # check r is correct shape
    if type(r) is list:
        if np.ndim(r) == 2:
            if np.shape(r)[1] == 2:
                isCum = True
                if any(z <= 0 for z in np.diff([y for x,y in r])) or any(z <= 0 for z in np.diff([x for x,y in r])):
                    raise ValueError('if r is a cumulative map [pos, r_cumulative], then both position and r_cumulative must be strictly increasing: r[i] > r[i-1] for r[][0] and r[][1].')
                if r[0][-1] != 0:
                    raise ValueError('r_cumulative must begin at 0')
            elif np.shape(r)[1] == 3:
                isCum = False
                if any([x != 0 for x in get_map_diffs(r)]):
                    raise ValueError("List r must be in increasing order by position with no gaps and no overlap: r[i][0] == r[i-1][1]") 
            else:
                raise ValueError('If r is list it must be of the form [[start, stop, r per bp]] or [[pos, r_cumulative]]')
        else: 
            raise ValueError('If r is list it must be of the form [[start, stop, r per bp]] or [[pos, r_cumulative]]')
    elif isinstance(r, (int, float)):
        if r < 0:
            raise ValueError('r must be greater than 0')
    else:
        raise ValueError("r must be a float or integer or list of the form [[start, stop, r per bp]] or [[pos, r_cumulative]]")
    
    # if maps provided make sure the mutation map is contained within rec map
    if type(r) is list and type(u) is list:
        if r[-1][-2] < u[-1][1]:
            raise ValueError("Final position of u is greater than r. r must span the entire chromosome")
        if r[0][0] > u[0][0]:
            raise ValueError('Start position of u is less than the start of r. r must span the entire chromosome')
    
    # check that chromosome size was specified by one of L, r, or u.
    # if not specified, create them
    # TODO L = r[-1][-2] for cum maps means type(L) is numpy.int64? does this matter?
    if L is None: 
        if type(u) is list:
            L = u[-1][1]
        # if both r and u are maps default to r for total size
        if type(r) is list:
            L = r[-1][-2] 
        # if neither of the two previous statements were triggered. L is still None
        if L is None:
            raise ValueError("If u and r are constant values than the chrom. size, L, must be specified.")
    # if L is specified make sure it matches the lengths given in r and/or u
    elif isinstance(L, (int, float)):
        if type(u) is list:
            if u[-1][1] > L:
                raise ValueError('There is positive mutation rate at locations greater than the provided L')
        if type(r) is list:
            if r[-1][-2] != L:
                raise ValueError("Final position or r does not match the provided L. r must span the entire chromosome")
    else:
        raise ValueError('L must be None, int or float')
                
    # repeat with minPos
    if minPos is None:
        if type(u) is list:
            minPos = u[0][0] 
        # if both r and u are maps default to r for total size
        if type(r) is list:
            minPos = r[0][0]
        # if neither of the two previous statements were triggered. L is still None
        if minPos is None:
            minPos = 0
    # if L is specified make sure it matches the lengths given in r and/or u
    elif isinstance(minPos, (int, float)):
        if type(u) is list:
            if u[0][0] < minPos:
                raise ValueError('There is positive mutation rate at locations less than the provided minPos')
        if type(r) is list:
            if r[0][0] != minPos:
                raise ValueError("start position or r does not match the provided minPos. r must span the entire chromosome")
    else:
        raise ValueError('minPos must be None, int or float')
                
    # simplify u map or convert constant rate to a map
    if type(u) is list:
        u = combine_and_split_regions(u)
    else:
        nregions = math.ceil((L - minPos)/targetSize)
        regionSize = (L - minPos) / nregions
        u = [[minPos + regionSize * i, minPos + regionSize * (i+1), u] for i in range(nregions)]      
                
    # check that ss is of the correct form
    if type(ss) is not list:
        raise ValueError("ss must a list of selection coefficients")
    else:
        if np.ndim(ss) != 1:
            raise ValueError("ss must be a 1D list")
        else:
            if any([not isinstance(s, (int, float)) for s in ss]):
                raise ValueError("ss must a list of int or float selection coefficients")

    # if ps specified make sure it is of right form otherwise create it
    if ps is None:
        ps = [1 / len(ss)] * len(ss) 
    elif type(ps) is list:
        if np.ndim(ps) != 1:
            raise ValueError('ps must be a 1D list')
        if any([not isinstance(p, (int, float)) for p in ps]):
            raise ValueError("ps must a list of probabilitie that a new mutation has selection coefficient given by ss")
        if sum(ps) != 1:
            raise ValueError('ps must sum to 1')
        if len(ps) != len(ss):
            raise ValueError('ps and ss must be the same length')
    else: 
        raise ValueError('ps must be None or a list of probabilities')
    
    # define position of point masses as center of each u region
    positions = [(x + y)/2 for x,y,z in u]

    # convert r to a cumulative interpolating function
    if type(r) is float: 
        r = [[minPos,L,r]]
        r = make_cum_map(r)
    elif type(r) is list:
        r = make_cum_map(r, isCum)
        
    # make sure only one of cs or g is specified
    if cs is not None and g is not None:
        raise ValueError('either a size function, cs, or a demes graph, g, should be provided but not both')
        
    # ensure focalPos and tol are specified and numbers
    if not isinstance(tol, (int, float)):
        raise ValueError('tol must be float or int')
    if not isinstance(focalPos, (int, float)):
        raise ValueError('focalPos must be int or float')
        
    # if cs is specified 
    if cs is not None:
        # check cs is the correct shape. Also define reference census size
        if type(cs) is list:  
            if np.ndim(cs) != 1 or len(cs) > 1:
                raise ValueError('size function only supported for single populations. if more than one is needed, use demes graphs')
            elif isinstance(cs[0], (int, float)):
                censusSize = cs[0]
                cs = lambda t: [censusSize]
                # equilibrium requires no integration forward in time
                totalT = 0
            else:
                raise ValueError('list cs must have a single float or int entry')
        elif isinstance(cs, types.FunctionType):
            if type(cs(0)) is list:
                if np.ndim(cs(0)) != 1 or len(cs(0)) > 1:
                    raise ValueError('size function only supported for single populations. if more than one is needed, use demes graphs')
                elif isinstance(cs(0)[0],(int, float)):
                    censusSize = cs(0)[0]
                else:
                    raise ValueError('size function must return a list of length one with float or integer population size as a function of time')
                
                if totalT is None:
                    raise ValueError('if cs is a function totalT must be specified')
        
        # check totalT is a number
        if not isinstance(totalT, (int, float)):
            raise ValueError('totalT must be float or int')

        # define scaling function
        f, ancTime, ancNe = get_scaling_fun(positions, u, ss, ps, r, focalPos, censusSize, totalT, tol)
        # rescale census size function to be in units of 2 * Ne_bgs generations 
        rescaledcs = rescale_cs(cs, totalT, ancTime, ancNe, censusSize) 
        
        # size function with BGS
        g = lambda t: [x * y for x,y in zip(f(t), rescaledcs(t))]
        
        # check sample size
        if type(sample_size) is list:
            if np.ndim(sample_size) != 1 or len(sample_size) > 1:
                raise ValueError('size function only supported for single populations. if more than one is needed, use demes graphs')
            elif type(sample_size[0]) is not int:
                raise ValueError('sample_size must be a list of integers')

        # use two different moments functions depending on whether selection acts on the focal allele
        if focal_s is None:
            fs = moments.Demographics1D.snm(sample_size)
            fs.integrate(g, ancTime / 2 / ancNe)
        elif isinstance(focal_s, (int, float)): 
            fs = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sample_size[0], gamma = - 2 * ancNe * focal_s)) 
            fs.integrate(g, ancTime / 2 / ancNe, gamma = - 2 * ancNe * focal_s, h = 1/2, adapt_dt=True)
        else:
            raise ValueError('focal_s must be None, int or float')
        return fs, ancNe
    # if a graph is specified
    elif g is not None:
        if not isinstance(g, demes.demes.Graph):
            raise ValueError('g must be a demes graph')
        if g.time_units != "generations":
            g = g.in_generations()
        # demesdraw.tubes(g);
        
        # define ancestral parameters from g
        totalgen, censusSize = getOldestEpoch(g)
        
        # if totalT is specified make sure it matches
        if totalT is not None:
            if totalgen / 2 / censusSize != totalT:
                raise ValueError('provided totalT does not match the first demographic event specified in g')
        else:
            totalT = totalgen / 2 / censusSize
            
        # define scaling function
        f, ancTime, ancNe = get_scaling_fun(positions, u, ss, ps, r, focalPos, censusSize, totalT, tol)    
        
        # define gamma 
        if isinstance(focal_s, (int, float)):
            gamma = 2 * ancNe * focal_s
        elif focal_s is None:
            gamma = None
        else:
            raise ValueError('focal_s must be None, int or float')
            
        if type(sampled_demes) is None:
            raise ValueError('sampled_demes must be specified for demes graph')
       
        # compute SFS
        fs = SFS_bgs(
           g,
           sampled_demes=sampled_demes,
           sample_sizes=sample_size,
           theta = 1,
           bgs_Ne = ancNe,
           anc_gen = ancTime,
           scaling_fun=f,
           gamma = gamma
       )
        
        return fs, ancNe
    else:
        raise ValueError('cs or g must be specified.')



#########################
# OOA three populations #
#########################

# hardcoding some parameters
# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# sample_size = 40
# sampled_demes = ["CEU"]

# fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)
# for i in range(3):
#     for j in range(1):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curdemo = ["ooa.yaml"][j]
        
#         os.chdir("/media/nathan/T7/BGSdemo/parsedooaThreepopData")

#         simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False) 
#         projData = simData.project([sample_size])
        
#         os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop")

#         # test
#         demo = demes.load(curdemo)
#         if demo.time_units != "generations":
#             demo = demo.in_generations()
#         # demesdraw.tubes(demo);
        

#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 sampled_demes=sampled_demes,
#                                 g = demo) 
        
#         fs_neu = moments.Spectrum.from_demes(
#             curdemo, 
#             sampled_demes=sampled_demes, 
#             sample_sizes=[sample_size]
#         )

#         oldestEpoch, ancCensusSize = getOldestEpoch(demo)
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize   
        
#         ax[i].plot(fs, ".-", ms=8, lw=1, label="BGS")
#         ax[i].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
#         ax[i].plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
#         ax[i].set_title("s = " + str(curs) + ", demo = " + curdemo)
#         ax[i].set_yscale('log')
#         if np.logical_and(i == 2, j == 0):
#             ax[i].legend();

#######################
# OOA two populations #
#######################

# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# sample_size = 40
# sampled_demes = ["OOA"]

# fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)

# for i in range(3):
#     for j in range(1):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curdemo = ["ooaTwoPop.yaml"][j]
        
#         os.chdir("/media/nathan/T7/BGSdemo/parsedooaTwoPopData")

#         simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False) 
#         projData = simData.project([sample_size])
        
#         os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/twopop")

#         # test
#         demo = demes.load(curdemo)
#         if demo.time_units != "generations":
#             demo = demo.in_generations()
#         # demesdraw.tubes(demo);
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                        r = r,
#                                        focalPos = focalPos,
#                                        sample_size = [sample_size],
#                                        ss = [curs],
#                                        L = L,
#                                        sampled_demes=sampled_demes,
#                                        g = demo) 
        
#         fs_neu = moments.Spectrum.from_demes(
#             curdemo, 
#             sampled_demes=["OOA"], 
#             sample_sizes=[sample_size]
#         )
        
#         oldestEpoch, ancCensusSize = getOldestEpoch(demo)
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize   
        
#         ax[i].plot(fs, ".-", ms=8, lw=1, label="BGS")
#         ax[i].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
#         ax[i].plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
#         ax[i].set_title("s = " + str(curs) + ", demo = " + curdemo)
#         ax[i].set_yscale('log')
#         if np.logical_and(i == 2, j == 0):
#             ax[i].legend();

# jdata = parseJointData()
# sampled_demes = ["OOA", "YRI"]
# fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)

# for i in range(3):
#     for j in range(1):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curdemo = ["ooaTwoPop.yaml"][j]
        
#         os.chdir("/media/nathan/T7/BGSdemo/ooaTwoPopData/joint")

#         # simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
#         # simData = simData[0].to_numpy()
#         # simData = moments.Spectrum(simData,data_folded=False) 
#         # projData = simData.project([proj_size])
        
#         simData = jdata[i]
#         simData = moments.Spectrum(simData, data_folded = False)
#         projData = simData.project([sample_size, sample_size])
        
#         os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/twopop")

#         # test
#         demo = demes.load(curdemo)
#         if demo.time_units != "generations":
#             demo = demo.in_generations()
#         # demesdraw.tubes(demo);
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size, sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 sampled_demes=sampled_demes,
#                                 g = demo)

#         # normalizing so singletons have freq 1, cause thats all I can think of right now
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         # fs_neu = fs_neu * 8 * 1e-8 * censusSize   
#         # ds = ds * 8 * 1e-8 * censusSize
        
#         # moments.Plotting.plot_single_2d_sfs(fs)
#         # moments.Plotting.plot_single_2d_sfs(projData)plot_3d_spectrum_mayavi
        
#         moments.Plotting.plot_2d_comp_Poisson(fs, projData)

# #########################
# # OOA single population #
# #########################

# # hardcoding some parameters
# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# sample_size = 40

# # for curs in [1e-3, 5e-3, 1e-2]:
# #     for curN in [1e3, 5e3, 1e4]:
# fig, ax = plt.subplots(3, 1, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)

# for i in range(3):
#     for j in range(1):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curdemo = ["ooaSinglePop.yaml"][j]
        
#         os.chdir("/media/nathan/T7/BGSdemo/parsedooaSinglePopData")

#         simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False) 
#         projData = simData.project([sample_size])
        
#         os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa")

#         demo = demes.load(curdemo)
#         demo = demo.in_generations()
#         # demesdraw.tubes(demo);
        
#         oldestEpoch, ancCensusSize = getOldestEpoch(demo)
#         ds = reversedCensusFun(demo, oldestEpoch, ancCensusSize) 
#         cs =  lambda t: [ds(t)[0] * ancCensusSize] 
        
#         # fig, ax = plt.subplots(1, 1, figsize=(8, 4))
#         # ax.plot(np.arange(0,oldestEpoch / 2 / ancCensusSize, 0.001),[cs(t) for t in np.arange(0,oldestEpoch / 2 / ancCensusSize, 0.001)], "-", ms=8, lw=1, label="cs(t)")
#         # ax.set_xlabel("Time in past")
#         # ax.set_ylabel("cs(t)")
#         # ax.legend();   
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 cs = cs,
#                                 totalT = oldestEpoch / 2 / ancCensusSize)
        
#         fs_neu = moments.Demographics1D.snm([sample_size])
#         fs_neu.integrate(ds, oldestEpoch / 2 / ancCensusSize)
        
#         sampled_demes = ["CEU"]

#         fs_demes = moments.Spectrum.from_demes(
#             curdemo, sampled_demes=sampled_demes, sample_sizes=[sample_size]
#         )
        
#         bgs_demes, ancNe2 = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 sampled_demes=sampled_demes,
#                                 g = demo) 
        

#         # normalizing so singletons have freq 1, cause thats all I can think of right now
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize   
#         fs_demes = fs_demes * 8 * 1e-8 * ancCensusSize
#         bgs_demes = bgs_demes * 8 * 1e-8 * ancNe2   
        
#         ax[i].plot(fs, ".-", ms=8, lw=1, label="BGS")
#         ax[i].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
#         ax[i].plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
#         ax[i].plot(fs_demes, "*-", ms=8, lw=1, label="demes")
#         ax[i].plot(bgs_demes, "*-", ms=8, lw=1, label="bgs_demes")        
#         ax[i].set_title("s = " + str(curs) + ", demo = " + curdemo)
#         ax[i].set_yscale('log')
#         if np.logical_and(i == 2, j == 0):
#             ax[i].legend();
            
#################################
# bottleneck neutral focal site #
#################################

# # hardcoding some parameters
# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# sample_size = 40

# fig, ax = plt.subplots(3, 2, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)
# for i in range(3):
#     for j in range(2):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curdemo = ["1k.yaml", "5k.yaml"][j]
        
        
#         os.chdir("/media/nathan/T7/BGSdemo/parsedbottleneckData")
        
#         simData = pd.read_csv(str(curs) + "_" + str(curdemo) + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False) 
#         projData = simData.project([sample_size])
        
#         os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims")

#         demo = demes.load(curdemo)
#         # demesdraw.tubes(demo);
        
#         oldestEpoch, ancCensusSize = getOldestEpoch(demo)
#         ds = reversedCensusFun(demo, oldestEpoch, ancCensusSize) 
#         cs =  lambda t: [ds(t)[0] * ancCensusSize] 
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 cs = cs,
#                                 totalT = oldestEpoch / 2 / ancCensusSize)
        
#         fs_neu = moments.Demographics1D.snm([sample_size])
#         fs_neu.integrate(ds, oldestEpoch / 2 / ancCensusSize)
        
#         sampled_demes = ["B"]

#         fs_demes = moments.Spectrum.from_demes(
#             curdemo, sampled_demes=sampled_demes, sample_sizes=[sample_size]
#         )
        
#         bgs_demes, ancNe2 = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 sampled_demes=sampled_demes,
#                                 g = demo) 
        

#         # normalizing so singletons have freq 1, cause thats all I can think of right now
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize   
#         fs_demes = fs_demes * 8 * 1e-8 * ancCensusSize
#         bgs_demes = bgs_demes * 8 * 1e-8 * ancNe2
        
#         ax[i,j].plot(fs, ".-", ms=8, lw=1, label="BGS")
#         ax[i,j].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
#         ax[i,j].plot(fs_neu, "+-", ms=8, lw=1, label="neutral")
#         ax[i,j].plot(fs_demes, "*-", ms=8, lw=1, label="demes")
#         ax[i,j].plot(bgs_demes, "*-", ms=8, lw=1, label="bgs_demes")
#         ax[i,j].set_title("s = " + str(curs) + ", demo = " + curdemo)
#         ax[i,j].set_yscale('log')
#         if np.logical_and(i == 2, j == 1):
#             ax[i,j].legend();
            
# ###################################
# # equilibrium selected focal site #
# ###################################

# os.chdir("/media/nathan/T7/BGSdemo/parsedMoreSelData")

# # hardcoding some parameters
# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# sample_size = 500
# proj_size = 40

# curwi = "wi5e4.csv"

# fig, ax = plt.subplots(3, 3, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.text(0.5, 0.96, curwi, ha='center')
# fig.subplots_adjust(hspace = .25)
# for i in range(3):
#     for j in range(3):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curN = [1e3, 5e3, 1e4][j]
#         simData = pd.read_csv(str(curs) + "_" + str(int(curN)) + "_" + curwi + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False) / 11 / 1000
#         projData = simData.project([proj_size])
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 cs = [curN],
#                                 focal_s = curs)
        
#         fs_indep = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sample_size, gamma = - 2 * curN * curs))

#         if curwi == "wi5e4.csv":
#             fs = fs *  4 * 2 * 5e4 * 1e-8 * ancNe
#             fs_indep = fs_indep * 4 * 2 * 5e4 * 1e-8 * curN
#             projData = projData 
#         if curwi == "wi1e5.csv":
#             fs = fs *  4 * 2 * 1e5 * 1e-8 * ancNe
#             fs_indep = fs_indep * 4 * 2 * 1e5 * 1e-8 * curN
#             projData = projData 
        
#         fs = fs.project([proj_size])
#         fs_indep = fs_indep.project([proj_size])
        
#         ax[i,j].plot(fs[0:21], ".-", ms=8, lw=1, label="BGS")
#         ax[i,j].plot(projData[0:21], "x-", ms=8, lw=1, label="fwdpy")
#         ax[i,j].plot(fs_indep[0:21], "+-", ms=8, lw=1, label="single locus")
#         ax[i,j].set_title("s = " + str(curs) + ", N = " + str(int(curN)))
#         ax[i,j].set_yscale('log')
#         if np.logical_and(i == 2, j == 2):
#             ax[i,j].legend();

# ###################
# # equilibrium dfe #
# ###################

# os.chdir("/media/nathan/T7/BGSdemo/dfeAFS")
# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# sample_size = 40


# fig, ax = plt.subplots(3, 3, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)
# for i in range(3):
#     for j in range(3):
#         p = [0.25, 0.5, 0.75][i]
#         curN = [1e3, 5e3, 1e4][j]
        
#         simData = pd.read_csv(str(p) + "_" + str(int(curN)) + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False)
#         projData = simData.project([sample_size])
                
#         ps = [p,1-p]
#         ss = [0.01, 0.005]
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = ss,
#                                 ps = ps,
#                                 L = L,
#                                 cs = [curN])
        
#         fs_neu = moments.Demographics1D.snm([sample_size])
        
#         # normalizing based on N and mu 
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         fs_neu = fs_neu * 8 * 1e-8 * curN
        
#         ax[i,j].plot(fs, ".-", ms=8, lw=1, label="BGS")
#         ax[i,j].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
#         ax[i,j].plot(fs_neu, "+-", ms=8, lw=1, label="SNM")
#         ax[i,j].set_title("p = " + str(p) + ", N = " + str(int(curN)))
#         ax[i,j].set_yscale('log')
#         if np.logical_and(i == 2, j == 2):
#             ax[i,j].legend();

##############################################
# equilibrium human maps neutral focal locus #
##############################################
os.chdir("/media/nathan/T7/BGSdemo/humanData")
recName = "YRI_recombination_map_hg38_chr_22.bed"
mutName = "roulette_tbl_chr22.csv"
exonName = "exons_chr22.bed"
recMap = read_rec_map(recName)
mutMap = read_mut_rates(mutName)
exonMap = read_exon_map(exonName)

recMap = simplify_rate_map(recMap)
mutMap = simplify_rate_map(mutMap)

exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap) 

os.chdir("/media/nathan/T7/BGSdemo/parsedequilHuman5027")
focalPos = 50270000
sample_size = 40

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
        
        fs, ancNe = bgs_wrapper(u = exonMutMap,
                                r = recMap,
                                focalPos = focalPos,
                                sample_size = [sample_size],
                                ss = [curs],
                                cs = [curN])
        
        fs_neu = moments.Demographics1D.snm([sample_size])
        
        # normalizing based on N and mu 
        # should it be ancNe
        fs = fs * 8 * 1e-8 * ancNe
        projData = projData * 1e-8
        fs_neu = fs_neu * 8 * 1e-8 * curN        
        
        ax[i,j].plot(fs, "-", ms=8, lw=1, label="BGS")
        ax[i,j].plot(projData, "-", ms=8, lw=1, label="fwdpy")
        ax[i,j].plot(fs_neu, "-", ms=8, lw=1, label="SNM")
        ax[i,j].set_title("s = " + str(curs) + ", N = " + str(int(curN)))
        ax[i,j].set_yscale('log')
        ax[i,j].annotate("B=" + str(round(fs.pi()/fs_neu.pi(),3)),xy=(sample_size / 10, max(fs[1:-1]) * 0.8),size="x-large")
        ax[i,j].annotate("B_sim=" + str(round(projData.pi()/fs_neu.pi(),3)),xy=(sample_size / 10, max(fs[1:-1]) * 0.55),size="x-large")
        if np.logical_and(i == 2, j == 2):
            ax[i,j].legend();
            
# ###################################
# # equilibrium neutral focal locus #
# ###################################
# os.chdir("/media/nathan/T7/BGSdemo/equilAFS")
# sample_size = 40
# u = 1e-8
# r = 1e-8
# L = 1e6
# focalPos = 5e5
# fig, ax = plt.subplots(3, 3, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)
# for i in range(3):
#     for j in range(3):
#         curs = [1e-3, 5e-3, 1e-2][i]
#         curN = [1e3, 5e3, 1e4][j]
#         simData = pd.read_csv(str(curs) + "_" + str(int(curN)) + ".csv", header = None)
#         simData = simData[0].to_numpy()
#         simData = moments.Spectrum(simData,data_folded=False)
#         projData = simData.project([sample_size])
        
#         fs, ancNe = bgs_wrapper(u = u,
#                                 r = r,
#                                 focalPos = focalPos,
#                                 sample_size = [sample_size],
#                                 ss = [curs],
#                                 L = L,
#                                 cs = [curN])
        
#         fs_neu = moments.Demographics1D.snm([sample_size])
        
#         # normalizing based on N and mu 
#         # should it be ancNe
#         fs = fs * 8 * 1e-8 * ancNe
#         projData = projData * 1e-8
#         fs_neu = fs_neu * 8 * 1e-8 * curN
        
#         ax[i,j].plot(fs, ".-", ms=8, lw=1, label="BGS")
#         ax[i,j].plot(projData, "x-", ms=8, lw=1, label="fwdpy")
#         ax[i,j].plot(fs_neu, "+-", ms=8, lw=1, label="SNM")
#         ax[i,j].set_title("s = " + str(curs) + ", N = " + str(int(curN)))
#         ax[i,j].set_yscale('log')
#         if np.logical_and(i == 2, j == 2):
#             ax[i,j].legend();
            
# ########
# # misc #
# ########
            
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[g(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="g(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("g(t)")
# ax.legend(); 
     
# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[rescaledcs(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="rescaledcs(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("rescaledcs(t)")
# ax.legend();     

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,totalT, 0.001),[cs(t) for t in np.arange(0,totalT, 0.001)], "-", ms=8, lw=1, label="cs(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("cs(t)")
# ax.legend();   

# [cs(0),rescaledcs(0)]

# [cs(totalT), rescaledcs(ancTime / 2/ ancNe)]

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime / 2 / ancNe, 0.001),[f(t) for t in np.arange(0,ancTime / 2 / ancNe, 0.001)], "-", ms=8, lw=1, label="B(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend();       

        