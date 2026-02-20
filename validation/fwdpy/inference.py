import moments
import numpy as np
import math
import matplotlib.pylab as plt
import os
import pandas as pd
import types
import demes
import demesdraw
from datetime import datetime
import copy
from scipy import interpolate
import warnings
# from __future__ import annotations
from typing import List
from scipy.special import gamma, gammaincc, exp1 # exp1(x) = E1(x) = Γ(0, x) (upper incomplete gamma at a=0)
from sklearn_extra.cluster import KMedoids
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score





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

def rescaledPointMassContribution(scaledu, s, t, r, ancestralNe, ancTime):
    return - scaledu / s * (s / (r + s) * (1 - math.exp(- r * (ancTime - t * 2 * ancestralNe) - s * (ancTime - t * 2 * ancestralNe))))**2

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

def discretize_deleterious_gamma_dfe_mean_shape(
    mean,
    shape,
    n_bins,
    representative = "conditional_mean"):
    """
    Discretize a deleterious Gamma DFE parameterized by (mean, shape)
    into a list [[p, s], ...], with s > 0.

    Gamma parameterization:
        S ~ Gamma(k=shape, theta=mean/shape)
        E[S] = mean

    Binning scheme:
        Equal-probability (quantile) bins, with a final [q, +inf) tail.
        No truncation.

    Parameters
    ----------
    mean : float
        Mean deleterious magnitude E[S] > 0.
    shape : float
        Gamma shape k > 0.
    n_bins : int
        Number of bins >= 1.
    representative : {"conditional_mean", "mid_quantile"}
        Representative value per bin.

    Returns
    -------
    dfe : list of [p, s]
        Discrete DFE with probabilities p and selection coefficients s < 0.

    Notes
    -----
    Requires SciPy.
    """
    if mean <= 0:
        raise ValueError("mean must be > 0")
    if shape <= 0:
        raise ValueError("shape must be > 0")
    if n_bins < 1 or int(n_bins) != n_bins:
        raise ValueError("n_bins must be an integer >= 1")

    representative = representative.lower()
    if representative not in {"conditional_mean", "mid_quantile"}:
        raise ValueError("representative must be 'conditional_mean' or 'mid_quantile'")

    try:
        from scipy.stats import gamma as gamma_dist
    except ImportError as e:
        raise ImportError("This function requires SciPy (pip install scipy).") from e

    k = float(shape)
    mean = float(mean)
    theta = mean / k
    n = int(n_bins)

    # Single-bin edge case
    if n == 1:
        return [[1.0, mean]]

    dist_k = gamma_dist(a=k, scale=theta)
    dist_k1 = gamma_dist(a=k + 1.0, scale=theta)

    p = 1.0 / n

    # Interior quantile cutpoints: F^{-1}(i/n)
    cutpoints = [dist_k.ppf(i / n) for i in range(1, n)]

    bounds = [(0.0, cutpoints[0])]
    for i in range(1, len(cutpoints)):
        bounds.append((cutpoints[i - 1], cutpoints[i]))
    bounds.append((cutpoints[-1], float("inf")))

    dfe = []

    if representative == "mid_quantile":
        for i in range(n):
            u = (i + 0.5) / n
            s_mag = dist_k.ppf(u)
            dfe.append([p, -float(s_mag)])
        return dfe

    # Conditional-mean representative
    for (a, b) in bounds:
        if b == float("inf"):
            denom = dist_k.sf(a)
            numer = dist_k1.sf(a)
        else:
            denom = dist_k.cdf(b) - dist_k.cdf(a)
            numer = dist_k1.cdf(b) - dist_k1.cdf(a)

        if denom <= 0:
            raise RuntimeError(f"Nonpositive bin mass for bin [{a}, {b}].")

        mean_mag = (theta * k) * (numer / denom)
        dfe.append([p, float(mean_mag)])

    return dfe
 
# u = exonMutMap
# r = recMap
# focalPos = focalPos
# sample_size = [sample_size]
# ss = ss
# sampled_demes=sampled_demes
# g = demo

# cs = None
# totalT = None
# L = None
# ps = None
# targetSize = 1e4
# tol = 1e-4
# minPos = None
# focal_s = None
# eval_thres = "step"
# R_cutoff = 1e-3

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
                eval_thres = "step",
                R_cutoff = 1e-3,
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
        # f, ancTime, ancNe = get_scaling_fun(positions, u, ss, ps, r, focalPos, censusSize, totalT, tol)
        f, ancTime, ancNe = get_scaling_fun_2(u, ss, ps, r, focalPos, censusSize, totalT, tol, eval_thres = eval_thres, R_cutoff = R_cutoff)
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
        # f, ancTime, ancNe = get_scaling_fun(positions, u, ss, ps, r, focalPos, censusSize, totalT, tol)    
        f, ancTime, ancNe = get_scaling_fun_2(u, ss, ps, r, focalPos, censusSize, totalT, tol, eval_thres = eval_thres, R_cutoff = R_cutoff)
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

# u = exonMutMap
# r = recMap
# focal_positions = focal_positions
# ss = ss

# L = None
# ps = None
# targetSize = 1e4
# tol = 1e-4
# minPos = None
# dt = 1e4
# R_cutoff = 1e-3
# n_clusters = 20
# equil_thres = "step"

def cluster_scale_funs(u,
                r,
                focal_positions,
                ss,
                L = None,
                ps = None,
                targetSize = 1e4,
                tol = 1e-4,
                minPos = None,
                dt = 1e4,
                R_cutoff = 1e-3,
                n_clusters = 20,
                equil_thres = "step"):
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
    # positions = [(x + y)/2 for x,y,z in u]

    # convert r to a cumulative interpolating function
    if type(r) is float: 
        r = [[minPos,L,r]]
        r = make_cum_map(r)
    elif type(r) is list:
        r = make_cum_map(r, isCum)
            
    # ensure focalPos and tol are specified and numbers
    if not isinstance(tol, (int, float)):
        raise ValueError('tol must be float or int') 
    if (type(focal_positions) is not list and type(focal_positions) is not np.ndarray) or not isinstance(focal_positions[0], (int, float)):
        raise ValueError('focal_positions must be a list of int or float')

    if equil_thres == "step":
        # compute grid of B(t) for each focal position, every dt generations until |B(t) - B(t-dt)| < tol
        B_grid = get_Bs(u, ss, ps, r, focal_positions, tol, dt, R_cutoff)
        mxt = max([len(x) for x in B_grid])
        t_vals = [dt * j for j in range(mxt)]
    elif equil_thres == "inf":
        B_grid, t_vals = get_Bs_2(u, ss, ps, r, focal_positions, tol, dt, R_cutoff)
        mxt = max([len(x) for x in t_vals])
        t_vals = [x for x in t_vals if len(x) == mxt][0]

    # pad ends of B functions that reached quil early
    B_rect = ragged_to_rect_repeat_last(B_grid)
    # compute distance matrix
    D = pairwise_ssd(B_rect)
    
    labels, medoids = cluster_B_kmedoids(D, K=n_clusters)
    
    return B_rect, labels, medoids, t_vals

    
    # directly computing features        
    startTime = datetime.now()
    B_inf, half_time, half_index = get_Binf_and_halftime(u, ss, ps, r, focal_positions, dt, R_cutoff)
    endTime = datetime.now()
    endTime - startTime
            
    features = np.column_stack([B_inf, half_time])
    labels, medoids, model, X = kmedoids_on_features(features, K=20)
    
    sil_feat = silhouette_score(X, labels, metric='euclidean')
    ari_feat = ari_stability_test(lambda X: cluster_features(X, K=20), features)
    ari_feat.mean()
    plot_kmedoids_clusters(B_rect, labels, medoids)

    features = np.column_stack([B_inf, np.log10(half_time)])
    labels, medoids, model, X = kmedoids_on_features(features, K=20)
    sil_feat_log = silhouette_score(X, labels, metric='euclidean')
    ari_feat_log = ari_stability_test(lambda X: cluster_features(X, K=20), features)
    ari_feat_log.mean()

    ###########
    
    startTime = datetime.now()
    B_grid_old = get_Bs(u, ss, ps, r, focal_positions, tol, dt, R_cutoff)
    endTime = datetime.now()
    endTime - startTime

    
    # from  full B matrix
    startTime = datetime.now()
    B_grid, t_vals = get_Bs_2(u, ss, ps, r, focal_positions, 5e-2, dt, R_cutoff)
    endTime = datetime.now()
    endTime - startTime
    B_rect = ragged_to_rect_repeat_last(B_grid)
    D = pairwise_ssd(B_rect)
    
    mxt = max([len(x) for x in t_vals])
    t_vals = [x for x in t_vals if len(x) == mxt][0]
    
    labels, medoids = cluster_B_kmedoids(D, K=n_clusters)

    plot_kmedoids_clusters(B_grid, labels, medoids, t_vals=t_vals)

    B_inf = get_Binf(u, ss, ps, r, focal_positions, dt, R_cutoff)
    B_equil = [x[-1] for x in B_grid]
    diff = B_equil - B_inf
    np.mean(diff)
    np.median(diff)
    max(diff)
    min(diff)
    size = [(y - x)  for x,y,z in u]
    
    ex = np.asarray(u, dtype=np.float64)   # (E,3)
    LL = ex[:,0]                               # (E,) start positions
    LU = ex[:,1]                               # (E,) end positions
    MU = ex[:,2]                               # (E,) local mu per bp
    
    positions = 0.5*(LL + LU)  
    
    recDist = getRecDist_3(positions, r, focal_positions)
    
    focal_idx=np.nonzero([x == max(diff) for x in diff])
    big_diffs = recDist[np.nonzero([x == max(diff) for x in diff])]
    nearby = [[y,z] for x in big_diffs for y,z in zip(x,size) if y < R_cutoff]
    
    plt.figure(figsize=(6,4))
    plt.plot(B_rect[focal_idx[0][0]], color='red', lw=3)
    plt.xlabel("time index")
    plt.ylabel("B(t)")
    plt.tight_layout()
    plt.show()
    
    
    
    max([x[1] for x in nearby])
    labels, medoids = cluster_B_kmedoids(D, K=20)

    sil_curve = silhouette_score(D, labels, metric='precomputed')
    ari_curve = ari_stability_precomputed(lambda Dsub: cluster_curves(Dsub, K=20), D)
    ari_curve.mean()
    plot_kmedoids_clusters(B_rect, labels, medoids)

    sil_vals = []
    ari_vals = []
    for kk in range(2,40):
        features = np.column_stack([B_inf, half_time])
        labels, medoids, model, X = kmedoids_on_features(features, K=kk)
        
        sil_feat = silhouette_score(X, labels, metric='euclidean')
        ari_feat = ari_stability_test(lambda X: cluster_features(X, K=kk), features)

        features = np.column_stack([B_inf, np.log10(half_time)])
        labels, medoids, model, X = kmedoids_on_features(features, K=kk)
        sil_feat_log = silhouette_score(X, labels, metric='euclidean')
        ari_feat_log = ari_stability_test(lambda X: cluster_features(X, K=kk), features)
        
        labels, medoids = cluster_B_kmedoids(D, K=kk)

        sil_curve = silhouette_score(D, labels, metric='precomputed')
        ari_curve = ari_stability_precomputed(lambda Dsub: cluster_curves(Dsub, K=kk), D)
        
        sil_vals.append([sil_feat, sil_feat_log, sil_curve])
        ari_vals.append([ari_feat.mean(), ari_feat_log.mean(), ari_curve.mean()])
    
    sil_feat = [x[0] for x in sil_vals]
    max_idx = np.nonzero(sil_feat == max(sil_feat))[0][0]
    print("feat: optimal K is " + str(max_idx + 2) + " with silhouette score " + str(max(sil_feat)))
    
    ari_feat = [x[0] for x in ari_vals]
    max_idx = np.nonzero(ari_feat == max(ari_feat))[0][0]
    print("feat: optimal K is " + str(max_idx + 2) + " with adjusted rand index " + str(max(ari_feat)))
    
    sil_curve = [x[2] for x in sil_vals]
    max_idx = np.nonzero(sil_curve == max(sil_curve))[0][0]
    print("curve: optimal K is " + str(max_idx + 2) + " with silhouette score " + str(max(sil_curve)))
    
    ari_curve = [x[2] for x in ari_vals]
    max_idx = np.nonzero(ari_curve == max(ari_curve))[0][0]
    print("curve: optimal K is " + str(max_idx + 2) + " with adjusted rand index " + str(max(ari_curve)))
    # B_inf, half_time, half_index = compute_Binf_and_halftime(B_grid, dt = 1e3)
    
    # features = np.column_stack([B_inf, half_time])

    # labels, medoids, model, X = kmedoids_on_features(features, K=20)
    

    
    # startTime = datetime.now()
    # best_K, labels, medoids, scores = choose_k_kmedoids_silhouette(D, k_min=2, k_max=25)
    # endTime = datetime.now()
    # endTime - startTime
    
    

# from sklearn.decomposition import PCA

# X = PCA(2).fit_transform(B_grid)

# plt.scatter(X[:,0], X[:,1], s=5)
# plt.xlabel("PC1")
# plt.ylabel("PC2")
# plt.show()

# plt.scatter(X[:,0], X[:,1], c=labels, s=5)
# plt.show()

def cluster_curves(D, K=2):
    model = KMedoids(n_clusters=K, metric="precomputed", random_state=0)
    return model.fit_predict(D)

def cluster_features(features, K=2):
    X = StandardScaler().fit_transform(features)
    model = KMedoids(n_clusters=K, metric="euclidean", random_state=0)
    return model.fit_predict(X)

def ari_stability_precomputed(cluster_func, D, n_runs=20, subsample_frac=0.8, random_state=0):
    """
    ARI stability for clustering functions that take a PRECOMPUTED distance matrix D (NxN).
    """
    rng = np.random.default_rng(random_state)
    N = D.shape[0]
    if D.shape[1] != N:
        raise ValueError("D must be square (N,N).")

    labels_full = cluster_func(D)
    ari_scores = np.empty(n_runs, dtype=np.float64)

    for r in range(n_runs):
        idx = rng.choice(N, size=int(subsample_frac * N), replace=False)
        idx.sort()

        D_sub = D[np.ix_(idx, idx)]          # <-- key fix: square submatrix
        labels_sub = cluster_func(D_sub)

        ari_scores[r] = adjusted_rand_score(labels_full[idx], labels_sub)

    return ari_scores

def ari_stability_features(cluster_func, X, n_runs=20, subsample_frac=0.8, random_state=0):
    """
    ARI stability for clustering functions that take a FEATURE matrix X (N,d).
    """
    rng = np.random.default_rng(random_state)
    N = X.shape[0]

    labels_full = cluster_func(X)
    ari_scores = np.empty(n_runs, dtype=np.float64)

    for r in range(n_runs):
        idx = rng.choice(N, size=int(subsample_frac * N), replace=False)
        idx.sort()

        X_sub = X[idx]
        labels_sub = cluster_func(X_sub)

        ari_scores[r] = adjusted_rand_score(labels_full[idx], labels_sub)

    return ari_scores

def ari_stability_test(cluster_func, data, n_runs=20, subsample_frac=0.8, random_state=0):
    """
    cluster_func: function that takes data_subset and returns labels
    data: full dataset (features or distance matrix)
    n_runs: number of subsampling trials
    subsample_frac: fraction of data to keep
    """
    rng = np.random.default_rng(random_state)
    N = len(data)

    # Full clustering
    labels_full = cluster_func(data)

    ari_scores = []

    for _ in range(n_runs):
        idx = rng.choice(N, size=int(subsample_frac * N), replace=False)

        data_sub = data[idx] if data.ndim == 2 else data[np.ix_(idx, idx)]

        labels_sub = cluster_func(data_sub)

        # Compare only overlapping points
        ari = adjusted_rand_score(labels_full[idx], labels_sub)
        ari_scores.append(ari)

    return np.array(ari_scores)

def kmedoids_on_features(features, K, random_state=0):
    X = StandardScaler().fit_transform(features)  # recommended
    model = KMedoids(
        n_clusters=K,
        metric="euclidean",
        init="k-medoids++",
        random_state=random_state
    )
    labels = model.fit_predict(X)
    medoids = model.medoid_indices_
    return labels, medoids, model, X

def compute_Binf_and_halftime(B_rect, dt = 1e3):
    """
    B_rect: (F,T) padded rectangular array, monotone decreasing to B_inf
    dt: time step size (in generations, or whatever units you're using)

    Returns:
      B_inf: (F,)
      half_time: (F,) in same time units as dt
      half_index: (F,) indices on the grid (0..T-1)
    """
    B = np.asarray(B_rect, dtype=np.float64)
    F, T = B.shape

    B_inf = B[:, -1]                         # (F,)
    thresh = 0.5 * (1.0 + B_inf)             # (F,)

    # first index where B <= threshold
    # If a curve never crosses (shouldn't if it decays), set to T-1.
    crossed = B <= thresh[:, None]           # (F,T) boolean
    half_index = np.where(crossed.any(axis=1), crossed.argmax(axis=1), T - 1)

    half_time = half_index.astype(np.float64) * float(dt)
    return B_inf, half_time, half_index

def choose_k_kmedoids_silhouette(D, k_min=2, k_max=20, random_state=0):
    """
    Returns (best_K, best_labels, best_medoids, scores_dict)
    D: (F,F) distance matrix (squared distances are OK; monotone transform preserves ordering reasonably)
    """
    best = (-np.inf, None, None, None)  # (score, K, labels, medoids)
    scores = {}

    for K in range(k_min, k_max + 1):
        kmed = KMedoids(
            n_clusters=K,
            metric="precomputed",
            init="k-medoids++",
            random_state=random_state
        )
        labels = kmed.fit_predict(D)

        # Silhouette requires at least 2 clusters and not all points in one cluster
        if len(np.unique(labels)) < 2:
            scores[K] = -np.inf
            continue

        score = silhouette_score(D, labels, metric="precomputed")
        scores[K] = score

        if score > best[0]:
            best = (score, K, labels, kmed.medoid_indices_)

    best_score, best_K, best_labels, best_medoids = best
    return best_K, best_labels, best_medoids, scores

def plot_kmedoids_clusters(B, labels, medoids, max_curves=50, t_vals = None, rect = False):
    """
    B        : (F,T) matrix of B(t) curves
    labels   : cluster labels for each curve
    medoids  : indices of medoid curves (length K)
    """

    K = len(medoids)
    if t_vals is None:
        t_vals = range(len(B[0]))

    for c in range(K):
        idx = np.where(labels == c)[0]
        m = medoids[c]

        if not rect:
            B_c = np.array(B, dtype=object)[idx[:max_curves]].tolist()
            B_c.append(B[m])
            B_c = ragged_to_rect_repeat_last(B_c)
            t_vals_c = t_vals[:len(B_c[0])]
        else:
            B_c = B[idx[:max_curves]]
            B_c.append(B[m])
        
        plt.figure(figsize=(6,4))

        # plot some member curves (light)
        for i in range(len(B_c)):
            plt.plot(t_vals_c, B_c[i], color='gray', alpha=0.25)

        # plot medoid (bold)
        plt.plot(t_vals_c, B_c[-1], color='red', lw=3)

        plt.title(f"Cluster {c} (n={len(idx)})")
        plt.xlabel("time index")
        plt.ylabel("B(t)")
        plt.tight_layout()
        plt.show()

def cluster_B_kmedoids(D, K, random_state=0):
    kmed = KMedoids(
        n_clusters=K,
        metric='precomputed',
        init='k-medoids++',
        random_state=random_state
    )

    labels = kmed.fit_predict(D)
    medoids = kmed.medoid_indices_

    return labels, medoids

def pairwise_ssd(B):
    """
    B: (F, T) array
    returns: (F, F) matrix D where D[i,k] = sum_j (B[i,j] - B[k,j])**2
    """

    # squared norms of each row: (F,)
    norms2 = np.einsum('ij,ij->i', B, B)

    # Gram matrix: (F,F) = B @ B.T
    G = B @ B.T

    # pairwise squared distances
    D = norms2[:, None] + norms2[None, :] - 2.0 * G

    # numerical cleanup: tiny negative values to 0
    np.maximum(D, 0.0, out=D)
    np.fill_diagonal(D, 0) # correct for small nonzer values due to numerical precision
    return D

def ragged_to_rect_repeat_last(all_bs, dtype=np.float64):
    F = len(all_bs)
    lengths = np.fromiter((len(row) for row in all_bs), count=F, dtype=np.int64)
    T_max = int(lengths.max())

    out = np.empty((F, T_max), dtype=dtype)

    for i, row in enumerate(all_bs):
        r = np.asarray(row)
        n = r.size
        out[i, :n] = r
        out[i, n:] = r[-1]   # repeat equilibrium value

    return out

def getRecDist_3(positions, r, focal_positions):
    Rpos = r(positions)   # shape (M,)
    Rfoc = r(focal_positions)   # shape (F,)
    
    bigR = np.abs(Rpos[None, :] - Rfoc[:, None])           # (F, M)
    return 0.5 * (1.0 - np.exp(-2.0 * bigR))            # (F, M)

# u = MU[j]
# ll = LL[j]
# lu = LU[j]
# focalPos = fp

def internal_containing_vec(s, u, l, t, r):
    """
    internalContribution(s,u,l,t,r)
          = (u/r) * [  (A(s) / (l*r + s))
                      - 2*s*t*Γ(0, s*t) + 2*s*t*Γ(0, 2*s*t)
                      + 2*s*t*(Γ(0, (l*r + s)*t) - Γ(0, 2*(l*r + s)*t))  ]

    where:
        A(s) = -(l*r)
               + s*exp(-2*(l*r + s)*t)
               - 2*s*exp(-(l*r + s)*t)
               - (l*r + s)*exp(-2*s*t)
               + 2*(l*r + s)*exp(-s*t)
    """
    if t <= 0:
        return np.zeros_like(s, dtype=np.float64)

    lr = l * r                     # scalar
    denom = lr + s                 # (K,)

    # ---- Exponential pieces (K,) ----
    e_lr_s  = np.exp(-denom * t)
    e2_lr_s = np.exp(-2.0 * denom * t)
    e_s     = np.exp(-s * t)
    e2_s    = np.exp(-2.0 * s * t)

    # ---- Rational/exponential bracket A(s)/(lr+s) ----
    frac_term = (-(lr)
                 + s * e2_lr_s
                 - 2.0 * s * e_lr_s
                 - denom * e2_s
                 + 2.0 * denom * e_s) / denom

    # ---- Γ(0, x) terms (upper incomplete gamma at a=0) ----
    #  inc_gamma(0, x) equals exp1(x).
    x1 = s * t
    x2 = 2.0 * s * t
    x3 = denom * t
    x4 = 2.0 * denom * t

    gamma_term = (-2.0 * s * t * exp1(x1)
                  + 2.0 * s * t * exp1(x2)
                  + 2.0 * s * t * (exp1(x3) - exp1(x4)))

    return (u * (frac_term + gamma_term)) / r

# Source - https://stackoverflow.com/a/52186952
# Posted by MuellerSeb
# Retrieved 2026-02-08, License - CC BY-SA 4.0

def inc_gamma(a, x):
    return exp1(x) if a == 0 else gamma(a)*gammaincc(a, x)

def make_external_logB_evaluator(recDist_i, scaledu, ss, ps, ext_mask):
    U = scaledu[ext_mask][:, None]   # (Mext,1) 
    R = recDist_i[ext_mask][:, None] # (Mext,1)
    S_row = ss[None, :]              # (1,K)
    P_row = ps[None, :]              # (1,K)

    D = R + S_row                                            # (Mext,k)
    # W = (-u*s / D^2) * p  (all t-independent, includes p)
    invD = np.divide(1.0, D, out=np.zeros_like(D), where=(D != 0.0))
    W = (-U) * S_row * (invD * invD) * P_row

    # reusable temp arrays to avoid reallocations each time
    E = np.empty_like(D)  # exp(-t*D)
    G = np.empty_like(D)  # 1 - exp(-t*D)

    def ext_logB(t):
        t = float(t)
        # E = exp(-t*D)
        np.multiply(D, -t, out=E)
        np.exp(E, out=E)
        # G = 1 - E
        np.subtract(1.0, E, out=G)
        # logB = sum(W * G^2)
        np.multiply(G, G, out=G)
        return float(np.sum(W * G))

    return ext_logB

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,ancTime, 1000),[B(scaledu, ss, ps, t, rd) for t in np.arange(0,ancTime,1000)], "-", ms=8, lw=1, label="B(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend(); 

# fig, ax = plt.subplots(1, 1, figsize=(8, 4))
# ax.plot(np.arange(0,totalgen, 100),[B(scaledu, ss, ps, t, rd) for t in np.arange(0,totalgen,100)], "-", ms=8, lw=1, label="B(t)")
# ax.set_xlabel("Time in past")
# ax.set_ylabel("B(t)")
# ax.legend(); 

def exon_contains_focal_log(ss, ps, u, ll, lu, focalPos, rbp, t):
    """
    Log-contribution for an exon [ll, lu] that CONTAINS the focal site.
    """
    if t <= 0:
        return 0.0
    
    # distance to endpoints
    lL = focalPos - ll
    lR = lu - focalPos
    
    contrib = internal_containing_vec(ss, u, lL, t, rbp) + internal_containing_vec(ss, u, lR, t, rbp)
    return float(np.sum(ps * contrib))

def get_scaling_fun_2(u, ss, ps, r, focalPos, censusSize, totalT, tol, grid_pts = 100, eval_thres = "step", R_cutoff = 1e-3):
    ex = np.asarray(u, dtype=np.float64)   # (E,3)
    LL = ex[:,0]                               # (E,) start positions
    LU = ex[:,1]                               # (E,) end positions
    MU = ex[:,2]                               # (E,) local mu per bp
    
    positions = 0.5*(LL + LU)                  # (E,) point mass positions
    scaledu   = (LU - LL)*MU                   # (E,) point-mass weight
    
    # recombination map at exon endpoints (compute once)
    RLL = r(LL)                  # (E,) 
    RLU = r(LU)                  # (E,) 
    RB  = (RLU - RLL) / (LU - LL)              # (E,) recomb per bp within exon
    
    ss = np.asarray(ss, np.float64) 
    ps = np.asarray(ps, np.float64)
    
    # startTime = datetime.now()
    recDist = getRecDist_3(positions, r, [focalPos]) # (F,E) recDist[i, j] is distance between focal_positions[i] and positions[j]
    # endTime = datetime.now()
    # endTime - startTime
    
    # startTime = datetime.now()
    
    rf = r(focalPos)
    B_log_at_t, idx_contains, idx_nearby, ext_mask = make_B_log_evaluator_3case(
        recDist_i=recDist[0], fp=focalPos, R_cutoff=R_cutoff,
        ss=ss, ps=ps, LL=LL, LU=LU, MU=MU, RLL=RLL, RLU=RLU, RB=RB, rf=rf,
        scaledu=scaledu
    )
    
    dt = censusSize / 10
    if eval_thres == "step":
        prev_B = 1 
        j = 0
        bvals = [1.0]
        t_vals = [0.0]
        
        equil = False
        while not equil:
            j += 1
            t = j * dt
            B = float(np.exp(B_log_at_t(t)))
            bvals.append(B) 
            t_vals.append(t)

            
            if abs(B - prev_B) <= tol: 
                equil = True
            prev_B = B
    elif eval_thres == "inf":
        B_inf = get_Binf(u, ss, ps, r, [focalPos], R_cutoff = R_cutoff)
        
        t_vals, bvals = eval_multires_until_equil_stride_tol(
            B_log_at_t,
            B_inf[0],
            dt,
            tol = 1e-2)
    
    ancTime = max(t_vals[-1], totalT * 2 * censusSize)
    if t_vals[-1] != ancTime or len(bvals) < grid_pts:
        t_vals = np.linspace(0, ancTime, max(grid_pts, len(bvals)))
        bvals = np.exp([B_log_at_t(t) for t in t_vals])
        
    ancB = bvals[-1]
    ancNe = ancB * censusSize
    
    rescaled_ts = [(ancTime - x) / 2 / ancNe for x in reversed(t_vals)]
    rescaled_bs = [x / ancB for x in reversed(bvals)]       
        
    tmp_fun = interpolate.interp1d(rescaled_ts, rescaled_bs)
    def q_fun(t):
        tt = t
        if tt < 0:
            tt = 0
        elif tt > ancTime / 2 / ancNe:
            tt = ancTime / 2 / ancNe
        return [tmp_fun(tt).tolist()]
    
    return(q_fun, ancTime, ancNe)

# plt.figure(figsize=(6,4))
# plt.plot(t_vals_inf, bvals_inf, color='blue', lw=3, label = "inf")
# plt.plot(t_vals_step, bvals_step, color='red', lw=3, label = "step")
# plt.xlabel("time (generations)")
# plt.ylabel("B(t)")
# plt.tight_layout()
# plt.legend()
# plt.show()

# plt.figure(figsize=(6,4))
# # plt.plot(t_vals_inf, B_rect_inf[focal_idx], color='blue', lw=3, label = "inf")
# plt.plot(np.linspace(0, ancTime / 2 / ancNe, 100), [q_fun(x) for x in np.linspace(0, ancTime / 2 / ancNe, 100)], color='red', lw=3, label = "step")
# plt.xlabel("time (generations)")
# plt.ylabel("B(t)")
# plt.tight_layout()
# plt.legend()
# plt.show()

def get_Binf(exons, ss, ps, r_func, focal_positions,
                          dt=1e3, R_cutoff=0.01):
    ex = np.asarray(exons, dtype=np.float64) # (E,3)
    LL = ex[:, 0] # (E,) start positions
    LU = ex[:, 1] # (E,) start positions
    MU = ex[:, 2] # (E,) start positions
    
    positions = 0.5 * (LL + LU) # (E,) point mass positions
    scaledu = (LU - LL) * MU # (E,) point mass weights
    
    # Precompute endpoint r-values and within-exon rbp
    RLL = r_func(LL)
    RLU = r_func(LU)
    RB = (RLU - RLL) / (LU - LL)
    
    # Precompute r at point-mass positions (for recDist per focal)
    recDist = getRecDist_3(positions, r_func, focal_positions) # (F,E) recDist[i, j] is distance between focal_positions[i] and positions[j]
    
    ss = np.asarray(ss, np.float64)
    ps = np.asarray(ps, np.float64)
    
    F = len(focal_positions)
    B_inf = np.empty(F, dtype=np.float64)
    
    for i, fp in enumerate(focal_positions):
        fp = float(fp)
        rf = float(r_func(fp))
        recDist_i = recDist[i]
    
        # Build 3-case evaluator 
        B_log_at_t, idx_contains, idx_nearby, ext_mask = make_B_log_evaluator_3case(
            recDist_i=recDist_i, fp=fp, R_cutoff=R_cutoff,
            ss=ss, ps=ps, LL=LL, LU=LU, MU=MU, RLL=RLL, RLU=RLU, RB=RB, rf=rf,
            scaledu=scaledu
        )
    
        # ---- B_inf using t->inf limits (3 cases) ----
        log_inf = 0.0
        log_inf += logB_inf_external(recDist_i, scaledu, ss, ps, ext_mask)
        log_inf += logB_inf_contains(fp, LL, LU, MU, RB, ss, ps, idx_contains)
        log_inf += logB_inf_nearby(fp, LL, LU, MU, RB, RLL, RLU, rf, ss, ps, idx_nearby)
    
        B_inf_i = float(np.exp(log_inf))
        B_inf[i] = B_inf_i
        
    return B_inf

def get_Binf_and_halftime(exons, ss, ps, r_func, focal_positions,
                          dt=1e3, R_cutoff=0.01, max_steps=200000):
    """
    exons: array-like (E,3) rows [start, stop, local_mu_per_bp]
    Returns:
      B_inf: (F,)
      half_time: (F,) (in same units as dt), np.nan if not reached within max_steps
      half_index: (F,) index j where crossing first occurs (so t = j*dt)
    """
    ex = np.asarray(exons, dtype=np.float64) # (E,3)
    LL = ex[:, 0] # (E,) start positions
    LU = ex[:, 1] # (E,) start positions
    MU = ex[:, 2] # (E,) start positions

    positions = 0.5 * (LL + LU) # (E,) point mass positions
    scaledu = (LU - LL) * MU # (E,) point mass weights

    # Precompute endpoint r-values and within-exon rbp
    RLL = r_func(LL)
    RLU = r_func(LU)
    RB = (RLU - RLL) / (LU - LL)

    # Precompute r at point-mass positions (for recDist per focal)
    recDist = getRecDist_3(positions, r_func, focal_positions) # (F,E) recDist[i, j] is distance between focal_positions[i] and positions[j]

    ss = np.asarray(ss, np.float64)
    ps = np.asarray(ps, np.float64)

    F = len(focal_positions)
    B_inf = np.empty(F, dtype=np.float64)
    half_time = np.full(F, np.nan, dtype=np.float64)
    half_index = np.full(F, -1, dtype=np.int64)

    for i, fp in enumerate(focal_positions):
        fp = float(fp)
        rf = float(r_func(fp))
        recDist_i = recDist[i]

        # Build 3-case evaluator 
        B_log_at_t, idx_contains, idx_nearby, ext_mask = make_B_log_evaluator_3case(
            recDist_i=recDist_i, fp=fp, R_cutoff=R_cutoff,
            ss=ss, ps=ps, LL=LL, LU=LU, MU=MU, RLL=RLL, RLU=RLU, RB=RB, rf=rf,
            scaledu=scaledu
        )

        # ---- B_inf using t->inf limits (3 cases) ----
        log_inf = 0.0
        log_inf += logB_inf_external(recDist_i, scaledu, ss, ps, ext_mask)
        log_inf += logB_inf_contains(fp, LL, LU, MU, RB, ss, ps, idx_contains)
        log_inf += logB_inf_nearby(fp, LL, LU, MU, RB, RLL, RLU, rf, ss, ps, idx_nearby)

        B_inf_i = float(np.exp(log_inf))
        B_inf[i] = B_inf_i

        # ---- half-time on uniform grid using full 3-case logB(t) ----
        thresh = 0.5 * (1.0 + B_inf_i)

        # If threshold is >= 1 (shouldn't), half-time is 0.
        if thresh >= 1.0:
            half_time[i] = 0.0
            half_index[i] = 0
            continue

        # evaluate B(j*dt) sequentially until first crossing
        for j in range(1, int(max_steps) + 1):
            t = j * float(dt)
            Bval = float(np.exp(B_log_at_t(t)))
            if Bval <= thresh:
                half_index[i] = j
                half_time[i] = t
                break

    return B_inf, half_time, half_index

def get_Bs_2(u, ss, ps, r_func, focal_positions, tol, dt = 1e3, R_cutoff = 0.01):
    # first define some variables that are reused often
    ex = np.asarray(u, dtype=np.float64)   # (E,3)
    LL = ex[:,0]                               # (E,) start positions
    LU = ex[:,1]                               # (E,) end positions
    MU = ex[:,2]                               # (E,) local mu per bp
    
    positions = 0.5*(LL + LU)                  # (E,) point mass positions
    scaledu   = (LU - LL)*MU                   # (E,) point-mass weight
    
    # recombination map at exon endpoints (compute once)
    RLL = r_func(LL)                  # (E,) 
    RLU = r_func(LU)                  # (E,) 
    RB  = (RLU - RLL) / (LU - LL)              # (E,) recomb per bp within exon
    
    ss = np.asarray(ss, np.float64) 
    ps = np.asarray(ps, np.float64)
    
    # startTime = datetime.now()
    recDist = getRecDist_3(positions, r_func, focal_positions) # (F,E) recDist[i, j] is distance between focal_positions[i] and positions[j]
    # endTime = datetime.now()
    # endTime - startTime
    
    # startTime = datetime.now()
    
    B_inf = get_Binf(u, ss, ps, r_func, focal_positions, dt, R_cutoff)
    
    all_bs = []
    all_ts = []
    
    for i,fp in enumerate(focal_positions):
        rf = r_func(fp)
        # dL = np.abs(fp - LL)   # (E,) distance to left endpoint for all exons
        # dU = np.abs(fp - LU)   # (E,) distance to right endpoint
        
        # lb_all = np.minimum(dL, dU)   # (E,) nearer boundary distance (bp)
        # ub_all = np.maximum(dL, dU)   # (E,) farther boundary distance (bp)
        
        B_log_at_t, idx_contains, idx_nearby, ext_mask = make_B_log_evaluator_3case(
            recDist_i=recDist[i], fp=fp, R_cutoff=R_cutoff,
            ss=ss, ps=ps, LL=LL, LU=LU, MU=MU, RLL=RLL, RLU=RLU, RB=RB, rf=rf,
            scaledu=scaledu
        )
        
        # B_log_at_t = make_B_log_evaluator_2case(recDist_i = recDist[i], 
        #                                         fp = fp,
        #                                         ss = ss,
        #                                         ps = ps,
        #                                         LL = LL,
        #                                         LU = LU,
        #                                         MU = MU,
        #                                         RB = RB,
        #                                         scaledu = scaledu)
        
        t_vals, B_vals = eval_multires_until_equil_stride_tol(
            B_log_at_t,
            B_inf[i],
            dt,
            tol)
        
        all_bs.append(B_vals)
        all_ts.append(t_vals)
    # endTime = datetime.now()
    # endTime - startTime
    return all_bs, all_ts
 
def eval_multires_until_equil_stride_tol(
    B_log_at_t,
    b_inf,
    dt,
    tol,
    phases=((1, 100), (100, 10100), (1e4, 1010100)),
    max_evals=1_000_000,
):
    """
    Returns:
        t_vals : array of sampled times
        B_vals : array of B(t) values

    t_vals are in the same units as dt (since t = j*dt)

    Stop when approximate per-dt change is below tol.
    """

    t_vals = [0.0]
    B_vals = [1.0]

    j = 0
    evals = 0

    for stride, j_max in phases:
        stride = int(stride)
        j_max = int(j_max)

        while j < j_max:
            j += stride
            t = j * float(dt)

            B = float(np.exp(B_log_at_t(t)))

            t_vals.append(t)
            B_vals.append(B)
            evals += 1

            if evals >= max_evals:
                return np.asarray(t_vals), np.asarray(B_vals)

            # stride-aware stopping condition
            if abs(B - b_inf) <= tol:
                return np.asarray(t_vals), np.asarray(B_vals)


    raise ValueError('didnt reach equil in phases provided')
    return np.asarray(t_vals), np.asarray(B_vals)
 
def get_Bs(u, ss, ps, r_func, focal_positions, tol, dt = 1e3, R_cutoff = 0.01):
    # first define some variables that are reused often
    ex = np.asarray(u, dtype=np.float64)   # (E,3)
    LL = ex[:,0]                               # (E,) start positions
    LU = ex[:,1]                               # (E,) end positions
    MU = ex[:,2]                               # (E,) local mu per bp
    
    positions = 0.5*(LL + LU)                  # (E,) point mass positions
    scaledu   = (LU - LL)*MU                   # (E,) point-mass weight
    
    # recombination map at exon endpoints (compute once)
    RLL = r_func(LL)                  # (E,) 
    RLU = r_func(LU)                  # (E,) 
    RB  = (RLU - RLL) / (LU - LL)              # (E,) recomb per bp within exon
    
    ss = np.asarray(ss, np.float64) 
    ps = np.asarray(ps, np.float64)
    
    # startTime = datetime.now()
    recDist = getRecDist_3(positions, r_func, focal_positions) # (F,E) recDist[i, j] is distance between focal_positions[i] and positions[j]
    # endTime = datetime.now()
    # endTime - startTime
    
    # startTime = datetime.now()
    
    all_bs = []
    
    for i,fp in enumerate(focal_positions):
        rf = r_func(fp)
        # dL = np.abs(fp - LL)   # (E,) distance to left endpoint for all exons
        # dU = np.abs(fp - LU)   # (E,) distance to right endpoint
        
        # lb_all = np.minimum(dL, dU)   # (E,) nearer boundary distance (bp)
        # ub_all = np.maximum(dL, dU)   # (E,) farther boundary distance (bp)
        
        B_log_at_t, idx_contains, idx_nearby, ext_mask = make_B_log_evaluator_3case(
            recDist_i=recDist[i], fp=fp, R_cutoff=R_cutoff,
            ss=ss, ps=ps, LL=LL, LU=LU, MU=MU, RLL=RLL, RLU=RLU, RB=RB, rf=rf,
            scaledu=scaledu
        )
        
        # B_log_at_t = make_B_log_evaluator_2case(recDist_i = recDist[i], 
        #                                         fp = fp,
        #                                         ss = ss,
        #                                         ps = ps,
        #                                         LL = LL,
        #                                         LU = LU,
        #                                         MU = MU,
        #                                         RB = RB,
        #                                         scaledu = scaledu)
        
        prev_B = 1 
        j = 0
        bvals = [1.0]
        
        equil = False
        while not equil:
            j += 1
            t = j * dt
            B = float(np.exp(B_log_at_t(t)))
            bvals.append(B) 
            
            if abs(B - prev_B) <= tol: 
                equil = True
            prev_B = B
            
        all_bs.append(bvals)
    # endTime = datetime.now()
    # endTime - startTime
    return all_bs
                                                
# t_equil = find_equil_time_by_doubling(B_log_at_t, t0=dt, tol=tol)
# bvals = eval_on_uniform_grid_with_padding(B_log_at_t, dt=dt, t_equil=t_equil)


        # B_log_at_t = make_B_log_evaluator_3case(recDist_i = recDist[i],  s
        #                                         fp = fp,
        #                                         R_cutoff = 0.01,
        #                                         ss = ss, 
        #                                         ps = ps, 
        #                                         LL = LL,
        #                                         LU = LU,
        #                                         MU = MU,
        #                                         RLL = RLL, 
        #                                         RLU = RLU,
        #                                         RB = RB,
        #                                         rf = rf,
        #                                         scaledu = scaledu,
        #                                         lb_all = lb_all,
        #                                         ub_all = ub_all)
                                              
def make_B_log_evaluator_2case(recDist_i, fp, ss, ps, LL, LU, MU, RB, scaledu):
    # contains mask
    contains = (LL <= fp) & (fp <= LU) # TODO implement searchsorted???

    # external is everything except the containing exon (or exons)
    external = ~contains # TODO copy and flip one element

    # external point-mass evaluator
    ext_logB = make_external_logB_evaluator(recDist_i, scaledu, ss, ps, external)

    idx_contains = np.nonzero(contains)[0]  # usually 0 or 1

    def logB(t):
        total = ext_logB(t)
        if len(idx_contains) > 1:
            raise ValueError("Multiple exons contain the focal site.")
        if len(idx_contains) == 1:
            j = idx_contains[0]
            total += exon_contains_focal_log(ss, ps, MU[j], LL[j], LU[j], fp, RB[j], t)
        return total

    return logB

def make_B_log_evaluator_3case(recDist_i, fp, R_cutoff,
                               ss, ps, LL, LU, MU, RLL, RLU, RB, rf, scaledu):
    """
    Builds a per-focal logB(t) evaluator for half-time search.
    """
    # creating masks
    #TODO deal with contains length >1 edge-case
    contains = (LL <= fp) & (fp <= LU)
    nearby = (~contains) & (recDist_i < R_cutoff)
    external = (~contains) & (~nearby)

    ext_logB = make_external_logB_evaluator(recDist_i, scaledu, ss, ps, external)

    idx_contains = np.nonzero(contains)[0]
    idx_nearby = np.nonzero(nearby)[0]

    near_logB = None
    if idx_nearby.size > 0:
        if R_cutoff == 0:
            raise ValueError("Should be no nearby with R_cutoff == 0.")
        LLn = LL[idx_nearby]; LUn = LU[idx_nearby]
        dL = np.abs(fp - LLn)
        dU = np.abs(fp - LUn)
        lb = np.minimum(dL, dU)
        ub = np.maximum(dL, dU)

        r0 = np.minimum(np.abs(RLL[idx_nearby] - rf), np.abs(RLU[idx_nearby] - rf))
        near_logB = make_nearby_logB_evaluator(MU[idx_nearby], lb, ub, RB[idx_nearby], r0, ss, ps)

    def logB(t):
        # Case A: far away
        total = ext_logB(t)

        # Case B: containing exon
        if idx_contains.size > 1:
            raise ValueError("Multiple exons contain focal site; handle overlaps explicitly.")
        if idx_contains.size == 1:
            j = int(idx_contains[0])
            total += exon_contains_focal_log(ss, ps, MU[j], LL[j], LU[j], fp, RB[j], t)

        # Case C: nearby
        if near_logB is not None:
            total += near_logB(t)

        return total

    return logB, idx_contains, idx_nearby, external

def make_B_log_evaluator_3case_dep(
    recDist_i, fp, R_cutoff, ss, ps, LL, LU, MU, RLL, RLU, RB, rf, scaledu, lb_all, ub_all
):
    # creating masks
    #TODO deal with contains length >1 edge-case
    contains = (LL <= fp) & (fp <= LU) # (E,)
    # nearby = (~contains) & (recDist_i <= R_cutoff)
    external = (~contains) # & (~nearby)
    
    # distance from fp to NEARER exon endpoint
    # r0 = np.minimum(np.abs(RLL - rf), np.abs(RLU - rf)) # (E,)
    
    # preparing lambda function for far away exons
    ext_logB = make_external_logB_evaluator(recDist_i, scaledu, ss, ps, external)

    # indices to evaluate over:
    idx_contains = np.nonzero(contains)[0]      # length 0 or 1
    # idx_nearby   = np.nonzero(nearby)[0]        # typically small

    def logB(t):
        # Case A: far away
        total = ext_logB(t)

        # Case B: containing exon
        if len(idx_contains) > 1:
            raise ValueError('Woah there, partner!')
        for j in idx_contains:
            total += exon_contains_focal_log(ss, ps, MU[j], LL[j], LU[j], fp, RB[j], t)

        # # Case C: nearby non-containing exons
        # if len(idx_nearby) > 0:
        #     lb_near  = lb_all[idx_nearby]
        #     ub_near  = ub_all[idx_nearby]
        #     rbp_near = RB[idx_nearby]
        #     r0_near  = r0[idx_nearby]
        #     mu_near  = MU[idx_nearby]
            
        #     total += exon_nearby_log_with_r0_batch(
        #                  ss, ps,
        #                  mu_near, lb_near, ub_near,
        #                  rbp_near, r0_near, t
        #              )


        return total

    return logB   

def make_nearby_logB_evaluator(mu_near, lb_near, ub_near, rbp_near, r0_near, ss, ps):
    """
    Time-dependent nearby evaluator with r0, batch over nearby exons and s-bins.
    Precomputes t-independent arrays; each call still does exp/exp1.
    """
    # Convert to arrays
    s = ss[None, :]        # (1,K)
    p = ps[None, :]        # (1,K)

    lb  = lb_near[:, None]  # (N,1)
    ub  = ub_near[:, None]
    rbp = rbp_near[:, None]
    r0  = r0_near[:, None]
    mu  = mu_near[:, None]

    # precompute common terms:
    lb_r = lb * rbp   # (N,1)
    ub_r = ub * rbp
    mu_over_r = mu / rbp

    a0 = lb_r + r0 + s # (N,K)
    b0 = r0 + s + ub_r
    denom = a0 * b0

    # reusable temp arrays to avoid reallocations each time
    tmp = np.empty_like(a0) # (N,K)
    e_tb  = np.empty_like(a0)
    e2_tb = np.empty_like(a0)
    e_ta  = np.empty_like(a0)
    e2_ta = np.empty_like(a0)

    E1_at  = np.empty_like(a0)
    E1_2at = np.empty_like(a0)
    E1_bt  = np.empty_like(a0)
    E1_2bt = np.empty_like(a0)

    numer = np.empty_like(a0)
    frac = np.empty_like(a0)
    inside = np.empty_like(a0)

    def near_logB(t):
        t = float(t)
        if t <= 0.0:
            return 0.0

        # Exponentials (N,K)
        np.multiply(b0, -t, out=tmp);      np.exp(tmp, out=e_tb) # exp(-t * b)
        np.multiply(b0, -2.0*t, out=tmp);  np.exp(tmp, out=e2_tb) # exp(-2 * t * b)
        np.multiply(a0, -t, out=tmp);      np.exp(tmp, out=e_ta) # exp(-t * a)
        np.multiply(a0, -2.0*t, out=tmp);  np.exp(tmp, out=e2_ta) # exp(-2 * t * a)

        # Gamma terms (N,K)
        E1_at[:]  = exp1(a0 * t)
        E1_2at[:] = exp1(2.0 * a0 * t)
        E1_bt[:]  = exp1(b0 * t)
        E1_2bt[:] = exp1(2.0 * b0 * t)

        # numer = -(lb*r) - a0*e2tb + 2*a0*etb + (ub*r) + b0*e2ta - 2*b0*eta
        numer[:] = -lb_r
        numer[:] += -a0 * e2_tb
        numer[:] +=  2.0 * a0 * e_tb
        numer[:] +=  ub_r
        numer[:] +=  b0 * e2_ta
        numer[:] += -2.0 * b0 * e_ta

        np.divide(numer, denom, out=frac, where=(denom != 0.0))

        inside[:] = frac
        inside[:] +=  2.0 * t * E1_at
        inside[:] += -2.0 * t * E1_2at
        inside[:] +=  2.0 * t * (-E1_bt + E1_2bt)

        # contrib = -(s * mu * inside) / rbp, with p weights
        return float(np.sum((-s) * mu_over_r * inside * p))

    return near_logB

# TODO rename
def exon_nearby_log_with_r0_batch(ss, ps,
                                   mu_near, lb_near, ub_near,
                                   rbp_near, r0_near, t):
    """
    Log-contribution for a NEARBY exon that does NOT contain the focal site
    """

    t = float(t)
    if t <= 0:
        return 0.0

    # Convert to arrays
    s = ss[None, :]        # (1,K)
    p = ps[None, :]        # (1,K)

    lb  = np.asarray(lb_near,  dtype=np.float64)[:, None]  # (N,1) # TODO check if as array is needed when we deal with nearbys
    ub  = np.asarray(ub_near,  dtype=np.float64)[:, None]
    rbp = np.asarray(rbp_near, dtype=np.float64)[:, None]
    r0  = np.asarray(r0_near,  dtype=np.float64)[:, None]
    mu  = np.asarray(mu_near,  dtype=np.float64)[:, None]

    # precompute common terms:
    a = lb * rbp + r0 + s        # (N,K)
    b = r0 + s + ub * rbp        # (N,K)

    # Exponentials
    e2tb = np.exp(-2.0 * t * b)  # (N,K)
    etb  = np.exp(-t * b)
    e2ta = np.exp(-2.0 * t * a)
    eta  = np.exp(-t * a)

    #  numerator
    # TODO some of this can be replaced with a and b
    numer = (-(lb * rbp)   # (N,K)
             - (lb * rbp + r0 + s) * e2tb
             + 2.0 * (lb * rbp + r0 + s) * etb
             + rbp * ub
             + (r0 + s + rbp * ub) * e2ta
             - 2.0 * (r0 + s + rbp * ub) * eta)

    denom = a * b   # (N,K)
    frac = np.divide(numer, denom,
                     out=np.zeros_like(numer),
                     where=(denom != 0.0))

    # Gamma terms
    gamma_part = (2.0 * t * exp1(a * t)
                  - 2.0 * t * exp1(2.0 * a * t)
                  + 2.0 * t * (-exp1(b * t) + exp1(2.0 * b * t)))

    inside = frac + gamma_part

    contrib = -(s * mu * inside) / rbp      # (N,K)

    # Sum over DFE bins and exons
    return float(np.sum(p * contrib))

def get_pred_loci(n_pos, r):
    loci = np.linspace(r[0][0], r[-1][-2], n_pos + 2)
    return loci[1:-1]
    
# ----------------------------
# t->infty limits (log-space)
# ----------------------------

def logB_inf_external(recDist_i, scaledu, ss, ps, ext_mask):
    """
    external limit:  -(s*u)/(s + r)^2  summed over exons and DFE bins.
    Here r is recDist_i for point-mass positions.
    """
    U = scaledu[ext_mask][:, None]  # (Mext,1)
    R = recDist_i[ext_mask][:, None]# (Mext,1)
    s = ss[None, :]                 # (1,K)
    p = ps[None, :]                 # (1,K)

    D = R + s                                               # (Mext,K)
    invD = np.divide(1.0, D, out=np.zeros_like(D), where=(D != 0.0))
    term = (-U) * s * (invD * invD)                          # -(u*s)/(r+s)^2
    return float(np.sum(term * p))

def logB_inf_contains(fp, LL, LU, MU, RB, ss, ps, idx_contains):
    """
    containing limit (evaluated twice, left and right):
        - (l * u) / (l*r + s)
    where:
      l = distance from focal to endpoint (bp)
      u = MU (local mu per bp for the exon)
      r = RB (recomb per bp within the exon)
    """
    if idx_contains.size == 0:
        return 0.0
    if idx_contains.size > 1:
        raise ValueError("Multiple exons contain focal site; handle overlaps explicitly.")

    j = int(idx_contains[0])
    ll, lu = float(LL[j]), float(LU[j])
    u = float(MU[j])
    r = float(RB[j])

    lL = fp - ll
    lR = lu - fp

    s = np.asarray(ss, np.float64)
    p = np.asarray(ps, np.float64)

    # left:  -(lL*u)/(lL*r + s), right: -(lR*u)/(lR*r + s)
    denomL = lL * r + s
    denomR = lR * r + s
    termL = - (lL * u) * np.divide(1.0, denomL, out=np.zeros_like(s), where=(denomL != 0.0))
    termR = - (lR * u) * np.divide(1.0, denomR, out=np.zeros_like(s), where=(denomR != 0.0))
    return float(np.sum(p * (termL + termR)))

def logB_inf_nearby(fp, LL, LU, MU, RB, RLL, RLU, rf, ss, ps, idx_nearby):
    """
    nearby limit:
        (s*u*(lb - ub)) / ((lb*r + r0 + s) * (r0 + s + r*ub))

    where:
      lb, ub are bp distances from focal to nearer/farther exon boundary
      r  = rbp within exon
      r0 = recomb distance from focal to nearer boundary (in your r-space)
      u  = MU (local mu per bp)
    Note lb-ub < 0, so term is negative.
    """
    if idx_nearby.size == 0:
        return 0.0

    s = ss[None, :]   # (1,K)
    p = ps[None, :]   # (1,K)

    j = idx_nearby
    LLn = LL[j]; LUn = LU[j]
    mu  = MU[j][:, None]     # (N,1)
    rbp = RB[j][:, None]     # (N,1)

    dL = np.abs(fp - LLn)
    dU = np.abs(fp - LUn)
    lb = np.minimum(dL, dU)[:, None]          # (N,1)
    ub = np.maximum(dL, dU)[:, None]          # (N,1)

    # r0 from precomputed endpoint r-values
    r0 = np.minimum(np.abs(RLL[j] - rf), np.abs(RLU[j] - rf))[:, None]  # (N,1)

    a = lb * rbp + r0 + s                      # (N,K)
    b = r0 + s + ub * rbp                      # (N,K)

    denom = a * b
    numer = s * mu * (lb - ub)                 # (N,K) via broadcast (lb-ub) (N,1)

    term = np.divide(numer, denom, out=np.zeros_like(denom), where=(denom != 0.0))
    return float(np.sum(term * p))







os.chdir("/home/nathan/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa/threepop/weakSelection")
recName = "YRI_recombination_map_hg38_chr_22.bed"
mutName = "roulette_tbl_chr22.csv"
exonName = "exons_chr22.bed"
recMap = read_rec_map(recName)
mutMap = read_mut_rates(mutName)
exonMap = read_exon_map(exonName)

recMap = simplify_rate_map(recMap)
mutMap = simplify_rate_map(mutMap)

exonMutMap, U = make_exon_only_mutmap(mutMap, exonMap) 

focal_positions = [(x+y)/2 for x,y in exonMap][1::4]
# focal_positions = get_pred_loci(int(1e3), recMap)

# 875 / 2 / 12378 * shape to get the 0.006 number
mean = 0.00657416 / 2
shape = 0.186
nbins = 10
dfe = discretize_deleterious_gamma_dfe_mean_shape(mean, shape, nbins)
ss = [y for x,y in dfe]

startTime = datetime.now()
B_rect, labels, medoids, t_vals = cluster_scale_funs(u = exonMutMap,
                    r = recMap,
                    focal_positions = focal_positions,
                    ss = ss,
                    tol = 1e-5,
                    dt = 730,
                    R_cutoff = 1e-3,
                    n_clusters = 20,
                    equil_thres = "step")
endTime = datetime.now()
endTime - startTime

B_rect_step = B_rect
t_vals_step = np.array(t_vals, dtype = 'float64')

startTime = datetime.now()
B_rect, labels, medoids, t_vals = cluster_scale_funs(u = exonMutMap,
                    r = recMap,
                    focal_positions = focal_positions,
                    ss = ss,
                    tol = 1e-2,
                    dt = 730,
                    R_cutoff = 1e-3,
                    n_clusters = 20,
                    equil_thres = "inf")
endTime = datetime.now()
endTime - startTime

        

B_rect_inf = B_rect
t_vals_inf = t_vals

focal_idx = 516
focalPos = focal_positions[focal_idx]

plt.figure(figsize=(6,4))
plt.plot(t_vals_inf, B_rect_inf[focal_idx], color='blue', lw=3, label = "inf")
plt.plot(t_vals_step, B_rect_step[focal_idx], color='red', lw=3, label = "step")
plt.xlabel("time (generations)")
plt.ylabel("B(t)")
plt.tight_layout()
plt.legend()
plt.show()

plt.figure(figsize=(6,4))
plt.plot(t_vals_inf[np.nonzero(t_vals_inf <= max(t_vals_step))], B_rect_inf[focal_idx][np.nonzero(t_vals_inf <= max(t_vals_step))], color='blue', lw=3, label = "inf")
plt.plot(t_vals_step, B_rect_step[focal_idx], color='red', lw=3, label = "step")
plt.xlabel("time (generations)")
plt.ylabel("B(t)")
plt.tight_layout()
plt.legend()
plt.show()

t_max = None
B_prev = B_rect_step[focal_idx][0]
for i in range(1, len(t_vals_step)):
    if B_prev == B_rect_step[focal_idx][i]:
        t_max = t_vals_step[i]
        break
    else:
        B_prev = B_rect_step[focal_idx][i]

plt.figure(figsize=(6,4))
plt.plot(t_vals_inf[np.nonzero(t_vals_inf <= t_max)], B_rect_inf[focal_idx][np.nonzero(t_vals_inf <= t_max)], color='blue', lw=3, label = "inf")
plt.plot(t_vals_step[np.nonzero(t_vals_step <= t_max)], B_rect_step[focal_idx][np.nonzero(t_vals_step <= t_max)], color='red', lw=3, label = "step")
plt.xlabel("time (generations)")
plt.ylabel("B(t)")
plt.tight_layout()
plt.legend()
plt.show()

sample_size = [40]
sampled_demes = ['CEU']
demo = demes.load('ooa.yaml')

def get_ceu_data(params, focalIndex, n_load, key="ceu", dtype=np.float64):
    seeds = [y for x, y in params][0:n_load]
    n = len(seeds)
    if n == 0:
        raise ValueError("No seeds")

    # init accumulator
    with np.load(f"{seeds[0]}.npz") as z:
        acc = np.array(z[key][focalIndex], dtype=dtype, copy=True)

    # in-place add
    for seed in seeds[1:]:
        with np.load(f"{seed}.npz") as z:
            acc += z[key][focalIndex]

    acc /= n
    return acc

def read_params():
    intervals = []
    with open('ooaGamma.txt') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split()
            demo   = str(fields[0])
            seed  = int(fields[1])
            intervals.append([demo, seed])
    return intervals

params = read_params()
focalIndex = [i for i,x in enumerate(exonMap) if x[0] < focalPos and x[1] > focalPos][0]
os.chdir('/media/nathan/T7/BGSdemo/ooaGammaExonWindows')
data = get_ceu_data(params, focalIndex,500, dtype=np.float32)
np.savez_compressed('site_' + str(focalPos), data=data)


fs, ancNe = bgs_wrapper(u = exonMutMap,
                        r = recMap,
                        focalPos = focalPos,
                        sample_size = sample_size,
                        ss = ss,
                        sampled_demes=sampled_demes,
                        g = demo,
                        tol = 1e-5,
                        eval_thres = "step")

fs_inf, ancNe_inf = bgs_wrapper(u = exonMutMap,
                        r = recMap,
                        focalPos = focalPos,
                        sample_size = sample_size,
                        ss = ss,
                        sampled_demes=sampled_demes,
                        g = demo,
                        tol = 1e-2,
                        eval_thres = "inf")

startTime = datetime.now()
fs, ancNe = bgs_wrapper(u = exonMutMap,
                        r = recMap,
                        focalPos = focalPos,
                        sample_size = sample_size,
                        ss = ss,
                        cs = [1e4],
                        tol = 1e-5,
                        eval_thres = "step",
                        totalT = 0)
endTime = datetime.now()
endTime - startTime


        
fs_neu = moments.Spectrum.from_demes(
    demo, 
    sampled_demes=sampled_demes, 
    sample_sizes=sample_size
)

oldestEpoch, ancCensusSize = getOldestEpoch(demo)
fs = fs * 4  * 1e-8 * ancNe
fs_inf = fs_inf * 4 * 1e-8 * ancNe_inf
fs_neu = fs_neu * 4  * 1e-8 * ancCensusSize    


fig, ax = plt.subplots(1, 1, figsize=(16, 8), sharex=True, sharey=False)
fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
fig.subplots_adjust(hspace = .25)
fig.suptitle(str(focalPos))
ax.plot(fs, "-", ms=8, lw=1, label="BGS step")
# ax.plot(fs_inf, "-", ms = 8, lw = 1, label = "BGS inf")
ax.plot(fs_neu, "-", ms=8, lw=1, label="SNM")
ax.set_title("nbins = " + str(nbins)) 
ax.set_yscale('log')
ax.annotate("B=" + str(round(fs.pi()/fs_neu.pi(),3)),xy=(sample_size[0] / 10, max(fs[1:-1]) * 0.8),size="x-large")
ax.annotate("B=" + str(round(fs_inf.pi()/fs_neu.pi(),3)),xy=(sample_size[0] / 10, max(fs[1:-1]) * 0.9),size="x-large")
ax.legend();

plot_kmedoids_clusters(B_rect, labels, medoids)



                                                         
# focalPos = 49885518.61638361
# demo = demes.load('ooa.yaml')
# sampled_demes = ['CEU']
# u = exonMutMap
# r = recMap
# focal_positions = focal_positions
# ss = ss
# sample_size = 40

# [[x,y] for x,y in exonMap if x < focalPos and y > focalPos]
# u = [[x,y,z] for x,y,z in u if x <= 49887179 and y >= 49883661]

# fs, ancNe = bgs_wrapper(u = u,
#                         r = recMap,
#                         focalPos = focalPos,
#                         sample_size = [sample_size],
#                         ss = ss,
#                         sampled_demes=sampled_demes,
#                         g = demo)

# fs_neu = moments.Spectrum.from_demes(
#     demo, 
#     sampled_demes=sampled_demes, 
#     sample_sizes=[sample_size]
# )

# oldestEpoch, ancCensusSize = getOldestEpoch(demo)
# fs = fs * 8 * 1e-8 * ancNe
# fs_neu = fs_neu * 8 * 1e-8 * ancCensusSize    


# fig, ax = plt.subplots(1, 1, figsize=(16, 8), sharex=True, sharey=False)
# fig.text(0.5, 0.04, 'Allele Frequency', ha='center')
# fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
# fig.subplots_adjust(hspace = .25)
# fig.suptitle(str(focalPos))
# ax.plot(fs, "-", ms=8, lw=1, label="BGS")
# ax.plot(fs_neu, "-", ms=8, lw=1, label="SNM")
# ax.set_title("nbins = " + str(nbins)) 
# ax.set_yscale('log')
# ax.annotate("B=" + str(round(fs.pi()/fs_neu.pi(),3)),xy=(sample_size / 10, max(fs[1:-1]) * 0.8),size="x-large")
# ax.legend();

#############
# OOA Gamma #
#############
# def get_ceu_data(params):
#     seeds = [y for x,y in params]
#     numSum = 0
#     data = None
#     for seed in seeds:
#         loaded = np.load(str(seed) + '.npz')
#         loaded = loaded['ceu']
#         loaded = loaded
        
#         if data is None:
#             data = loaded
#         else:
#             data = np.add(data, loaded)
        
#         numSum += 1
    
#     data /= numSum
#     return data

# def read_params():
#     intervals = []
#     with open('ooaGamma.txt') as f:
#         for line in f:
#             if line.startswith("#") or line.strip() == "":
#                 continue
#             fields = line.strip().split()
#             demo   = str(fields[0])
#             seed  = int(fields[1])
#             intervals.append([demo, seed])
#     return intervals

# def get_windows(recMap, mutMap):
#     chrom_start = recMap[0][0]
#     chrom_end = get_max_positions(recMap, mutMap)

    
#     pos = [chrom_start + 5e5]
#     done = False
#     while not done:
#         new_pos = pos[-1] + 1e6
#         if new_pos <= chrom_end:
#             pos.append(new_pos)
#         else:
#             done = True
            
#     windows = [[p - 1, p + 1] for p in pos]
#     windows = [x for w in windows for x in w]
#     windows.insert(0, 0)
#     windows.append(chrom_end)
    
#     return windows

# def get_focal_pos(windows):
#     tmp = []
#     for i in range(len(windows) - 1):
#         if windows[i+1] - windows[i] == 2:
#             tmp.append([windows[i], windows[i+1]])
    
#     return [(x+y)/2 for x,y in tmp]

# def get_max_positions(recMap, mutMap):
#     return max(recMap[-1][1], mutMap[-1][1])

# 875 / 2 / 12378 * shape to get the 0.006 number
# mean = 0.00657416 / 2
# shape = 0.186
# n_bins = 10
# representative = 'conditional_mean'


# windows = get_windows(recMap, mutMap)
# focal_positions = get_focal_pos(windows)
# # import random
# # focalPos = random.sample(focal_positions, 1)[0]
# chosen = [pos for i, pos in enumerate(focal_positions) if (i + 1) % 10 == 0 ]
# sample_size = 40
# sampled_demes = ['CEU']
# demo = demes.load('ooa.yaml')
# os.chdir('/media/nathan/T7/BGSdemo/gammaOOAData')

# focalPos = chosen[0]
# focalIndex = [i for i,x in enumerate(focal_positions) if x == focalPos][0]
# # 875 / 2 / 12378 * shape to get the 0.006 number
# mean = 0.00657416 / 2
# shape = 0.186

# nbins = 10

# dfe = discretize_deleterious_gamma_dfe_mean_shape(mean, shape, nbins)

# ss = [y for x,y in dfe]

# fs, ancNe = bgs_wrapper(u = exonMutMap,
#                         r = recMap,
#                         focalPos = focalPos,
#                         sample_size = [sample_size],
#                         ss = ss,
#                         sampled_demes=sampled_demes,
#                         g = demo)


                           # all_bs = []
                           # dt = 1e3
                           
                           # for i,fp in enumerate(len(focal_positions)):
                               
                               
                               
                           #     B_log_at_t = make_B_evaluator_for_focal(recDist[i], S_row, U_col, P_row) 
                           
                           #     prev_logB = 0.0
                         
                           #     j = 0
                           #     bvals = [1.0]
                               
                           #     found = False
                           
                           #     while True:
                           #         j += 1
                           #         t = j * dt
                           #         logB = B_log_at_t(t)
                           #         bvals.append(float(np.exp(logB)))
                                   
                           #         if float(np.exp(logB)) < 1e-5:
                           #             found = True
                           #             break
                           
                           #         if abs(np.exp(logB)-np.exp(prev_logB)) <= tol:   
                           #             break
                           #         prev_logB = logB
                           #     if found:
                           #         break
                           #     all_bs.append(bvals)
                           # endTime = datetime.now()
                           # endTime - startTime