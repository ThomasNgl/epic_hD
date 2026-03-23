import numpy as np
from itertools import product
from . import c_informME
from . import utils_informME

def get_n_CpGs(df_informME):
    return len(df_informME)

def get_ME_pattern_prob(ME_pattern, 
                        margprob_x1,
                        transprob_ME_arr, 
                        trans_indices = None,
                        transprob_unME_arr = None):
    # If there is only one CpG return the marg prob of the CpG
    init_prob = margprob_x1 if ME_pattern[0] == 1 else 1 - margprob_x1
    if len(ME_pattern) == 1:
        return init_prob

    if transprob_unME_arr is None:
        transprob_unME_arr = 1 - transprob_ME_arr        # P(next=0|current)
    if trans_indices is None:
        n_Cpgs = len(ME_pattern)
        trans_indices = np.arange(n_Cpgs-1)
    # Prepare bits
    n_bits = ME_pattern[:-1]
    nP1_bits = ME_pattern[1:]
    # Get the corresponding transition probs
    ME_pattern_transprob = np.where(
        nP1_bits == 1,
        transprob_ME_arr[trans_indices, n_bits],         # P(next=1|current)
        transprob_unME_arr[trans_indices, n_bits]        # P(next=0|current)
    )
    # Get the final ME pattern prob
    ME_pattern_prob = init_prob * np.prod(ME_pattern_transprob)

    return ME_pattern_prob

def sample_ME_pattern(margprob_x1, transprob_ME_arr):
    n_CpGs = len(transprob_ME_arr) + 1
    ME_pattern_sample = np.zeros(n_CpGs, dtype=np.uint8)
    ME_pattern_sample[0] = np.random.rand() <= margprob_x1

    for n in range(n_CpGs - 1):
        if ME_pattern_sample[n] == 1:
            # Draw next bit: P(next=1 | current=1)
            p = transprob_ME_arr[n, 1]
        else:
            # Draw next bit: P(next=1 | current=0)
            p = transprob_ME_arr[n, 0]
        ME_pattern_sample[n + 1] = np.random.rand() <= p

    return ME_pattern_sample

####

def get_ME_level_distribution(margprob_arr, transprob_ME_arr, max_CpGs=18):
    margprob_x1 = margprob_arr[0]
    transprob_unME_arr = 1 - transprob_ME_arr  # P(next=0|current)

    n_Cpgs = get_n_CpGs(df_informME=margprob_arr)
    trans_indices = np.arange(n_Cpgs-1)

    ml_distribution = np.zeros((n_Cpgs + 1))
    ME_level_list = np.arange(n_Cpgs + 1) / n_Cpgs
    if n_Cpgs <= max_CpGs:  # exact
        ME_pattern_list = product([0, 1], repeat=n_Cpgs)
        for ME_pattern in ME_pattern_list:
            ME_amplitude = np.sum(ME_pattern)
            ME_pattern_prob = get_ME_pattern_prob(
                ME_pattern=ME_pattern,
                transprob_ME_arr=transprob_ME_arr,
                transprob_unME_arr=transprob_unME_arr,
                trans_indices=trans_indices,
                margprob_x1=margprob_x1,
            )
            ml_distribution[ME_amplitude] += ME_pattern_prob

    else:  # Monte Carlo + maxent
        n_moments = 4
        moments = np.zeros(n_moments)
        powers = np.arange(1, n_moments+1)
        n_monte_carlo_samples = 2**(max_CpGs-1)

        samples = []
        for _ in range(n_monte_carlo_samples):
            ME_pattern_sample = sample_ME_pattern(margprob_x1, transprob_ME_arr)
            ME_level = np.mean(ME_pattern_sample)
            samples.append(ME_level)
            moments += ME_level ** powers
        moments = moments / n_monte_carlo_samples

        # Try maxent
        _, ml_distribution, _ = utils_informME.maxent(
            mu=moments,
            x=ME_level_list,
            lambda0=None,
            max_iters=6000,
            tol=1e-5,
        )

        # Fallback if maxent failed
        if ml_distribution is None:
            print("⚠️ Falling back to histogram approximation")
            ml_distribution, _ = np.histogram(
                samples, bins=len(ME_level_list), range=(0,1), density=True
            )

    # Cleanup
    ml_distribution[ml_distribution < 0] = 0
    ml_distribution /= ml_distribution.sum()
    return ml_distribution, ME_level_list


def get_mml(margprob_arr):
    mml = np.mean(margprob_arr)
    return mml

def get_nme(ml_distribution):
    ml_distribution = ml_distribution[ml_distribution > 0]
    ME_entropy = np.dot(ml_distribution, -np.log2(ml_distribution))

    n_ME_levels = len(ml_distribution) 
    nme = ME_entropy / np.log2(n_ME_levels)

    # Clip to [0, 1]
    nme = max(0, min(nme, 1))    
    return nme 

def get_ME_level_classification(ml_distribution, 
                                ME_level_list, 
                                thresh=[0, 0.25, 0.5, 0.75, 1], 
                                threshold = 0.4):
    """
    Classify methylation state from a methylation level distribution.

    Returns:
        meth_class (int):
            -2 = highly unmethylated
            -1 = partially unmethylated
             0 = variable
             1 = partially methylated
             2 = highly methylated
    """
    ml_distribution = np.array(ml_distribution)
    ME_level_list = np.array(ME_level_list)
    thresh = np.array(thresh)

    # Find index for L=0.5
    idx_05 = np.where(ME_level_list == thresh[2])[0]
    halph_p_correction = ml_distribution[idx_05[0]]/2 if idx_05.size > 0 else 0

    # Compute probabilities in each region
    p_low_ME = np.sum(ml_distribution[(thresh[0] <= ME_level_list) & (ME_level_list <= thresh[1])])
    p_inter_low_ME = np.sum(ml_distribution[(thresh[1] < ME_level_list) & (ME_level_list < thresh[2])]) + halph_p_correction
    p_inter_high_ME = np.sum(ml_distribution[(thresh[2] < ME_level_list) & (ME_level_list < thresh[3])]) + halph_p_correction
    p_high_ME = np.sum(ml_distribution[(thresh[3] <= ME_level_list) & (ME_level_list <= thresh[4])])
    p_unME = p_low_ME + p_inter_low_ME
    p_ME = p_inter_high_ME + p_high_ME
    # Assign class
    # Unmeth
    var_class = 0
    if 1-threshold <= p_unME <= 1:
        if p_low_ME <= 1-threshold:
            meth_class = -1
        else:
            meth_class = -2

    # Meth
    elif 0 <= p_unME <= threshold:
        if p_high_ME <= 1-threshold:
            meth_class = 1
        else:
            meth_class = 2
    
    # Var
    else:
        meth_class = 0
        pratio1 = p_low_ME / p_unME if p_unME > 0 else 0
        pratio2 = p_high_ME / p_ME if p_ME > 0 else 0

        if pratio1 <= threshold and pratio2 <= threshold:
            var_class = 1
        elif threshold <= pratio1 <= 1-threshold and threshold <= pratio2 <= 1-threshold:
            var_class = 2
        elif 1-threshold <= pratio1 <= 1 and 1-threshold <= pratio2 <= 1:
            var_class = 3

    return meth_class, var_class

def get_ME_entropy_classification(nme):
    if nme >= 0.99:
        me_class = 2
    elif nme >= 0.92:
        me_class = 1
    elif nme > 0.44:
        me_class = 0
    elif nme > 0.28:
        me_class = -1
    elif nme >= 0:
        me_class = -2
    else:
        me_class = float('nan')  # or raise ValueError("Invalid NME value")
    return me_class

def get_mml_nme(margprob_arr, transprob_ME_arr, 
                max_CpGs = 18, 
                thresh = [0, 0.25, 0.5, 0.75, 1],
                threshold = 0.4):
    mml = get_mml(margprob_arr)
    ml_distribution, ME_level_list = get_ME_level_distribution(margprob_arr, transprob_ME_arr, max_CpGs)
    nme = get_nme(ml_distribution)
    meth_class, var_class = get_ME_level_classification(ml_distribution, 
                                                        ME_level_list, 
                                                        thresh = thresh, 
                                                        threshold = threshold)
    me_class = get_ME_entropy_classification(nme)

    return ml_distribution, mml, nme, meth_class, var_class, me_class

####

def eval_informME_epsilon(distance_arr, density_arr, theta, epsilon,
                        ):
    theta_eps = theta * (1 + epsilon)
    informME_epsilon = c_informME.eval_informME(distance_arr, 
                                    density_arr, 
                                    theta = theta_eps,
                                    )
    informME_epsilon["distance_arr"] = distance_arr
    informME_epsilon["density_arr"] = density_arr
    informME_epsilon["informME_params"] = theta_eps
    return informME_epsilon

# Entropic sensistivity index
def get_esi(ml_distribution, ml_distribution_epsilon, epsilon):
    nme_epsilon = get_nme(ml_distribution_epsilon)
    nme = get_nme(ml_distribution)
    n_CpGs = len(ml_distribution)-1
    da  = np.abs(nme_epsilon-nme)/np.log2(n_CpGs+1)
    esi = da / epsilon
    return esi

# # methylation sensitivity indices
def get_msi(ml_distribution, ml_distribution_epsilon, epsilon):
    ml_distribution_mask = ml_distribution > 0
    ml_distribution_epsilon_mask = ml_distribution_epsilon > 0

    term1 = np.sum(
        ml_distribution[ml_distribution_mask] *
        np.log2(2 * ml_distribution[ml_distribution_mask] /
                (ml_distribution[ml_distribution_mask] + ml_distribution_epsilon[ml_distribution_mask]))
    ) / 2

    term2 = np.sum(
        ml_distribution_epsilon[ml_distribution_epsilon_mask] *
        np.log2(2 * ml_distribution_epsilon[ml_distribution_epsilon_mask] /
                (ml_distribution[ml_distribution_epsilon_mask] + ml_distribution_epsilon[ml_distribution_epsilon_mask]))
    ) / 2

    # Clip to avoid sqrt of tiny negatives
    total = term1 + term2
    if total < 0 and total > -1e-12:  # tolerance threshold
        total = 0

    Dpqa = np.sqrt(total)
    msi = Dpqa / epsilon
    return msi

def get_esi_msi(epsilon, 
                margprob_eps_arr,
                transprob_eps_arr,
               ml_distribution,
               max_CpGs,
               ):
    ml_distribution_epsilon, _ = get_ME_level_distribution(margprob_eps_arr, 
                                                           transprob_ME_arr = transprob_eps_arr,
                                                            max_CpGs = max_CpGs)
    esi = get_esi(ml_distribution, ml_distribution_epsilon, epsilon)
    msi = get_msi(ml_distribution, ml_distribution_epsilon, epsilon)
    return ml_distribution_epsilon, esi, msi

####

def get_odds_ratio(margprob_arr):
    odds_ratio = margprob_arr / (1-margprob_arr)
    return odds_ratio

# turnover ratios
def get_turn(odds_ratio):
    turn = np.mean(np.log2(odds_ratio))
    return turn 

# channel capacities
def get_cap(odds_ratio):
    capvals = np.zeros_like(odds_ratio)
    mask_ge1 = odds_ratio >= 1
    mask_lt1 = odds_ratio < 1
    # For lambda >= 1
    x_ge1 = odds_ratio[mask_ge1] / (1 + odds_ratio[mask_ge1])
    capvals[mask_ge1] = (
        1 - 0.52 * utils_informME.h_func(x_ge1) / (1 + odds_ratio[mask_ge1])
    )
    # For lambda < 1
    x_lt1 = odds_ratio[mask_lt1] / (1 + odds_ratio[mask_lt1])
    capvals[mask_lt1] = (
        1 - 0.52 * utils_informME.h_func(x_lt1) * odds_ratio[mask_lt1] / (1 + odds_ratio[mask_lt1])
    )
    return np.mean(capvals)

# relative dissipated energies
def get_rde(odds_ratio):
    rdevals = np.zeros_like(odds_ratio)
    mask_ge1 = odds_ratio >= 1
    mask_lt1 = odds_ratio < 1
    rdevals[mask_lt1] = np.log2((1 + odds_ratio[mask_lt1]) / (2 * odds_ratio[mask_lt1])) + 4.76
    rdevals[mask_ge1] = np.log2((1 + odds_ratio[mask_ge1]) / 2) + 4.76
    return np.mean(rdevals)

def get_turn_cap_rde(margprob_arr):
    odds_ratio = get_odds_ratio(margprob_arr)
    turn = get_turn(odds_ratio)
    cap = get_cap(odds_ratio)
    rde = get_rde(odds_ratio)
    return turn, cap, rde

####

# Diff analysis

def get_diff(val1, val2):
    diff_vals = val1-val2
    return diff_vals

def _norm(p):
    p = np.asarray(p, float); s = p.sum()
    return p/s if s > 0 else p

from scipy.spatial.distance import jensenshannon as _scipy_js

def _jsd(p1, p2):
    p1 = _norm(p1); p2 = _norm(p2)
    s1, s2 = p1.sum(), p2.sum()
    if s1 == 0.0 and s2 == 0.0:
        return 0.0
    # SciPy expects proper distributions; our _norm ensures that.
    try:
        d = float(_scipy_js(p1, p2, base=2.0))
    except Exception:
        d = float("nan")
    if not np.isfinite(d) or d < 0.0:
        # Fallback: the theoretical max JSD distance (base-2) is 1.
        # If one vector is all-zero and the other nonzero, the true value is sqrt(0.5) ≈ 0.7071.
        if (s1 == 0.0) ^ (s2 == 0.0):
            return float(np.sqrt(0.5))
        return 0.0
    return d


def _dmml(p1, p2):
    p1, p2 = _norm(p1), _norm(p2)
    n = len(p1) - 1
    pD = np.convolve(p1, p2[::-1])
    Dvals = np.linspace(-1, 1, 2*n + 1)
    return float(np.dot(pD, Dvals))

def _dmu(p1, p2, thresh_dmu, thresh, min_num_cpg):
    p1, p2 = _norm(p1), _norm(p2)
    n = len(p1) - 1
    pD = np.convolve(p1, p2[::-1])
    D = np.linspace(-1, 1, 2*n + 1)
    t = thresh_dmu
    q = np.array([
        pD[(D>=t[0]) & (D<=t[1])].sum(),
        pD[(D> t[1]) & (D<=t[2])].sum(),
        pD[(D> t[2]) & (D< t[3])].sum(),
        pD[(D>=t[3]) & (D< t[4])].sum(),
        pD[(D>=t[4]) & (D<=t[5])].sum()
    ])
    if n >= min_num_cpg:
        if q[2] > thresh: return 0
        denom = max(1e-12, 1 - q[2])
        if (q[0]+q[1]) / denom > thresh:
            if q[0] > thresh: return -3
            elif (q[0]+q[1]) > thresh: return -2
            else: return -1
        if (q[3]+q[4]) / denom > thresh:
            if q[4] > thresh: return 3
            elif (q[3]+q[4]) > thresh: return 2
            else: return 1
    else:
        if q[2] > thresh: return 0
        if (q[3]+q[4]) > thresh: return 2
        if (q[0]+q[1]) > thresh: return -2
        denom = max(1e-12, 1 - q[2])
        if (q[3]+q[4]) / denom > thresh: return 1
        if (q[0]+q[1]) / denom > thresh: return -1
    return np.nan