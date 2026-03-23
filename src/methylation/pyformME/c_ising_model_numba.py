# c_ising_model_numba.py
# ----------------------------------------------------------------------
# Accelerated Ising model utilities:
#   - (A) Fuse contig loop and JIT with Numba (parallel) for avg log-likelihood
#   - (C) Deduplicate identical segments before evaluation
#
# This file mirrors the original numpy implementation but replaces
# get_ising_aveloglikelihood_fast(...) with a dedup + JIT-compiled version.
# All other functions keep the original behavior/signatures so you can
# drop this in as a replacement module.
# ----------------------------------------------------------------------

import numpy as np

from numba import njit, prange

# ---------- Coefficients ----------
def get_ising_coef(distance_arr, density_arr, theta,
                   is_first_CpG=None, is_last_CpG=None):
    """
    Ar has length R (sites); Cr has length R-1 (bonds).
    distance_arr must be length R-1, density_arr length R.
    """
    Ar = theta[0] + theta[1] * np.asarray(density_arr, dtype=np.float64)
    Cr = (theta[2] / np.asarray(distance_arr, dtype=np.float64)).astype(np.float64)

    if is_first_CpG or is_first_CpG is None:
        Ar[0] = theta[3]
    if is_last_CpG or is_last_CpG is None:
        Ar[-1] = theta[4]
    return Ar.astype(np.float64, copy=False), Cr.astype(np.float64, copy=False)[:-1]


# ---------- Forward messages (Z) ----------
def get_logZ_np(Ar, Cr):
    """
    Fast/stable forward recursion using numpy.float64 and logaddexp.
    Returns: logZ1, logZ0, logZ (all numpy float64)
    """
    Ar = np.asarray(Ar, dtype=np.float64)
    Cr = np.asarray(Cr, dtype=np.float64)
    R = Ar.shape[0]
    assert Cr.shape[0] == R - 1, "Cr must have length R-1"

    logZ1 = np.zeros(R, dtype=np.float64)  # base: log(1)
    logZ0 = np.zeros(R, dtype=np.float64)

    # r = R-2 down to 1
    for r in range(R - 2, 0, -1):
        a_r1 = Ar[r + 1]
        c_r = Cr[r]
        logZ1[r] = np.logaddexp(-a_r1 - c_r + logZ0[r + 1],
                                 +a_r1 + c_r + logZ1[r + 1])
        logZ0[r] = np.logaddexp(-a_r1 + c_r + logZ0[r + 1],
                                 +a_r1 - c_r + logZ1[r + 1])

    # boundary (n = 0)
    a1, a0, c0 = Ar[1], Ar[0], Cr[0]
    logZ1[0] = np.logaddexp(a0 - a1 - c0 + logZ0[1],
                            a0 + a1 + c0 + logZ1[1])
    logZ0[0] = np.logaddexp(-a0 - a1 + c0 + logZ0[1],
                            -a0 + a1 - c0 + logZ1[1])

    logZ = np.logaddexp(logZ0[0], logZ1[0])
    return logZ1, logZ0, logZ


# ---------- Backward messages (Z~) ----------
def get_logZtilde_np(Ar, Cr):
    """
    Fast/stable backward recursion using numpy.float64 and logaddexp.
    Returns: logZ1tilde, logZ0tilde, logZtilde
    """
    Ar = np.asarray(Ar, dtype=np.float64)
    Cr = np.asarray(Cr, dtype=np.float64)
    R = Ar.shape[0]
    assert Cr.shape[0] == R - 1, "Cr must have length R-1"

    logZ1t = np.zeros(R, dtype=np.float64)
    logZ0t = np.zeros(R, dtype=np.float64)

    # first two boundary values (index 1)
    a0, a1, c0 = Ar[0], Ar[1], Cr[0]
    logZ1t[1] = np.logaddexp(-a0 + a1 - c0, +a0 + a1 + c0)
    logZ0t[1] = np.logaddexp(-a0 - a1 + c0, +a0 - a1 - c0)

    # r = 2..R-1
    for r in range(2, R):
        a_r = Ar[r]
        c_rm1 = Cr[r - 1]
        logZ1t[r] = np.logaddexp(a_r - c_rm1 + logZ0t[r - 1],
                                 a_r + c_rm1 + logZ1t[r - 1])
        logZ0t[r] = np.logaddexp(-a_r + c_rm1 + logZ0t[r - 1],
                                 -a_r - c_rm1 + logZ1t[r - 1])

    logZt = np.logaddexp(logZ0t[-1], logZ1t[-1])
    return logZ1t, logZ0t, logZt


# ---------- Single-site marginals, vectorized ----------
def get_margprob_arr_np(logZ1, logZ0, logZ, logZ1tilde, logZ0tilde):
    """
    P(x_r = 1) for all r, vectorized:
      log P(x_r=1) = -logZ + logZ1tilde[r] + logZ1[r]
    """
    logP1 = -logZ + logZ1tilde + logZ1
    return np.exp(logP1)


# ---------- Transition probabilities, vectorized ----------
def get_transprob_arr_np(logZ1, logZ0, Ar, Cr):
    """
    For r=0..R-2:
      P(x_{r+1}=1 | x_r=0) and P(x_{r+1}=1 | x_r=1)
    """
    Ar = np.asarray(Ar, dtype=np.float64)
    Cr = np.asarray(Cr, dtype=np.float64)
    R = Ar.shape[0]
    trans = np.zeros((R, 2), dtype=np.float64)

    a_next = Ar[1:]           # length R-1
    c = Cr                    # length R-1

    # log phi terms (note boundary add at r==0)
    logphi_0to1 = a_next - c
    logphi_1to1 = a_next + c
    if R > 1:
        logphi_0to1[0] += -Ar[0]
        logphi_1to1[0] += +Ar[0]

    log_trans_0to1 = logphi_0to1 + logZ1[1:] - logZ0[:-1]
    log_trans_1to1 = logphi_1to1 + logZ1[1:] - logZ1[:-1]

    trans[:-1, 0] = np.exp(log_trans_0to1)
    trans[:-1, 1] = np.exp(log_trans_1to1)
    return trans


# ---------- General segment probability (log) ----------
def get_log_margpmf_fast(q, s, x_q_qPLUSs,
                         logZ1, logZ0, logZ,
                         logZ1tilde, logZ0tilde,
                         Ar, Cr):
    """
    Vectorized compute of log P(x_q...x_{q+s}) using messages + pair terms.
    x_q_qPLUSs: array-like of 0/1, length s+1
    """
    x = np.asarray(x_q_qPLUSs, dtype=np.int8)
    xs = 2 * x - 1  # in {-1, +1}

    # -log Z
    logp = -logZ

    # endpoint messages
    logp += (logZ1tilde[q] if x[0] else logZ0tilde[q])
    logp += (logZ1[q + s]  if x[-1] else logZ0[q + s])

    # internal pair factors sum
    if s > 0:
        a_next = Ar[q + 1:q + s + 1]   # length s
        c_seg  = Cr[q:q + s]           # length s
        logp += np.sum((a_next + c_seg * xs[:-1]) * xs[1:])

        # boundary phi add only if the segment includes r=0
        if q == 0:
            logp += Ar[0] * xs[0]

    return float(logp)


# ---------- Helpers for dedup + packing ----------
def _pack_segments_bool_to_flat(list_of_segments):
    """
    Flatten list of 0/1 arrays and build offsets (so we can avoid Python slicing in JIT).
    Returns:
      flat (uint8), offsets (int64, len = n_unique+1)
    """
    lengths = np.fromiter((len(x) for x in list_of_segments), dtype=np.int64)
    offsets = np.empty(len(lengths) + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(lengths, out=offsets[1:])
    flat = np.empty(offsets[-1], dtype=np.uint8)
    pos = 0
    for x in list_of_segments:
        n = len(x)
        flat[pos:pos+n] = np.asarray(x, dtype=np.uint8, order="C")
        pos += n
    return flat, offsets


def _dedup_segments(MEmat_contig_list, CpGs_start, CpGs_end):
    """
    Deduplicate identical (q, end, x) triples.
    Returns unique_x_list, unique_q (int64), unique_end (int64), counts (int64)
    """
    from collections import OrderedDict

    uniq = OrderedDict()
    for x, q, end in zip(MEmat_contig_list, CpGs_start, CpGs_end):
        x_u8 = np.asarray(x, dtype=np.uint8)
        key = (int(q), int(end), np.packbits(x_u8).tobytes(), len(x_u8))
        if key in uniq:
            uniq[key][2] += 1
        else:
            # store a representative x alongside q, end, and count
            uniq[key] = [int(q), int(end), 1, x_u8]

    unique_q   = []
    unique_end = []
    counts     = []
    unique_x   = []
    for (_, _ , _, L), (q, end, c, xrep) in zip(uniq.keys(), uniq.values()):
        unique_q.append(q)
        unique_end.append(end)
        counts.append(c)
        unique_x.append(xrep)

    return unique_x, np.asarray(unique_q, dtype=np.int64), np.asarray(unique_end, dtype=np.int64), np.asarray(counts, dtype=np.int64)


# ---------- JIT kernel (fused, weighted) ----------
@njit(fastmath=True, parallel=True, cache=True)
def _sum_loglik_numba(flat_x, offsets, q_arr, end_arr, counts,
                        Ar, Cr, logZ1, logZ0, logZ, logZ1t, logZ0t):
    """
    Compute weighted sum of log P(x_q...x_{q+s}) over unique segments.
    """
    n = q_arr.shape[0]
    total = 0.0

    for i in prange(n):
        q = int(q_arr[i])
        end = int(end_arr[i])
        s = end - q

        start = offsets[i]
        stop  = offsets[i+1]

        # -log Z
        logp = -logZ

        # endpoints
        x0 = flat_x[start]           # 0/1
        xN = flat_x[stop - 1]        # 0/1
        logp += (logZ1t[q]  if x0 == 1 else logZ0t[q])
        logp += (logZ1[q+s] if xN == 1 else logZ0[q+s])

        # boundary phi add only if the segment includes r=0
        if q == 0:
            xs0 = 1 if x0 == 1 else -1
            logp += Ar[0] * xs0

        # internal pair/factor sum
        if s > 0:
            prev = 1 if x0 == 1 else -1
            k = start
            for j in range(s):
                cur = 1 if flat_x[k + 1] == 1 else -1
                logp += (Ar[q + j + 1] + Cr[q + j] * prev) * cur
                prev = cur
                k += 1

        # weight by multiplicity
        total += logp * counts[i]

    return total

# ---------- Average log-likelihood (DEDUP + JIT) ----------
def get_ising_aveloglikelihood_fast(MEmat_contig_list, CpGs_start, CpGs_end,
                                    Ar, Cr):
    """
    Compute average log-likelihood over contigs, with:
      - (C) deduplication of identical segments (q, end, x)
      - (A) fused + JIT-parallel evaluation of the unique set

    Returns average over the ORIGINAL number of contigs.
    """
    # Precompute messages once
    logZ1, logZ0, logZ = get_logZ_np(Ar, Cr)
    logZ1t, logZ0t, _  = get_logZtilde_np(Ar, Cr)

    # Deduplicate identical segments
    uniq_x_list, uniq_q, uniq_end, counts = _dedup_segments(
        MEmat_contig_list, CpGs_start, CpGs_end
    )

    # Pack ragged unique segments into flat buffers
    flat_x, offsets = _pack_segments_bool_to_flat(uniq_x_list)

    # Ensure contiguous dtypes for JIT
    Ar   = np.ascontiguousarray(Ar, dtype=np.float64)
    Cr   = np.ascontiguousarray(Cr, dtype=np.float64)
    logZ1  = np.ascontiguousarray(logZ1, dtype=np.float64)
    logZ0  = np.ascontiguousarray(logZ0, dtype=np.float64)
    logZ1t = np.ascontiguousarray(logZ1t, dtype=np.float64)
    logZ0t = np.ascontiguousarray(logZ0t, dtype=np.float64)
    uniq_q   = np.ascontiguousarray(uniq_q, dtype=np.int64)
    uniq_end = np.ascontiguousarray(uniq_end, dtype=np.int64)
    counts   = np.ascontiguousarray(counts, dtype=np.int64)
    flat_x   = np.ascontiguousarray(flat_x, dtype=np.uint8)
    offsets  = np.ascontiguousarray(offsets, dtype=np.int64)

    # Weighted sum over unique segments
    total_weighted = _sum_loglik_numba(flat_x, offsets, uniq_q, uniq_end, counts,
                                       Ar, Cr, logZ1, logZ0, logZ, logZ1t, logZ0t)

    # Average over the original number of contigs
    n_total = int(len(CpGs_start))
    return total_weighted / n_total


# ---------- Objective for optimizer ----------
def get_ising_func2minimize(theta,
                            MEmat_contig_list, CpGs_start, CpGs_end,
                            distance_arr, density_arr, precision=None):
    """
    Returns NEGATIVE average log-likelihood (suitable for minimizers).
    """
    Ar, Cr = get_ising_coef(distance_arr, density_arr, theta)
    aveloglik = get_ising_aveloglikelihood_fast(MEmat_contig_list, CpGs_start, CpGs_end, Ar, Cr)
    return -aveloglik
