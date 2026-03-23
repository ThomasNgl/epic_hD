import numpy as np

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


# ---------- Average log-likelihood ----------
def get_ising_aveloglikelihood_fast(MEmat_contig_list, CpGs_start, CpGs_end,
                                    Ar, Cr):
    """
    Compute average log-likelihood over contigs.
    """
    logZ1, logZ0, logZ = get_logZ_np(Ar, Cr)
    logZ1t, logZ0t, _  = get_logZtilde_np(Ar, Cr)

    total = 0.0
    for x, q, end in zip(MEmat_contig_list, CpGs_start, CpGs_end):
        s = int(end) - int(q)
        total += get_log_margpmf_fast(q, s, x,
                                      logZ1, logZ0, logZ,
                                      logZ1t, logZ0t,
                                      Ar, Cr)
    return total / len(CpGs_start)


# ---------- Objective for optimizer ----------
def get_ising_func2minimize(theta,
                            MEmat_contig_list, CpGs_start, CpGs_end,
                            distance_arr, density_arr, precision = None):
    """
    Returns NEGATIVE average log-likelihood (suitable for minimizers).
    """
    Ar, Cr = get_ising_coef(distance_arr, density_arr, theta)
    aveloglik = get_ising_aveloglikelihood_fast(MEmat_contig_list, CpGs_start, CpGs_end, Ar, Cr)
    return -aveloglik
