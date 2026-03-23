from mpmath import mp, mpf, log, exp, zeros
import numpy as np

def get_ising_coef(distance_arr, density_arr, theta,
                   is_first_CpG = None,
                    is_last_CpG = None):
    # theta = [a , b ,g, a', a'']
    Ar = theta[0] + theta[1] * density_arr

    if is_first_CpG or is_first_CpG is None:
        Ar[0] = theta[3]

    if is_last_CpG or is_last_CpG is None:
        Ar[-1] = theta[4]

    Cr = theta[2]/distance_arr
    return Ar, Cr


def get_logphi_r(Ar, Cr, r, x_r, x_rPLUS1):
    a_rPlus1 = Ar[r+1]
    c_r2rPlus1 = Cr[r]
    logphi_r = (a_rPlus1 + c_r2rPlus1 * (2*x_r - 1) ) * (2*x_rPLUS1 - 1) 
    if r == 0:
        logphi_r += Ar[0] * (2*x_r - 1) 
    return logphi_r


def get_logZ(Ar, Cr, precision = 200):
    """
    Recursive Computation of Partition Function
    usage:
        logZ1, logZ0, logZ = computeZ(An, Cn)
    An, Cn should be Python lists or arrays of floats.
    precision: decimal digits of precision (default 200)
    """
    # Set precision (decimal digits)
    mp.dps = precision
    
    R = len(Ar)
    
    # Convert inputs to mpmath.mpf for precision
    Ar = [mpf(a) for a in Ar]
    Cr = [mpf(c) for c in Cr]

    # Initialize output arrays
    logZ1 = zeros(R, 1)
    logZ0 = zeros(R, 1)
    
    # Recurse through non-boundary values
    for r in range(R-2, 0, -1):
        a_r_1 = Ar[r+1]
        c_r = Cr[r]
        logZ0_next = logZ0[r+1, 0]
        logZ1_next = logZ1[r+1, 0]

        # For Z1[r]:
        logZ1[r] = log(
            exp(-a_r_1 - c_r + logZ0_next) + # x_{r+1} = 0
              exp(a_r_1 + c_r + logZ1_next)   # x_{r+1} = 1
        )
        
        # For Z0[r]:
        logZ0[r] = log(
            exp(-a_r_1 + c_r + logZ0_next) +  # x_{r+1} = 0
            exp(a_r_1 - c_r + logZ1_next)   # x_{r+1} = 1
        )
        
    # Calculate last boundary values (n = 0)
    a_1 = Ar[1]
    a_0 = Ar[0]
    c_0 = Cr[0]
    logZ0_1 = logZ0[1, 0]
    logZ1_1 = logZ1[1, 0]

    logZ1[0] = log(
        exp( a_0 - a_1 - c_0 + logZ0_1) +
        exp( a_0 + a_1 + c_0 + logZ1_1)
    )
    logZ0[0] = log(
        exp(-a_0 - a_1 + c_0 + logZ0_1) +
        exp(-a_0 + a_1 - c_0 + logZ1_1)
    )

    # Compute log partition function
    logZ = log(exp(logZ0[0, 0]) + exp(logZ1[0, 0]))

    # Return as python lists for compatibility
    return logZ1, logZ0, logZ


def get_logZtilde(Ar, Cr, precision=200):
    """
    Backwards Recursive Computation of Partition Function

    usage:
        logZ1tilde, logZ0tilde, logZtilde = computeZtilde(An, Cn)

    An, Cn should be Python lists or arrays of floats.
    precision: decimal digits of precision (default 200)
    """
    # Set precision (decimal digits)
    mp.dps = precision

    R = len(Ar)

    # Convert to mpf for high-precision
    Ar = [mpf(a) for a in Ar]
    Cr = [mpf(c) for c in Cr]

    # Initialize outputs
    logZ1tilde = zeros(R, 1)
    logZ0tilde = zeros(R, 1)

    # Calculate first two boundary values (indexing: Python is 0-based)
    a_0 = Ar[0]
    a_1 = Ar[1]
    c_0 = Cr[0]
    logZ1tilde[1] = log(
        exp(-a_0 + a_1 - c_0) +
        exp( a_0 + a_1 + c_0)
    )
    logZ0tilde[1] = log(
        exp(-a_0 - a_1 + c_0) +
        exp( a_0 - a_1 - c_0)
    )

    # Recurse through non-boundary values
    for r in range(2, R):
        a_r = Ar[r]
        c_rM1 = Cr[r-1]
        logZ0tilde_prev = logZ0tilde[r-1, 0]
        logZ1tilde_prev = logZ1tilde[r-1, 0]
        
        logZ1tilde[r] = log(
            exp( a_r - c_rM1 + logZ0tilde_prev) +
            exp( a_r + c_rM1 + logZ1tilde_prev)
        )
        logZ0tilde[r] = log(
            exp(-a_r + c_rM1 + logZ0tilde_prev) +
            exp(-a_r - c_rM1 + logZ1tilde_prev)
        )

    # Compute log partition function
    logZtilde = log(exp(logZ0tilde[R-1, 0]) + exp(logZ1tilde[R-1, 0]))

    return logZ1tilde, logZ0tilde, logZtilde


def get_log_margpmf(
    q, s, x_q_qPLUSs,
    logZ1, logZ0, logZ,
    logZ1tilde, logZ0tilde,
    Ar, Cr
):
    """
    Compute the Marginal Probability (in log space).
    
    r, s: integers (indices, Python 0-based)
    x_r_rPLUSs: array of ints (length s+1, each entry 0 or 1)
    logZ1, logZ0: arrays of logZ1/logZ0 (output from computeZ, length >= r+s+1)
    logZ: log partition function (scalar)
    logZ1tilde, logZ0tilde: arrays from computeZtilde (length >= r+1)
    Ar, Cr: arrays of floats (lengths at least as long as needed)

    Follow equation S19 from paper.
    """
    # Set to log(1/Z)
    log_margprob = -logZ
    # Add log[\tilde{Z}_{q}(x_q)]
    if x_q_qPLUSs[0] > 0:  
        log_margprob += logZ1tilde[q, 0]
    else:
        log_margprob += logZ0tilde[q, 0]

    # Add log[Z_{r+s}(x_{r+s})]
    if x_q_qPLUSs[-1] > 0:  
        log_margprob += logZ1[q + s, 0]
    else:
        log_margprob += logZ0[q + s, 0]

    # Add sum_{n=r}^{r+s-1} log[phi(x_n, x_{n+1})]
    if s == 0:
        return log_margprob
    
    for r in range(q, q+s): # from q to q+s-1
        log_margprob += get_logphi_r(Ar, Cr, r, 
                                        x_r = x_q_qPLUSs[r-q], 
                                        x_rPLUS1 = x_q_qPLUSs[r+1-q])

    return log_margprob

def get_margpmf(
    q, s, x_q_qPLUSs,
    logZ1, logZ0, logZ,
    logZ1tilde, logZ0tilde,
    Ar, Cr
):
    log_margpmf = get_log_margpmf(
                    q, s, x_q_qPLUSs,
                    logZ1, logZ0, logZ,
                    logZ1tilde, logZ0tilde,
                    Ar, Cr
                    )
    return float(exp(log_margpmf))

def get_margprob_arr(
    logZ1, logZ0, logZ,
    logZ1tilde, logZ0tilde,
    Ar, Cr
):
    R = len(Ar)
    margprob_arr = np.zeros((R))
    for r in range(0, R):
        margprob_xr_1 = get_margpmf(q = r, 
                                    s = 0, 
                                    x_q_qPLUSs = [1],
                                    logZ1=logZ1, logZ0=logZ0, logZ=logZ,
                                    logZ1tilde=logZ1tilde, logZ0tilde=logZ0tilde,
                                    Ar=Ar, Cr = Cr)
        margprob_arr[r] = margprob_xr_1
    return margprob_arr

def get_transprob_arr(
    logZ1, logZ0, logZ,
    Ar, Cr
):
    R = len(Ar)
    transprob_arr = np.zeros((R, 2))

    for r in range(0, R-1):
        logphi_xr_0to1 = get_logphi_r(Ar, Cr, r, x_r = 0, x_rPLUS1 = 1)
        log_transp_xr_0to1 = logphi_xr_0to1 + logZ1[r+1, 0] - logZ0[r, 0]
        transp_xr_0to1 = float(exp(log_transp_xr_0to1))

        logphi_xr_1to1 = get_logphi_r(Ar, Cr, r, x_r = 1, x_rPLUS1 = 1)
        log_transp_xr_1to1 = logphi_xr_1to1 + logZ1[r+1, 0] - logZ1[r, 0]
        transp_xr_1to1 = float(exp(log_transp_xr_1to1))

        transprob_arr[r] = np.array([transp_xr_0to1, transp_xr_1to1])

    return transprob_arr

def get_ising_aveloglikelihood(MEmat_contig_list, CpGs_start, CpGs_end, 
                               An, Cn,
                               precision = 200):
    """
    Compute the average log likelihood for a dataset.
    
    An, Cn: lists/arrays of floats, length N
    dataMatrix: 2D numpy array (N, K), each column an observation (can be np.int8)
    CpGs_start, CpGs_end: lists/arrays of integers, length K, 0-based indices
    """
    n_contig_reads = len(CpGs_start)
 
    # Compute Z values
    logZ1, logZ0, logZ = get_logZ(An, Cn, precision = precision)
    # Compute Ztilde values
    logZ1tilde, logZ0tilde, _ = get_logZtilde(An, Cn, precision = precision)

    sum_log_margprob = 0.0

    for contig_read_idx in range(n_contig_reads):  # loop over all observations
        q = CpGs_start[contig_read_idx]
        s = CpGs_end[contig_read_idx] - CpGs_start[contig_read_idx]
        x_q_qPLUSs = MEmat_contig_list[contig_read_idx]
        # Fill in observation vector for current read

        sum_log_margprob += get_log_margpmf(
            q, s, x_q_qPLUSs, 
            logZ1, logZ0, logZ, 
            logZ1tilde, logZ0tilde,
            An, Cn
        )
    return float(sum_log_margprob / n_contig_reads)


def get_ising_func2minimize(theta, 
                          MEmat_contig_list, CpGs_start, CpGs_end,
                          distance_arr, density_arr, precision):
    Ar, Cr = get_ising_coef(distance_arr, density_arr, theta)

    ising_aveloglikelihood = get_ising_aveloglikelihood(MEmat_contig_list, CpGs_start, CpGs_end, Ar, Cr, precision=precision)

    ising_negaveloglikelihood = -ising_aveloglikelihood

    return ising_negaveloglikelihood

