
import numpy as np
from scipy.optimize import root

meth_dict = {-2: "highly unmethylated",
            -1: "partially unmethylated",
            0: "variably methylated",
            1: "partially methylated",
            2: "highly methylated",}

meth_var_dict = {1: "mixed",
                2: "highly mixed",
                3: "bistable",}

entr_dict = {-2: "highly ordered",
            -1: "moderatly ordered",
            0: "weakly ordered/disordered",
            1: "moderately disordered",
            2: "highly disordered",}

# Differential methylation (DMU) classes: -3..3
d_ml_dict = {
    -3: "strongly hypomethylated",
    -2: "moderately hypomethylated",
    -1: "weakly hypomethylated",
     0: "isomethylated",
     1: "weakly hypermethylated",
     2: "moderately hypermethylated",
     3: "strongly hypermethylated",
}

# Entropy-based (DEU) classes: -3..3
d_me_dict = {
    -3: "strongly hypoentropic",
    -2: "moderately hypoentropic",
    -1: "weakly hypoentropic",
     0: "isoentropic",
     1: "weakly hyperentropic",
     2: "moderately hyperentropic",
     3: "strongly hyperentropic",
}


chr_map = {
        '1': 'NC_000067.6',
        '2': 'NC_000068.7',
        '3': 'NC_000069.6',
        '4': 'NC_000070.6',
        '5': 'NC_000071.6',
        '6': 'NC_000072.6',
        '7': 'NC_000073.6',
        '8': 'NC_000074.6',
        '9': 'NC_000075.6',
        '10': 'NC_000076.6',
        '11': 'NC_000077.6',
        '12': 'NC_000078.6',
        '13': 'NC_000079.6',
        '14': 'NC_000080.6',
        '15': 'NC_000081.6',
        '16': 'NC_000082.6',
        '17': 'NC_000083.6',
        '18': 'NC_000084.6',
        '19': 'NC_000085.6',
        'X': 'NC_000086.7',
        'Y': 'NC_000087.7',
        'M': 'NC_005089.1'
    }

def chr2refseq(chr_num):
    """
    Maps mouse chromosome names/numbers (GRCm39) to their RefSeq accession numbers.
    Accepts input like '1', 'chr1', 'X', 'chrX', 1, etc.
    """

    # Handle integer input
    if isinstance(chr_num, int):
        chr_key = str(chr_num)
    elif isinstance(chr_num, str):
        s = chr_num.strip().upper()
        # Accept "CHR1", "chrX", etc.
        if s.startswith('CHR'):
            chr_key = s[3:]
        else:
            chr_key = s
        # Handle "MT" as an alias for "M"
        if chr_key == 'MT':
            chr_key = 'M'
    else:
        raise ValueError("Input chromosome must be a string or integer.")

    if chr_key in chr_map:
        return chr_map[chr_key]
    else:
        raise ValueError(f"Chromosome '{chr_num}' not recognized. Valid inputs: {list(chr_map.keys()) + ['chr'+k for k in chr_map.keys()] + ['MT']}")

def check_chr_name(chr_name, bam_chromosome_names):
    " Check if the chromosome is availbale in the bam file"
    tried_chr_names = [chr_name]
    # Remove 'chr' prefix if present
    if chr_name.startswith('chr'):
        tried_chr_names.append(chr_name[3:])
    # Try ref conversion 
    refchr_name = chr2refseq(chr_num=chr_name)
    tried_chr_names.append(refchr_name)

    # Find a matching chromosome name
    found_chr = None
    for cname in tried_chr_names:
        if cname in bam_chromosome_names:
            found_chr = cname
            break

    return found_chr, tried_chr_names

def get_refseq(fasta_file, chr_name, ref_pos):
    ref_seq = fasta_file[chr_name][ref_pos[0]:ref_pos[-1]]
    return ref_seq


def maxent(mu, x, lambda0=None, max_iters=5000, tol=1e-5, debug=True):
    """
    Find maximum entropy distribution with moment constraints.

    Parameters
    ----------
    mu : array-like, shape (N,)
        Moments E[x^k], for k=1..N (mu[0] = E[x])
    x : array-like
        Points where to evaluate p(x)
    lambda0 : array-like, optional
        Initial guess for Lagrange multipliers (shape N)
    max_iters : int
        Maximum iterations for solver
    tol : float
        Tolerance for solver
    debug : bool
        If True, print diagnostics when something goes wrong

    Returns
    -------
    lambda_ : ndarray
        Optimal Lagrange multipliers
    p : ndarray
        Maximum entropy distribution over x
    entropy : float
        Shannon entropy of the distribution
    """
    mu = np.asarray(mu).flatten()
    mu = np.concatenate([[1], mu])  # add mu_0 = 1 (normalization)
    x = np.asarray(x).flatten()
    N = len(mu)
    lx = len(x)

    # Build power matrix
    fin = np.ones((lx, 2*N-1))
    for n in range(1, 2*N-1):
        fin[:, n] = x * fin[:, n-1]

    if lambda0 is None:
        lambda_ = np.zeros(N)
        lambda_[0] = np.log(x[-1] - x[0]) if x[-1] > x[0] else 0.0
    else:
        lambda_ = np.array(lambda0).flatten()

    def safe_normalize(exponents, context=""):
        """Exponentiate and normalize with stability shift."""
        exponents = exponents - np.max(exponents)  # stability trick
        p = np.exp(exponents)
        Z = np.sum(p)

        if Z == 0 or not np.isfinite(Z):
            if debug:
                print(f"⚠️ Normalization failed in {context}")
                print(f"   mu = {mu}")
                print(f"   x = {x}")
                print(f"   lambda = {lambda_}")
                print(f"   exponents range = [{np.min(exponents)}, {np.max(exponents)}]")
            return None
        return p / Z

    def equations(lambda_):
        exponents = -(fin[:, :N] @ lambda_)
        p = safe_normalize(exponents, context="equations()")
        if p is None:
            return np.ones(N) * np.inf  # force solver failure

        moments = np.zeros(N)
        for k in range(N):
            moments[k] = np.sum(fin[:, k] * p)
        return moments - mu

    # Solve for lambda
    sol = root(equations, lambda_, method="lm", tol=tol, options={"maxiter": max_iters})

    if not sol.success:
        if debug:
            print("⚠️ Warning: Maximum Entropy did not converge:", sol.message)
        return None, None, None  # signal failure

    lambda_ = sol.x
    exponents = -(fin[:, :N] @ lambda_)
    p = safe_normalize(exponents, context="final normalization")
    if p is None:
        return None, None, None

    # Shannon entropy
    entropy = -np.sum(p * np.log(p + 1e-20))
    return lambda_, p, entropy


def h_func(P):
    """
    Computes log2-based binary entropy for each element of P.

    Parameters
    ----------
    P : array_like
        QxR matrix of probabilities.

    Returns
    -------
    ENTR : ndarray
        QxR matrix of entropies for each (q, r) element.
    """
    P = np.asarray(P)
    ENTR = np.zeros_like(P)
    mask = (P > 0) & (P < 1)
    ENTR[mask] = -P[mask] * np.log2(P[mask]) - (1 - P[mask]) * np.log2(1 - P[mask])
    return ENTR


