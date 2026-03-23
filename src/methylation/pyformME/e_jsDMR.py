import re
import math
import numpy as np
import pandas as pd

# ----------------- parsing helpers -----------------

def _parse_subregion(sr):
    """
    Accepts pandas.Interval, '[start, end)' string, or (start,end) tuple/list.
    Returns (start, end) as ints.
    """
    # pandas.Interval?
    if hasattr(sr, "left") and hasattr(sr, "right"):
        return int(sr.left), int(sr.right)
    # string "[a, b)"
    if isinstance(sr, str):
        m = re.match(r"\s*\[(\d+)\s*,\s*(\d+)\s*\)\s*$", sr)
        if not m:
            raise ValueError(f"Cannot parse subregion string: {sr!r}")
        return int(m.group(1)), int(m.group(2))
    # 2-tuple/list
    if isinstance(sr, (tuple, list)) and len(sr) == 2:
        return int(sr[0]), int(sr[1])
    raise TypeError(f"Unsupported subregion type: {type(sr)}")

def flatten_concat_df(df_chr, chrom_name):
    """
    From one 'concat_subregions' per-chromosome DataFrame -> flat table
    with columns: chrom,start,end,jsd.
    - Requires a 'subregion' column and a 'jsd' column.
    """
    starts, ends = zip(*df_chr["subregion"].map(_parse_subregion))
    flat = pd.DataFrame({
        "chrom": chrom_name,
        "start": np.asarray(starts, dtype=int),
        "end":   np.asarray(ends,   dtype=int),
        "jsd":   pd.to_numeric(df_chr["jsd"], errors="coerce"),
    })
    # clean & sort
    flat = flat[np.isfinite(flat["jsd"])].copy()
    flat["jsd"] = flat["jsd"].clip(lower=0.0, upper=1.0)  # JSD distance ∈ [0,1]
    flat.sort_values(["chrom","start","end"], inplace=True, ignore_index=True)
    return flat

# ----------------- minimal stats -----------------

def _ecdf(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x.sort()
    n = x.size
    if n == 0:
        return lambda t: np.zeros_like(np.asarray(t, dtype=float))
    def F(t):
        t = np.asarray(t, dtype=float)
        return np.searchsorted(x, t, side="right") / n
    return F

def _logit_safe(p):
    p = np.asarray(p, dtype=float)
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return np.log(p) - np.log1p(-p)

from scipy.special import erf
def _norm_cdf(z):
    return 0.5 * (1.0 + erf(np.asarray(z, dtype=float) / np.sqrt(2.0)))


def pvals_from_null_ecdf(jsd_alt, jsd_null):
    F = _ecdf(jsd_null)
    return 1.0 - F(jsd_alt)

def pvals_from_two_gaussian_mixture_logit(jsd, fit_subsample=50000, random_state=0):
    """
    Fit 2-Gaussian mixture on logit(JSD) using a (possibly) subsampled set,
    take lower-mean component as null; then score all points.
    """
    import numpy as np
    rng = np.random.default_rng(random_state)

    x_all = _logit_safe(jsd).astype(float)
    x_fit = x_all[np.isfinite(x_all)]
    n = x_fit.size
    if n == 0:
        return np.ones_like(jsd, dtype=float)

    if fit_subsample and n > fit_subsample:
        idx = rng.choice(n, size=fit_subsample, replace=False)
        x = x_fit[idx]
    else:
        x = x_fit

    # tiny EM
    mu = np.array([np.percentile(x, 25), np.percentile(x, 75)], dtype=float)
    sig = np.array([np.std(x), np.std(x)], dtype=float) + 1e-6
    lam = np.array([0.5, 0.5], dtype=float)

    prev_ll = -np.inf
    for _ in range(50):  # fewer iters with early stop
        # responsibilities
        r0 = lam[0] * np.exp(-0.5*((x - mu[0])/sig[0])**2) / (sig[0] + 1e-12)
        r1 = lam[1] * np.exp(-0.5*((x - mu[1])/sig[1])**2) / (sig[1] + 1e-12)
        s  = r0 + r1 + 1e-16
        w0, w1 = r0/s, r1/s

        # M-step
        N0, N1 = w0.sum(), w1.sum()
        lam = np.array([N0, N1]) / x.size
        mu  = np.array([(w0 @ x)/N0, (w1 @ x)/N1])
        sig = np.sqrt(np.array([(w0 @ (x - mu[0])**2)/N0,
                                (w1 @ (x - mu[1])**2)/N1]) + 1e-6)

        # early stop on loglik
        ll = np.sum(np.log(s))
        if ll - prev_ll < 1e-6:
            break
        prev_ll = ll

    null_idx = int(mu[0] < mu[1])
    mu0, sd0 = mu[null_idx], sig[null_idx]

    z = (x_all - mu0) / (sd0 + 1e-12)
    p = 1.0 - _norm_cdf(z)
    p[~np.isfinite(p)] = 1.0
    return np.clip(p, 0, 1)



def fdr_adjust(pvals, method="BY"):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    v = p[order]
    c_m = (np.sum(1.0/np.arange(1, m+1)) if method.upper()=="BY" else 1.0)
    q = np.empty_like(v)
    prev = 1.0
    for i in range(m-1, -1, -1):
        rank = i + 1
        prev = min(prev, (v[i]*m*c_m)/rank)
        q[i] = prev
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out

def smooth_jsd(flat, bandwidth_bp=50_000, sigma_bp=None, step_bp=200):
    """
    Fast Nadaraya–Watson-style smoothing:
    1) bin midpoints to a fixed grid (weighted by interval length),
    2) 1D Gaussian smooth on the grid,
    3) interpolate back to original positions.

    step_bp: grid step (bp). 150–500 bp works well; 200 is a good default.
    """
    import numpy as np
    try:
        from scipy.ndimage import gaussian_filter1d
        _have_scipy = True
    except Exception:
        _have_scipy = False

    start = flat["start"].to_numpy()
    end   = flat["end"].to_numpy()
    mids  = 0.5 * (start + end)
    lens  = np.maximum(1, end - start).astype(float)  # weights
    jsd   = flat["jsd"].to_numpy().astype(float)

    if mids.size <= 2:
        return jsd.copy()

    lo = int(mids.min())
    hi = int(mids.max()) + 1
    edges = np.arange(lo, hi + step_bp, step_bp, dtype=int)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Bin with length weights: average JSD per grid bin
    jsd_sum = np.histogram(mids, bins=edges, weights=jsd * lens)[0]
    w_sum   = np.histogram(mids, bins=edges, weights=lens)[0]
    grid = np.divide(jsd_sum, np.maximum(w_sum, 1e-12))

    # Fill empty bins by linear interpolation to avoid flat zeros
    empty = w_sum <= 0
    if empty.any():
        idx = np.flatnonzero(~empty)
        if idx.size >= 2:
            grid[empty] = np.interp(np.flatnonzero(empty), idx, grid[idx])
        else:
            grid[empty] = grid[~empty][0]

    # Choose sigma in bins (paper: bw=50kb, sd≈18.5kb)
    if sigma_bp is None:
        sigma_bp = 18500.0 if bandwidth_bp == 50_000 else max(1.0, bandwidth_bp / 2.7)
    sigma_bins = max(1e-6, sigma_bp / float(step_bp))

    if _have_scipy:
        sm = gaussian_filter1d(grid, sigma=sigma_bins, mode="nearest")
    else:
        # Fallback: manual Gaussian kernel + conv (truncate at ±3σ)
        radius = int(np.ceil(3.0 * sigma_bins))
        x = np.arange(-radius, radius + 1)
        k = np.exp(-0.5 * (x / sigma_bins) ** 2)
        k /= k.sum()
        sm = np.convolve(grid, k, mode="same")

    # Interpolate smoothed grid back to original midpoints
    return np.interp(mids, centers, sm).astype(float)


def call_sites(jsd, method="BY", alpha=0.01, null_jsd=None, max_sqs=250.0,
               mixture_subsample=50000, random_state=0):
    jsd = np.asarray(jsd, dtype=float)
    if null_jsd is not None and len(null_jsd) > 0:
        p = pvals_from_null_ecdf(jsd, np.asarray(null_jsd, dtype=float))
    else:
        p = pvals_from_two_gaussian_mixture_logit(
            jsd, fit_subsample=mixture_subsample, random_state=random_state
        )
    q = fdr_adjust(p, method=method)
    sqs = -10.0 * np.log10(np.clip(q, 1e-300, 1.0))
    sqs = np.minimum(sqs, max_sqs)
    sig = q <= alpha
    return q, sqs, sig


# ----------------- merging (simple "closing") -----------------

def merge_significant_regions(df, sig_mask, bandwidth_bp=50_000, gu_size=150.0):
    """
    Merge consecutive significant subregions if the gap <= bandwidth_bp.
    Score ~ sum(SQS * length) / GUsize (similar spirit to R's binnedSumms/GUsize).
    """
    chrom = df["chrom"].to_numpy()
    start = df["start"].to_numpy()
    end   = df["end"].to_numpy()
    sqs   = df["SQS"].to_numpy()

    out = []
    i, n = 0, len(df)
    while i < n:
        if not sig_mask[i]:
            i += 1
            continue
        cchr = chrom[i]
        lo   = start[i]
        hi   = end[i]
        wsum = sqs[i] * max(1, end[i]-start[i])
        cnt  = 1
        j = i + 1
        while j < n and sig_mask[j] and chrom[j] == cchr and (start[j] - hi) <= bandwidth_bp:
            wsum += sqs[j] * max(1, end[j]-start[j])
            hi = max(hi, end[j])
            cnt += 1
            j += 1
        score = wsum / float(gu_size)
        out.append((cchr, int(lo), int(hi), float(score), int(cnt)))
        i = j

    return pd.DataFrame(out, columns=["chrom", "start", "end", "score", "n_subregions"])

# ----------------- driver for your concat_subregions() output -----------------

def call_dmrs_from_concat_list(concat_list, chr_names=None, alpha=0.01, method="BY",
                               bandwidth_bp=50_000, gu_size=150.0, null_jsd_by_chr=None):
    """
    concat_list: list of per-chrom DataFrames as returned by your concat_subregions(df_list).
                 Each must have columns: 'subregion' and 'jsd'.
    chr_names:   list of chromosome names, same length/order as concat_list.
                 If None, uses 'chr1','chr2',... by index.
    null_jsd_by_chr: optional dict {chrom: 1D array of null JSDs} for ref/ref ECDF p-values.
    Returns dict chrom -> (flat_df_with_q_sqs, dmrs_df)
    """
    out = {}
    if chr_names is None:
        chr_names = [f"chr{k+1}" for k in range(len(concat_list))]

    for df_chr, chrom in zip(concat_list, chr_names):
        flat = flatten_concat_df(df_chr, chrom)
        null_jsd = None
        if null_jsd_by_chr is not None:
            null_jsd = null_jsd_by_chr.get(chrom, None)

        q, sqs, sig = call_sites(flat["jsd"].to_numpy(), method=method, alpha=alpha,
                                 null_jsd=null_jsd, max_sqs=250.0)
        flat = flat.assign(qval=q, SQS=sqs)
        dmrs = merge_significant_regions(flat, sig, bandwidth_bp=bandwidth_bp, gu_size=gu_size)
        out[chrom] = (flat, dmrs)
    return out

# ----------------- example usage -----------------
# Given:
#   df_list_per_chr = [...]  # one original df per chr (each has 'region' and 'diffME_analysis')
#   concat_list = concat_subregions(df_list_per_chr)  # -> one df per chr, each row is a subregion and has 'subregion' + 'jsd'
#
# Then:
# results = call_dmrs_from_concat_list(concat_list,
#                                      chr_names=["chr1","chr2",...],
#                                      alpha=0.01, method="BY",
#                                      bandwidth_bp=50_000, gu_size=150.0)
#
# # Save BEDs:
# for chrom, (flat, dmrs) in results.items():
#     dmrs[["chrom","start","end","score"]].to_csv(f"{chrom}_dmrs.bed", sep="\t", header=False, index=False)
