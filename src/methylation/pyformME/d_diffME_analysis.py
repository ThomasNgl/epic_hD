import os
import pandas as pd
from . import io_informME
from . import a_refCpGs_process
from . import d_ME_ana_features
from . import utils_informME
import re
import numpy as np
import pandas as pd

def diffME_analysis_region(
    df1, df2,
    on=("subregion",),
    opt_metrics=("esi","msi","turn","cap","rde"),
    thresh_diff_nme_class=(-1,-0.5,-0.3,-0.05,0.05,0.3,0.5,1),
    thresh_dmu_class=(-1,-0.55,-0.1,0.1,0.55,1),
    dmu_thresh=0.55,
    min_num_cpg=2,
):
    """
    Merge two per-subregion DataFrames and compute differential metrics.
    Adds informME-style DMU classification via convolution over methylation
    level distributions, coarse binning into q1..q5, and the ratio test.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Must share keys in `on`. If present, will use:
        - ml_distribution_[1|2] : probability vector over L (length Ncg+1)
        - Ncg_[1|2] : CpG count per subregion (optional; checked if present)
        - isModeled_[1|2] : boolean flags (optional; respected if present)
        - mml_[1|2], nme_[1|2], and any in `opt_metrics` for simple diffs.
    on : tuple[str]
        Merge keys.
    opt_metrics : tuple[str]
        Extra scalar metrics to diff if both sides exist.
    thresh_diff_nme_class : sequence[float]
        Bin edges for dnme_class (DEU-style class). Uses pd.cut (right-closed).
    thresh_dmu_class : sequence[float] (len=6)
        Edges [a,b,c,d,e,f] for the five coarse bins used by DMU classification.
    dmu_thresh : float
        Threshold used in the ratio test (default 0.55).
    min_num_cpg : int
        CpG threshold separating the high- vs low-CpG DMU branches.

    Returns
    -------
    pandas.DataFrame
        Contains `on` keys, d* columns, dnme_class, jsd, dmml, dmu, dmu_class
        when available. Columns are de-duplicated while preserving order.
    """
    df_ME_analysis_1 = df1["ME_analysis"]
    df_ME_analysis_2 = df2["ME_analysis"]
    metrics = ("mml","nme") + opt_metrics
    out = df_ME_analysis_1.merge(df_ME_analysis_2, on=list(on), suffixes=("_1","_2"))
    out["position"] = df_ME_analysis_1["position"]
    # Simple diffs for scalar metrics
    for m in metrics:
        c1, c2 = f"{m}_1", f"{m}_2"
        if c1 in out and c2 in out:
            out[f"d{m}"] = out[c2] - out[c1]

    # dNME + coarse class (note: MATLAB uses NME1 - NME2; here it's 2-1)
    if "dnme" not in out and "nme_1" in out and "nme_2" in out:
        out["dnme"] = out["nme_2"] - out["nme_1"]
    if "dnme" in out:
        out["dnme_class"] = pd.cut(
            out["dnme"],
            bins=thresh_diff_nme_class,
            labels=[-3,-2,-1,0,1,2,3],
            include_lowest=True
        )

    # Distribution-based metrics and DMU classification
    if "ml_distribution_1" in out and "ml_distribution_2" in out:
        # Optional helper functions if provided by your module
        out["jsd"] = out.apply(
            lambda r: d_ME_ana_features._jsd(r["ml_distribution_1"], r["ml_distribution_2"]),
            axis=1
        )
        out["dmml"] = out.apply(
            lambda r: d_ME_ana_features._dmml(r["ml_distribution_1"], r["ml_distribution_2"]),
            axis=1
        )
        out["dmu"] = out.apply(
            lambda r: d_ME_ana_features._dmu(r["ml_distribution_1"], r["ml_distribution_2"],
                                             thresh_dmu_class, dmu_thresh, min_num_cpg),
            axis=1
        )

        # Strict check like informME if CpG counts are provided
        if "Ncg_1" in out.columns and "Ncg_2" in out.columns:
            if not (out["Ncg_1"].fillna(-1).values == out["Ncg_2"].fillna(-1).values).all():
                raise ValueError("Inconsistent CpG counts between phenotypes (Ncg_1 != Ncg_2).")

        a, b, c, d, e, f = thresh_dmu_class

        def _dmu_class_row(row):
            # Respect modeling flags if present
            if ("isModeled_1" in row and "isModeled_2" in row) and (not bool(row["isModeled_1"]) or not bool(row["isModeled_2"])):
                return np.nan

            pL1 = np.asarray(row["ml_distribution_1"], dtype=float)
            pL2 = np.asarray(row["ml_distribution_2"], dtype=float)
            if pL1.ndim != 1 or pL2.ndim != 1 or len(pL1) != len(pL2) or len(pL1) == 0:
                return np.nan

            # Infer CpG count
            num_cpg = len(pL1) - 1
            if num_cpg <= 0:
                return np.nan

            # Normalize safely
            s1, s2 = pL1.sum(), pL2.sum()
            if s1 <= 0 or s2 <= 0:
                return np.nan
            pL1 = pL1 / s1
            pL2 = pL2 / s2

            # pD for D = L1 - L2, Dvals = -1 .. 1 with step 1/num_cpg
            pD = np.convolve(pL1, pL2[::-1])
            pD = np.maximum(pD, 0)
            sD = pD.sum()
            if sD <= 0:
                return np.nan
            pD /= sD

            Dvals = np.linspace(-1.0, 1.0, 2 * num_cpg + 1)

            # Coarse bins with MATLAB's inequalities:
            # q1: [a,b]; q2: (b,c]; q3: (c,d); q4: [d,e); q5: [e,f]
            q1 = pD[(Dvals >= a) & (Dvals <= b)].sum()
            q2 = pD[(Dvals >  b) & (Dvals <= c)].sum()
            q3 = pD[(Dvals >  c) & (Dvals <  d)].sum()
            q4 = pD[(Dvals >= d) & (Dvals <  e)].sum()
            q5 = pD[(Dvals >= e) & (Dvals <= f)].sum()

            denom = 1.0 - q3

            # High-CpG branch
            if num_cpg >= min_num_cpg:
                if q3 > dmu_thresh:
                    return 0
                if denom > 0 and ((q1 + q2) / denom) > dmu_thresh:
                    if q1 > dmu_thresh:
                        return -3
                    elif (q1 + q2) > dmu_thresh:
                        return -2
                    else:
                        return -1
                if denom > 0 and ((q4 + q5) / denom) > dmu_thresh:
                    if q5 > dmu_thresh:
                        return 3
                    elif (q4 + q5) > dmu_thresh:
                        return 2
                    else:
                        return 1
                return np.nan
            # Low-CpG branch
            else:
                if q3 > dmu_thresh:
                    return 0
                if (q4 + q5) > dmu_thresh:
                    return 2
                if (q1 + q2) > dmu_thresh:
                    return -2
                if denom > 0 and ((q4 + q5) / denom) > dmu_thresh:
                    return 1
                if denom > 0 and ((q1 + q2) / denom) > dmu_thresh:
                    return -1
                return np.nan

        out["dmu_class"] = out.apply(_dmu_class_row, axis=1)

    # ------- Build output columns with your preferred de-dup --------
    candidate_keep = (
        ["subregion", "position"]
        + [c for c in out.columns if c.startswith("d")]
        + ["dnme_class", "jsd", "dmml", "dmu", "dmu_class"]
    )
    # de-dup while preserving order, then filter to existing columns
    keep = [c for c in dict.fromkeys(candidate_keep) if c in out.columns]

    return df1["region"], out[keep]







