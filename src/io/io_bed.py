import re
import pandas as pd

def df_to_bed(
    df: pd.DataFrame,
    out_path: str,
    chrom_col: str | None = None,
    start_col: str | None = None,
    end_col: str | None = None,
    add_chr_prefix: bool = False,
    drop_na: bool = True,
    sort_by: str | None = "genome",   # "genome" (chr/start), "position" (start/end only), or None
):
    """
    Convert a DataFrame with chromosome, start (closed), end (open) to a BED file (0-based, half-open).
    Writes 3-column BED: chrom  start  end  (no header).

    Parameters:
      df              : input DataFrame
      out_path        : where to save the .bed
      chrom_col       : name of chromosome column (auto-detected if None)
      start_col       : name of start column (auto-detected if None)
      end_col         : name of end column (auto-detected if None)
      add_chr_prefix  : if True, ensure 'chr' prefix on chrom names
      drop_na         : drop rows with NA in any of the three columns
      sort_by         : "genome" -> sort by chrom (natural order) then start;
                        "position" -> sort by start then end; None -> keep order
    """
    # --- 1) Auto-detect column names if needed
    def _guess(cols, pats):
        pats = [re.compile(p, re.I) for p in pats]
        for p in pats:
            for c in cols:
                if p.fullmatch(str(c)) or p.search(str(c)):
                    return c
        return None

    cols = list(df.columns)
    chrom_col = chrom_col or _guess(cols, [r"^chr(om|omosome)?$", r"^seq(name)?$", r"^chrom$", r"^contig$"])
    start_col = start_col or _guess(cols, [r"^start$", r"^chromStart$", r"^begin$", r"^pos(ition)?_?start?$"])
    end_col   = end_col   or _guess(cols, [r"^end$", r"^chromEnd$", r"^stop$", r"^pos(ition)?_?end?$"])

    if not all([chrom_col, start_col, end_col]):
        raise ValueError(
            f"Could not auto-detect columns. "
            f"Detected -> chrom:{chrom_col}, start:{start_col}, end:{end_col}. "
            f"Pass them explicitly."
        )

    # --- 2) Select and clean
    bed = df[[chrom_col, start_col, end_col]].copy()
    if drop_na:
        bed = bed.dropna()
    # coerce to correct types
    bed.iloc[:, 1] = pd.to_numeric(bed.iloc[:, 1], errors="coerce")
    bed.iloc[:, 2] = pd.to_numeric(bed.iloc[:, 2], errors="coerce")
    bed = bed.dropna(subset=[bed.columns[1], bed.columns[2]])
    bed[bed.columns[1]] = bed.iloc[:, 1].astype("int64")
    bed[bed.columns[2]] = bed.iloc[:, 2].astype("int64")

    # Ensure start < end and non-negative
    bed = bed[bed.iloc[:, 1] < bed.iloc[:, 2]]
    bed = bed[bed.iloc[:, 1] >= 0]

    # normalize chrom names
    if add_chr_prefix:
        bed.iloc[:, 0] = bed.iloc[:, 0].astype(str).map(lambda x: x if x.startswith("chr") else f"chr{x}")

    # --- 3) Sort (optional)
    if sort_by == "position":
        bed = bed.sort_values([bed.columns[1], bed.columns[2], bed.columns[0]]).reset_index(drop=True)
    elif sort_by == "genome":
        # natural chrom order chr1..chr19, chrX, chrY, chrM/chrMT, then others
        chroms = bed.iloc[:, 0].astype(str)
        def _key(c):
            m = re.match(r"^chr?(\d+)$", c)
            if m: return (0, int(m.group(1)))
            if c in {"chrX","X"}: return (1, 23)
            if c in {"chrY","Y"}: return (1, 24)
            if c in {"chrM","chrMT","M","MT"}: return (1, 25)
            return (2, c)
        order = chroms.map(_key)
        bed = bed.assign(_ord=order).sort_values(by=["_ord", bed.columns[1], bed.columns[2]]).drop(columns="_ord")

    # --- 4) Write BED (no header, tab-separated)
    bed.to_csv(out_path, sep="\t", header=False, index=False)

    return out_path
