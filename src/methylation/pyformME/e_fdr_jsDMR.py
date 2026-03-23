import os
import pickle as pkl
import pandas as pd
import re
from io_utils import io_bed, utils
from pyformME import io_informME

def get_significant_diffME_df(fdr_jsd_sites_df, diffME_df):
    # --- 1) Parse region into numeric start/end
    r = diffME_df['region'].str.extract(r'^chr[^:]+:(\d+)-(\d+)$').astype(int)
    diffME_df = diffME_df.assign(reg_start=r[0], reg_end=r[1])

    # --- 2) Build an IntervalIndex over the windows
    # Build IntervalIndex for [start, end)
    iv = pd.IntervalIndex.from_arrays(
        diffME_df['reg_start'], diffME_df['reg_end'], closed='left'
    )

    # Map each fdr_jsd_sites_df interval [s,e) into bins: use e-1 so both ends fall in same bin
    start_bins = iv.get_indexer(fdr_jsd_sites_df['start'])
    end_bins   = iv.get_indexer((fdr_jsd_sites_df['end'] - 1).clip(lower=fdr_jsd_sites_df['start']))

    mask_same_bin = (start_bins == end_bins) & (start_bins != -1)
    covering_idx  = pd.unique(start_bins[mask_same_bin])

    fdr_diffME_df= diffME_df.iloc[covering_idx].copy()
    return fdr_diffME_df


def get_significant_diffME_subregion_df(fdr_jsd_sites_df: pd.DataFrame,
                                        fdr_diffME_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each row in fdr_diffME_df, take the per-row DataFrame in `diffME_analysis`,
    parse its `subregion` into [sub_start, sub_end], and inner-join with the set
    of [start, end] pairs from fdr_jsd_sites_df. Concatenate all matches and
    add a `chr` column derived from the row's `region` string (e.g., 'chr15', 'chrX').

    Returns: DataFrame with the same columns as each `diffME_analysis` plus `chr`.
    """

    def parse_pair(x):
        # Accept [start, end], (start, end), or strings like "[141, 498]"
        if isinstance(x, (list, tuple)) and len(x) >= 2:
            return int(x[0]), int(x[1])
        m = re.search(r'(\d+)\D+(\d+)', str(x))
        return (int(m.group(1)), int(m.group(2))) if m else (None, None)

    # Build the key set once
    keys = (fdr_jsd_sites_df[['start', 'end']]
            .dropna()
            .astype({'start': 'int64', 'end': 'int64'})
            .drop_duplicates()
            .rename(columns={'start': 'sub_start', 'end': 'sub_end'}))

    out_frames = []

    for _, row in fdr_diffME_df.iterrows():
        region = row.get('region', None)

        # extract 'chr...' from the region string (e.g. "chr15:123-456" -> "chr15")
        chrom = None
        if pd.notna(region):
            m = re.match(r'^(chr[^:]+):', str(region))
            if m:
                chrom = m.group(1)

        region_res = row.get('diffME_analysis', None)

        # we expect region_res to be a DataFrame with a 'subregion' column
        if not isinstance(region_res, pd.DataFrame) or 'subregion' not in region_res.columns:
            continue

        tmp = region_res.copy()

        # Parse subregion into integers
        sub_pairs = tmp['subregion'].map(parse_pair)
        tmp[['sub_start', 'sub_end']] = pd.DataFrame(sub_pairs.tolist(), index=tmp.index)

        # Only keep rows with valid ints for the join
        tmp = tmp.dropna(subset=['sub_start', 'sub_end']).astype({'sub_start': 'int64', 'sub_end': 'int64'})

        # Join with keys and keep only the original region_res columns
        matched = tmp.merge(keys, on=['sub_start', 'sub_end'], how='inner')
        matched = matched[region_res.columns].copy()

        # Add chromosome column
        matched['chr'] = chrom

        out_frames.append(matched)

    if out_frames:
        return pd.concat(out_frames, ignore_index=True)

    # Fallback empty frame with expected columns if nothing matched
    # (keep the union of possible columns where available)
    try:
        sample_cols = list(fdr_diffME_df.iloc[0]['diffME_analysis'].columns)
    except Exception:
        sample_cols = []
    return pd.DataFrame(columns=sample_cols + ['chr'])


def save_significant_diffME_subregions_chr(jsd_sites_path, 
                                           diffME_path, 
                                           threshold, 
                                           chr, 
                                           save_dir,
                                           fdr_jsd_sites_name = 'fdr_jsd_sites',
                                           fdr_diffME_name = 'fdr_diffME'):
    jsd_sites_df = io_informME.load_pkl(jsd_sites_path)
    # Get significant
    fdr_jsd_sites_df = jsd_sites_df.loc[jsd_sites_df['qval'] < threshold].copy()
    fdr_jsd_sites_df = fdr_jsd_sites_df.sort_values('start', ascending=True)
    n_dmrs = len(fdr_jsd_sites_df)
    print(n_dmrs, " regions significant in ", chr)
    res_matched_df = None
    if n_dmrs > 0:
        out_path =  os.path.join(save_dir, f"{chr}/{chr}_{fdr_jsd_sites_name}.bed")

        io_bed.df_to_bed(fdr_jsd_sites_df, out_path = out_path) 
        diffME_df = io_informME.load_pkl(diffME_path)

        fdr_diffME_df = get_significant_diffME_df(fdr_jsd_sites_df, diffME_df)
        res_matched_df = get_significant_diffME_subregion_df(fdr_jsd_sites_df, fdr_diffME_df)
        out_path =  os.path.join(save_dir, f"{chr}/{chr}_{fdr_diffME_name}.pkl")
        io_informME.save_pkl(res_matched_df,out_path)
    return res_matched_df, fdr_jsd_sites_df


# def save_significant_diffME_subregions(base_dir, threshold, save_dir):
#     chr_list = utils.get_chr_names()
#     for chr in chr_list:
#         chr_jsd_sites_path =  os.path.join(base_dir, "06_jsDMR", f"{chr}/{chr}_jsd_sites.pkl")
#         chr_diffME_path =  os.path.join(base_dir,"05_diff_analysed_informME", f"{chr}/{chr}_oocyte_diffME_analysis.pkl")

#         save_significant_diffME_subregions_chr(jsd_sites_path = chr_jsd_sites_path,
#                                                 diffME_path = chr_diffME_path, 
#                                                 threshold = threshold,
#                                                 save_dir = save_dir)
#     return