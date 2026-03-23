import pandas as pd
from functools import reduce

def region_intersection(df_list):
    """
    Restrict all DataFrames to the intersection of regions.
    """
    if not df_list:
        return []
    common_regions = reduce(lambda x, y: x.intersection(y),
                            [set(df["region"]) for df in df_list])
    return [df[df["region"].isin(common_regions)].copy() for df in df_list]


def concat_subregions(df_list, col_candidates=("diffME_analysis", "ME_analysis")):
    out = []
    for df in df_list:
        col = next((c for c in col_candidates if c in df.columns), None)
        if col is None:
            raise KeyError(f"None of {col_candidates} found. Columns: {list(df.columns)}")
        out.append(pd.concat(df[col].tolist(), ignore_index=True))
    return out


def prep_features(df_list, feature):
    """
    Extract a given column from each DataFrame in df_list.

    Parameters
    ----------
    df_list : list of pd.DataFrame
        List of concatenated DataFrames (output of concat_subregions).
    feature : str
        Column name to extract.

    Returns
    -------
    list of pd.Series
        One Series per DataFrame containing the requested feature.
    """
    return [df[feature].reset_index(drop=True) for df in df_list]
