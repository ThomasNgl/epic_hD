import numpy as np
import pandas as pd

def get_Cs_chr(fasta_file, chr_name):
    seq = str(fasta_file[chr_name][:])

    positions = []
    next_bases = []
    
    for i in range(len(seq) - 1):  # Avoid indexing beyond last base
        if seq[i].upper() == "C":
            positions.append(i + 1)  # 1-based position
            next_bases.append(seq[i + 1].upper())


    df = pd.DataFrame({
        "position": positions,
        "next_base": next_bases,
    })
    return df


def get_CpGs_chr(fasta_file, chr_name):
    df_Cs = get_Cs_chr(fasta_file, chr_name)
    df_CpGs = df_Cs[df_Cs['next_base'] == 'G'] 
    df_CpGs.reset_index(drop=True, inplace=True)
    return df_CpGs


def get_CpGs_distance(df_CpGs):
    # Distance to next CpG (last one is NaN or a big number, as in MATLAB)
    pos = df_CpGs['position'].values
    dist = np.diff(pos, append=1e9)  
    df_CpGs['distance_to_next'] = dist
    return df_CpGs


def get_CpGs_density(df_CpGs, window=1000):
    positions = df_CpGs['position'].values
    positions.sort()  # just in case
    half_w = window // 2

    # Compute window boundaries for each CpG
    left_bounds = positions - half_w
    right_bounds = positions + half_w

    # For each CpG, count how many CpGs are within the window using searchsorted
    left_idcs = np.searchsorted(positions, left_bounds, side='left')
    right_idcs = np.searchsorted(positions, right_bounds, side='right')
    CpGs_counts = right_idcs - left_idcs

    # Optionally, normalize by window size
    CpGs_densities = CpGs_counts / window

    df_CpGs = df_CpGs.copy()
    df_CpGs['density'] = CpGs_densities
    return df_CpGs

def _get_chr_regions(df_CpGs, region_length):
    # Returns a list of boundaries. 
    # The lower is included in the interval and the upper is not.
    # The last value of the list is the upper boundary of the last interval.
    # This values is the position of the bp after the last G.
    first_CpG_pos = df_CpGs.iloc[0]['position']
    final_CpG_pos = df_CpGs.iloc[-1]['position']

    region_lims_pos = np.arange(first_CpG_pos, final_CpG_pos + 1, region_length)
    # To include the last CpG with the complement GpC
    region_lims_pos = np.append(region_lims_pos, final_CpG_pos+2)
    return region_lims_pos 


def get_CpGs_per_region(df_CpGs,
                        region_length,
                        min_n_CpGs_per_region):
    region_lims_pos = _get_chr_regions(df_CpGs = df_CpGs, 
                                    region_length = region_length)
    
    # right=False: As said in get_chr_regions the lower value is included and could be a C position.
    # And the upper is not included in the interval as it is the included lower boundary in the next interval 
    df_CpGs = df_CpGs.copy()  # Avoid modifying original DataFrame
    df_CpGs['region'] = pd.cut(df_CpGs['position'], bins=region_lims_pos, right=False)
    
    # Get number of CpGs per region
    n_CpGs_per_region = df_CpGs['region'].value_counts().sort_index()
    regions_to_keep = n_CpGs_per_region[n_CpGs_per_region >= min_n_CpGs_per_region].index

    # Filter DataFrame to only those regions
    df_CpGs_filtered = df_CpGs[df_CpGs['region'].isin(regions_to_keep)].copy()

    # Return grouped DataFrame

    return df_CpGs_filtered.groupby('region', observed=True)
