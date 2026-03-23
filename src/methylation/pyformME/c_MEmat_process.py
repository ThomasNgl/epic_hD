import numpy as np
import pandas as pd
from scipy.sparse import vstack

def merge_MEmats(df_MEmats):
    # Expand MEmats into a big dataframe with columns 'region', 'MEmat'
    all_MEmat_df = pd.concat(df_MEmats['MEmats'].values, ignore_index=True)

    # Group by region and vertically stack the MEmat matrices
    records = []
    for region, group in all_MEmat_df.groupby('region'):
        matrices = list(group['MEmat'])  # or group['MEmat'] if that's your column name
        stacked = vstack(matrices)

        CpG_positions = np.array(group['position'].iloc[0])
        distance_arr = np.array(group['distance_to_next'].iloc[0])
        density_arr = np.array(group['density'].iloc[0])

        records.append({'region': region, 'MEmat': stacked,
                         "position": CpG_positions,
                        "distance_to_next": distance_arr,
                        "density": density_arr,})

    return pd.DataFrame(records)

def vect_MEmat(MEmat_dense):
    MEmat_contig_list = []
    CpGs_start = []
    CpGs_end = []

    for read in MEmat_dense:
        obs_indices = np.where(read > -1)[0]  # observed CpG indices
        if len(obs_indices) == 0:
            continue
        # Find contiguous blocks
        splits = np.split(obs_indices, np.where(np.diff(obs_indices) != 1)[0]+1)
        for block in splits:
            CpGs_start.append(block[0]) 
            CpGs_end.append(block[-1])
            MEmat_contig_list.append(read[block].tolist())

    return MEmat_contig_list, CpGs_start, CpGs_end


def sparse2dense_MEmat(MEmat_sparse):
    # Recreate observed matrix from sparse storage structure
    MEmat_dense = MEmat_sparse.toarray()
    MEmat_dense =- 1
    return MEmat_dense