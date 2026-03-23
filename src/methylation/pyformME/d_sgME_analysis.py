import os
import pandas as pd
from . import io_informME
from . import a_refCpGs_process
from . import d_ME_ana_features
import re
import numpy as np

def sgME_analysis_region(dict_inforME, 
                         sub_region_length, 
                         max_CpGs = 18, 
                         thresh = [0, 0.25, 0.5, 0.75, 1],
                         threshold = 0.4,
                         epsilon = 0.01):
    
    dict_inforME_region = dict_inforME['informME']
    theta = dict_inforME_region['informME_params']['x']

    df_CpGs_informME_region = pd.DataFrame(dict_inforME_region['position'], columns=["position"])

    df_CpGs_informME_region['marginal_probs'] = dict_inforME_region['marginal_probs']
    transition_probs = dict_inforME_region['transition_probs']  # shape (n cpg -1, 2) 
    transition_probs_list = [list(row) for row in transition_probs]
    df_CpGs_informME_region['transition_probs'] = transition_probs_list

    distance_arr = dict_inforME_region['distance_arr']
    df_CpGs_informME_region['distance_arr'] = distance_arr
    density_arr = dict_inforME_region['density_arr']
    df_CpGs_informME_region['density_arr'] = density_arr

    # epsilon
    dict_inforME_epsilon_region = d_ME_ana_features.eval_informME_epsilon(distance_arr, density_arr, theta, epsilon,)
    df_CpGs_informME_region['marginal_eps_probs'] = dict_inforME_epsilon_region['marginal_probs']
    transition_probs = dict_inforME_epsilon_region['transition_probs']  
    transition_probs_list = [list(row) for row in transition_probs]
    df_CpGs_informME_region['transition_eps_probs'] = transition_probs_list

    # Group per subregion
    df_CpGs_informME_subregion = a_refCpGs_process.get_CpGs_per_region(
        df_CpGs_informME_region,
        region_length=sub_region_length,
        min_n_CpGs_per_region=1
    )

    ME_analysis_results_list = []
    for sub_region_lim, df_CpGs_informME_subregion in df_CpGs_informME_subregion:
        margprob_arr_subregion = df_CpGs_informME_subregion['marginal_probs'].tolist()
        transprob_arr_subregion = df_CpGs_informME_subregion['transition_probs'].tolist()
        margprob_arr_subregion = np.array(margprob_arr_subregion)
        transprob_arr_subregion = np.array(transprob_arr_subregion)
    
        CpG_pos_subregion =  df_CpGs_informME_subregion['position'].tolist()
        CpG_pos_subregion = np.array(CpG_pos_subregion)

        ml_distribution, mml, nme, meth_class, var_class, me_class = d_ME_ana_features.get_mml_nme(
            margprob_arr=margprob_arr_subregion, 
            transprob_ME_arr=transprob_arr_subregion, 
            max_CpGs=max_CpGs, 
            thresh=thresh,
            threshold=threshold
        )

        margprob_eps_arr_subregion = df_CpGs_informME_subregion['marginal_eps_probs'].tolist()
        transprob_eps_arr_subregion = df_CpGs_informME_subregion['transition_eps_probs'].tolist()
        margprob_eps_arr_subregion = np.array(margprob_eps_arr_subregion)
        transprob_eps_arr_subregion = np.array(transprob_eps_arr_subregion)

        ml_distribution_epsilon, esi, msi = d_ME_ana_features.get_esi_msi(
            epsilon = epsilon,
            margprob_eps_arr = margprob_eps_arr_subregion,
            transprob_eps_arr = transprob_eps_arr_subregion,
            ml_distribution = ml_distribution,
            max_CpGs=max_CpGs,
        )
        turn, cap, rde = d_ME_ana_features.get_turn_cap_rde(margprob_arr_subregion)

        # Flattened result:
        result_row = {
            "subregion": sub_region_lim,
            "position": CpG_pos_subregion,
            "ml_distribution": ml_distribution,
            "mml": mml,
            "nme": nme,
            "meth_class": meth_class,
            "var_class": var_class,
            "me_class": me_class,
            "informME_epsilon": dict_inforME_epsilon_region,
            "ml_distribution_epsilon": ml_distribution_epsilon,
            "esi": esi,
            "msi": msi,
            "turn": turn,
            "cap": cap,
            "rde": rde
        }
        ME_analysis_results_list.append(result_row)

    df_ME_analysis = pd.DataFrame(ME_analysis_results_list)
    return dict_inforME['region'], df_ME_analysis

def sgME_analysis_chr(chr_name,
                      sub_region_length,
                      informME_file_path,
                      pheno_label_value_dict,
                      save_path,
                      max_CpGs = 18, 
                      thresh = [0, 0.25, 0.5, 0.75, 1],
                      threshold = 0.4,
                      epsilon = 0.01):

    df_informME = io_informME.load_informME(informME_file_path)
    ME_analysis_chr = []
    for i, (_, row_informME) in enumerate(df_informME.iterrows(), start=1):
        # your processing here ...

        if i == 100 or i == 500 or i == 1000 or i % 1000 == 0:
            print(f"Processed {i} rows")
        region = row_informME['region']
        dict_inforME_region = row_informME['informME']
        if dict_inforME_region['informME_params'] is None:
            continue
        region, df_ME_analysis_region = sgME_analysis_region(dict_inforME_region, 
                                                        sub_region_length, 
                                                        max_CpGs = max_CpGs, 
                                                        thresh = thresh,
                                                        threshold = threshold,
                                                        epsilon = epsilon)
        ME_analysis_chr.append({"region": region, "ME_analysis": df_ME_analysis_region})

    df_ME_analysis_chr = pd.DataFrame(ME_analysis_chr, columns=["region", "ME_analysis"])

    # Build save path using phenotypes
    pheno_vals = [str(pheno_label_value_dict[label]) for label in list(pheno_label_value_dict.keys())]
    # Safe join for folder structure
    subdir = os.path.join(save_path, *pheno_vals)
    os.makedirs(subdir, exist_ok=True)
    # Build filename: group_condition_replicate_MEmats.pkl (adapt as needed)
    fname = f"{'_'.join(pheno_vals)}_ref_{chr_name}_ME_analysis.pkl"
    out_path = os.path.join(subdir, fname)
    # Save
    df_ME_analysis_chr.to_pickle(out_path)
    print(f"Saved ME analysis results in {out_path}")
    return df_ME_analysis_chr


