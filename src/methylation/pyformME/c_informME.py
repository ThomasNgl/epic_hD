import os
import numpy as np
import scipy
import pandas as pd
import re
from scipy.optimize import differential_evolution
from scipy.optimize import dual_annealing, minimize
from scipy.optimize import shgo
import secrets

from . import io_informME
from . import a_refCpGs_process
from . import c_ising_model, c_ising_model_numpy, c_ising_model_numba
from . import c_MEmat_process

# def fit_informME_region(MEmat_sparse, distance_arr, density_arr, theta_bounds, precision, optimizer = None):
#     if isinstance(MEmat_sparse, scipy.sparse.csr_matrix):
#         MEmat_dense = MEmat_sparse.toarray()-1
#     else:  # If already dense
#         MEmat_dense = MEmat_sparse

#     _, n_CpGs = MEmat_dense.shape  # Number of columns
#     # Calculate the percent of CpGs covered at least by one read
#     mask = MEmat_dense > -1  # Boolean matrix
#     n_covered_columns = np.sum(np.sum(mask, axis=0) > 0)
#     percent_covered_columns = n_covered_columns / n_CpGs
#     # Calculate the average number of reads per CpG
#     n_CpG_info = np.sum(mask)
#     depth_cov = n_CpG_info / n_CpGs
#     # Check data sufficiency
#     if (depth_cov < 2.5) or (percent_covered_columns < (2/3)):
#         # Insufficient data - do not build statistical model.
#         fitted_ising = None
#         return fitted_ising

#     # Filter out reads wit full -1
#     MEmat_dense = MEmat_dense[np.sum(MEmat_dense > -1, axis=1) > 0, :]

#     # Split reads into contiguous CpG reads
#     MEmat_contig_list, CpGs_start, CpGs_end = c_MEmat_process.vect_MEmat(MEmat_dense)
#     def my_callback(xk, convergence=None):
#         print(f"Current best parameters: {xk}")
#         if convergence is not None:
#             print(f"Convergence metric: {convergence}")

#     if optimizer is None:
#         fitted_ising = differential_evolution(func = c_ising_model.get_ising_func2minimize,
#                                     bounds = theta_bounds,
#                                     args = (MEmat_contig_list, CpGs_start, CpGs_end, distance_arr, density_arr, precision),
#                                         strategy='best1bin',    # Standard, robust strategy
#                                         maxiter=1000,           # You can adjust for more thorough search
#                                         popsize=15,             # Default is fine, can increase for more exploration
#                                         tol=0.01,               # Convergence tolerance
#                                         polish=True,            # Refine result with local search
#                                         # disp=True,               # Prints progress
#                                             # callback=my_callback  # <--- Here!

#                                     )     

#     elif optimizer == "da":
#         # Global search (simulated annealing style)
#         fitted_ising = dual_annealing(
#             func=c_ising_model.get_ising_func2minimize,
#             bounds=theta_bounds,
#             args=(MEmat_contig_list, CpGs_start, CpGs_end, distance_arr, density_arr, precision),
#             maxiter=1000,      # keep similar budget as DE for fair comparison
#             seed=42            # reproducibility
#         )

#     elif optimizer == "shgo":
#         fitted_ising = shgo(
#             c_ising_model.get_ising_func2minimize,
#             bounds=theta_bounds,
#             args=(MEmat_contig_list, CpGs_start, CpGs_end, distance_arr, density_arr, precision),
#             sampling_method="sobol", iters=3
#         )

#     return fitted_ising


# # --- put this at module import time so child processes inherit ---
# os.environ.setdefault("OMP_NUM_THREADS", "1")
# os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
# os.environ.setdefault("MKL_NUM_THREADS", "1")
import time
def fit_informME_region(
    MEmat_sparse, distance_arr, density_arr, theta_bounds, precision,
    optimizer=None, do_polish=True, seed=None, is_numpy=True,
):
    if seed is None:
        seed = secrets.randbits(32)
    # Convert (keep this top-level for pickling)
    # Sparse matrices are stored with +1 applied so that NA values are 0 and not -1
    if isinstance(MEmat_sparse, scipy.sparse.csr_matrix):
        MEmat_dense = MEmat_sparse.toarray() - 1
    else:
        MEmat_dense = MEmat_sparse

    _, n_CpGs = MEmat_dense.shape
    mask = MEmat_dense > -1
    n_covered_columns = np.sum(np.sum(mask, axis=0) > 0)
    percent_covered_columns = n_covered_columns / n_CpGs
    n_CpG_info = np.sum(mask)
    depth_cov = n_CpG_info / n_CpGs
    if (depth_cov < 2.5) or (percent_covered_columns < (2/3)):
        return None  # insufficient data

    # Filter out reads with full -1
    MEmat_dense = MEmat_dense[np.sum(MEmat_dense > -1, axis=1) > 0, :]

    # Split reads into contiguous CpG reads
    MEmat_contig_list, CpGs_start, CpGs_end = c_MEmat_process.vect_MEmat(MEmat_dense)
    if is_numpy:
        f = c_ising_model_numpy.get_ising_func2minimize
    else:
        f = c_ising_model_numba.get_ising_func2minimize

    common_args = (MEmat_contig_list, CpGs_start, CpGs_end, distance_arr, density_arr, precision)

    # --- Choose optimizer (all serial inside the process) ---
    if optimizer is None or optimizer == "de":
        res = differential_evolution(
            func=f, bounds=theta_bounds, args=common_args,
            strategy='best1bin', maxiter=1000, popsize=15, tol=0.01,
            polish=False, seed=seed  # serial: no workers/updating
        )
        if do_polish:
            res = minimize(f, res.x, args=common_args, method="L-BFGS-B", bounds=theta_bounds)

    elif optimizer == "da":
        res = dual_annealing(
            func=f, bounds=theta_bounds, args=common_args,
            maxiter=1000, seed=seed
        )
        if do_polish:
            res = minimize(f, res.x, args=common_args, method="L-BFGS-B", bounds=theta_bounds)

    elif optimizer == "shgo":
        res = shgo(
            f, bounds=theta_bounds, args=common_args,
            sampling_method="sobol", iters=3,
            minimizer_kwargs=dict(method="L-BFGS-B")
        )
    else:
        raise ValueError("optimizer must be None/'de', 'da', or 'shgo'")

    return res



def eval_informME(distance_arr, 
                  density_arr, 
                  theta,
                  is_numpy = True,
                    ):
    if is_numpy:
        Ar, Cr = c_ising_model_numpy.get_ising_coef(distance_arr, density_arr, theta,
                                            )
        logZ1, logZ0, logZ = c_ising_model_numpy.get_logZ_np(Ar, Cr)
        logZ1tilde, logZ0tilde, _ = c_ising_model_numpy.get_logZtilde_np(Ar, Cr)
        # Compute transitions probabilities
        transprob_arr = c_ising_model_numpy.get_transprob_arr_np(logZ1=logZ1, 
                                                        logZ0=logZ0, 
                                                        Ar = Ar, 
                                                        Cr = Cr)
        # compute marginal probabilities
        margprob_arr = c_ising_model_numpy.get_margprob_arr_np(logZ1=logZ1, 
                                                        logZ0=logZ0, 
                                                        logZ = logZ,
                                                        logZ1tilde = logZ1tilde,
                                                        logZ0tilde = logZ0tilde,
                                                        )
    else:
        Ar, Cr = c_ising_model_numba.get_ising_coef(distance_arr, density_arr, theta,
                                            )
        logZ1, logZ0, logZ = c_ising_model_numba.get_logZ_np(Ar, Cr)
        logZ1tilde, logZ0tilde, _ = c_ising_model_numba.get_logZtilde_np(Ar, Cr)
        # Compute transitions probabilities
        transprob_arr = c_ising_model_numba.get_transprob_arr_np(logZ1=logZ1, 
                                                        logZ0=logZ0, 
                                                        Ar = Ar, 
                                                        Cr = Cr)
        # compute marginal probabilities
        margprob_arr = c_ising_model_numba.get_margprob_arr_np(logZ1=logZ1, 
                                                        logZ0=logZ0, 
                                                        logZ = logZ,
                                                        logZ1tilde = logZ1tilde,
                                                        logZ0tilde = logZ0tilde,
                                                        )

    # Check p = 0 or 1 because potential issue with log
    eps = np.finfo(float).eps
    margprob_arr = np.clip(margprob_arr, eps, 1 - eps)
    transprob_arr = np.clip(transprob_arr, eps, 1 - eps)
    # Return as dict
    result = {
        "marginal_probs": margprob_arr,
        "transition_probs": transprob_arr,
        "logZ1": np.array(logZ1),
        "logZ0": np.array(logZ0),
        "logZ": float(logZ),
        "logZ1tilde": np.array(logZ1tilde),
        "logZ0tilde": np.array(logZ0tilde),
        "Ar": Ar,
        "Cr": Cr,
    }
    return result


def informME_chr(CpGs_ref_folder_path,
                 chr_name,
                 region_length,
                 min_n_CpGs_per_region,
                 MEmats_folder_path,
                 pheno_label_value_dict,
                 optimizer,
                 theta_bounds,
                 precision,
                 save_path):
    # 1. Load reference CpG positions for the given chromosome
    df_CpGs = io_informME.load_CpGs_chr(CpGs_ref_folder_path, chr_name)
    df_CpGs_filtered_per_region = a_refCpGs_process.get_CpGs_per_region(df_CpGs=df_CpGs,
                                                                    region_length=region_length,
                                                                    min_n_CpGs_per_region=min_n_CpGs_per_region
                                                                    )
    
    # 2. Load and merge MEmats
    df_MEmats_pheno_chr  = io_informME.load_MEmats_pheno(base_folder_path = MEmats_folder_path,
                                                        pheno_label_list = list(pheno_label_value_dict.keys()),
                                                        pheno_values_dict = pheno_label_value_dict,
                                                        )
    df_MEmat_merged_per_region = c_MEmat_process.merge_MEmats(df_MEmats = df_MEmats_pheno_chr)

    informME_result_list = []
    for _, MEmat_region in df_MEmat_merged_per_region.iterrows():
        # ME mat information
        region = MEmat_region['region']
        MEmat_sparse = MEmat_region['MEmat']
        # Ref information
        # Extract start and end
        match = re.search(r'(\d+)-(\d+)', region)
        start, end = map(int, match.groups())
        # Create interval
        region_interval = pd.Interval(left=start, right=end, closed='left')
        df_CpGs_region = df_CpGs_filtered_per_region.get_group(region_interval)
        CpG_positions = np.array(df_CpGs_region['position'])
        distance_arr = np.array(df_CpGs_region['distance_to_next'])
        density_arr = np.array(df_CpGs_region['density'])
        print("Train on:", region)
        # Fit the parameters
        informME_params_region = fit_informME_region(
            MEmat_sparse=MEmat_sparse, 
            distance_arr=distance_arr,
            density_arr=density_arr,
            optimizer=optimizer, 
            theta_bounds=theta_bounds,
            precision = precision
        )

        print("\nFinished Ising model parameters fitting\n")
        if informME_params_region is not None:
            best_theta = informME_params_region['x']
            # Use the fitted parameters to evaluate the model
            inforME_region_result_dict = eval_informME(
                distance_arr=distance_arr, 
                density_arr=density_arr, 
                theta=best_theta,
            )
        else:
            continue
        inforME_region_result_dict["distance_arr"] = distance_arr
        inforME_region_result_dict["density_arr"] = density_arr
        inforME_region_result_dict["position"] = CpG_positions
        inforME_region_result_dict["informME_params"] = informME_params_region

        # Append result as a tuple (region, informME_region_result_dict)
        informME_result_list.append({"region": region, "informME": inforME_region_result_dict})

    # Convert to DataFrame
    df_informME = pd.DataFrame(informME_result_list, columns=["region", "informME"])

    # Build save path using phenotypes
    pheno_vals = [str(pheno_label_value_dict[label]) for label in list(pheno_label_value_dict.keys())]
    # Safe join for folder structure
    subdir = os.path.join(save_path, *pheno_vals)
    os.makedirs(subdir, exist_ok=True)
    # Build filename: group_condition_replicate_MEmats.pkl (adapt as needed)
    fname = f"{'_'.join(pheno_vals)}_ref_{chr_name}_informME.pkl"
    out_path = os.path.join(subdir, fname)
    # Save

    df_informME.to_pickle(out_path)
    print(f"Saved fitting InformME results in {out_path}")
    
    me_merged_name = f"{'_'.join(pheno_vals)}_ref_{chr_name}_MEmats.pkl"
    me_merged_path = os.path.join(subdir, me_merged_name)
    df_MEmat_merged_per_region.to_pickle(me_merged_path)
    print(f"Merged ME matrices saved in {me_merged_path}")

    return df_informME


# Helper function to process a single region
def process_region(MEmat_region_dict, theta_bounds, precision, optimizer, is_numpy = True):
    region = MEmat_region_dict['region']
    MEmat_sparse = MEmat_region_dict['MEmat']
    CpG_positions = np.array(MEmat_region_dict['position'])
    distance_arr = np.array(MEmat_region_dict['distance_to_next'])
    density_arr = np.array(MEmat_region_dict['density'])

    informME_params_region = fit_informME_region(
        MEmat_sparse=MEmat_sparse, 
        distance_arr=distance_arr,
        density_arr=density_arr,
        optimizer=optimizer, 
        theta_bounds=theta_bounds,
        precision=precision,
        do_polish=True, 
        seed=None,
        is_numpy= is_numpy
    )

    if informME_params_region is not None:
        best_theta = informME_params_region['x']
        inforME_region_result_dict = eval_informME(
            distance_arr=distance_arr, 
            density_arr=density_arr, 
            theta=best_theta,
            is_numpy= is_numpy

        )
        inforME_region_result_dict["distance_arr"] = distance_arr
        inforME_region_result_dict["density_arr"] = density_arr
        inforME_region_result_dict["position"] = CpG_positions
        inforME_region_result_dict["informME_params"] = informME_params_region
        return {"region": region, "informME": inforME_region_result_dict}
    else:
        return None