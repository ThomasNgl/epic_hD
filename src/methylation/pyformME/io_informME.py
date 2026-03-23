
import os
import numpy as np
import pickle as pkl
import pandas as pd
import pysam
import glob
import shutil
import sys


def load_fasta(fasta_file_path):
    fa = pysam.FastaFile(fasta_file_path)
    return fa

def load_CpGs_chr(CpG_folder_path, chr_name):
    fname = f"CpGslocation_{chr_name}.csv"
    file_path = os.path.join(CpG_folder_path, fname)
    df_CpGs = pd.read_csv(file_path)
    return df_CpGs

def load_CpGs(CpG_folder_path):
    """
    Load all per-chromosome CpG location CSVs into a dictionary of DataFrames.
    Returns:
        results: dict {chr_name: df}
    """
    prefix = "CpGslocation_"
    results = {}
    for fname in os.listdir(CpG_folder_path):
        if fname.startswith(prefix) and fname.endswith(".csv"):
            chr_name = fname[len(prefix):-len(".csv")]
            file_path = os.path.join(CpG_folder_path, fname)
            df_CpGs = pd.read_csv(file_path)
            results[chr_name] = df_CpGs
    return results


def load_bam(bam_file_path):
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    return bam


def load_MEmat(MEmat_file_path):
    with open(MEmat_file_path, "rb") as f:
        df_MEmat = pkl.load(f)
    return df_MEmat

def store_rename_files(input_folder_path, 
                       pheno_path, 
                       pheno_label_list, 
                       suffix, 
                       output_folder_path):
    """
    input_folder_path: str, path to folder with files
    pheno_path: str, path to csv file
    pheno_label_list: list of str, column names in pheno.csv to build new path
    suffix: str, file extension/suffix (e.g., '.txt')
    output_folder_path: str, base path for output files
    """
    # Step 1: Find all files with the given suffix
    files = glob.glob(os.path.join(input_folder_path, f"*{suffix}"))
    old_file_paths = files
    file_names = [os.path.basename(f) for f in files]

    # Step 2: Load pheno (register) CSV
    pheno_df = pd.read_csv(pheno_path)
    print(output_folder_path)
    # Step 3: Get new file paths
    new_file_paths = []
    for filename in file_names:
        # Strip suffix for matching if needed
        base_name = os.path.splitext(filename)[0]
        # Find matching row
        row_match = pheno_df[pheno_df['Sample'].apply(lambda s: s in base_name)]
        if not row_match.empty:
            row = row_match.iloc[0]
            row_new_path_parts = [output_folder_path]
            for pheno_label in pheno_label_list:
                value = str(row[pheno_label]) if pd.notna(row[pheno_label]) else "None"
                row_new_path_parts.append(value)
        else:
            # If no match, use "None" for all pheno labels
            row_new_path_parts = [output_folder_path] + ["None"] * len(pheno_label_list)
        
        # Build new directory path and ensure it exists
        new_dir = os.path.join(*row_new_path_parts)
        os.makedirs(new_dir, exist_ok=True)
        # Compose the full file path
        new_file_path = os.path.join(new_dir, filename)
        new_file_paths.append(new_file_path)

    # Step 4: Build DataFrame
    df = pd.DataFrame({
        "file_name": file_names,
        "old_file_path": old_file_paths,
        "new_file_path": new_file_paths
    })

    # Step 5: Copy files and record success
    copied_list = []
    for idx, row in df.iterrows():
        try:
            shutil.copy2(row['old_file_path'], row['new_file_path'])
            copied_list.append(True)
            print("Copied ", row['new_file_path'])
        except Exception as e:
            print(f"Failed to copy {row['filename']}")
            copied_list.append(False)
    df['copied'] = copied_list

    return df


def get_df_files(base_folder_path, pheno_label_list, suffix):
    """
    Scan for files in base_folder_path at a depth corresponding to len(pheno_label_list),
    extract subfolder names as phenotypes, and return a DataFrame.
    """
    depth = len(pheno_label_list)
    # Construct a glob pattern: one * per folder level
    pattern = os.path.join(
        base_folder_path,
        *['*' for _ in range(depth)],
        f'*{suffix}'
    )
    files = glob.glob(pattern, recursive=False)
    
    records = []
    for file_path in files:
        rel_path = os.path.relpath(file_path, base_folder_path)
        parts = rel_path.split(os.sep)
        # Ensure there are enough parts (subfolders + filename)
        if len(parts) == depth + 1:
            record = {
                'file_name': parts[-1],
                'file_path': file_path,
            }
            for i, label in enumerate(pheno_label_list):
                record[label] = parts[i]
            records.append(record)
    df = pd.DataFrame(records)
    return df

     
def load_MEmats_pheno(base_folder_path,
                        pheno_label_list,
                        pheno_values_dict,
                        ):
    df_info = get_df_files(base_folder_path = base_folder_path, 
                                            pheno_label_list = pheno_label_list, 
                                            suffix = "MEmats.pkl")
    
    # Build a boolean mask for all conditions
    mask = np.ones(len(df_info), dtype=bool)
    for key, value in pheno_values_dict.items():
        mask &= df_info[key] == value

    # Apply the mask
    df_info_MEmats_pheno = df_info[mask].copy()

    df_MEmats = [pd.read_pickle(fp) for fp in df_info_MEmats_pheno['file_path']]

    df_info_MEmats_pheno['MEmats'] = df_MEmats

    return df_info_MEmats_pheno

def load_informME(informME_file_path):
    with open(informME_file_path, "rb") as f:
        df_informME = pkl.load(f)
    return df_informME

def load_pkl(file_path):
    """Load pickle safely handling numpy core/_core renaming issues."""
    try:
        with open(file_path, "rb") as f:
            return pkl.load(f)
    except ModuleNotFoundError as e:
        # Patch depending on which module is missing
        if "numpy.core.numeric" in str(e):
            sys.modules["numpy.core.numeric"] = np.core.numeric
        elif "numpy._core.numeric" in str(e):
            sys.modules["numpy._core.numeric"] = np.core.numeric
        else:
            raise  # re-raise if it's a different error

        # Retry once after patching
        with open(file_path, "rb") as f:
            return pkl.load(f)


def save_pkl(data, file_path):
    os.makedirs(os.path.dirname(file_path) or ".", exist_ok=True)

    with open(file_path, "wb") as f:
        pkl.dump(data, f)

    return 