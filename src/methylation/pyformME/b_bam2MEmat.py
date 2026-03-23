import os
import numpy as np
import pandas as pd
import pysam
from . import utils_informME
from . import a_refCpGs_process
from . import b_bam_process
from . import io_informME
from scipy.sparse import lil_matrix, csr_matrix

def reads2MEmat(bam_reads,
                 df_CpGs_region,
                 min_base_q = 1,
                 min_phred_score = 1):
    n_reads = len(bam_reads)
    n_CpGs = len(df_CpGs_region)
    MEmat = -1 * np.ones((n_reads, n_CpGs))

    # Build lookup: CpG position (1-based) --> column in MEmat
    positions = list(df_CpGs_region['position'])
    position_to_column = {pos: idx for idx, pos in enumerate(positions)}

    n_has_idnship = 0
    n_no_xm = 0
    n_low_mq = 0
    for read_index, read in enumerate(bam_reads):  # reads = iterable/list of pysam AlignedSegment
        if read.mapping_quality <= min_phred_score:
            n_low_mq += 1
            continue
        has_idnshp = bool(set('IDNSHP') & set(read.cigarstring))
        if has_idnshp:
            n_has_idnship+=1
            continue
        try:
            XM = read.get_tag('XM')
        except KeyError:
            n_no_xm += 1
            continue  # No XM tag, skip this read
        ref_pos = read.get_reference_positions(full_length=False)  # list of 0-based positions
        qualities = read.query_qualities

        try:
            GA_strand = read.get_tag('XG') == "GA"
        except KeyError:
            GA_strand = False  # or handle as needed
        
        # Iterate through a read sequence
        for meth_call, pos0based, base_q in zip(XM, ref_pos, qualities):
            # Check quality of the base
            if base_q < min_base_q:
                continue
            if meth_call in ("Z", "z"):
                pos1based = pos0based + 1  # convert to 1-based
                pos1based -= 1 if GA_strand else 0 # If reverse strand, the C is one bp above the forward ref CpG
                col = position_to_column.get(pos1based)
                if col is not None:
                    MEmat[read_index, col] = 1 if meth_call == "Z" else 0
                
    # Find rows where not all elements are -1
    mask = ~(MEmat == -1).all(axis=1)

    # Filter out those rows
    MEmat_filtered = MEmat[mask]

    # Shift values so: -1 → 0, 0 → 1, 1 → 2
    MEmat_plus1 = MEmat_filtered + 1

    # Convert to sparse matrix
    ME_sparse = csr_matrix(MEmat_plus1)
    return ME_sparse, {"n_low_mq":n_low_mq, "n_has_idnship":n_has_idnship, "n_no_xm":n_no_xm}

def reads2MEmat2(bam_reads,
                 df_CpGs_region):
    n_reads = len(bam_reads)
    n_CpGs = len(df_CpGs_region)
    MEmat = lil_matrix((n_reads, n_CpGs), dtype=int)  

    # Build lookup: CpG position (1-based) --> column in MEmat
    # CpGs_positions = np.array(df_CpGs_region['position'], dtype=int) - 1 #from 1-based to 0-based
    CpGs_positions = list(df_CpGs_region['position'])
    for read_index, read in enumerate(bam_reads):  # reads = iterable/list of pysam AlignedSegment
        # Check matching
        has_idnshp = bool(set('IDNSHP') & set(read.cigarstring))
        if has_idnshp:
            continue
        # Check methylation data
        try:
            xm_tag = read.get_tag('XM')
        except KeyError:
            print(f"No XM tag for methylation info in read {read}")
            continue  # No XM tag, skip this read
        
        ref_0_based_positions = read.get_reference_positions(full_length=False)  # Only aligned positions

        # Map reference positions to their methylation state
        refpos_to_meth = dict(zip(ref_0_based_positions, xm_tag))
        start_read, end_read = ref_0_based_positions[0], ref_0_based_positions[-1]
        for col, CpG_pos in enumerate(CpGs_positions):
            if CpG_pos < start_read:
                continue
            if end_read < CpG_pos:
                break
            meth_state = refpos_to_meth.get(CpG_pos-1, None) # From 1-base to 0-based
            if meth_state == "Z":
                MEmat[read_index, col] = 2
            elif meth_state == "z":
                MEmat[read_index, col] = 1
            # else leave as -1


    # Convert to sparse matrix
    MEmat_csr = MEmat.tocsr()
    # Find rows where not all elements are -1
    mask = MEmat_csr.getnnz(axis=1) > 0
    ME_sparse_filtered = MEmat_csr[mask, :]
    return ME_sparse_filtered


def bam2MEmats(bam_file_path,
                paired_ends,
                df_CpGs_per_region,
                chr_name,
                min_base_q = 1,
                min_phred_score = 1):
    # Get the reads of chr_name and create a matrix for each region in the chr_name
    bam = io_informME.load_bam(bam_file_path)
    bam_chromosome_names = bam.references

    chr_name, tried_chr_names = utils_informME.check_chr_name(chr_name = chr_name,
                    bam_chromosome_names = bam_chromosome_names)
    if chr_name is None:
        print(f"No reads from chromosome {chr_name} (tried: {tried_chr_names}) in the BAM file {bam_file_path}")
        bam.close()
        return pd.DataFrame(columns=["region", "MEmat", "dup_qc", "MEmat_qc", "position", "distance_to_next", "density"])

    results = []
    # Loop through your regions
    for region_lim, df_CpGs_region in df_CpGs_per_region:
        start = int(region_lim.left)
        end = int(region_lim.right)
        region_str = f"{chr_name}:{start}-{end}"

        # PySAM fetch: note end is exclusive
        # -1 because fetch is 0 based
        reads = bam.fetch(chr_name, start-1, end)
        reads = list(reads)
        if len(reads) != 0:
            sg_reads, dup_qc = b_bam_process.deduplicate_reads(reads)
        else:
            # print(f"Error in pySAMtools read at {region_str}:")
            # print(f"No reads found.\n")
            MEmat = None
            dup_qc = None
            continue

        n_reads = len(sg_reads)
        if n_reads == 0 or n_reads is None:
            # print(f"All reads were filtered out in {region_str}")
            MEmat = None
            continue

        if n_reads >= 5000:
            # print(f"Too many reads mapped to region {region_str}")
            MEmat = None
            continue

        else:
            MEmat, MEmat_qc = reads2MEmat(bam_reads=sg_reads, 
                                   df_CpGs_region=df_CpGs_region,
                                   min_phred_score = min_phred_score,
                                   min_base_q = min_base_q)
            if MEmat.nnz == 0:
                MEmat = None
                MEmat_qc = None
                continue

            CpG_positions = np.array(df_CpGs_region['position'])
            distance_arr = np.array(df_CpGs_region['distance_to_next'])
            density_arr = np.array(df_CpGs_region['density'])
            results.append({
                "region": region_str,
                "MEmat": MEmat,
                "dup_qc": dup_qc,
                "MEmat_qc": MEmat_qc,
                "position": CpG_positions,
                "distance_to_next": distance_arr,
                "density": density_arr,
            })
    
    bam.close()

    # Convert results to DataFrame
    df_MEmats = pd.DataFrame(results, columns=["region", "MEmat", "dup_qc", "MEmat_qc", "position", "distance_to_next", "density"])
    num_rows_no_none_MEmat = df_MEmats["MEmat"].notnull().sum()

    print(f"{num_rows_no_none_MEmat} regions processed over {len(df_CpGs_per_region)}")
    return df_MEmats


def get_MEmats_chr(
    bam_folder_path,
    pheno_label_list,
    suffix,
    paired_ends,
    min_phred_score,
    CpGs_ref_folder_path,
    chr_name,
    region_length,
    save_path,
    min_n_CpGs_per_region=10,
):
    # 1. Load reference CpG positions for the given chromosome
    df_CpGs = io_informME.load_CpGs_chr(CpGs_ref_folder_path, chr_name)
    df_CpGs_filtered_per_region = a_refCpGs_process.get_CpGs_per_region(
        df_CpGs=df_CpGs,
        region_length=region_length,
        min_n_CpGs_per_region=min_n_CpGs_per_region
    )

    # 2. Get BAM files metadata as DataFrame
    df_bam = io_informME.get_df_files(
        base_folder_path=bam_folder_path,
        pheno_label_list=pheno_label_list,
        suffix=suffix
    )

    # 3. Compute ME mats and save
    for _, row in df_bam.iterrows():
        bam_file_path = row['file_path']
        bam_file_name = row['file_name']
        bam_file_name_wo_ext = os.path.splitext(os.path.basename(bam_file_name))[0]

        # Compute ME matrix
        df_MEmats = bam2MEmats(
            bam_file_path=bam_file_path,
            paired_ends=paired_ends,
            min_phred_score=min_phred_score,
            df_CpGs_per_region=df_CpGs_filtered_per_region,
            chr_name=chr_name
        )
        
        # Build save path using phenotypes
        pheno_vals = [str(row[label]) for label in pheno_label_list]
        # Safe join for folder structure
        subdir = os.path.join(save_path, *pheno_vals)
        os.makedirs(subdir, exist_ok=True)
        # Build filename: group_condition_replicate_MEmats.pkl (adapt as needed)
        fname = f"{'_'.join(pheno_vals)}_ref_{chr_name}_{bam_file_name_wo_ext}_MEmats.pkl"
        out_path = os.path.join(subdir, fname)
        
        # Save
        df_MEmats.to_pickle(out_path)
        print(f"Saved MEmat in {out_path}")

    return 

