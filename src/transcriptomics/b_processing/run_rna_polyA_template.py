#!/usr/bin/env python3
"""
Runner template for RNA-seq poly-A pipeline.

Copy this file to scripts/transcriptomics/ and rename it for your dataset,
e.g. run_sg_oocyte_MSUS37_DK.py

Edit the sections marked with ← to configure your experiment.
"""
import os
import argparse
from snakemake import snakemake
from config import EnvConfig  # repo-level config: sets REPO_DIR, paths, etc.

# ── Environment paths (from EnvConfig) ────────────────────────────────────
env_config  = EnvConfig()
results_dir = env_config.MAIN_DIRS["results"]
repo_dir    = env_config.REPO_DIR
env_dir     = env_config.EXP_DIRS["envs"]

# ── Experiment settings  ← edit these ─────────────────────────────────────
experiment_folder = "sg_oocyte_MSUS37_DK"   # subfolder under results_dir
modality          = "SMARTseq"               # e.g. SMARTseq, NovaSeq, TotalRNA
genome_folder     = "GRCm39"                 # genome subfolder under results_dir

# ── Step parameters  ← set run_* to True/False to enable/disable each step ─
# suffix: append a string to the output directory name for this step
# (useful to compare results with different parameters without overwriting)

fqc_params = {
    "run_fqc"  : True,
    "threads"  : 2,
    "suffix"   : ""
}

fq_screen_params = {
    "run_fq_screen"    : True,
    "threads"          : 2,
    "fastq_screen_conf": "/mnt/auxiliary/fastq_screen/fastq_screen.conf",
    "suffix"           : ""
}

trimming_params = {
    "run_trimming": True,
    "quality"     : 30,
    "length"      : 30,
    "stringency"  : 3,
    "threads"     : 2,
    "suffix"      : ""
}

fqc_trimmed_params = {
    "run_fqc" : True,
    "threads" : 2,
    "suffix"  : ""
}

star_align_params = {
    "run_star"       : True,
    "threads"        : 6,
    "star_genome_dir": "/mnt/auxiliary/star/mm39_star_125",
    "gencode_folder" : "gencode_primary",
    "gtf_file"       : "gencode.vM37.primary_assembly.annotation.gtf",
    "suffix"         : ""
}

feature_count_params = {
    "run_feature_count": True,
    "threads"          : 6,
    "suffix"           : ""
}

rseq_params = {
    "run_rseqc"     : True,
    "threads"       : 2,
    "gencode_folder": "gencode_primary",
    "bed_file"      : "gencode.vM37.mm39.bed12",
    "suffix"        : ""
}

salmon_quant_params = {
    "run_salmon_quant": True,
    "threads"         : 6,
    "salmon_index_dir": "/mnt/auxiliary/salmon/mm39_salmon_index",
    "suffix"          : ""
}

# ── Assemble full config ───────────────────────────────────────────────────
full_config = {
    "results_dir"         : results_dir,
    "repo_dir"            : repo_dir,
    "env_dir"             : env_dir,
    "experiment_folder"   : experiment_folder,
    "modality"            : modality,
    "genome_folder"       : genome_folder,
    "fqc_params"          : fqc_params,
    "fq_screen_params"    : fq_screen_params,
    "trimming_params"     : trimming_params,
    "fqc_trimmed_params"  : fqc_trimmed_params,
    "star_align_params"   : star_align_params,
    "feature_count_params": feature_count_params,
    "rseq_params"         : rseq_params,
    "salmon_quant_params" : salmon_quant_params,
}

# ← Change to rna_polyA_pe.smk for paired-end data
SNAKEFILE = os.path.join(repo_dir, "src", "transcriptomics", "pipeline", "rna_polyA_se.smk")

# ── CLI ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run RNA-seq poly-A pipeline")
    parser.add_argument(
        "--n_samples", default="all",
        help="Samples to process: 'all', an integer (first N), "
             "or comma-separated substrings (e.g. '3_2,4_5'). Default: all"
    )
    parser.add_argument(
        "--cores", type=int, default=6,
        help="Number of cores (default: 6)"
    )
    parser.add_argument(
        "--dry-run", "-n", action="store_true",
        help="Dry run — print jobs without executing"
    )
    args = parser.parse_args()

    # Parse n_samples
    if args.n_samples == "all":
        full_config["n_samples"] = "all"
    else:
        try:
            full_config["n_samples"] = int(args.n_samples)
        except ValueError:
            full_config["n_samples"] = [x.strip() for x in args.n_samples.split(",")]

    print(f"\nPipeline  : {SNAKEFILE}")
    print(f"Experiment: {full_config['experiment_folder']} / {full_config['modality']}")
    print(f"Genome    : {full_config['genome_folder']}")
    print(f"Samples   : {full_config['n_samples']}")
    print(f"Cores     : {args.cores}")
    print(f"Dry run   : {args.dry_run}\n")

    success = snakemake(
        snakefile       = SNAKEFILE,
        config          = full_config,
        cores           = args.cores,
        use_conda       = True,
        conda_frontend  = "conda",
        dryrun          = args.dry_run,
        printshellcmds  = True,
        keepgoing       = False,
        force_incomplete= True
    )

    exit(0 if success else 1)
