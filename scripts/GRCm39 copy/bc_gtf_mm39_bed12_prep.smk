import os

# Load variables from config file
configfile: "_config.yml"

# --- Config & paths ----------------------------------------------------------
ENV_DIR = os.path.abspath(config["env_dir"])

# Path to the directory where genome ref will be stored
LOCAL_PROJECT_DIR = os.path.abspath(config["analysis_dir"])
LOCAL_GENOME_DIR    = os.path.join(LOCAL_PROJECT_DIR, config["genome_name"])
LOCAL_GENCODE_DIR = os.path.join(LOCAL_GENOME_DIR, config["gencode_name"])

local_gtf_file = os.path.join(LOCAL_GENCODE_DIR, config["annotation_gtf"])
print("input file: ", local_gtf_file)

rule all:
    input:
        os.path.join(LOCAL_GENCODE_DIR, "gencode.vM37.mm39.bed12")

rule gtf_to_bed12:
    input:
        gtf = local_gtf_file  # e.g., .../gencode.vM37.annotation.gtf
    output:
        bed = os.path.join(LOCAL_GENCODE_DIR, "gencode.vM37.mm39.bed12")
    log:
        os.path.join(LOCAL_GENCODE_DIR, "gtf2bed.log")
    conda:
        os.path.join(ENV_DIR, "gtf2bed.yml")
    shell:
        r"""
        set -euo pipefail
        export TMPDIR="{LOCAL_GENCODE_DIR}"

        # Convert GTF -> genePred (extended) -> BED12 and sort
        gtfToGenePred -genePredExt -allErrors "{input.gtf}" stdout \
        | genePredToBed stdin stdout \
        | sort -k1,1 -k2,2n > "{output.bed}" 2> "{log}"
        """
