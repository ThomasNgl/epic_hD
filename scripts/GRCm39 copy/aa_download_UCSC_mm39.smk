import os

# Load variables from config file
configfile: "_config.yml"

# Path to the directory where genome ref will be stored
MATRILINE_DIR = os.path.abspath(config["project_root"])
GENOME_DIR    = os.path.join(MATRILINE_DIR, config["genome_name"])

LOCAL_ROOT    = os.path.expanduser(config["local_root"])
LOCAL_GENOME_DIR  = os.path.join(LOCAL_ROOT, config["genome_name"])

print("")
print(f"server path: {MATRILINE_DIR}\n")
print(f"local path: {LOCAL_ROOT}\n")
print("")

# Flag to say it is sync finished
sync_flag = os.path.join(LOCAL_GENOME_DIR, ".fa_sync2serv")

# The main workflow target: obtain the uncompressed reference genome FASTA file
rule all:
    input:
        sync_flag

# Download the compressed genome FASTA and the global md5sum.txt from UCSC
rule download_genome:
    output:
        f"{LOCAL_GENOME_DIR}/mm39.fa.gz",      # The gzipped FASTA file (temporary)
        f"{LOCAL_GENOME_DIR}/md5sum.txt"             # The checksum file for all UCSC files
    conda:
        "envs/download_tools.yml"
    shell:
        """
        mkdir -p {LOCAL_GENOME_DIR}
        wget -O {output[0]} http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
        wget -O {output[1]} http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/md5sum.txt
        """

# Verify that the downloaded genome matches its md5 checksum
rule verify_genome:
    input:
        fa_zip = f"{LOCAL_GENOME_DIR}/mm39.fa.gz",
        md5s = f"{LOCAL_GENOME_DIR}/md5sum.txt"
    output:
        f"{LOCAL_GENOME_DIR}/mm39.checked"           # Simple file to indicate successful check
    conda:
        "envs/download_tools.yml"
    shell:
        """
        cd {LOCAL_GENOME_DIR}
        grep 'mm39.fa.gz' md5sum.txt > mm39.fa.gz.md5
        md5sum -c mm39.fa.gz.md5
        touch mm39.checked
        """

# Uncompress the genome FASTA file for downstream use
rule unzip_genome:
    input:
        fa_zip = f"{LOCAL_GENOME_DIR}/mm39.fa.gz",
        checked = f"{LOCAL_GENOME_DIR}/mm39.checked"
    output:
        fa = f"{LOCAL_GENOME_DIR}/mm39.fa"
    conda:
        "envs/download_tools.yml"
    shell:
        """
        gunzip -c {input.fa_zip} > {output.fa}
        """

rule index_genome:
    input:
        fa = f"{LOCAL_GENOME_DIR}/mm39.fa",
    output:
        fai = f"{LOCAL_GENOME_DIR}/mm39.fa.fai"
    conda:
        "envs/download_tools.yml"
    shell:
        """
        samtools faidx {input.fa}
        """

rule sync_to_server:
    input:
        fai = f"{LOCAL_GENOME_DIR}/mm39.fa.fai"
    output:
        # Use a flag file to indicate completion (not required but useful)
        flag = temp(sync_flag)
    params:
        local_dir = LOCAL_GENOME_DIR,
        server_dir = GENOME_DIR
    shell:
        r"""
        # Sync all files from local genome directory to the server genome 
        rsync -av "{params.local_dir}/" "{params.server_dir}/"
        # Create a flag file to signal completion
        touch {output.flag}
        """