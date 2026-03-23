import os

# Load variables from config file
configfile: "_config.yml"

# Path to the directory where genome ref will be stored
MATRILINE_DIR = os.path.abspath(config["project_root"])
GENOME_DIR    = os.path.join(MATRILINE_DIR, config["genome_name"])
GENCODE_DIR = os.path.join(GENOME_DIR, config["gencode_name"])

LOCAL_ROOT    = os.path.expanduser(config["local_root"])
LOCAL_GENOME_DIR  = os.path.join(LOCAL_ROOT, config["genome_name"])
LOCAL_GENCODE_DIR = os.path.join(LOCAL_GENOME_DIR, config["gencode_name"])

print("")
print(f"server path: {GENCODE_DIR}\n")
print(f"local path: {LOCAL_GENCODE_DIR}\n")
print("")

# Flag to say it is sync finished
gencode_flag = os.path.join(LOCAL_GENCODE_DIR, ".gencode_sync2serv")

genome_zip_file = f"{LOCAL_GENCODE_DIR}/GRCm39.primary_assembly.genome.fa.gz"
genome_file     = f"{LOCAL_GENCODE_DIR}/GRCm39.primary_assembly.genome.fa"

transcripts_zip_file = f"{LOCAL_GENCODE_DIR}/gencode.vM37.transcripts.fa.gz"
transcripts_file     = f"{LOCAL_GENCODE_DIR}/gencode.vM37.transcripts.fa"

gtf_zip_file = f"{LOCAL_GENCODE_DIR}/gencode.vM37.primary_assembly.annotation.gtf.gz"
gtf_file     = f"{LOCAL_GENCODE_DIR}/gencode.vM37.primary_assembly.annotation.gtf"  # Uncompressed, optional


# The main workflow target: obtain the uncompressed reference genome FASTA file
rule all:
    input:
        gencode_flag

# Download the compressed genome FASTA and the global md5sum.txt from GENCODE
rule download_gencode:
    output:
        genome_zip = genome_zip_file,
        gtf_zip = gtf_zip_file,
        transcripts_zip = transcripts_zip_file
    conda:
        "envs/download_tools.yml"
    shell:
        """
        mkdir -p {LOCAL_GENCODE_DIR}
        wget -O {output.genome_zip} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.primary_assembly.genome.fa.gz
        wget -O {output.gtf_zip} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.primary_assembly.annotation.gtf.gz
        wget -O {output.transcripts_zip} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.transcripts.fa.gz
        """


rule unzip_gencode:
    input:
        genome_zip = genome_zip_file,
        transcripts_zip = transcripts_zip_file,
        gtf_zip = gtf_zip_file
    output:
        genome = genome_file,
        transcripts = transcripts_file,
        gtf = gtf_file
    conda:
        "envs/download_tools.yml"
    shell:
        """
        gunzip -c {input.genome_zip} > {output.genome}
        gunzip -c {input.transcripts_zip} > {output.transcripts}
        gunzip -c {input.gtf_zip} > {output.gtf}
        """


rule make_decoys:
    input:
        genome = genome_file,
    output:
        decoys = f"{LOCAL_GENCODE_DIR}/decoys.txt"
    conda:
        "envs/download_tools.yml"
    shell:
        """
        grep '^>' {input.genome} | cut -d " " -f 1 | sed 's/^>//' > {output.decoys}
        """

rule make_gentrome:
    input:
        transcripts = transcripts_file,
        genome = genome_file  # Uncompressed genome FASTA
    output:
        gentrome = f"{LOCAL_GENCODE_DIR}/gentrome.fa"
    conda:
        "envs/download_tools.yml"
    shell:
        """
        cat {input.transcripts} {input.genome} > {output.gentrome}
        """

# rule update_config:
#     input:
#         renamed_dir=RENAMED_FOLDER
#     output:
#         touch(".config_updated.log")
#     run:
#         import yaml
#         config_path = "config.yml"
#         with open(config_path, "r") as f:
#             config_data = yaml.safe_load(f)
#         config_data["raw_data_renamed_dir"] = renamed_dir
#         with open(config_path, "w") as f:
#             yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)
#         print(f"Updated {config_path} with raw_data_renamed_dir: {renamed_dir}")


rule sync_to_server:
    input:
        decoys = f"{LOCAL_GENCODE_DIR}/decoys.txt",
        gentrome = f"{LOCAL_GENCODE_DIR}/gentrome.fa"
    output:
        flag = temp(gencode_flag)
    params:
        local_dir = LOCAL_GENCODE_DIR,
        server_dir = GENCODE_DIR
    conda:
        "envs/download_tools.yml"
    shell:
        r"""
        # Sync all files from local genome directory to the server genome 
        rsync -av "{params.local_dir}/" "{params.server_dir}/"
        # Create a flag file to signal completion
        touch {output.flag}
        """