import os

# Load variables from config file
configfile: "_config.yml"

MATRILINE_DIR = os.path.abspath(config["project_root"])
GENOME_DIR    = os.path.join(MATRILINE_DIR, config["genome_name"])

# Path to the input genome fasta file (from config)
GENOME_FA = os.path.join(GENOME_DIR, config["genome_fa"])
GENOME_FAI = os.path.join(GENOME_DIR, config["genome_fai"])
# Directory where Bismark will prepare the genome (from config)
BISMARK_DIR = os.path.join(GENOME_DIR, config["bismark_genome_name"])

LOCAL_ROOT    = os.path.expanduser(config["local_root"])
LOCAL_GENOME_DIR  = os.path.join(LOCAL_ROOT, config["genome_name"])

LOCAL_GENOME_FA = os.path.join(LOCAL_GENOME_DIR, config["genome_fa"])
LOCAL_GENOME_FAI = os.path.join(LOCAL_GENOME_DIR, config["genome_fai"])
LOCAL_BISMARK_DIR =  os.path.join(LOCAL_GENOME_DIR, config["bismark_genome_name"])

LOCAL_BISMARK_GENOME_FA = os.path.join(LOCAL_BISMARK_DIR, config["genome_fa"])

CHR_NAMES = config["chr_names_file"]
LOCAL_CHR_NAMES_TXT = os.path.join(LOCAL_BISMARK_DIR, CHR_NAMES)
LOCAL_BISULFITE_DIR = os.path.join(LOCAL_BISMARK_DIR, config["bisulfite_conversion_name"])

LOCAL_CT_CONV_DIR    = os.path.join(LOCAL_BISULFITE_DIR, config["ct_conv_name"])
LOCAL_GA_CONV_DIR    = os.path.join(LOCAL_BISULFITE_DIR, config["ga_conv_name"])

# Flag to say it is sync finished
sync_flag = os.path.join(LOCAL_BISMARK_DIR, ".bismarkconv_sync2serv")

print("")
print(f"server path: {MATRILINE_DIR}\n")
print(f"local path: {LOCAL_BISMARK_DIR}\n")
print("")

# Final rule: this workflow only runs the genome preparation step
rule all:
    input:
        sync_flag

rule copy_fa:
    input:
        src_fa=GENOME_FA,
        src_fai=GENOME_FAI
    output:
        dest_fa=LOCAL_GENOME_FA,
        dest_fai=LOCAL_GENOME_FAI
    shell:
        """
        mkdir -p $(dirname {output.dest_fa})
        cp {input.src_fa} {output.dest_fa}
        cp {input.src_fai} {output.dest_fai}
        """

# Rule to run bismark_genome_preparation on the reference genome
# USES 2 x THREADS (one per strand)
rule bismark_genome_prep:
    input:
        fa=LOCAL_GENOME_FA,
        fai=LOCAL_GENOME_FAI
    output:
        # Dummy output: just the genome directory to trigger the rule
        ct_conv_dir = directory(LOCAL_CT_CONV_DIR),
        ga_conv_dir = directory(LOCAL_GA_CONV_DIR),
        fasta_file_path = LOCAL_BISMARK_GENOME_FA
    params:
        local_dir = LOCAL_BISMARK_DIR,
    conda:
        "envs/bismark_env.yml"
    threads: 6
    shell:
        """
        echo "Preparing Bismark genome in: {params.local_dir}"
        export TMPDIR={params.local_dir}
        # Make the output directory if it doesn't exist
        mkdir -p {params.local_dir}
        # Copy the fasta file into the Bismark genome directory
        cp {input.fa} {params.local_dir}/
        cp {input.fai} {params.local_dir}/
        # Run Bismark genome preparation (builds all indexes and conversions)
        bismark_genome_preparation --bowtie2 --parallel {threads} {params.local_dir}
        """


rule get_chr_names:
    input:
        fasta_file_path = LOCAL_BISMARK_GENOME_FA
    output:
        chr_names = LOCAL_CHR_NAMES_TXT
    run:
        from pyformME import io_informME
        fasta_file = io_informME.load_fasta(input.fasta_file_path)
        chr_name_list = list(fasta_file.references)
        with open(output.chr_names, "w") as f:
            for chr_name in chr_name_list:
                f.write(chr_name + "\n")


rule sync_to_server:
    input:
        ct_conv_dir = LOCAL_CT_CONV_DIR,
        ga_conv_dir = LOCAL_GA_CONV_DIR,
        chr_names = LOCAL_CHR_NAMES_TXT
    output:
        # Use a flag file to indicate completion (not required but useful)
        flag = temp(sync_flag)
    params:
        local_dir = LOCAL_BISMARK_DIR,
        server_dir = BISMARK_DIR
    shell:
        r"""
        # Sync all files from local experiment directory to the server project/experiment
        rsync -av "{params.local_dir}/" "{params.server_dir}/"
        # Create a flag file to signal completion
        touch {output.flag}
        """