import os

configfile: "_config.yml"

MATRILINE_DIR = os.path.abspath(config["project_root"])
GENOME_DIR    = os.path.join(MATRILINE_DIR, config["genome_name"])
BISMARK_DIR = os.path.join(GENOME_DIR, config["bismark_genome_name"])
BISULFITE_DIR = os.path.join(BISMARK_DIR, config["bisulfite_conversion_name"])
CT_CONV_DIR    = os.path.join(BISULFITE_DIR, config["ct_conv_name"])
CPG = "CpGs_location"
CHR_NAMES_TXT = os.path.join(BISMARK_DIR, CHR_NAMES)

GA_CONV_DIR    = os.path.join(BISULFITE_DIR, config["ga_conv_name"])
LOCAL_CPG_DIR = os.path.join(BISMARK_DIR, CPG)

LOCAL_ROOT    = os.path.expanduser(config["local_root"])
LOCAL_GENOME_DIR  = os.path.join(LOCAL_ROOT, config["genome_name"])
LOCAL_BISMARK_DIR =  os.path.join(LOCAL_GENOME_DIR, config["bismark_genome_name"])
LOCAL_BISMARK_GENOME_FA = os.path.join(LOCAL_BISMARK_DIR, config["genome_fa"])
LOCAL_CPG_DIR = os.path.join(LOCAL_BISMARK_DIR, CPG)
CHR_NAMES = config["chr_names_file"]

# Flag to say it is sync finished
sync2loc_flag = os.path.join(BISMARK_DIR, ".cpgloc_sync2loc")
sync2serv_flag = os.path.join(LOCAL_BISMARK_DIR, ".cpgloc_sync2serv")

print("")
print(f"server path: {BISMARK_DIR}\n")
print(f"local path: {LOCAL_BISMARK_DIR}\n")
print("")

def read_chr_names(wildcards):
    with open(CHR_NAMES_TXT) as f:
        return [line.strip() for line in f]

rule all:
    input:
        sync2serv_flag,
        expand(os.path.join(LOCAL_CPG_DIR, "CpGslocation_{chr}.csv"), chr=read_chr_names)

rule sync_to_local:
    input:
        # To be sure they are there as they are necessary
        ct_conv_dir = CT_CONV_DIR,
        ga_conv_dir = GA_CONV_DIR
    output:
        flag = temp(sync2loc_flag)
    params:
        server_dir = BISMARK_DIR,
        local_dir = LOCAL_BISMARK_DIR,
        cpg_dir = LOCAL_CPG_DIR
    shell:
        r"""
        # Sync all files from server bismarkl genome directory to the local 
        mkdir -p "{params.local_dir}"
        rsync -av "{params.server_dir}/" "{params.local_dir}/"
        # Create a flag file to signal completion
        mkdir -p "{params.cpg_dir}"
        touch {output.flag}
        """


rule extract_CpGs:
    input:
        fasta=LOCAL_BISMARK_GENOME_FA
    output:
       os.path.join(LOCAL_CPG_DIR, "CpGslocation_{chr}.csv")
    params:
        cpg_dir = LOCAL_CPG_DIR,
        window = config["density_window"]
    run:
        from pyformME import a_fasta2CpG
        a_fasta2CpG.fasta2CpGs_chr(
            fasta_file=input.fasta,
            out_dir=params.cpg_dir,
            chr_name=wildcards.chr,
            window=params.window
        )

rule sync_to_server:
    input:
        # To be sure they are there as they are necessary
        expand(os.path.join(LOCAL_CPG_DIR, "CpGslocation_{chr}.csv"), chr=read_chr_names)
    output:
        flag = temp(sync2serv_flag)
    params:
        local_dir = LOCAL_CPG_DIR,
        server_dir = CPG_DIR
    shell:
        r"""
        # Sync all files from server bismarkl genome directory to the local 
        mkdir -p "{params.server_dir}"
        rsync -av "{params.local_dir}/" "{params.server_dir}/"
        # Create a flag file to signal completion
        touch {output.flag}
        """