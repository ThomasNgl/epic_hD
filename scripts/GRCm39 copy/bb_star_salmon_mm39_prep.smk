import os

# Load variables from config file
configfile: "_config.yml"

# Path to the directory where genome ref will be stored
MATRILINE_DIR = os.path.abspath(config["project_root"])
GENOME_DIR    = os.path.join(MATRILINE_DIR, config["genome_name"])
GENCODE_DIR = os.path.join(GENOME_DIR, config["gencode_name"])
STAR_DIR = os.path.join(GENOME_DIR, config["star_name"])
SALMON_DIR = os.path.join(GENOME_DIR, config["salmon_name"])

LOCAL_ROOT    = os.path.expanduser(config["local_root"])
LOCAL_GENOME_DIR  = os.path.join(LOCAL_ROOT, config["genome_name"])

LOCAL_GENCODE_DIR = os.path.join(LOCAL_GENOME_DIR, config["gencode_name"])
LOCAL_STAR_DIR = os.path.join(LOCAL_GENOME_DIR, config["star_name"])
LOCAL_SALMON_DIR = os.path.join(LOCAL_GENOME_DIR, config["salmon_name"])

# Flag to say it is sync finished
sync2loc_flag = os.path.join(GENCODE_DIR, ".gencode_sync2loc")
sync2serv_flag = os.path.join(LOCAL_GENCODE_DIR, ".gencode_sync2serv")

local_star_tag = os.path.join(LOCAL_STAR_DIR, ".star_indexed")
local_salmon_tag = os.path.join(LOCAL_SALMON_DIR, ".salmon_indexed")

genome_file = os.path.join(GENCODE_DIR, config["gencode_genome_fa"])
annotation_file = os.path.join(GENCODE_DIR, config["annotation_gtf"])

gentrom_file = os.path.join(GENCODE_DIR, config["gentrome_fa"])
decoys_file = os.path.join(GENCODE_DIR, config["decoys_txt"])

local_genome_file = os.path.join(LOCAL_GENCODE_DIR, config["gencode_genome_fa"])
local_annotation_file = os.path.join(LOCAL_GENCODE_DIR, config["annotation_gtf"])
local_gentrom_file = os.path.join(LOCAL_GENCODE_DIR, config["gentrome_fa"])
local_decoys_file = os.path.join(LOCAL_GENCODE_DIR, config["decoys_txt"])

print("")
print(f"server path: {GENCODE_DIR}\n")
print(f"local path: {LOCAL_STAR_DIR}\n")
print("")

rule all:
    input:
        sync2serv_flag

rule sync_to_local:
    input:
        # To be sure they are there as they are necessary
        gentrom = gentrom_file
    output:
        flag = temp(sync2loc_flag),
        genome = local_genome_file,
        annotation = local_annotation_file,
        gentrome = local_gentrom_file,
        decoys = local_decoys_file
    params:
        server_dir = GENCODE_DIR,
        local_dir = LOCAL_GENCODE_DIR,
        local_star_dir = LOCAL_STAR_DIR,
        local_salmon_dir = LOCAL_SALMON_DIR

    shell:
        r"""
        # Sync all files from server bismarkl genome directory to the local 
        mkdir -p "{params.local_dir}"
        rsync -av "{params.server_dir}/" "{params.local_dir}/"
        mkdir -p "{params.local_star_dir}"
        mkdir -p "{params.local_salmon_dir}"
        # Create a flag file to signal completion
        touch {output.flag}
        """


rule index_star:
    input:
        flag = sync2loc_flag,
        genome = local_genome_file,
        annotation = local_annotation_file
    output:
        star_dir = directory(LOCAL_STAR_DIR),
        star_tag = temp(local_star_tag)
    params:
        overhang = config["overhang"]
    threads: 1
    conda:
        "envs/star.yml"
    shell:
        """
        cd {LOCAL_STAR_DIR}
        
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output.star_dir} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.annotation} \
             --sjdbOverhang {params.overhang}
        
        touch {output.star_tag}
        """


rule index_salmon:
    input:
        flag = sync2loc_flag,
        gentrome = local_gentrom_file,
        decoys = local_decoys_file
    output:
        salmon_dir = directory(LOCAL_SALMON_DIR),
        salmon_tag = temp(local_salmon_tag)
    threads: 6
    conda:
        "envs/salmon.yml"
    shell:
        """
        cd {LOCAL_SALMON_DIR}

        salmon index -p {threads} \
            -t {input.gentrome} \
            -d {input.decoys} \
            -i {output.salmon_dir} \
            --gencode
        
        touch {output.salmon_tag}
        """


rule sync_to_server:
    input:
        star_tag = local_star_tag,
        salmon_tag = local_salmon_tag
    output:
        flag = temp(sync2serv_flag)
    params:
        local_star_dir = LOCAL_STAR_DIR,
        local_salmon_dir = LOCAL_SALMON_DIR,
        star_dir = STAR_DIR,
        salmon_dir = SALMON_DIR

    shell:
        r"""
        mkdir -p "{params.star_dir}"
        mkdir -p "{params.salmon_dir}"
        
        rsync -av --progress "{params.local_star_dir}/" "{params.star_dir}/"
        rsync -av --progress "{params.local_salmon_dir}/" "{params.salmon_dir}/"
        
        # Create a flag file to signal completion
        touch {output.flag}
        """