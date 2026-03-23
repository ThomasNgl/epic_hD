"""
RNA-seq poly-A paired-end pipeline (PE)
=========================================
Runs: FastQC → TrimGalore (paired) → FastQC trimmed → FastQ Screen →
      STAR align → featureCounts → RSeQC (TIN + genebody) → Salmon
      → MultiQC per step

Each step is independently toggled via config["<step>_params"]["run_<step>"].
Output directories accept a suffix via config["<step>_params"]["suffix"] to
allow parameter comparisons without overwriting existing results.

Called from a dataset-specific runner script (scripts/transcriptomics/).
All paths are injected by the runner — nothing is hardcoded here.
"""
import os

# ── Paths ──────────────────────────────────────────────────────────────────
ENV_DIR          = os.path.abspath(config["env_dir"])
RESULTS_DIR      = os.path.abspath(config["results_dir"])
REPO_DIR         = os.path.abspath(config["repo_dir"])

EXPERIMENT_DIR   = os.path.join(RESULTS_DIR, config["experiment_folder"], config["modality"])
RAW_DIR          = os.path.join(EXPERIMENT_DIR, "raw_data_renamed")
GENOME_DIR       = os.path.join(RESULTS_DIR, config["genome_folder"])
GENCODE_DIR      = os.path.join(GENOME_DIR, config["star_align_params"]["gencode_folder"])

GTF_FILE         = os.path.join(GENCODE_DIR, config["star_align_params"]["gtf_file"])
BED_FILE         = os.path.join(GENCODE_DIR, config["rseq_params"]["bed_file"])

# ── Output directory names (with optional suffix) ──────────────────────────
PREPROCESS_DIR   = os.path.join(EXPERIMENT_DIR, "01_preprocess")
LOG_DIR          = os.path.join(EXPERIMENT_DIR, "log")
QC_DIR           = os.path.join(EXPERIMENT_DIR, "02_QC_plot_FastQ2Count")

FASTQC_RAW       = "01_FastQC_raw_data"         + config["fqc_params"]["suffix"]
TRIMMED          = "02_trimmed_data"             + config["trimming_params"]["suffix"]
FASTQC_TRIM      = "03_FastQC_trimmed_data"      + config["fqc_trimmed_params"]["suffix"]
FQ_SCREEN        = "04_FastQScreen_trimmed_data" + config["fq_screen_params"]["suffix"]
STAR_ALIGN       = "05_STAR_aligned_sorted_indexed_data" + config["star_align_params"]["suffix"]
FEATURE_COUNT    = "06_STAR_feature_counts"      + config["feature_count_params"]["suffix"]
COVERAGE_TIN     = "07_STAR_coverage_tin_score"  + config["rseq_params"]["suffix"]
SALMON_QUANT     = "05_Salmon_quant"             + config["salmon_quant_params"]["suffix"]

# Convenience path builders
def pre(step):  return os.path.join(PREPROCESS_DIR, step)
def log(step):  return os.path.join(LOG_DIR, "01_preprocess", step)
def qc(step):   return os.path.join(QC_DIR, step)

# ── Sample selection ───────────────────────────────────────────────────────
# Paired-end: discover samples from R1 files, assume R2 exists
all_r1 = set(glob_wildcards(os.path.join(RAW_DIR, "{sample}_R1.fq.gz")).sample)
all_r2 = set(glob_wildcards(os.path.join(RAW_DIR, "{sample}_R2.fq.gz")).sample)
all_samples = sorted(all_r1 & all_r2)

r1_only = all_r1 - all_r2
r2_only = all_r2 - all_r1
if r1_only:
    print(f"WARNING: R1 without R2: {r1_only}")
if r2_only:
    print(f"WARNING: R2 without R1: {r2_only}")

print(f"\nExperiment : {EXPERIMENT_DIR}")
print(f"Paired samples found ({len(all_samples)}): {all_samples}\n")

n = config["n_samples"]
if n == "all":
    samples = all_samples
elif isinstance(n, int):
    samples = all_samples[:n]
    print(f"Using first {n} samples: {samples}")
elif isinstance(n, list):
    samples = [s for s in all_samples if any(sub in s for sub in n)]
    print(f"Using samples matching {n}: {samples}")
else:
    raise ValueError(f"Unsupported n_samples type: {type(n)}")

# ── MultiQC target list (only enabled steps) ───────────────────────────────
active_steps = []
if config["fqc_params"]["run_fqc"]:            active_steps.append(FASTQC_RAW)
if config["trimming_params"]["run_trimming"]:   active_steps.append(TRIMMED)
if config["fqc_trimmed_params"]["run_fqc"]:     active_steps.append(FASTQC_TRIM)
if config["fq_screen_params"]["run_fq_screen"]: active_steps.append(FQ_SCREEN)
if config["star_align_params"]["run_star"]:     active_steps.append(STAR_ALIGN)
if config["feature_count_params"]["run_feature_count"]: active_steps.append(FEATURE_COUNT)
if config["rseq_params"]["run_rseqc"]:          active_steps.append(COVERAGE_TIN)
if config["salmon_quant_params"]["run_salmon_quant"]: active_steps.append(SALMON_QUANT)

# ── Rule all ──────────────────────────────────────────────────────────────
rule all:
    input:
        expand(os.path.join(QC_DIR, "{step}/multiqc_report.html"), step=active_steps)


# ── Step 1: FastQC raw (both reads) ───────────────────────────────────────
if config["fqc_params"]["run_fqc"]:
    rule fqc_raw:
        input:
            fastq = os.path.join(RAW_DIR, "{sample}_R{read}.fq.gz")
        output:
            html = os.path.join(pre(FASTQC_RAW), "{sample}_R{read}_fastqc.html"),
            zip  = os.path.join(pre(FASTQC_RAW), "{sample}_R{read}_fastqc.zip")
        log:
            os.path.join(log(FASTQC_RAW), "{sample}_R{read}.log")
        conda:
            os.path.join(ENV_DIR, "fastqc.yml")
        threads: config["fqc_params"]["threads"]
        shell:
            """
            export TMPDIR={pre(FASTQC_RAW)}
            fastqc -t {threads} -o {pre(FASTQC_RAW)} {input.fastq} > {log} 2>&1
            """


# ── Step 2: Trimming (paired-end) ─────────────────────────────────────────
if config["trimming_params"]["run_trimming"]:
    rule trimming:
        input:
            R1 = os.path.join(RAW_DIR, "{sample}_R1.fq.gz"),
            R2 = os.path.join(RAW_DIR, "{sample}_R2.fq.gz")
        output:
            trimmed1  = os.path.join(pre(TRIMMED), "{sample}_R1_val_1.fq.gz"),
            trimmed2  = os.path.join(pre(TRIMMED), "{sample}_R2_val_2.fq.gz"),
            report1   = os.path.join(pre(TRIMMED), "{sample}_R1.fq.gz_trimming_report.txt"),
            report2   = os.path.join(pre(TRIMMED), "{sample}_R2.fq.gz_trimming_report.txt"),
            unpaired1 = os.path.join(pre(TRIMMED), "{sample}_R1_unpaired_1.fq.gz"),
            unpaired2 = os.path.join(pre(TRIMMED), "{sample}_R2_unpaired_2.fq.gz")
        log:
            os.path.join(log(TRIMMED), "{sample}.log")
        params:
            quality    = config["trimming_params"]["quality"],
            length     = config["trimming_params"]["length"],
            stringency = config["trimming_params"]["stringency"],
            outdir     = pre(TRIMMED)
        threads: config["trimming_params"]["threads"]
        conda:
            os.path.join(ENV_DIR, "trimg.yml")
        shell:
            r"""
            set -euo pipefail
            export TMPDIR="{params.outdir}"
            trim_galore --paired \
                -q {params.quality} \
                --length {params.length} \
                --trim-n \
                --stringency {params.stringency} \
                --output_dir "{params.outdir}" \
                --retain_unpaired \
                --cores {threads} \
                {input.R1} {input.R2} &> {log}
            """


# ── Step 3: FastQC trimmed (both reads) ───────────────────────────────────
if config["fqc_trimmed_params"]["run_fqc"]:
    rule fqc_trimmed:
        input:
            fastq = os.path.join(pre(TRIMMED), "{sample}_R{read}_val_{read}.fq.gz")
        output:
            html = os.path.join(pre(FASTQC_TRIM), "{sample}_R{read}_val_{read}_fastqc.html"),
            zip  = os.path.join(pre(FASTQC_TRIM), "{sample}_R{read}_val_{read}_fastqc.zip")
        log:
            os.path.join(log(FASTQC_TRIM), "{sample}_R{read}.log")
        conda:
            os.path.join(ENV_DIR, "fastqc.yml")
        threads: config["fqc_trimmed_params"]["threads"]
        shell:
            """
            export TMPDIR={pre(FASTQC_TRIM)}
            fastqc -t {threads} -o {pre(FASTQC_TRIM)} {input.fastq} > {log} 2>&1
            """


# ── Step 4: FastQ Screen ──────────────────────────────────────────────────
if config["fq_screen_params"]["run_fq_screen"]:
    rule fastq_screen:
        input:
            fastq = os.path.join(pre(TRIMMED), "{sample}_R{read}_val_{read}.fq.gz")
        output:
            txt  = os.path.join(pre(FQ_SCREEN), "{sample}_R{read}_val_{read}_screen.txt"),
            html = os.path.join(pre(FQ_SCREEN), "{sample}_R{read}_val_{read}_screen.html")
        params:
            conf = config["fq_screen_params"]["fastq_screen_conf"]
        log:
            os.path.join(log(FQ_SCREEN), "{sample}_R{read}.log")
        conda:
            os.path.join(ENV_DIR, "fqscreen.yml")
        threads: config["fq_screen_params"]["threads"]
        shell:
            """
            fastq_screen --aligner bowtie2 \
                         --threads {threads} \
                         --conf {params.conf} \
                         --outdir {pre(FQ_SCREEN)} \
                         {input.fastq} > {log} 2>&1
            """


# ── Step 5a: STAR align (paired-end) ──────────────────────────────────────
if config["star_align_params"]["run_star"]:
    rule star_align:
        input:
            R1 = os.path.join(pre(TRIMMED), "{sample}_R1_val_1.fq.gz"),
            R2 = os.path.join(pre(TRIMMED), "{sample}_R2_val_2.fq.gz")
        output:
            bam      = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam"),
            bai      = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.bai"),
            idxstat  = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.idxstats.txt"),
            stats    = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.stats.txt"),
            flagstat = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.flagstat.txt")
        params:
            star_dir   = config["star_align_params"]["star_genome_dir"],
            gtf        = GTF_FILE,
            out_prefix = os.path.join(pre(STAR_ALIGN), "{sample}_star_aligned_"),
            star_bam   = os.path.join(pre(STAR_ALIGN), "{sample}_star_aligned_Aligned.sortedByCoord.out.bam")
        log:
            os.path.join(log(STAR_ALIGN), "{sample}.log")
        conda:
            os.path.join(ENV_DIR, "star.yml")
        threads: config["star_align_params"]["threads"]
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {params.out_prefix})"
            STAR --runThreadN {threads} \
                 --genomeDir {params.star_dir} \
                 --sjdbGTFfile {params.gtf} \
                 --readFilesCommand zcat \
                 --readFilesIn {input.R1} {input.R2} \
                 --outFileNamePrefix {params.out_prefix} \
                 --outSAMtype BAM SortedByCoordinate \
                 > {log} 2>&1
            mv {params.star_bam} {output.bam}
            samtools index    -@ {threads} {output.bam}           2>> {log}
            samtools idxstats {output.bam} > {output.idxstat}     2>> {log}
            samtools stats    {output.bam} > {output.stats}       2>> {log}
            samtools flagstat {output.bam} > {output.flagstat}    2>> {log}
            """


# ── Step 6: featureCounts ─────────────────────────────────────────────────
if config["feature_count_params"]["run_feature_count"]:
    rule feature_counts:
        input:
            bams = expand(os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam"), sample=samples)
        output:
            counts  = os.path.join(pre(FEATURE_COUNT), "gcounts_exon_s0.txt"),
            summary = os.path.join(pre(FEATURE_COUNT), "gcounts_exon_s0.txt.summary")
        params:
            gtf = GTF_FILE
        log:
            os.path.join(log(FEATURE_COUNT), "counts.log")
        conda:
            os.path.join(ENV_DIR, "featurecounts.yml")
        threads: config["feature_count_params"]["threads"]
        shell:
            """
            featureCounts -T {threads} \
                          -s 0 \
                          -t exon \
                          -g gene_id \
                          -a {params.gtf} \
                          -o {output.counts} \
                          {input.bams} > {log} 2>&1
            """


# ── Step 7: RSeQC — TIN score ─────────────────────────────────────────────
if config["rseq_params"]["run_rseqc"]:
    rule calculate_tin:
        input:
            bam = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam")
        output:
            tsv = os.path.join(pre(COVERAGE_TIN), "{sample}_tin.tsv")
        params:
            bed = BED_FILE
        conda:
            os.path.join(ENV_DIR, "rseqc.yml")
        threads: config["rseq_params"]["threads"]
        shell:
            """
            calculate-tin.py \
                -r {params.bed} \
                -i {input.bam} \
                --names={wildcards.sample} \
                -p {threads} > {output.tsv}
            """

    rule merge_tin:
        input:
            expand(os.path.join(pre(COVERAGE_TIN), "{sample}_tin.tsv"), sample=samples)
        output:
            os.path.join(pre(COVERAGE_TIN), "merged_tin.tsv")
        conda:
            os.path.join(ENV_DIR, "rseqc.yml")
        shell:
            "merge-tin.py --input-files {input} --output-file {output}"

    rule plot_tin:
        input:
            os.path.join(pre(COVERAGE_TIN), "merged_tin.tsv")
        output:
            os.path.join(qc(COVERAGE_TIN), "tin_scores.png")
        params:
            script = os.path.join(REPO_DIR, "src/transcriptomics/pipeline/plot_tin_scores.R")
        conda:
            os.path.join(ENV_DIR, "r_ggplot2.yml")
        shell:
            "Rscript {params.script} {input} {output}"

    rule genebody_coverage:
        input:
            bam = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam"),
            bai = os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.bai")
        output:
            curves = temp(os.path.join(pre(COVERAGE_TIN), "{sample}.geneBodyCoverage.curves.pdf")),
            r      = temp(os.path.join(pre(COVERAGE_TIN), "{sample}.geneBodyCoverage.r")),
            txt    = os.path.join(pre(COVERAGE_TIN), "{sample}.geneBodyCoverage.txt")
        params:
            prefix = os.path.join(pre(COVERAGE_TIN), "{sample}"),
            bed    = BED_FILE
        log:
            os.path.join(log(COVERAGE_TIN), "{sample}_genebody.log")
        conda:
            os.path.join(ENV_DIR, "rseqc.yml")
        shell:
            """
            touch {input.bai}
            geneBody_coverage.py \
                -r {params.bed} \
                -i {input.bam} \
                -o {params.prefix} > {log} 2>&1
            """


# ── Step 5b: Salmon quant (paired-end, runs in parallel with STAR) ─────────
if config["salmon_quant_params"]["run_salmon_quant"]:
    rule salmon_quant:
        input:
            R1 = os.path.join(pre(TRIMMED), "{sample}_R1_val_1.fq.gz"),
            R2 = os.path.join(pre(TRIMMED), "{sample}_R2_val_2.fq.gz")
        output:
            quant_dir = directory(os.path.join(pre(SALMON_QUANT), "{sample}_quant")),
            quant     = os.path.join(pre(SALMON_QUANT), "{sample}_quant/quant.sf")
        params:
            index = config["salmon_quant_params"]["salmon_index_dir"]
        log:
            os.path.join(log(SALMON_QUANT), "{sample}.log")
        conda:
            os.path.join(ENV_DIR, "salmon.yml")
        threads: config["salmon_quant_params"]["threads"]
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {log})"
            salmon quant \
                -i {params.index} \
                -l A \
                -1 {input.R1} \
                -2 {input.R2} \
                -p {threads} \
                --validateMappings \
                --output {output.quant_dir} \
                &> {log}
            """


# ── MultiQC inputs per step ────────────────────────────────────────────────
READS = ["1", "2"]
MULTIQC_INPUTS = {}

if config["fqc_params"]["run_fqc"]:
    MULTIQC_INPUTS[FASTQC_RAW] = expand(
        os.path.join(pre(FASTQC_RAW), "{sample}_R{read}_fastqc.html"),
        sample=samples, read=READS)

if config["trimming_params"]["run_trimming"]:
    MULTIQC_INPUTS[TRIMMED] = (
        expand(os.path.join(pre(TRIMMED), "{sample}_R1.fq.gz_trimming_report.txt"), sample=samples) +
        expand(os.path.join(pre(TRIMMED), "{sample}_R2.fq.gz_trimming_report.txt"), sample=samples))
    if config["fqc_trimmed_params"]["run_fqc"]:
        MULTIQC_INPUTS[FASTQC_TRIM] = expand(
            os.path.join(pre(FASTQC_TRIM), "{sample}_R{read}_val_{read}_fastqc.html"),
            sample=samples, read=READS)
    if config["fq_screen_params"]["run_fq_screen"]:
        MULTIQC_INPUTS[FQ_SCREEN] = (
            expand(os.path.join(pre(FQ_SCREEN), "{sample}_R{read}_val_{read}_screen.txt"),
                   sample=samples, read=READS) +
            expand(os.path.join(pre(FQ_SCREEN), "{sample}_R{read}_val_{read}_screen.html"),
                   sample=samples, read=READS))

if config["star_align_params"]["run_star"]:
    MULTIQC_INPUTS[STAR_ALIGN] = (
        expand(os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam"), sample=samples) +
        expand(os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.bai"), sample=samples) +
        expand(os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.flagstat.txt"), sample=samples) +
        expand(os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.idxstats.txt"), sample=samples) +
        expand(os.path.join(pre(STAR_ALIGN), "{sample}_star_sorted.bam.stats.txt"), sample=samples))

if config["feature_count_params"]["run_feature_count"]:
    MULTIQC_INPUTS[FEATURE_COUNT] = [
        os.path.join(pre(FEATURE_COUNT), "gcounts_exon_s0.txt"),
        os.path.join(pre(FEATURE_COUNT), "gcounts_exon_s0.txt.summary")]

if config["rseq_params"]["run_rseqc"]:
    MULTIQC_INPUTS[COVERAGE_TIN] = (
        expand(os.path.join(pre(COVERAGE_TIN), "{sample}.geneBodyCoverage.txt"), sample=samples) +
        expand(os.path.join(pre(COVERAGE_TIN), "{sample}.geneBodyCoverage.curves.pdf"), sample=samples) +
        expand(os.path.join(pre(COVERAGE_TIN), "{sample}.geneBodyCoverage.r"), sample=samples) +
        [os.path.join(qc(COVERAGE_TIN), "tin_scores.png")])

if config["salmon_quant_params"]["run_salmon_quant"]:
    MULTIQC_INPUTS[SALMON_QUANT] = expand(
        os.path.join(pre(SALMON_QUANT), "{sample}_quant/quant.sf"), sample=samples)


# ── MultiQC rule ──────────────────────────────────────────────────────────
rule multiqc:
    input:
        lambda wc: MULTIQC_INPUTS[wc.step]
    output:
        html = os.path.join(QC_DIR, "{step}/multiqc_report.html")
    params:
        inputdir   = lambda wc: os.path.join(PREPROCESS_DIR, wc.step),
        outdir     = lambda wc: os.path.join(QC_DIR, wc.step),
        configfile = os.path.join(ENV_DIR, "multiqc_config.yml")
    conda:
        os.path.join(ENV_DIR, "multiqc124.yml")
    shell:
        r"""
        export TMPDIR={params.outdir}
        mkdir -p {params.outdir}
        multiqc \
            -c {params.configfile} \
            -o {params.outdir} {params.inputdir} \
            -s -d --interactive
        """
