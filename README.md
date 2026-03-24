# epic_hD

**Source code for my PhD — Epigenetics and Intergenerational Chromatin**

Investigation of exposure-sensitive parental RNA–protein complexes as mechanistic links between environmental experience and chromatin states in the early embryo.

---

## Project overview

A parent's environment can modify phenotypes in offspring across generations. The abundance of several RNA classes in mature sperm and oocyte is altered by parental exposures, and injection of sperm RNA content into zygotes can recapitulate offspring phenotypes. Yet how these exposure-sensitive gamete RNAs influence the embryo remains unclear.

This project tests the hypothesis that **gamete-delivered RNAs modify the embryonic epigenome by forming ribonucleoprotein complexes (RNPs) with RNA-binding proteins (RBPs) that recruit chromatin-modifying enzymes to specific genomic loci** — providing a mechanistic route by which diverse environmental exposures converge on a shared chromatin-based logic to produce heritable phenotypic variation.

### Aims

| | Aim | Core approach |
|---|---|---|
| **1** | Establish candidate exposure-sensitive gamete RNAs and their predicted RBP partners | Integration of public and in-house sperm/oocyte RNA-seq and proteomics datasets; computational RNA–RBP interaction prediction |
| **2** | Determine whether RNA–RBP complexes form in the early embryo and associate with chromatin | smFISH in mouse zygotes; chromatin fractionation; ChIRP-MS |
| **3** | Test whether RNA–RBP complexes causally influence histone/DNA modifications in the preimplantation embryo | RNA/RNP microinjection; RNA depletion; low-input CUT&Tag |

---

## Repository structure

```
epic_hD/
├── src/                  Python package (import epic_hd) and R functions
├── scripts/
│   └── GRCm39/           Reference genome preparation scripts
│       ├── 00_download_gencode_m38.sh
│       ├── 01_build_annotation.sh
│       ├── 02_build_bed12.sh
│       ├── 10_prepare_salmon_inputs.sh
│       ├── 11_build_salmon_index.sh
│       └── 12_build_star_index.sh
├── environment.yml       Conda environment specification
└── pyproject.toml        Python package metadata
```

Data and results are stored separately on the cluster at `/mnt/thomas/epic_hD/` and are not tracked by git.

---

## Environment setup

```bash
# Clone the repo
git clone https://github.com/ThomasNgl/epic_hD.git
cd epic_hD

# Create the conda environment
conda env create --prefix /path/to/envs/epic_hd_env \
    --file environment.yml

# Activate
conda activate epic_hd_env
```

To update an existing environment after changes to `environment.yml`:

```bash
conda env update --prefix /path/to/envs/epic_hd_env \
    --file environment.yml --prune
```

---

## Reference genome setup (GRCm39 / GENCODE M38)

Scripts are numbered by stage. Run them in order from the `scripts/GRCm39/` directory.

```bash
cd scripts/GRCm39

# Download + verify + decompress + faidx
bash 00_download_gencode_m38.sh

# Gene annotation table (gene_id, symbol, biotype, coordinates)
bash 01_build_annotation.sh

# BED12 transcript structures (for RSeQC, deeptools)
bash 02_build_bed12.sh

# Genomic annotation BEDs (gene structure, CpG context, repeats, cCREs)
bash 03_build_genomic_annotations.sh

# Salmon: decoys + gentrome
bash 10_prepare_salmon_inputs.sh

# Salmon index (~15 min, ~8 GB RAM)
bash 11_build_salmon_index.sh --threads 10

# Salmon digest 
Rscript 12_build_salmon_digest.R

# STAR index (~45 min, ~32 GB RAM)
bash 20_build_star_index.sh --threads 10 --overhang 99
# wait for it to finish, then:
bash 20_build_star_index.sh --threads 10 --overhang 149
```

All indices are stored in `/mnt/auxiliary/` on the cluster and are shared across projects.

---

## Contact

Thomas Negrello — [negrello@hifo.uzh.ch](mailto:negrello@hifo.uzh.ch)  
Mansuy Lab, Brain Research Institute, UZH / ETH Zürich
