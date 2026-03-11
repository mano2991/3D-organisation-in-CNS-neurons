# 3D Chromatin Organization in CNS Neurons
### Hi-C Data Analysis Across Developmental and Injury States

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Genome: mm10](https://img.shields.io/badge/Genome-mm10-blue.svg)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/)
[![Tools: HOMER](https://img.shields.io/badge/Tools-HOMER-green.svg)](http://homer.ucsd.edu/homer/)

---

## Overview

This repository contains a comprehensive computational pipeline for analyzing **three-dimensional chromatin organization** in central nervous system (CNS) neurons across key developmental timepoints and pathological conditions. Using Hi-C sequencing data, the project investigates how the spatial architecture of the genome is established, maintained, and altered in neurons at:

| Stage | Description |
|---|---|
| **E12.5** | Embryonic day 12.5 — early neurogenesis |
| **P0** | Postnatal day 0 — perinatal transition |
| **Adult (Uninjured)** | Mature neuronal state |
| **Adult (Injured / 7dpi)** | 7 days post-injury — injury response |

The mammalian genome is non-randomly organized within the nucleus, and this three-dimensional architecture plays a fundamental role in regulating gene expression. In CNS neurons, dynamic changes in chromatin conformation are essential for orchestrating the transcriptional programs that govern **neuronal differentiation, maturation, and responses to damage**. This pipeline provides a multi-scale, systematic framework to investigate these changes.

---

## Scientific Background

Higher-order chromatin organization is analyzed at three hierarchical scales:

```
Megabase scale      →   A/B Compartments   (active vs. inactive chromatin)
Sub-megabase scale  →   TADs               (topologically associating domains)
Kilobase scale      →   Chromatin Loops    (enhancer–promoter contacts)
```

### A/B Compartments
At the megabase scale, **A/B compartment analysis** distinguishes transcriptionally active (A) and inactive (B) chromatin domains, enabling genome-wide tracking of **compartment switching events** across developmental transitions and injury.

### Topologically Associating Domains (TADs)
At the sub-megabase scale, **TAD analysis** identifies self-interacting chromatin domains and characterizes TAD boundary dynamics — revealing how **insulation landscapes are remodeled** during neuronal maturation.

### Chromatin Loops
At the kilobase scale, **chromatin loop analysis** maps specific long-range enhancer–promoter and promoter–promoter interactions, providing mechanistic insight into **gene regulatory rewiring** across conditions.

---

## Repository Structure

```
.
├── run_hic_pipeline.sh       # Main automated pipeline script
├── Tad_ann_CLI.R             # TAD & Loop annotation and comparison (R script)
├── data.txt                  # Sample sheet (name | comparison group | R1 | R2)
├── gene_split.csv            # Gene annotation reference for R script
├── README.md                 # This file
│
├── HiC_Results/              # Generated output directory
│   ├── SAM/                  # Alignment files
│   ├── TagDirs/              # HOMER tag directories (one per sample)
│   │   └── <SampleName>/
│   │       ├── *.tad.2D.bed
│   │       ├── *.loop.2D.bed
│   │       ├── *.PC1.bedGraph
│   │       └── *.hic          # (if Juicer export enabled)
│   └── Comparisons/          # Pairwise comparison outputs
│       └── <SampleA_vs_SampleB>/
│           ├── merged_TADs.bed
│           ├── merged_Loops.bed
│           ├── TADLoop_score.*
│           ├── Chromatin_compartment_*.txt
│           └── R_annotation/
│               ├── Loop/
│               │   ├── Loop_comparison_final.csv
│               │   └── Loop_Summary_Plot.png
│               └── TAD/
│                   ├── TAD_comparison_final.csv
│                   └── TAD_Summary_Plot.png
```

---

## Dependencies

### Required Tools

| Tool | Version | Purpose |
|---|---|---|
| [HOMER](http://homer.ucsd.edu/homer/) | ≥ 4.11 | Hi-C processing, TAD/Loop calling, PCA |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | ≥ 2.4 | Read alignment |
| [R](https://www.r-project.org/) | ≥ 4.1 | TAD/Loop annotation and comparison |
| [Java](https://www.java.com/) | ≥ 8 | Juicer tools (.hic export) |

### R Packages

```r
install.packages(c("dplyr", "ggplot2", "tidyr"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```

### Optional

| Tool | Purpose |
|---|---|
| [Juicer Tools](https://github.com/aidenlab/juicer) | Export `.hic` files for Juicebox visualization |

---

## Installation

```bash
# 1. Clone this repository
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>

# 2. Install HOMER (if not already installed)
# See: http://homer.ucsd.edu/homer/introduction/install.html

# 3. Download mm10 genome and Bowtie2 index
# HOMER genome:
perl /path/to/homer/configureHomer.pl -install mm10

# Bowtie2 index (UCSC):
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
unzip mm10.zip

# 4. Install R packages (see above)
```

---

## Quick Start

### 1. Prepare your sample sheet (`data.txt`)

Tab-separated file with **no header**. Samples sharing the same **Comparison Group** will be compared pairwise.

```
# SampleName    ComparisonGroup   R1_fastq                     R2_fastq
HiC-E12.5       E12vsP0           /data/E12_R1.fastq.gz        /data/E12_R2.fastq.gz
HiC-P0          E12vsP0           /data/P0_R1.fastq.gz         /data/P0_R2.fastq.gz
HiC-Uninj       AdultvsInjured    /data/Uninj_R1.fastq.gz      /data/Uninj_R2.fastq.gz
HiC-7dpi        AdultvsInjured    /data/7dpi_R1.fastq.gz       /data/7dpi_R2.fastq.gz
```

> Any samples sharing a comparison group are automatically compared against each other. Three samples in one group will produce 3 pairwise comparisons.

### 2. Run the full pipeline

```bash
bash run_hic_pipeline.sh \
  -d data.txt \
  -b /path/to/mm10_Bowtie2_index \
  -g mm10 \
  -t 15 \
  -r gene_split.csv \
  -o HiC_Results \
  -j /path/to/juicer_tools.jar    # optional
```

### 3. Run annotation only (if tag directories already exist)

```bash
bash run_hic_pipeline.sh \
  -d data.txt \
  -b /path/to/mm10_Bowtie2_index \
  -s \          # skip trimming and alignment
  -o HiC_Results
```

---

## Pipeline Options

| Flag | Description | Default |
|---|---|---|
| `-d` | Sample sheet (required) | — |
| `-b` | Bowtie2 index prefix (required) | — |
| `-g` | Genome build | `mm10` |
| `-t` | CPU threads | `15` |
| `-r` | Gene annotation CSV for R script | `gene_split.csv` |
| `-o` | Output root directory | `HiC_Results` |
| `-j` | Juicer tools JAR path (enables `.hic` export) | disabled |
| `-s` | Skip trimming + alignment | `false` |
| `-h` | Show help | — |

---

## Pipeline Steps in Detail

```
Step 1  ──  Adapter Trimming        homerTools trim  (-3 GATC)
Step 2  ──  Read Alignment          bowtie2          (single-end, per read)
Step 3  ──  Tag Directory           makeTagDirectory (HiC-specific filters)
Step 4  ──  PCA / Compartments      runHiCpca.pl     (50kb res, 100kb window)
Step 5  ──  Compaction Stats        analyzeHiC
Step 6  ──  TAD Calling             findTADsAndLoops.pl find
Step 7  ──  Loop Calling            findTADsAndLoops.pl find
Step 8  ──  [Optional] .hic Export  tagDir2hicFile.pl + juicer_tools
Step 9  ──  Merge TAD/Loop beds     merge2Dbed.pl    (across comparison group)
Step 10 ──  Score TADs & Loops      findTADsAndLoops.pl score
Step 11 ──  Compartment Compare     annotatePeaks.pl (PC1 bedGraph comparison)
Step 12 ──  R Annotation            Tad_ann_CLI.R    (TAD + Loop, both modes)
```

---

## `Tad_ann_CLI.R` — TAD & Loop Annotation Script

This custom R script performs pairwise comparison and gene annotation of TADs and chromatin loops.

### Features
- Classifies TADs/Loops as **Conserved**, **Shifted**, or **Condition-Specific**
- Uses a **greedy 1-to-1 matching** algorithm based on genomic overlap (≥40% → candidate; ≥80% → conserved)
- Annotates anchors with overlapping gene names using `GenomicRanges`
- Generates summary bar plots and annotated CSV output

### Direct usage

```bash
Rscript Tad_ann_CLI.R \
  --mode  both \
  --n1    HiC-P0 \
  --n2    HiC-7dpi \
  --t1    HiC-P0/HiC-P0.tad.2D.bed \
  --t2    HiC-7dpi/HiC-7dpi.tad.2D.bed \
  --l1    HiC-P0/HiC-P0.loop.2D.bed \
  --l2    HiC-7dpi/HiC-7dpi.loop.2D.bed \
  --gene  gene_split.csv \
  --out   P0_vs_7dpi_comparison
```

### Arguments

| Flag | Description |
|---|---|
| `--mode` | `tad`, `loop`, or `both` |
| `--n1` / `--n2` | Sample names (used in output labels) |
| `--t1` / `--t2` | TAD `.2D.bed` files |
| `--l1` / `--l2` | Loop `.2D.bed` files |
| `--gene` | Gene annotation CSV (`Chromosome`, `Start`, `End`, `Gene_name`) |
| `--out` | Output directory |

### Output files

```
<out>/
├── Loop/
│   ├── Loop_comparison_final.csv    # All loops with status, overlap %, gene annotations
│   └── Loop_Summary_Plot.png        # Bar chart: Conserved / Shifted / Specific
└── TAD/
    ├── TAD_comparison_final.csv
    └── TAD_Summary_Plot.png
```

---

## `gene_split.csv` Format

The gene annotation file must contain these columns:

```
Chromosome,Start,End,Gene_name
chr1,3205901,3671498,Xkr4
chr1,4807788,4848410,Gm37586
...
```

> You can generate this file from Ensembl BioMart or UCSC Table Browser for mm10.

---

## Key Parameters

| Parameter | Value | Description |
|---|---|---|
| Restriction enzyme | GATC (MboI) | Hi-C ligation site |
| PCA resolution | 50 kb | A/B compartment calling |
| PCA window | 100 kb | Smoothing window |
| TAD resolution | 25 kb | TAD boundary calling |
| Loop resolution | 5 kb | Loop anchor resolution |
| TAD overlap threshold | 40% / 80% | Shifted / Conserved cutoff |

---

## Citation

If you use this pipeline in your research, please cite:

> [Your manuscript title]. *Journal*. Year. doi:XXX

And the underlying tools:

- **HOMER**: Heinz S et al. (2010) *Molecular Cell*. doi:10.1016/j.molcel.2010.05.004
- **Bowtie2**: Langmead B & Salzberg SL (2012) *Nature Methods*. doi:10.1038/nmeth.1923
- **Juicer**: Durand NC et al. (2016) *Cell Systems*. doi:10.1016/j.cels.2016.07.002

---

## Contact

For questions or issues, please open a [GitHub Issue](../../issues) or contact:

**[Your Name]** — [your.email@institution.edu]  
[Your Lab / Institution]

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
