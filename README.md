# Haplotype Methylation Snakemake Pipeline

![Snakemake](https://img.shields.io/badge/snakemake-≥7.x-brightgreen)
![Conda](https://img.shields.io/badge/conda-mambaforge-blue)
![Status](https://img.shields.io/badge/status-active-success)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

A Snakemake workflow for **haplotype-resolved methylation analysis**.

This workflow:

1. Splits methylation calls into **maternal** and **paternal** haplotypes and sorts them (rule split_haplotypes_and_sort)
2. Extract gene coordinates from files maternal_genes.bed and paternal_genes.bed (rule extract_genes_coordinates)
4. Calculates mean methylations for all genes. (rule map_methylation)
5. Normalization of the
4. Generates a **summary table**  
5. Uses **conda environments** automatically (no Docker required)

---

## 📋 Requirements

To run this workflow, you need:

- Linux or macOS (Windows via WSL2 works)
- **Conda or Mamba**
- **Snakemake ≥ 7.x**
- Internet access to download conda packages

We strongly recommend **MambaForge** for fast environment solves.

The inputs annotation beds for this pipeline (maternal.bed, paternal.bed) are available in [IS folder](https://is.muni.cz/auth/www/bulantova.l/pv269_project/). We need also 5mC methylation data from ONT: Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC, which are available in this [methylation foldes](https://public.gi.ucsc.edu/~mcechova/HG002/). You can download it manually, or first should do it for you. 

## 🟦 1. Install Conda / Mamba

### **Recommended: MambaForge**

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

Activate:

```bash
source ~/mambaforge/etc/profile.d/conda.sh
conda activate
```

---

## 🟩 2. Install Snakemake

### Option A — Mamba (recommended)
```bash
mamba create -n snakemake -c conda-forge snakemake
conda activate snakemake
```

### Option B — Conda
```bash
conda create -n snakemake -c conda-forge snakemake
conda activate snakemake
```

Check installation:

```bash
snakemake --version
```

---

## 📁 3. Project Structure

```
Snakefile.smk
envs/
    bedtools.yml
    python.yml
scripts/
    compute_mean_methylation.py
data/
    Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed
results/
    split/
    intersect/
    summarize/
logs/
```

---

## ▶️ 4. Running the Workflow
Snakemake doesn´t run whole pipeline, if outputs are already there, to rerun all add --forcerun all, or delete the results folder. 
### Dry-run (recommended)

```bash
snakemake -n --snakefile Snakefile.smk
```

### Run with 5 cores and conda support

```bash
snakemake --snakefile Snakefile.smk --cores 5 --use-conda
```

---

## 📜 5. Saving Logs

Create a log folder:

```bash
mkdir -p logs
```

Run and save a timestamped log:

```bash
snakemake --snakefile Snakefile.smk --cores 20 --use-conda &> logs/run_$(date +%Y%m%d_%H%M).log
```

---

## 📦 6. Conda Environments

Each Snakemake rule automatically creates and activates its own conda environment from:

```
envs/bedtools.yml
envs/python.yml
envs/stats.yml
```

To force a rebuild of all environments:

```bash
rm -r .snakemake/conda
```

---

## 🛠 Troubleshooting

### Snakemake cannot find the Snakefile
Use:

```bash
snakemake --snakefile Snakefile.smk
```

### Conda environment solving slow?
Install mamba:

```bash
conda install -n base -c conda-forge mamba
```

### Permission issues
Run the workflow inside a directory where you have full write access (for example, `$HOME`).

### Rebuild everything
```bash
snakemake --snakefile Snakefile.smk --cores 20 --use-conda --forcerun all
```

---
