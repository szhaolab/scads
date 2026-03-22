# scads (Single Cell ATAC-seq Disease Score)  <img src="man/figures/scads_logo_gemini.png" align="right" height="139" style="float:right; height:139px;" />

The goal of `scads` is to provide a statistical framework for integrating GWAS and scATAC-seq data for quantifying genetic disease risk in single cells.

## Prerequisites

- **R** (>= 4.0)
- **Python 3** (>= 3.8)
- **Conda** (Miniconda or Anaconda) — needed to manage Python dependencies
- **Git**

## Installation

### Step 1: Install Miniconda (if not already installed)

Skip this step if you already have conda available.

**macOS (Apple Silicon):**
```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init zsh   # or bash, depending on your shell
```

**macOS (Intel):**
```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init zsh
```

**Linux:**
```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
```

Restart your shell after running `conda init`.

### Step 2: Install polyFUN and LDSC

Clone both repositories:

```bash
git clone https://github.com/omerwe/polyfun.git ~/polyfun
git clone https://github.com/bulik/ldsc.git ~/ldsc
```

### Step 3: Create the Python environment

Create a conda environment with the required Python dependencies:

```bash
conda create -n polyfun python=3.11 -y
conda activate polyfun
pip install numpy scipy pandas scikit-learn pyarrow tqdm bitarray networkx pandas-plink packaging pybedtools
```

Verify the installation:

```bash
python ~/polyfun/ldsc.py -h
python ~/ldsc/make_annot.py -h
```

Both commands should print usage information without errors.

### Step 4: Install Bioconductor dependencies

`scads` depends on several Bioconductor packages that are not available on CRAN. Install them first from within R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "SummarizedExperiment",
  "rtracklayer",
  "plyranges",
  "IRanges",
  "GenomicRanges",
  "BSgenome",
  "Biostrings"
))
```

> **Note:** If `plyranges` fails to install due to a `vctrs` version conflict, install `vctrs` from source first:
> ```r
> install.packages("vctrs", type = "source")
> ```
> Then retry `BiocManager::install("plyranges")`.

### Step 5: Install scads

**Option A:** Install directly from GitHub (recommended):

```r
install.packages("remotes")
remotes::install_github("szhaolab/scads")
```

**Option B:** Install from a local clone:

```bash
git clone https://github.com/szhaolab/scads.git
```

```r
remotes::install_local("path/to/scads")
```

### Step 6: Download S-LDSC reference data

`scads` requires reference files from the S-LDSC baseline LD v2.2 model (hg19). Download and extract them into a local directory (e.g., `~/LDSCORE/`).

> **Disk usage:** The compressed downloads total ~1 GB. After extraction, the reference data requires ~4 GB.

```bash
mkdir -p ~/LDSCORE && cd ~/LDSCORE

# Download from Zenodo (~1 GB total)
curl -fSL -o 1000G_Phase3_plinkfiles.tgz \
  "https://zenodo.org/records/10515792/files/1000G_Phase3_plinkfiles.tgz?download=1"
curl -fSL -o 1000G_Phase3_baselineLD_v2.2_ldscores.tgz \
  "https://zenodo.org/records/10515792/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz?download=1"
curl -fSL -o 1000G_Phase3_frq.tgz \
  "https://zenodo.org/records/10515792/files/1000G_Phase3_frq.tgz?download=1"
curl -fSL -o 1000G_Phase3_weights_hm3_no_MHC.tgz \
  "https://zenodo.org/records/10515792/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1"
curl -fSL -o hm3_no_MHC.list.txt \
  "https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt?download=1"

# Extract
tar -xzf 1000G_Phase3_plinkfiles.tgz
tar -xzf 1000G_Phase3_frq.tgz
tar -xzf 1000G_Phase3_weights_hm3_no_MHC.tgz

# Baseline LD extracts flat files, so create a subdirectory first
mkdir -p 1000G_Phase3_baselineLD_v2.2_ldscores
tar -xzf 1000G_Phase3_baselineLD_v2.2_ldscores.tgz -C 1000G_Phase3_baselineLD_v2.2_ldscores
```

### Step 7: Verify the installation

```r
library(scads)
data("demo_count_matrix")
dim(demo_count_matrix)
```

### Running scads

When running `scads`, you must activate the `polyfun` conda environment **before** starting R, so that `python3` resolves to the environment with the required dependencies:

```bash
conda activate polyfun
R
```

In your R session, point `scads` to the polyFUN and LDSC directories:

```r
scads_res <- scads(
  ...,
  polyfun_code_dir = "~/polyfun/",
  ldsc_code_dir    = "~/ldsc/",
  ...
)
```

See the [tutorial](https://szhaolab.github.io/scads/articles/Introduction.html) for a complete example.

# Reference

 LDSC/polyFUN: See [software](https://github.com/omerwe/polyfun) and [paper](https://www.nature.com/articles/s41588-020-00735-5)
 
 fastTopics: See [method](https://github.com/stephenslab/fastTopics) and [paper](https://doi.org/10.1186/s13059-023-03067-9)
 
 
