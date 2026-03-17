# scads (Single Cell ATAC-seq Disease Score)  <img src="man/figures/scads_logo_gemini.png" align="right" height="139" style="float:right; height:139px;" />

The goal of `scads` is to provide a statistical framework for integrating GWAS and scATAC-seq data for quantifying genetic disease risk in single cells.

## Installing instructions

1. Install LDSC (via polyFUN)

  - step 1: link download to folder : `~/polyfun`. 
  
  - step 2: test run polyfun (ldsc)
 
2. Install scads

  - Clone repository 
  
```
git clone https://github.com/szhaolab/scads
```
  - OR use below to install via R `remotes`

```
  # Install remotes package if not already installed
  install.packages("remotes")
  
  # Install scads package from your GitHub repository
  remotes::install_github("szhaolab/scads", ref = "v0.0.1")
```
 
# How to run scads: 

  - see [tutorial](https://github.com/szhaolab/scads/articles/Introduction.html).
 
# Reference

 LDSC/polyFUN: See [software](https://github.com/omerwe/polyfun) and [paper](https://www.nature.com/articles/s41588-020-00735-5)
 
 fastTopics: See [method](https://github.com/stephenslab/fastTopics) and [paper](https://doi.org/10.1186/s13059-023-03067-9)
 
 
