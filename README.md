# scads
Single Cell ATAC-seq Disease Score

# Installing instructions
1. Install LDSC (via polyFUN)
 step 1: link download to folder : `~/temp/polyfun`. 
 step 2: test run polyfun (ldsc) 
 
2. Install scads

  # Install remotes package if not already installed
  install.packages("remotes")
  
  # Install scads package from your GitHub repository
  remotes::install_github("szhaolab/scads", ref = "v0.0.1")
 
# Run scads: see [tutorial](docs/articles/Introduction.html).
 
# Reference
 LDSC/polyFUN: See [software](https://github.com/omerwe/polyfun) and [paper](https://www.nature.com/articles/s41588-020-00735-5)
 
 fastTopics: See [method](https://github.com/stephenslab/fastTopics) and [paper](https://doi.org/10.1186/s13059-023-03067-9)
 
 