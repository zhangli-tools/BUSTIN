# InSBut
InSBut is an R package, which identifies cell subpopulations that express the genes of your interest, e.g. the drug resistant genes, by integrating bulk and single-cell transcriptome data.

# Installation
devtools::install_github("zhangli-tools/InSBut")<br>
install.packages(c("Seurat","scSorter"))
# Usage
library(Seurat)<br>
library(scSorter)<br>
library(InSBut)<br>

## Calculate the contribution of the cell types to the overall gene expression.
cell.prop= calculate.contribution(seurat.obj)
## Identify the cell types primarily expressing the genes of interest.<br>
celltype= test.contribution(seurat.obj, resistant.genes) <br>
## Identify the cell subpopulations primarily expressing the genes of interest.<br>
subpop=identify.subpop(seurat.obj, resistant.genes) <br>
