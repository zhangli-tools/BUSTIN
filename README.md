# BUSTIN
BUSTIN (Bulk and Single-cell transcriptome data integration ) is an R package, which identifies cell subpopulations that primarily express the genes of your interest, e.g. the drug resistant genes, by integrating bulk and single-cell transcriptome data.

# Installation
devtools::install_github("zhangli-tools/BUSTIN")<br>
install.packages(c("Seurat","scSorter"))
# Usage
library(Seurat)<br>
library(scSorter)<br>
library(InSBut)<br>
## Idnetify the drug resistant genes from bulk data
resistant.genes=run.limma(bulk.data, pdata, resistant=T, padj=0.05, log2fc=0.5)
## if you use the bulk RNA-seq data
resistant.genes=run.DESeq(bulk.data, pdata, resistant=T, padj=0.05, log2fc=0.5)
## Calculate the contribution of the cell types to the overall gene expression.
cell.prop= calculate.contribution(seurat.obj)
## Identify the cell types primarily expressing the genes of interest.<br>
celltype= test.contribution(seurat.obj, resistant.genes) <br>
## Identify the cell subpopulations primarily expressing the genes of interest.<br>
subpop=identify.subpop(seurat.obj, resistant.genes) <br>
