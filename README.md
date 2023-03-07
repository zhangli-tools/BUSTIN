# BUSTIN
BUSTIN (BUlk and Single-cell Transcriptome data INtegration) is an R package, which identifies cell subpopulations that primarily express the genes of your interest, e.g. the drug resistant genes, by integrating bulk and single-cell transcriptome data.

# Installation
```Rscript
BiocManager::install(c("Seurat","DESeq2","limma","cluster")) <br>
devtools::install_github("zhangli-tools/BUSTIN")<br>
```
# Usage
```Rscript
library(Seurat)<br>
library(BUSTIN)<br>
```
## Idnetify the drug resistant genes from bulk data <br>
```Rscript
resistant.genes=run.limma(bulk.data, pdata, resistant=T, padj=0.05, log2fc=0.5)
```
## if you use the bulk RNA-seq data <br>
```Rscript
resistant.genes=run.DESeq(bulk.data, pdata, resistant=T, padj=0.05, log2fc=0.5)
```
## Build a hierachical clustering tree for resistant genes based on scRNA-seq data <br>
```Rscript
hc.resistant.genes=hclust.geneset(seurat.obj = seurat.obj,maxK = 10,geneset = resistant.genes) <br>
```
## Predict the phenotype-associated cells<br>
```Rscript
BUSTIN.out=predict.PAC(seurat.obj = seurat.obj,k.out.list = list(hc.resistant.genes),minModuleSize = 15) <br>
```
## Identify the cell subpopulations associated with phenotype (drug resistance).<br>
```Rscript
cell.type=ORA.celltype(BUSTIN.out$seurat.obj,features = grep("m\[0-9\]\+$", colnames(BUSTIN.out$seurat.obj@meta.data), value=T),group.by = "label") <br>
```
