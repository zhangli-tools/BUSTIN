#' @export
hclust.geneset=function(seurat.obj,geneset,maxK=10,p=0.05,percent=0.01){
  get_ave_sil_width <- function(d, cluster){
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("cluster package needed for this function to work. Please install it.")
    }
    ss <- cluster::silhouette(cluster, d)
    ss[,3]
  }
  geneset=intersect(geneset,rownames(seurat.obj@assays$RNA@data))
  expr.percent=sign(seurat.obj@assays$RNA@data[geneset,]) %*% matrix(1,nrow=ncol(seurat.obj@assays$RNA@data),ncol=1)
  expr.percent=expr.percent/ncol(seurat.obj@assays$RNA@data)
  geneset=geneset[as.matrix(expr.percent)[,1]>percent]
  geneset.expr=as.matrix(seurat.obj@assays$RNA@data[geneset,])
  scale.data=scale2zscore(geneset.expr)
  binary.data=t(apply(scale.data,1,function(x){ifelse(x>qnorm(1-p),1,0)}))
  dist.genes=binary.dist(binary.data)
  hc=hclust(dist(1-dist.genes),method = "ward.D")
  silhouette.index=lapply(2:maxK,function(i){
    cls=cutree(hc,k = i)
    y=get_ave_sil_width(dist.genes,cluster = cls)
  })
  names(silhouette.index)=2:maxK
  res=list(silhouette=silhouette.index,tree=hc,binary.data=binary.data,distance.mat=dist.genes)
  return(res)
}

binary.dist=function(mat){
  d1 = mat %*% t(mat)
  d2 = ncol(mat) - (1-mat) %*% t(1-mat)
  b.dist=1 - d1/d2
}

scale2zscore=function(mat){
  sum.bygene=mat %*% matrix(1,nrow=ncol(mat),ncol=1)
  mean.bygene=sum.bygene/ncol(mat)
  mean.mat=mean.bygene %*% matrix(1,nrow=1,ncol=ncol(mat))
  sd.mat=matrix(apply(mat,1,sd),ncol=1) %*% matrix(1,nrow=1,ncol=ncol(mat))
  scaledata=(mat-mean.mat)/sd.mat
  scaledata
}

get.gene.modules=function(k.out,minModuleSize=30,median.sihouette.cutoff=0.01){
  maxK=length(k.out$silhouette)+1
  class.mat=sapply(2:maxK,function(i){
    cls=cutree(k.out$tree,k=i)
    pval=tapply(k.out$silhouette[[i-1]],cls,function(x){wilcox.test(x,alternative="greater")$p.value})
    median.sihouette=tapply(k.out$silhouette[[i-1]],cls,median)
    cls[cls %in% names(which(pval>0.05 | median.sihouette<median.sihouette.cutoff)) | cls %in% names(which(table(cls)<minModuleSize))]=0
    cls
  })

  best.cls=class.mat[,1]
  for(i in 3:maxK){
    if(length(unique(class.mat[,i-1]))>length(unique(best.cls)))
    {
      best.cls=class.mat[,i-1]
    }else if(length(unique(class.mat[class.mat[,i-1]>0,i-1]))==length(unique(best.cls[best.cls>0]))) {
      next.cls=class.mat[,i-1]
      exclude.cls=unique(best.cls[next.cls==0 & best.cls>0])
      if(sum(best.cls %in% exclude.cls & next.cls>0)>0){
        best.cls=next.cls
      }
    }
  }
  if(sum(best.cls)==0){
    class.mat=sapply(2:maxK,function(i){
      cls=cutree(k.out$tree,k=i)
      pval=tapply(k.out$silhouette[[i-1]],cls,function(x){wilcox.test(x,alternative="greater")$p.value})
      median.sihouette=tapply(k.out$silhouette[[i-1]],cls,median)
      cls[cls %in% names(which(pval>0.05 | median.sihouette<median.sihouette.cutoff)) ]=0
      cls
    })
    max.cluster.num=apply(class.mat,2,function(cls){
      max(table(cls[cls>0]))
    })
    message(paste("warning: the maximal module size =",max(max.cluster.num),"< minModuleSize you provided!\nPACS returns the module with the largest size"))
    best.cls=class.mat[,which(max.cluster.num==max(max.cluster.num))[1]]
  }

  best.cls
}

predict.PAC=function(seurat.obj,k.out.list,p=0.05,prob=0.95,minModuleSize=30,weight=NULL,median.sihouette.cutoff=0.01){

  gene.modules.list=lapply(k.out.list,function(k.out){
    gene.modules=get.gene.modules(k.out,minModuleSize = minModuleSize,median.sihouette.cutoff  = median.sihouette.cutoff)
    gene.modules=gene.modules[gene.modules>0]
    genes=names(gene.modules)
    gene.modules=match(gene.modules,sort(unique(gene.modules)))
    names(gene.modules)=genes
    gene.modules
  })
  cell.prob.list=list()
  cell.score.list=list()
  weight=abs(weight)/sum(abs(weight))*length(abs(weight))
  for(i in 1:length(gene.modules.list)){
    gene.modules=gene.modules.list[[i]]
    geneset=names(gene.modules)
    binary.data=k.out.list[[i]]$binary.data
    gene.modules.weight.scale=tapply(geneset,gene.modules,function(x){weight[x]})
    genecount.list=lapply(gene.modules.weight.scale,function(x){
      tmp.binary=binary.data[names(x),]
      tmp.binary=matrix(x,nrow=1) %*% tmp.binary
    })
    genecount=eval(as.call(c(rbind,genecount.list)))
    rand.genecount.list=lapply(gene.modules.weight.scale,function(x){
      tmp.binary=t(apply(binary.data[names(x),],1,function(x){y=x[sample(1:length(x),length(x),replace = F)]}))
      tmp.binary=matrix(x,nrow=1) %*% tmp.binary
    })
    rand.genecount=eval(as.call(c(rbind,rand.genecount.list)))


    if(is.null(nrow(genecount))){genecount=t(as.matrix(genecount))}
    cell.prob=sapply(1:nrow(genecount),function(j){
      x=genecount[j,]
      x1=rand.genecount[j,]
      zscore=(x-mean(x1))/sd(x1)
      p=pnorm(zscore)
    })
    cell.score=sapply(1:nrow(genecount),function(j){
      x=genecount[j,]
      x1=rand.genecount[j,]
      zscore=(x-mean(x1))/sd(x1)
    })
    if(is.null(names(k.out.list))){
      colnames(cell.prob)=paste0("gs",i,".","m",1:ncol(cell.prob))}else{
        colnames(cell.prob)=paste0(names(k.out.list)[i],".","m",1:ncol(cell.prob))
      }
    colnames(cell.score)=colnames(cell.prob)
    cell.prob.list[i]=list(cell.prob)
    cell.score.list[i]=list(cell.score)
  }
  cell.prob.mat=eval(as.call(c(cbind,cell.prob.list)))
  cell.score.mat=eval(as.call(c(cbind,cell.score.list)))
  if(length(prob)==1){
    prob=rep(prob,ncol(cell.score.mat))
  }else{
    prob=rep(prob,sapply(gene.modules.list,function(x){length(unique(x))}))
  }
  cell.status=t(sapply(1:nrow(cell.score.mat),function(i){
    y=rep(0,ncol(cell.score.mat))
    prob.x=cell.prob.mat[i,]
    score.x=cell.score.mat[i,]
    max.idx=which.max(score.x)
    test=score.x==max(score.x) & prob.x>prob[max.idx]
    if(sum(test)>1){
      y[which.max(score.x)]=1
    }else{
      y[score.x==max(score.x) & prob.x>prob[max.idx]]=1
    }

    y
  }))
  colnames(cell.status)=colnames(cell.score.mat)
  colnames(cell.score.mat)=paste0(colnames(cell.score.mat),".score")
  seurat.obj@meta.data=cbind(seurat.obj@meta.data,cell.status,cell.score.mat)
  return(list(seurat.obj=seurat.obj,gene.modules.list=gene.modules.list))
}

ORA.celltype=function(seurat.obj,features=c("M1"),group.by=NULL)
{
  if(is.null(group.by)){
    celltype=as.matrix(Seurat::Idents(seurat.obj))[,1]
  }else{
    celltype=as.matrix(seurat.obj[[group.by]])[,1]
  }
  if(max(table(celltype)/length(celltype))>0.5)
  {
    stop("The maximal cluster accounted for 50% of total cells, please cluster with a higher resoultion!")
  }

  db.data=data.frame(celltype=celltype,cellid=names(Seurat::Idents(seurat.obj)))
  df.list=lapply(features,function(x){
    y=data.frame(cellid=names(which(as.matrix(seurat.obj[[x]])[,1]>0)),group=x)
  })
  enrich.tab=lapply(features,function(x){
    Module=x
    tab=table(as.matrix(seurat.obj[[x]])[,1],celltype)
    total=apply(tab,1,sum)
    pvalue=apply(tab,2,function(x){
      prop.test(c(x[2],sum(x)),c(sum(x),sum(total)),alternative = "greater")$p.value
    })
    CellRatio=apply(tab,2,function(x){
      y=paste(x[2],sum(x),sep="/")
    })
    BgRatio=apply(tab,2,function(x){
      y=paste(sum(x),sum(tab),sep="/")
    })
    Count=apply(tab,2,function(x){
      y=x[2]
    })
    res=data.frame(Module,Cell=colnames(tab),CellRatio,BgRatio,pvalue,Count)
  })
  enrich.tab=eval(as.call(c(rbind,enrich.tab)))
  enrich.tab$padj=p.adjust(enrich.tab$pvalue,method="fdr")
  rownames(enrich.tab)=NULL
  return(enrich.tab)
}

run.limma=function(bulk.data, pdata, resistant=T, padj=0.05, log2fc=0.5)
{
  design.group=cbind(yes=pdata,no=1-pdata)
  group.fit = limma::lmFit(bulk.data, design.group)
  group.fit = limma::eBayes(group.fit)
  group.contrast.matrix = limma::makeContrasts(CancervNormal =yes -
                                            no, levels = design.group)
  group.fit2 = limma::contrasts.fit(group.fit, group.contrast.matrix)
  group.fit2 = limma::eBayes(group.fit2)
  group.results = limma::topTable(group.fit2, number = nrow(bulk.data), sort.by = "p",
                             adjust.method = "BH")
  if(resistant)
  {
    res=rownames(group.results)[which(group.results$adj.P.Val<padj & group.results$logFC>log2fc)]
  }else{
    res=rownames(group.results)[which(group.results$adj.P.Val<padj & group.results$logFC< -log2fc)]
  }
  res
}

run.DESeq=function(bulk.data, pdata, resistant=T, padj=0.05, log2fc=0.5)
{
  dds=DESeq2::DESeqDataSetFromMatrix(bulk.data,colData = data.frame(g=pdata),design = ~g)
  dds=DESeq2::DESeq(dds)
  group.results = DESeq2::results(dds,contrast = c("g",1,0))
  if(resistant)
  {
    res=rownames(group.results)[which(group.results$padj<padj & group.results$log2FoldChange>log2fc)]
  }else{
    res=rownames(group.results)[which(group.results$padj<padj & group.results$log2FoldChange< -log2fc)]
  }
  res
}

get.BUSTIN.obj=function(seurat.obj,cell.type)
{
  labels=seurat.obj$label
  for(i in 1:nrow(cell.type)){
    id=as.matrix(cell.type$Module[i])[1]
    labels[seurat.obj[[id]]==1 & labels==cell.type$Cell[i]]=paste(cell.type$Cell[i],cell.type$Module[i],sep=".")
  }
  seurat.obj$BUSTIN.label=labels
  seurat.obj
}
