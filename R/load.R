#' @export

calculate.contribution=function(seurat.obj)
{
  cls=as.matrix(Seurat::Idents(seurat.obj))[,1]
  cls.uniq=unique(cls)
  count.mat=as.matrix(seurat.obj@assays$RNA@data)
  obs.prop=NULL
  for(i in cls.uniq)
  {
    cls.mat=matrix(0,nrow=length(cls),ncol=1)
    cls.mat[cls==i,1]=1
    obs.prop=cbind(obs.prop,count.mat%*%cls.mat)
  }
  colnames(obs.prop)=cls.uniq
  obs.prop=obs.prop/(obs.prop%*%matrix(1,nrow=ncol(obs.prop),ncol=ncol(obs.prop)))
  return(obs.prop)
}

enrich.test=function(geneset1,geneset2,bg)
{
  mat=matrix(c(length(intersect(geneset1,geneset2)),
               length(intersect(geneset1,setdiff(bg,geneset2))),
               length(intersect(setdiff(bg,geneset1),geneset2)),
               length(intersect(setdiff(bg,geneset2),setdiff(bg,geneset2)))),ncol=2)
  p=fisher.test(mat)$p.value
  return(p)
}



test.contribution=function(seurat.obj,geneset,bg=NULL)
{
  if(is.null(bg))
  {
    bg=rownames(seurat.obj@assays$RNA@data)
  }
  seurat.obj.rand=seurat.obj
  Seurat::Idents(seurat.obj.rand)=factor(as.matrix(Seurat::Idents(seurat.obj.rand))[sample(1:ncol(seurat.obj.rand),size = ncol(seurat.obj.rand),replace = F),1],levels=levels(Seurat::Idents(seurat.obj)))
  obs.prop=calculate.contribution(seurat.obj)
  rand.prop=calculate.contribution(seurat.obj.rand)
  rand.prop=rand.prop[apply(rand.prop,1,function(x){sum(is.na(x))})==0,colnames(obs.prop)]
  z=sapply(1:ncol(obs.prop),function(i){
    (obs.prop[,i]-median(rand.prop[,i]))/mad(rand.prop[,i])*0.6745
  })
  colnames(z)=colnames(obs.prop)
  p=1-pnorm(z)
  OR.genes=lapply(1:ncol(p),function(i){
    y=rownames(obs.prop)[which(p.adjust(p[,i],method="fdr")<0.05)]
    y1=names(which(apply(obs.prop[y,],1,which.max)==i))
  })
  names(OR.genes)=colnames(p)
  res=list(pvalue=sapply(OR.genes,function(x){enrich.test(x,geneset,bg)}),OR.genes=sapply(OR.genes,intersect,geneset))
}



identify.subpop=function(seurat.obj,geneset)
{
  seurat.obj=Seurat::FindVariableFeatures(seurat.obj,nfeatures = 2000,verbose =F)
  expr.normalized=as.matrix(seurat.obj@assays$RNA@data)
  geneset=intersect(rownames(expr.normalized),geneset)
  score=matrix(1,nrow=1,ncol=length(geneset))%*%expr.normalized[geneset,]
  cor.mat=cor(t(score),t(expr.normalized[setdiff(seurat.obj@assays$RNA@var.features,geneset),]),method="spearman")
  mk.list=list(pos=geneset,neg=names(head(sort(cor.mat[1,],decreasing = F),n=length(geneset)))) 
  anno=data.frame(Type=rep(c("pos","neg"),c(length(mk.list$pos),length(mk.list$neg))),
                    Marker=c(mk.list$pos,mk.list$neg),Weight=2)
  expr=seurat.obj@assays$RNA@counts
  topgenes = scSorter::xfindvariable_genes(expr, ngenes = 2000)
  expr = scSorter::xnormalize_scData(expr)
  topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
  topgenes = topgenes[topgene_filter]
  picked_genes = unique(c(anno$Marker, topgenes))
  expr = expr[rownames(expr) %in% picked_genes, ]
  rts <- scSorter::scSorter(expr, anno)
  Seurat::Idents(seurat.obj)=rts$Pred_Type
  return(seurat.obj)
}

test.contribution.m=function(seurat.obj,geneset,bg=NULL,sample="orig.ident")
{
  seurat.obj.list=Seurat::SplitObject(seurat.obj,split.by = sample)
  seurat.obj.out=lapply(seurat.obj.list,function(x){
    y=test.contribution(seurat.obj = x,geneset = geneset,bg = bg)
  })
  return(seurat.obj.out)
}

identify.subpop.m=function(seurat.obj,geneset,CellType,sample="orig.ident")
{
  seurat.obj=subset(seurat.obj,idents=CellType)
  seurat.obj.list=Seurat::SplitObject(seurat.obj,split.by = sample)
  seurat.obj.list=lapply(seurat.obj.list,function(x){
    identify.subpop(x,geneset = geneset)
  })
  cell.class=unlist(lapply(seurat.obj.list,function(x){clas=as.matrix(Seurat::Idents(x))[,1]}),use.names = F)
  names(cell.class)=unlist(lapply(seurat.obj.list,function(x){names(Seurat::Idents(x))}))
  cell.class=cell.class[names(Seurat::Idents(seurat.obj))]
  Idents(seurat.obj)=cell.class
  return(seurat.obj)
}



