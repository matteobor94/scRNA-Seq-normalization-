library(DESeq2)
library(scran)
library(edgeR)
library(Linnorm)
library(Seurat)
NUM_OF_THREAD=8

raw_count = function(sce){
  tp = system.time({
    logcounts(sce) = counts(sce)
  })
  
  method_name = "raw_count"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

scran_norm = function(sce){
  tp = system.time({
    sce = computeSumFactors(sce)
    scra = normalizeData(sce) # create a normalized matrix with the sum factors
  })
  assay(sce,"Scran")<-scra
  
  method_name = "scran"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#PARETO
#assay(sce,"Pareto")<-0*counts(sce)
pareto.MLE <- function(X)
{
  n <- length(X)
  m <- min(X)
  a <- n/sum(log(X)-log(m))
  return(c(m,a)) 
}
pareto_norm<-function(sce){
  m<-counts(sce)*0
  c<-counts(sce)
  tp = system.time({
    alfa<-apply(c+1, 2, pareto.MLE)[2,]
    for (i in 1:length(alfa)) {
      m[,i]<-log2(c[,i]*alfa[i]+1)
    }
    assay(sce,"Pareto")<-m
  })
  method_name = "pareto"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#CLR
clr_norm<-function(sce){
  c<-counts(sce)
  clr<-c*0
  tp=system.time({
    clr<-NormalizeData(c, margin<-2, normalization:method<-"CLR")
  })
  assay(sce, "CLR")<-clr  
  method_name = "CLR"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#DESeq
DESeq2_norm = function(sce){
c<-counts(sce)
  tp = system.time({
    sizeFactors(sce) <- estimateSizeFactorsForMatrix(c)
    des <- normalizeData(sce)
    assay(sce,"DESeq2")<-des
  })
  
  method_name = "DESeq2"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#TMM
TMM_norm = function(sce){
  c<-counts(sce)
  tp = system.time({
    sizeFactors(sce) <- calcNormFactors(c, method = "TMM")
    tmm <- normalizeData(sce)
    assay(sce, "TMM")<-tmm
  })
  
  method_name = "TMM"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#logCPM
logCPM_norm = function(sce){
c<-counts(sce)
  tp = system.time({
    assay(sce, "logCPM") = log2(edgeR::cpm(c) + 1)
  })
  
  method_name = "logCPM"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#LINNORM
linnorm_norm = function(sce){
c<-counts(sce)
  tp = system.time({
    assay(sce, "Linnorm") = Linnorm(c)
  })
  
  method_name = "Linnorm"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

#SCT
SCT_norm<-function(sce){
  c<-counts(sce)
  tp = system.time({
    pbmc <- CreateSeuratObject(counts = c)
    pbmc <- SCTransform(object = pbmc, verbose = FALSE)
    assay(sce,"SCT")<-as.matrix(pbmc@assays$SCT@data)
  })
  method_name = "SCT"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time,
                                       data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,
                                            method_type=method_type, time=unname(tp)[1])
  }
  return(sce)
}

norm_method <- list(
  scran = scran_norm,
  pareto=pareto_norm,
  DESeq2=DESeq2_norm,
  TMM=TMM_norm,
  logCPM=logCPM_norm,
  Linnorm=linnorm_norm,
  CLR=clr_norm,
  SCT=SCT_norm
)
