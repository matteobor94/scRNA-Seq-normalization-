library(scran)
library(scater)
library(CellBench)
library(DrImpute)
library(SAVER)
set_cellbench_threads(nthreads = 1)
load("C:/Users/Matteo/Desktop/Borsa studio/sincell_with_class.rdata") #my personal directory
load("C:/Users/Matteo/Desktop/Borsa studio/sincell_with_class_5cl.rdata")
#the datasets are available at the link:
# https://github.com/LuyiTian/sc_mixology/tree/master/data
datasets <- list(
  sc_CELseq2=sce_sc_CELseq2_qc,
  sc_10x=sce_sc_10x_qc#, #trying with only these 2 datasets
  # sc_Droseq=sce_sc_Dropseq_qc,
  # sc_10x_5cl=sce_sc_10x_5cl_qc,
  # sc_Celseq2_5cl_p1=sc_Celseq2_5cl_p1,
  # sc_Celseq2_5cl_p2=sc_Celseq2_5cl_p2,
  # sc_Celseq2_5cl_p3=sc_Celseq2_5cl_p3
)
remove(sc_CELseq2=sce_sc_CELseq2_qc,
       sc_10x=sce_sc_10x_qc,
       sc_Droseq=sce_sc_Dropseq_qc,
       sc_10x_5cl=sce_sc_10x_5cl_qc,
       sc_Celseq2_5cl_p1=sc_Celseq2_5cl_p1,
       sc_Celseq2_5cl_p2=sc_Celseq2_5cl_p2,
       sc_Celseq2_5cl_p3=sc_Celseq2_5cl_p3)

library(DESeq2)
library(scran)
library(edgeR)
library(Linnorm)
library(SCnorm)
library(scone)
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
    sce = normalize(sce) # goes to `logcounts` by default
    assay(sce,"Scran")<-logcounts(sce)
  })
  
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

DESeq2_norm = function(sce){
  tp = system.time({
    sizeFactors(sce) <- estimateSizeFactorsForMatrix(counts(sce))
    sce <- normalize(sce)
    assay(sce,"DESeq2")<-logcounts(sce)
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

TMM_norm = function(sce){
  tp = system.time({
    sizeFactors(sce) <- calcNormFactors(counts(sce), method = "TMM")
    sce <- normalize(sce)
    assay(sce, "TMM")<-logcounts(sce)
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

logCPM_norm = function(sce){
  tp = system.time({
    assay(sce, "logCPM") = log2(edgeR::cpm(counts(sce)) + 1)
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

linnorm_norm = function(sce){
  tp = system.time({
    assay(sce, "Linnorm") = Linnorm(counts(sce))
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

SCnorm_norm = function(sce){
  tp = system.time({
    SCnorm_out = SCnorm(Data=counts(sce),Conditions = rep(1,ncol(sce)),FilterCellNum = 10, NCores=NUM_OF_THREAD)
    assay(sce, "SCnorm") = log2(normcounts(SCnorm_out)+1)
  })
  
  method_name = "SCnorm"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}

scone_norm = function(sce){
  tp = system.time({
    scaling=list(none=identity, # Identity - do nothing
                 sum = SUM_FN, # SCONE library wrappers...
                 tmm = TMM_FN, 
                 uq = UQ_FN,
                 fq = FQT_FN,
                 deseq = DESEQ_FN)
    results = scone(SconeExperiment(counts(sce)), 
                    scaling=scaling,
                    run=TRUE, k_qc=0, k_ruv=0,
                    return_norm = "in_memory",
                    zero = "postadjust",
                    bpparam = BiocParallel::SerialParam())
    out_norm = get_normalized(results,
                              method = rownames(get_params(results))[1])
    assay(sce,"Scone") = log2(out_norm + 1)
  })
  
  method_name = "scone"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}
norm_method <- list(
  none = raw_count,
  scran = scran_norm,
  pareto=pareto_norm,
  DESeq2=DESeq2_norm,
  TMM=TMM_norm,
  logCPM=logCPM_norm,
  Linnorm=linnorm_norm#, for the moment scone and SCnorm will be excluded cause they are very time consuming
  #scone=scone_norm,
  #SCnorm=SCnorm_norm
)

res1 <- datasets %>%
  apply_methods(norm_method)

r_time<-c(metadata(res1$result[[1]])$running_time[3],
          metadata(res1$result[[2]])$running_time[3],
          metadata(res1$result[[3]])$running_time[3],
          metadata(res1$result[[4]])$running_time[3],
          metadata(res1$result[[5]])$running_time[3],
          metadata(res1$result[[6]])$running_time[3],
          metadata(res1$result[[7]])$running_time[3],
          metadata(res1$result[[8]])$running_time[3],
          metadata(res1$result[[9]])$running_time[3],
          metadata(res1$result[[10]])$running_time[3],
          metadata(res1$result[[11]])$running_time[3],
          metadata(res1$result[[12]])$running_time[3],
          metadata(res1$result[[13]])$running_time[3],
          metadata(res1$result[[14]])$running_time[3])
r_time<-as.numeric(r_time)
data<-cbind(res1[,1:2], r_time)
data #the running time for each normalization and dataset

#PCAplot (excluding none normalization)
#-----------celseq2
sce<-res1$result[[2]]
sce<-runPCA(sce, exprs_values="Scran")
Scran_cs2<-plotPCA(sce, colour_by="cell_line")+ggtitle("Scran", subtitle = "CELSeq2")

sce<-res1$result[[3]]
sce<-runPCA(sce, exprs_values="Pareto")
Pareto_cs2<-plotPCA(sce, colour_by="cell_line")+ggtitle("Pareto")

sce<-res1$result[[4]]
sce<-runPCA(sce, exprs_values="DESeq2")
DESeq2_cs2<-plotPCA(sce, colour_by="cell_line")+ggtitle("DESeq2")

sce<-res1$result[[5]]
sce<-runPCA(sce, exprs_values="TMM")
TMM_cs2<-plotPCA(sce, colour_by="cell_line")+ggtitle("TMM")

sce<-res1$result[[6]]
sce<-runPCA(sce, exprs_values="logCPM")
CPM_cs2<-plotPCA(sce, colour_by="cell_line")+ggtitle("logCPM")

sce<-res1$result[[7]]
sce<-runPCA(sce, exprs_values="Linnorm")
Linnorm_cs2<-plotPCA(sce, colour_by="cell_line")+ggtitle("Linnorm")

multiplot(Scran_cs2, Pareto_cs2, DESeq2_cs2,
          TMM_cs2, CPM_cs2, Linnorm_cs2, cols = 3)

#---------10X
sce<-res1$result[[9]]
sce<-runPCA(sce, exprs_values="Scran")
Scran_10X<-plotPCA(sce, colour_by="cell_line")+ggtitle("Scran", subtitle = "10X")

sce<-res1$result[[10]]
sce<-runPCA(sce, exprs_values="Pareto")
Pareto_10X<-plotPCA(sce, colour_by="cell_line")+ggtitle("Pareto")

sce<-res1$result[[11]]
sce<-runPCA(sce, exprs_values="DESeq2")
DESeq2_10X<-plotPCA(sce, colour_by="cell_line")+ggtitle("DESeq2")

sce<-res1$result[[12]]
sce<-runPCA(sce, exprs_values="TMM")
TMM_10X<-plotPCA(sce, colour_by="cell_line")+ggtitle("TMM")

sce<-res1$result[[13]]
sce<-runPCA(sce, exprs_values="logCPM")
CPM_10X<-plotPCA(sce, colour_by="cell_line")+ggtitle("logCPM")

sce<-res1$result[[14]]
sce<-runPCA(sce, exprs_values="Linnorm")
Linnorm_10X<-plotPCA(sce, colour_by="cell_line")+ggtitle("Linnorm")

multiplot(Scran_10X, Pareto_10X, DESeq2_10X,
          TMM_10X, CPM_10X, Linnorm_10X, cols = 3)


