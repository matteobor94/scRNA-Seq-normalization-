library(scran)
library(scater)
library(ggplot2)
library(DT)
library(Seurat)
library(DESeq)

load("C:/Users/Matteo/Desktop/Borsa studio/sincell_with_class.rdata") #my personal directory
#datasets at https://github.com/LuyiTian/sc_mixology

#rename datasets
sce_10x<-sce_sc_10x_qc
sce_cseq2<-sce_sc_CELseq2_qc
sce_dropseq<-sce_sc_Dropseq_qc
remove(sce_sc_10x_qc, sce_sc_CELseq2_qc, sce_sc_Dropseq_qc)

#Scran-norm 
logcounts(sce_10x)[1:10, 1:5] #scran matrix already exists

#Pareto-norm
pareto.MLE <- function(X)
{
  n <- length(X)
  m <- min(X)
  a <- n/sum(log(X)-log(m))
  return(c(m,a)) 
}
assay(sce_10x,"pareto")<-matrix(nrow = nrow(assay(sce_10x)), 
                                ncol = ncol(assay(sce_10x)))
alfa<-apply(counts(sce_10x)+1, 2, pareto.MLE)[2,] 
for (i in 1:length(alfa)) {
  assay(sce_10x,"pareto")[,i]<-log2(assay(sce_10x,"counts")[,i]*alfa[i]+1)
}
assay(sce_10x, "pareto")[1:10, 1:5] #values normalized with the Pareto parameter

#PCAplot
sce_10x<-runPCA(sce_10x, exprs_values="counts")
raw_counts<-plotPCA(sce_10x, colour_by="cell_line", #colour by different cell-type
                    shape_by="demuxlet_cls")+ggtitle("Raw counts") #different shapes whether the cell is a doublet or not

sce_10x<-runPCA(sce_10x, exprs_values="logcounts")
scran<-plotPCA(sce_10x, colour_by="cell_line", 
               shape_by="demuxlet_cls")+ggtitle("Scran norm") 

sce_10x<-runPCA(sce_10x, exprs_values="pareto")
pareto<-plotPCA(sce_10x, colour_by="cell_line",
                shape_by="demuxlet_cls")+ggtitle("Pareto norm")
multiplot(raw_counts, scran, pareto, cols=3)


