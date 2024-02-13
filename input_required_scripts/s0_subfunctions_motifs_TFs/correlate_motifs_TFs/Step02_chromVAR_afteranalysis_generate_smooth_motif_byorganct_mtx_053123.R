###################################################################################################
###                             chromVAR analysis of scATAC data                                ###
###################################################################################################

##updating 122222 get the average of motifs information
##updating 012921 generate a sparse for the output of deviation score
##add note information

# load libraries
library(chromVAR)
library(motifmatchr)
library(BiocParallel)
#library(BSgenome.Gmax.a4.v1)
#library(BSgenome.Zmays.AGP.v4)
library(Matrix)
library(SummarizedExperiment)
library(GenomicAlignments)
library(dplyr)
library(TFBSTools)
#library(JASPAR2018)
#library(JASPAR2020)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(Matrix)
library(gtools)
library(edgeR)
library(preprocessCore)
library(matrixStats)

# arguments
args <- commandArgs(TRUE)
input_sparse_fl <- as.character(args[1])
meta   <- as.character(args[2])
targetClust <- as.character(args[3])

input_output_dir <- as.character(args[4])


a <- read.table(input_sparse_fl,stringsAsFactors = T)

b <- read.delim(meta,header=T,row.names = 1)
rownames(b) <- b$cellID


# format
a <- sparseMatrix(i=as.numeric(a$V1),j=as.numeric(a$V2),x=as.numeric(a$V3),dimnames=list(levels(a$V1),levels(a$V2)))

##allow them to be same dim
share_id <- intersect(rownames(b),colnames(a))
b <- b[share_id,]

# iterate over clusters
clusts <- mixedsort(unique(b[,targetClust]))

mat <- matrix(0,nrow=nrow(a),ncol=length(clusts))

it <- 0
for (i in clusts){
  it <- it+1
  message(" - iterate over all genes for cluster ", i)
  
  #b <- read.table('opt_tfidf_NMF_30PCs_win30000_res2_metadata.txt')
  #head(b)
  #targetClust <- 'LouvainClusters'
  #ids <- rownames(subset(b,b[,targetClust]==i))
  #b$LouvainClusters
  
  ids <- rownames(subset(b, b[,targetClust]==i))
  #ids <- rownames(subset(b, b$Final_annotation_clust==i))
  aa <- a[,ids]
  #mat[,it] <- Matrix::rowSums(aa)
  mat[,it] <- Matrix::rowMeans(aa)
  

}

colnames(mat) <- clusts
rownames(mat) <- rownames(a)

write.table(mat, file=paste0(input_output_dir,'/opt_motif_average_smooth_accessible_clusters.txt'), quote=F, row.names=T, col.names=T, sep="\t")




