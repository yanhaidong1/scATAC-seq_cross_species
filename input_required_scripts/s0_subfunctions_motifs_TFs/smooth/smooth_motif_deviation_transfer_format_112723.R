###################################################################################################
###                             chromVAR analysis of scATAC data                                ###
###################################################################################################

##updating 012921 generate a sparse for the output of deviation score
##add note information

# load libraries
#library(chromVAR)
#library(motifmatchr)
#library(BiocParallel)
#library(BSgenome.Gmax.a4.v1)
#library(BSgenome.Zmays.AGP.v4)
library(Matrix)
#library(SummarizedExperiment)
#library(GenomicAlignments)
#library(dplyr)
#library(TFBSTools)
#library(JASPAR2018)
#library(JASPAR2020)
#library(pheatmap)
#library(ComplexHeatmap)
#library(circlize)
library(Seurat)

# arguments
args <- commandArgs(TRUE)
input_ds_fl <- as.character(args[1])
input_output_dir <- as.character(args[2])

score_dt <- read.table(input_ds_fl)

tpm.mat <- as.sparse(as.matrix(t(score_dt)))

##transfer to sparse
ia <- as.data.frame(summary(tpm.mat))
ia$i <- rownames(tpm.mat)[as.numeric(ia$i)]
ia$j <- colnames(tpm.mat)[as.numeric(ia$j)]
write.table(ia, file=paste0(input_output_dir,"/opt_motif_smooth_deviations.sparse"), quote=F, row.names=F, col.names=F, sep="\t")





