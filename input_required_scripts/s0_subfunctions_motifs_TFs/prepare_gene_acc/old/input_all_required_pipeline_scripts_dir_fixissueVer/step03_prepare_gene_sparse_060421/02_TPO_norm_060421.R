##this script we will normalize TPO per one instead of million for the gene sparse
##this script we will conduct TPM normalization for the gene sparse
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(sctransform)
library(ggplot2)
library(Seurat)

#the main idea of the g is to check the length of genes and conduct a filtration

# load arguments
args <- commandArgs(T)
if(length(args)!=3){stop("Rscript norm_TPO_gene_sparse <gene.sparse> <gene_length_fl> <output_dir>")}
gene_sparse_fl <- as.character(args[1])
gene_length_fl <- as.character(args[2])
output_dir <- as.character(args[3])


a <- read.table(gene_sparse_fl,stringsAsFactors = T)
a <- sparseMatrix(i=as.numeric(a$V1),
                  j=as.numeric(a$V2),
                  x=as.numeric(a$V3),
                  dimnames=list(levels(a$V1),levels(a$V2)))
#a <- head(a)[1:6,1:50]

gene_len_dt <- read.delim(gene_length_fl,header = FALSE)
rownames(gene_len_dt) <- gene_len_dt$V1
order_gene_len_dt <- gene_len_dt[rownames(a),]
len_vector <- as.numeric(order_gene_len_dt$V2)

a_div <- a/len_vector

tpm.mat <- t(t(a_div) / colSums(a_div)) 
#tpm.mat <- t(t(a_div) * 1e6 / colSums(a_div)) 
#tpm.mat[1:6,1:50]
tpm.mat[is.na(tpm.mat)] <- 0

tpm.mat <- as.sparse(tpm.mat)

##transfer to sparse
ia <- as.data.frame(summary(tpm.mat))
ia$i <- rownames(tpm.mat)[as.numeric(ia$i)]
ia$j <- colnames(tpm.mat)[as.numeric(ia$j)]
write.table(ia, file=paste0(output_dir,"/opt_gene_sparse_TPO_norm.sparse"), quote=F, row.names=F, col.names=F, sep="\t")

