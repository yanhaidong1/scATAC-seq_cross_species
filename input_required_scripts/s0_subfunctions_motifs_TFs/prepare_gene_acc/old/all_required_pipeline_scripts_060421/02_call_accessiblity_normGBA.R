##updating 070921 set message to check the error
##mute the gene_attr
##updating 061621 use the raw tn5 to be input and open the CTnorm
##this script is to generate accessibility from gene body based on tn5 mapping to the body regions

#the main idea of the g is to check the length of genes and conduct a filtration

# load arguments
args <- commandArgs(T)
if(length(args)!=5){stop("Rscript calculate_normGBA.R <gene.sparse> <meta> <TAIR10_genes_500bpTSS.bed> <prefix> <output_dir>")}
input <- as.character(args[1])
meta <- as.character(args[2])
gene <- as.character(args[3])
prefix <- as.character(args[4])
output_dir <- as.character(args[5])

##create sub output to store the files generated from this script
genebody_accessibility_dir <- paste0(output_dir,'/genebody_accessibility_dir')
dir.create(genebody_accessibility_dir)

##load libraries
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(sctransform)
library(ggplot2)

##load data
message(" - loading data ...")
if (file.exists(paste0(output_dir,'/temp_mtx.rds'))){
  
  a <- readRDS(paste0(output_dir,'/temp_mtx.rds'))
  
}else{
  
  a <- read.table(input,stringsAsFactors = T)
  b <- read.table(meta)
  g <- read.table(gene)
  
  ##process
  ##filter out genes without meeting the requirements of length
  a <- sparseMatrix(i=as.numeric(a$V1),
                    j=as.numeric(a$V2),
                    x=as.numeric(a$V3),
                    dimnames=list(levels(a$V1),levels(a$V2)))
  a <- a[,colnames(a) %in% rownames(b)]
  b <- b[colnames(a),]
  rownames(g) <- g$V4
  g <- g[rownames(a),]
  g$gene.len <- g$V3 - g$V2
  g <- subset(g, g$gene.len > 100)
  a <- a[rownames(g),]
  
  ##save the sparseMatrix
  saveRDS(a,paste0(output_dir,'/temp_mtx.rds'))

}
  
# gene attribute
##rownames of a is the gene information
##generate a gene_attr dataframe
#message(" - generate gene_attr data ...")
#gene_attr <- data.frame(mean = Matrix::rowMeans(a), 
#                        detection_rate = Matrix::rowMeans(a > 0),
#                        var = apply(a, 1, var))
#saveRDS(gene_attr,paste0(output_dir,'/temp_gene_attr.rds'))
#message(' - add the log_mean and log_var')
#gene_attr$log_mean <- log10(gene_attr$mean)
#gene_attr$log_var <- log10(gene_attr$var)
#rownames(gene_attr) <- rownames(a)
##generate a cell_attr dataframe
message(' - generate cell_attr data ...')
cell_attr <- data.frame(n_umi = Matrix::colSums(a),
                        n_gene = Matrix::colSums(a > 0),
                        log_umi = log10(Matrix::colSums(a)))
rownames(cell_attr) <- colnames(a)

##Poisson regression is used to model count variables.
##this figure is to check whether the gene_attr fit the poisson_model
#x = seq(from = -4, to = 1.5, length.out = 1000)
#poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
#ggplot(gene_attr, aes(log_mean, detection_rate)) + 
#    geom_point(alpha=0.3, shape=16) + 
#    geom_line(data=poisson_model, color='red') +
#    theme_gray(base_size = 8)

# cell - gene accessibility relationships
#ggplot(cell_attr, aes(n_umi, n_gene)) + 
#    geom_point(alpha=0.3, shape=16) + 
#    geom_density_2d(size = 0.3)
#ggsave(file=paste0(genebody_accessibility_dir,'/',prefix,".normalizationMetrics.pdf"), width=6, height=6, device="pdf")

#a_flt <- a[,Matrix::colSums(a>0)>300]
#a_flt <- a_flt[Matrix::rowSums(a_flt>0)>100,]
#norm <- a
##fail to do the sctransform 
# normalize accessibility by gene length
message(" - normalizing gene accessibility ...")
norm <- sctransform::vst(a, cell_attr = cell_attr,
                         latent_var = c('log_umi'), 
                         return_gene_attr = TRUE, 
                         return_cell_attr = TRUE, 
                         show_progress = TRUE)

message(' - save the norm information')
saveRDS(norm, paste0(output_dir,'/temp_norm.rds'))

message(' - correct and set matrix for the norm')
norm <- correct(norm, do_round = T)
norm <- Matrix(norm, sparse=T)

# filter rows/columns - 100 genes accessible per cell & at least 10 cells accessible per gene
message(' - filter norm')
norm <- norm[,Matrix::colSums(norm>0)>0]
norm <- norm[Matrix::rowSums(norm>0)>0,]

# normalize to sum gene accessibility to 1
message(' - normalize the gene acc to be 1')
ids <- colnames(norm)
norm <- norm %*% Diagonal(x=1/Matrix::colSums(norm))
colnames(norm) <- ids

# write output
message(' - write out sparse file')
out <- as.data.frame(summary(norm))
out$i <- rownames(norm)[out$i]
out$j <- colnames(norm)[out$j]
write.table(out, file=paste0(genebody_accessibility_dir,'/',prefix,".GBaccessibility.sparse"),
            quote=F, row.names=F, col.names=F, sep="\t")



