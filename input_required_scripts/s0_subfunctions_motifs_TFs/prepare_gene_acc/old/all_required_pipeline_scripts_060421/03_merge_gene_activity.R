###################################
##merge gene body and gene activity
##updating 071221 add the middle temp file

##load libraries
library(Matrix)

args <- commandArgs(T)
geneact <- as.character(args[1]) ##geneActivity
genebod <- as.character(args[2]) ##genebody
output  <- as.character(args[3]) ##set prefix
output_dir <- as.character(args[4]) ##output_dir

##load data
act <- read.table(geneact)
bod <- read.table(genebod)

##keep the same cells
message(' - check shared cells')
shared.cells <- intersect(unique(as.character(act$V2)), unique(as.character(bod$V2)))
bod <- bod[as.character(bod$V2) %in% shared.cells,]
act <- act[as.character(act$V2) %in% shared.cells,]

##merge data
colnames(act) <- c("gene","cell","activity")
colnames(bod) <- c("gene","cell","geneaccess")

##set the weight to merge
message (' - combine gene ciciro acc and gene body acc')
cboth <- merge(act, bod, by=c("gene","cell"), all=T)
cboth[is.na(cboth)] <- 0
message (' - transform to the dataframe')
cboth <- as.data.frame(cboth)
message (' - combine act and geneacc')
cboth$combined <- cboth$activity + (3*cboth$geneaccess)
cboth$activity <- NULL
cboth$geneaccess <- NULL
message (' - finish combination')
cboth$gene <- as.factor(cboth$gene)
cboth$cell <- as.factor(cboth$cell)

##write the temp file
message (' - save temp cboth file')
saveRDS(cboth,paste0(output_dir,'/cboth.rds'))

##normalize the columns to be 1
message (' normalize the new gene acc to be 1')
combine_mat <- sparseMatrix(i=as.numeric(cboth$gene),
			 j=as.numeric(cboth$cell),
			 x=as.numeric(cboth$combined),
			 dimnames=list(levels(cboth$gene),levels(cboth$cell)))
mat_id <- colnames(combine_mat)
combine_mat <- combine_mat %*% Diagonal(x=1/Matrix::colSums(combine_mat))
colnames(combine_mat) <- mat_id

##save the temp file
message (' - save temp combine_mat file')
saveRDS(combine_mat,paste0(output_dir,'/combine_mat.rds'))

##transfer to the sparse
message (' transfer to sparse')
cboth <- as.data.frame(summary(combine_mat))
cboth$i <- rownames(combine_mat)[cboth$i]
cboth$j <- colnames(combine_mat)[cboth$j]
write.table(cboth, file=paste0(output_dir,'/',output,".GAadjusted.sparse"),
            row.names=F, col.names=F, quote=F, sep="\t")



