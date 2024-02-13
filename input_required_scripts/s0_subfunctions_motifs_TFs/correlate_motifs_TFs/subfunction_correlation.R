##this script is to do the correlation 
##updating 060223 we will also return a pvalue mtx dataframe

library(Hmisc)

# arguments
args <- commandArgs(TRUE)

##col name is organct rowname is gene ID
input_tf_matrix_fl <- as.character(args[1])
##col name is organct rowname is motif
input_motif_matrix_fl <- as.character(args[2])

input_output_dir <- as.character(args[3])


gene_dt <- read.delim(input_tf_matrix_fl)
motif_dt <- read.delim(input_motif_matrix_fl)

##after reverse colname is IDs and row name is organct
#corr_dt <- cor(t(gene_dt),t(motif_dt),method = 'spearman')

##we will do the correlation and returns the corr
outs <- rcorr(as.matrix(t(gene_dt)),as.matrix(t(motif_dt)),type="spearman")
gene_dt_ncol <- ncol(t(gene_dt))
motif_dt_ncol <- ncol(t(motif_dt))
corr_p <- outs$P
corr_p_right <- corr_p[1:gene_dt_ncol,(gene_dt_ncol+1):(gene_dt_ncol+motif_dt_ncol)]
corr_cor <- outs$r
corr_cor_right <- corr_cor[1:gene_dt_ncol,(gene_dt_ncol+1):(gene_dt_ncol+motif_dt_ncol)]


##write out the corr_dt
write.table(corr_cor_right,paste0(input_output_dir,'/opt_corr_dt.txt'),quote = F, sep = '\t')
write.table(t(corr_cor_right),paste0(input_output_dir,'/opt_corr_dt_rev.txt'),quote = F, sep = '\t')


##write out the pval_dt
write.table(corr_p_right,paste0(input_output_dir,'/opt_pval_dt.txt'),quote = F, sep = '\t')
write.table(t(corr_p_right),paste0(input_output_dir,'/opt_pval_dt_rev.txt'),quote = F, sep = '\t')



