##updating 020123 we will also consider the tn5 for each cell to remove this effect
##updating 121622 we will check how many motifs we already anlayzed
##updating 080922 we need to add an option to check which we have provide the whole motif x peak matrix
##we need to provide input from Step01_s6_opt2
##updating 080722 create dir
##updating 080422 this script is to do a binomial regression test
##prefdic total motif counts using two input variables:
##1) target_cluster_annotation
##2) log(total number of nonzero entries in input peak matrix) ?? log(peak number) per cell?


library(Matrix)
library(parallel)
library(doSNOW)
library(VGAM)
library(methods)
library(tcltk)
library(iterators)
library(itertools)
library(dplyr)
library(MASS)

##examples
##outs <- diffaccess_sub1(d, factor(df.sub$Cluster2), df.sub$log10nSites)
##diffaccess_sub1 <- function(x, y, z){
#df.s <- data.frame(open=as.numeric(x), cluster=as.factor(y), readdp=as.numeric(z))
#fit1 <- suppressWarnings(VGAM::vglm(formula=open~cluster+readdp, family=binomialff(),data=df.s))
#fit2 <- suppressWarnings(VGAM::vglm(formula=open~readdp, family=binomialff(), data=df.s))
#lrt <- VGAM::lrtest(fit1,fit2)
#pval=lrt@Body["Pr(>Chisq)"][2,]
#beta <- coef(fit1)["cluster1"]
#out <- c(pval, beta)
#unname(out, force=T)
#return(out)
#}



args <- commandArgs(TRUE)


ipt_all_rds_dir <- as.character(args[1])
#ipt_peak_by_cell_mtx_rds_fl <- as.character(args[1])
#ipt_peak_by_motif_mtx_rds_fl <- as.character(args[2])
ipt_meta_fl <- as.character(args[2]) ##the meta file contains the number ACR for each cell
##we do not want to filer cells just try to use all of them
##we just use the opt_v4.5_addACRnum.txt
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/11_TFmotif_analysis_040522/output_dir_072522/step01_celltype_enrich_motif_dir/s5_prepare_regr_enrich_dir/opt_downsample_meta.txt
input_output_dir <- as.character(args[3])
target_cluster_colnm <- as.character(args[4])
threads <- as.character(args[5])


##updating 080722 create a dir
##updating 121622
finish_motif_names <- 'na'
if (!dir.exists(paste0(input_output_dir,'/store_motif_results_dir'))){
  dir.create((paste0(input_output_dir,'/store_motif_results_dir')))
  
  ##it shows the dir is not existing it is the first time to do the analysis, we will not check already analysed motifs
  finish_motif_names <- 'na'
  
} else {
  print("Dir already exists!")
  
  ##we will check how many files are 
  motif_files <- list.files(path = paste0(input_output_dir,'/store_motif_results_dir'))
  
  ##check whether motif_files is empty
  if (length(motif_files) != 0){
    message ('some motifs already analyzed, we will store them in a list')
    
    motif_files_a <- gsub('opt_','',motif_files)
    motif_files_b <- gsub('\\.txt','',motif_files_a)
    finish_motif_names <- motif_files_b
    
  }else{
    
    message('no motif files found, we still need to analyze motif from begining')
    finish_motif_names <- 'na'
  }

}



##step01 load files
ipt_meta_dt <- read.delim(ipt_meta_fl,row.names = 1)
peak_by_cell_rds_fl_list <- list.files(path = ipt_all_rds_dir,pattern = '*peak_cell*')
peak_by_motif_rds_fl <- paste0(ipt_all_rds_dir,'/','opt_all_motif_matches_mtx.rds')
#peak_by_motif_rds_fl_list <- list.files(path = ipt_all_rds_dir,pattern = '*motif_matches_mtx.rds')

if (file.exists(paste0(input_output_dir,'/opt_final_motif_by_cell_mtx.rds'))){
  
  ##updating 121622
  message (' - opt_motif_by_cell_mtx.rds is existed, and we will directly read it')
  motif_by_cell_mtx <- readRDS(paste0(input_output_dir,'/opt_final_motif_by_cell_mtx.rds'))
  
  ##filter motifs based on the finish_motif_names
  intersect_motifs <- intersect(rownames(motif_by_cell_mtx),finish_motif_names)
  outs <- lapply(rownames(motif_by_cell_mtx),function(x){
    if ((x %in% intersect_motifs) == FALSE){
      return(x)
    }
  })
  combine <- as.data.frame(do.call(rbind,outs))
  
  motif_by_cell_mtx <- motif_by_cell_mtx[combine$V1,]

  
}else{
  message(' - load rds fl')
  ##we need to read each rds and then combine them together
  outs <- lapply(peak_by_cell_rds_fl_list,function(x){
    
    ipt_peak_by_cell_mtx = readRDS(paste0(ipt_all_rds_dir,'/',x))
    organnm_1 <- gsub('opt_','',x)
    organnm_2 <- gsub('_peak_cell_mtx.rds','',organnm_1)
    ipt_peak_by_motif_mtx = readRDS(peak_by_motif_rds_fl)
    
    message(' - multiply two matrix')
    ipt_motif_by_peak_mtx <- t(ipt_peak_by_motif_mtx)
    ##make two matrix has same dimension
    shared_peaks <- intersect(colnames(ipt_motif_by_peak_mtx),rownames(ipt_peak_by_cell_mtx))
    ipt_motif_by_peak_mtx <- ipt_motif_by_peak_mtx[,shared_peaks]
    ##intersect cells
    shared_cells <- intersect(rownames(ipt_meta_dt),colnames(ipt_peak_by_cell_mtx))
    ipt_peak_by_cell_mtx <- ipt_peak_by_cell_mtx[shared_peaks,shared_cells]
    
    motif_by_cell_mtx <- ipt_motif_by_peak_mtx %*% ipt_peak_by_cell_mtx
    saveRDS(motif_by_cell_mtx,paste0(input_output_dir,'/opt_',organnm_2,'_motif_by_cell_mtx.rds'))
    
    return(motif_by_cell_mtx)
  })
  
  outs_final_combine <- do.call(cbind, outs)
  saveRDS(outs_final_combine,paste0(input_output_dir,'/opt_final_motif_by_cell_mtx.rds'))
  motif_by_cell_mtx <-  readRDS(paste0(input_output_dir,'/opt_final_motif_by_cell_mtx.rds'))
  
  
  ##updating 121622
  message (' - opt_motif_by_cell_mtx.rds is existed, and we will directly read it')

  ##filter motifs based on the finish_motif_names
  intersect_motifs <- intersect(rownames(motif_by_cell_mtx),finish_motif_names)
  outs <- lapply(rownames(motif_by_cell_mtx),function(x){
    if ((x %in% intersect_motifs) == FALSE){
      return(x)
    }
  })
  combine <- as.data.frame(do.call(rbind,outs))
  
  motif_by_cell_mtx <- motif_by_cell_mtx[combine$V1,]
  
  
}



message ('- input motif_by_cell_mtx')
message(paste0('-row number is ', nrow(motif_by_cell_mtx)))
message(paste0('-col number is ', ncol(motif_by_cell_mtx)))


#ipt_peak_by_cell_mtx <- readRDS(ipt_peak_by_cell_mtx_rds_fl)
#ipt_peak_by_motif_mtx <- readRDS(ipt_peak_by_motif_mtx_rds_fl)

#message(' - multiply two matrix')
#ipt_motif_by_peak_mtx <- t(ipt_peak_by_motif_mtx)
##make two matrix has same dimension
#shared_peaks <- intersect(colnames(ipt_motif_by_peak_mtx),rownames(ipt_peak_by_cell_mtx))
#ipt_motif_by_peak_mtx <- ipt_motif_by_peak_mtx[,shared_peaks]
##intersect cells
#shared_cells <- intersect(rownames(ipt_meta_dt),colnames(ipt_peak_by_cell_mtx))
#ipt_peak_by_cell_mtx <- ipt_peak_by_cell_mtx[shared_peaks,shared_cells]
#motif_by_cell_mtx <- ipt_motif_by_peak_mtx %*% ipt_peak_by_cell_mtx
#saveRDS(motif_by_cell_mtx,paste0(input_output_dir,'/opt_motif_by_cell_mtx.rds'))


##step02 merge meta to prepare the regression analysis
#message(' - merge to a new meta')
##updating 020123 we will 
ipt_meta_dt$log10tn5 <- log10(ipt_meta_dt$total)

target_meta_dt <- ipt_meta_dt[c('ACRnum','log10tn5','library',target_cluster_colnm)]
cell_by_motif_mtx <- as.data.frame(as.matrix(t(motif_by_cell_mtx)))
merged_meta_dt <- merge(target_meta_dt,cell_by_motif_mtx,by = 'row.names')


##step03 conduct regression analysis
##start from 7 since we add the log10tn5 and library the previous is 5 as we have no log10tn5 item
outs <- mclapply(colnames(merged_meta_dt)[c(7:ncol(merged_meta_dt)-1)], function(z){
  
  outs2 <- lapply(unique(merged_meta_dt[[target_cluster_colnm]]), function(x){
    
    # verbose
    message(" - check motif ", z , ' and check cluster ', x)
    
    merged_meta_dt$Cluster2 <- ifelse(merged_meta_dt[,target_cluster_colnm]==x,1,0)
    
    #df.sub$Cluster2 <- ifelse(df.sub[,var]==snnclusts[j],1,0)
    
    ##we need to remove the ACR number to influence the model like if one cell has more ACR number it will have more motif count
    df.s <- data.frame(motifc=as.numeric(merged_meta_dt[[z]]), cluster=as.factor(merged_meta_dt$Cluster2), ACRnum=log(as.numeric(merged_meta_dt$ACRnum)), tn5num =as.numeric(merged_meta_dt$log10tn5),
                       library = as.factor(merged_meta_dt$library))
    
    message(' - fit model')
    M1 <- glm.nb(motifc ~ cluster + ACRnum + tn5num + library,
              data = df.s)
  
    coefficient <- M1$coefficients[2]
    intercept <- M1$coefficients[1]
    message(" - check motif ", z , ' and check cluster ', x, " | intercept = ", intercept, " | coefficient = ", coefficient)
    
    fc <- as.numeric(exp(intercept + coefficient)/exp(intercept))
    pval <- as.numeric(summary(M1)$coefficients[,4][2])
    
    res <- c()
    res$pval <- pval
    res$fc <- fc
    res$beta <- coefficient
    res$motif <- z
    
    return(as.data.frame(res))
    
  })
  
  outs2_combine <- data.frame(do.call(rbind, outs2))
  outs2_combine$organct <- unique(merged_meta_dt[[target_cluster_colnm]])
  #rownames(outs2_combine) <- unique(merged_meta_dt[[target_cluster_colnm]])[1:2]
  
  outs2_combine <- as.matrix(outs2_combine)
  write.table(outs2_combine,paste0(input_output_dir,'/store_motif_results_dir/','opt_',z,'.txt'),sep = '\t',quote = F)
  
  return(outs2_combine)
  
}, mc.cores=threads)

saveRDS(outs,paste0(input_output_dir,'/opt_outs.rds'))

##it is important to have a as.data.frame
outs_final_combine <- as.data.frame(do.call(rbind, outs))

##we need to generate the qvalue
outs_final_combine$FDR <- p.adjust(outs_final_combine$pval, method = 'fdr')

outs_final_combine <- as.matrix(outs_final_combine)

write.table(outs_final_combine,paste0(input_output_dir,'/opt_motif_enrichment_cluster.txt'),sep = '\t',quote = F)





#res$table$FDR <- p.adjust(res$table[,4], method="fdr")
#qvals <- apply(pvals, 2, function(x){p.adjust(x, method="fdr")})
























