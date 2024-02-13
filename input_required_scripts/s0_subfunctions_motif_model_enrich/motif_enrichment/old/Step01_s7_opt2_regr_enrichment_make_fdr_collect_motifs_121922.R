##updating 121922 this script is to do the fdr per file and combine all the results
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


ipt_motif_results_dir <- as.character(args[1])

input_output_dir <- as.character(args[2])


all_fl_list <- list.files(path = ipt_motif_results_dir)


outs <- lapply(all_fl_list, function(x){
  
  target_fl <- paste0(ipt_motif_results_dir,'/',x)
  
  target_dt <- read.delim(target_fl)
  
  target_dt$FDR <- p.adjust(target_dt$pval, method = 'fdr')
  
  return(target_dt)
  
})

combine_dt <- do.call(rbind, outs)

write.table(combine_dt,paste0(input_output_dir,'/opt_motif_enrichment_cluster.txt'),sep = '\t',quote = F)







