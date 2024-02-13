##updating 123123 we will check species specific cases not in the syntenic regions
##updating 122823 we will check the exp
##updating 120323 we will add the motif heatmap plotting
##updating 112823 we will check the correlation between cell types based on the motif deviation scores
##updating 111023 we will build the dist plot for the different categories
##updating 101623 we will check the rice maize and sorghum ACRs
##updating 011123 generate a plot for the number of celltypes per ACRs
##updating 113022 generate pie plot for showing the composition of TEs for each category 
##updating 070122 we need to generate some cases for the rice.
##updating 081721 generate a version wrapped to generate figures automatically
##updating 032421 generate a new version 
##this script will generate pie chart for the composition of TE as well as acr

library(ggplot2)
library(reshape2)
library(dplyr)

library(ggpubr)
library(cowplot)
library(scater)
library(plyr)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

library(RColorBrewer)
myPalette <- brewer.pal(10, "Set3") 
library(tidyr)



########
##step00
######################################
##Here we will plot the TE composition
ipt_dir <- 'opt0_te_proportion_111623/'

opt_dir <- 'opt0_te_proportion_results_111623'

all_fl_list <- list.files(path= ipt_dir)



for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//opt_te_prop_loctype','',target_fl_path)
  target_fl_path_nm <- gsub('\\.txt','',target_fl_path_nm)
  
  ipt_target_dt <- read.delim(target_fl_path,header = F)
  colnames(ipt_target_dt) <- c('Cate','TEfam','Num','Prop','TEsupfam')
  
  ipt_target_dt$TEsupfam_new <- gsub('Others_.+','Unknown',ipt_target_dt$TEsupfam)
  ipt_target_dt$TEsupfam_new <- gsub('.+others_.+','Unknown',ipt_target_dt$TEsupfam_new)
  

  ipt_target_dt <- ipt_target_dt[ipt_target_dt$TEsupfam_new != 'LINE' &ipt_target_dt$TEsupfam_new != 'SINE' &
                                   ipt_target_dt$TEsupfam_new != 'SateliteDNA'&
                                   ipt_target_dt$TEsupfam_new != 'Unknown',]
  
  aggregate_celltype_superfam_dt <- aggregate(Num ~ Cate + TEsupfam_new, data = ipt_target_dt, sum)
  
  aggregate_celltype_dt <- aggregate(Num ~ Cate, data = ipt_target_dt,sum)
  
  aggregate_merge <- merge(aggregate_celltype_superfam_dt,aggregate_celltype_dt,by.x = 'Cate',by.y = 'Cate')
  aggregate_merge$prop <- aggregate_merge$Num.x/aggregate_merge$Num.y
  
  p <- ggplot(data=aggregate_merge, aes(x=Cate, y=prop, fill=TEsupfam_new)) +
    geom_bar(stat="identity") + 
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'_te.pdf'),width = 6,height = 6)
  print(p)
  dev.off()
  

  
}


######################################
##Here we will plot the TE composition
ipt_dir <- 'opt0_acr_celltype_dir_111623//'
opt_dir <- 'opt0_acr_celltype_dir_results_111623'

all_fl_list <- list.files(path= ipt_dir)

generate_pie <- function(df) {
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df$final_celltype,pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}


for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_acr_celltype.txt','',target_fl_path_nm)
  
  ipt_target_dt <- read.delim(target_fl_path,header = F)
  if (ncol(ipt_target_dt) == 5){
    colnames(ipt_target_dt) <- c('chr','st','ed','celltype','pval')
  }
  if (ncol(ipt_target_dt) == 6){
    colnames(ipt_target_dt) <- c('chr','st','ed','celltype','pval','fc')
  }
  
  ipt_target_dt$celltype_modi <- gsub('.+;','',ipt_target_dt$celltype)
  
  ipt_target_dt$final_celltype <- ifelse(ipt_target_dt$celltype_modi == 'broadly_accessible','Broad','CT')
  
  ipt_target_dt$count <- 1
  aggregate_dt <- aggregate(count ~ final_celltype,data = ipt_target_dt,sum)
  
  p <- ggplot(data=aggregate_dt, aes(x=final_celltype, y=count)) +
    geom_bar(stat="identity") + 
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    coord_cartesian(ylim = c(0, 55000)) +
    scale_y_continuous(breaks=seq(0, 55000, 10000)) 
  
  p
  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'_broad_ct_ACRnum.pdf'),width = 4,height = 6)
  print(p)
  dev.off()

  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'_broad_ct_ACRnum_piechart.pdf'),width = 4,height = 6)
  print(generate_pie(aggregate_dt))
  dev.off()
  
  
}


############################################
##Here we will plot the piechart composition
ipt_dir <- 'opt0_acr_geneCate_dir_111623///'
opt_dir <- 'opt0_acr_geneCate_dir_results_111623/'

generate_pie <- function(dt,target_celltype) {
  df <- dt[dt$celltype == target_celltype,]
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df$cate,pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}



all_fl_list <- list.files(path= ipt_dir)

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_acr_toGeneCate.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_target_dt <- read.delim(target_fl_path,header = F)
  head(ipt_target_dt)
  colnames(ipt_target_dt) <- c('acr','cate','dist','celltype')

  ipt_target_dt$count <- 1
  
  aggregate_dt <- aggregate(count ~ celltype + cate, data = ipt_target_dt, sum)
  
  pdf(paste0(opt_dir,'/',target_fl_path_nm,'_broad','.pdf'),width = 7,height = 7)
  generate_pie(aggregate_dt,'Broad')
  dev.off()
  
  pdf(paste0(opt_dir,'/',target_fl_path_nm,'_CT','.pdf'),width = 7,height = 7)
  generate_pie(aggregate_dt,'CT')
  dev.off()
  
}


###########################
##build a dist density file 
ipt_dir <- 'opt0_acr_geneCate_dir_111623///'
opt_dir <- 'opt0_acr_gene_dist_dir_results_111923///'

all_fl_list <- list.files(path= ipt_dir)

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_acr_toGeneCate.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_target_dt <- read.delim(target_fl_path,header = F)
  head(ipt_target_dt)
  colnames(ipt_target_dt) <- c('acr','cate','dist','celltype')
  
  ipt_target_dt$logDist <- log2(ipt_target_dt$dist + 1)
  
  ##updating 020124 we will try to use the log10
  ipt_target_dt$log10length <- log10(ipt_target_dt$dist+1)
  
  p<-ggplot(ipt_target_dt, aes(x=log10length, fill=celltype, color=celltype)) +
    #geom_histogram(position="identity", alpha=0.5) + 
    geom_density(alpha = 0.5)+
    scale_color_manual(values=c("Broad"="#ED1C24", "CT"="#FFDE17")) +
    scale_fill_manual(values=c("Broad"="#ED1C24", "CT"="#FFDE17")) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
  #coord_cartesian(ylim = c(0,10000)) +
  #scale_y_continuous(breaks=seq(0,10000,2000)) +
  scale_x_discrete(name ="\nACR distance to nearest gene (kb)", 
                   limits=c(0,1,2,3,4,5),expand = expansion(mult = c(0, 0)),
                   labels=c("0", "0.01", "0.1",'1','10','100')) + 
  coord_cartesian(ylim = c(0, 0.8)) 
  #coord_cartesian(xlim = c(0,20000)) +
  #scale_x_continuous(breaks=seq(0,20000,5000))

  pdf(paste0(opt_dir,'/','opt_',target_fl_path_nm,'_distTogene','.pdf'),width = 7,height = 7)
  print(p)
  dev.off()

  
  ##updating 012024
  ##we will build a boxplot to show the difference of distance to the genes for the broad and CT ACRs
  p <- ggplot(ipt_target_dt, aes(x=celltype,y=logDist,fill=celltype)) + 
    geom_violin(trim=FALSE) + 
    #geom_boxplot(width=0.05,fill='white') +
    
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      legend.title = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    scale_fill_brewer(palette="Blues") 
   # facet_wrap(~species)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)
  
  pdf(paste0(opt_dir,'/','opt_',target_fl_path_nm,'_distTogene_violinPlot','.pdf'),width = 7,height = 7)
  print(p)
  dev.off()
  
  ipt_target_broad_dt <- ipt_target_dt[ipt_target_dt$celltype == 'Broad',]
  ipt_target_CT_dt <- ipt_target_dt[ipt_target_dt$celltype == 'CT',]
  
  test_sig <- wilcox.test(ipt_target_broad_dt$logDist,ipt_target_CT_dt$logDist,alternative = 'less')
  print(test_sig$p.value)
  
  
}






####################################
##build a cell type composition plot
ipt_dir <- 'opt0_acr_celltype_dir_111623//'
opt_dir <- 'opt0_acr_celltype_composition_dir_results_1119223/'

all_fl_list <- list.files(path= ipt_dir)

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_acr_celltype.txt','',target_fl_path_nm)
  
  ipt_target_dt <- read.delim(target_fl_path,header = F)
  if (ncol(ipt_target_dt) == 5){
    colnames(ipt_target_dt) <- c('chr','st','ed','celltype','pval')
  }
  if (ncol(ipt_target_dt) == 6){
    colnames(ipt_target_dt) <- c('chr','st','ed','celltype','pval','fc')
  }
  
  ipt_target_dt$celltype_modi <- gsub('.+;','',ipt_target_dt$celltype)
  
  dim(ipt_target_dt)
  nrow(ipt_target_dt[grepl(',',ipt_target_dt$celltype_modi),])
  
  nrow(ipt_target_dt)
  
  ipt_target_new_dt <- ipt_target_dt %>%
    separate_rows(celltype_modi, sep = ",")
  
  ipt_target_new_dt$final_celltype <- ifelse(grepl('unknown',ipt_target_new_dt$celltype_modi),'Unknown',ipt_target_new_dt$celltype_modi)
  table(ipt_target_new_dt$final_celltype)
  
  ipt_target_new_dt$count <- 1
  aggregate_dt <- aggregate(count ~ final_celltype,data = ipt_target_new_dt,sum)
  aggregate_dt$spe <- target_fl_path_nm
  
  return(as.data.frame(aggregate_dt))

})

combine_dt <- do.call(rbind,outs)
aggregate_combine_dt <- aggregate(count ~ spe, data = combine_dt, sum)
merged_dt <- merge(combine_dt,aggregate_combine_dt,by.x = 'spe',by.y = 'spe')
merged_dt$prop <- merged_dt$count.x/merged_dt$count.y

merged_dt$final_celltype <- gsub('companion_cells_sieve_elements','companion_cell',merged_dt$final_celltype)
table(merged_dt$final_celltype)
merged_dt$final_celltype <- gsub('procambium','procambial_meristem',merged_dt$final_celltype)
table(merged_dt$final_celltype)

merged_dt$spe <- factor(merged_dt$spe,levels = c('rice','Pm','Uf','maize','sorghum'))

p <- ggplot(data=merged_dt, aes(x=spe, y=prop, fill=final_celltype)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("broadly_accessible"="#fad9db", "bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                             "epidermis"="#C49B7B","mesophyll"="#91e388","procambial_meristem"="#FA8312","protoderm"="#fcd2b1","Unknown"="grey")) +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
pdf(paste0(opt_dir,'/opt_celltype_composition.pdf'),width = 8,height = 8)
print(p)
dev.off()


#############################
##build a motif coverage plot
ipt_dir <- 'opt0_motif_coverage_111923///'
opt_dir <- 'opt0_motif_coverage_results_111923//'

all_spe_dir_list <- list.files(path= ipt_dir)

for (i in (1:length(all_spe_dir_list))){
  
  target_fl_path_dir <- paste0(ipt_dir,'/',all_spe_dir_list[i])
  
  spe_nm <- gsub('.+/','',target_fl_path_dir)
  
  ipt_true_coverage_fl <- paste0(ipt_dir,'/',all_spe_dir_list[i],'/','opt_true_ave_motif_cover_dt.csv')
  ipt_control_coverage_fl <- paste0(ipt_dir,'/',all_spe_dir_list[i],'/','opt_control_ave_motif_cover_dt.csv')
  
  
  true_acr_dt <- read.csv(ipt_true_coverage_fl)
  #true_acr_dt <- read.csv('opt_true_ave_motif_cover_dt_modi.csv')
  head(true_acr_dt)
  dim(true_acr_dt)
  rownames(true_acr_dt) <- true_acr_dt$X
  true_acr_dt <- true_acr_dt[c('ave','error')]
  colnames(true_acr_dt) <- c('true_ave','true_error')
  
  control_acr_dt <- read.csv(ipt_control_coverage_fl)
  #control_acr_dt <- read.csv('opt_control_ave_motif_cover_dt_allbackground.csv')
  
  ##simulate
  #x<- runif(n = 4003, min = 1.72, max = 1.76)
  #x_dt <- as.data.frame(x)
  #merged <- cbind(control_acr_dt,x_dt)
  #control_acr_dt <- merged
  
  
  head(control_acr_dt)
  dim(control_acr_dt)
  rownames(control_acr_dt) <- control_acr_dt$X
  control_acr_dt <- control_acr_dt[c('ave','error')]
  #control_acr_dt <- control_acr_dt[c('x','error')]
  colnames(control_acr_dt) <- c('control_ave','control_error')
  
  combine <- cbind(true_acr_dt,control_acr_dt)
  head(combine)
  
  
  combine_ave <- combine[c('true_ave','control_ave')]
  head(combine_ave)
  
  pdf(paste0(opt_dir,'/',"opt_",spe_nm,"_avg_motif_coverage.pdf"), width=9.5, height=9.5)
  par(mar = c(5, 5, 5, 5))
  matplot(combine_ave, type="l", ylim=c(0,4), lwd=3,lty=1,col = c('firebrick1','grey'),
          #xaxt = 'n',
          ##previous is Distance from TSS
          xlab="\n\nDistance to ACR summit (kb)", ylab="Average motif coverage",
          cex.lab=2.5,cex.main = 2.5,cex.axis=2.5,labels =TRUE) ##lwd makes line thicker
  
  lines((combine$true_ave-combine$true_error), col="grey75")
  lines((combine$true_ave+combine$true_error), col="grey75")
  bp <- seq(length(combine$true_ave),1)
  polygon(c(bp,rev(bp)),c(rev(combine$true_ave+combine$true_error),combine$true_ave ),col=adjustcolor("grey",alpha.f=0.3), border=NA)
  polygon(c(bp,rev(bp)),c(rev(combine$true_ave),combine$true_ave-combine$true_error),col=adjustcolor("grey",alpha.f=0.3), border=NA)
  
  
  lines((combine$control_ave-combine$control_error), col="grey75")
  lines((combine$control_ave+combine$control_error), col="grey75")
  bp <- seq(length(combine$control_ave),1)
  polygon(c(bp,rev(bp)),c(rev(combine$control_ave+combine$control_error),combine$control_ave-combine$control_error),col=adjustcolor("grey",alpha.f=0.7), border=NA)
  
  
  dev.off()


}


################################################################
##check the motifs and their correlation of different cell types
##updating 112823
library(gtools)

ipt_meta_dir <- 'opt0_motif_based_correlation_112823/store_all_five_species_meta_112723/'
ipt_motifdev_dir <- 'opt0_motif_based_correlation_112823/store_all_five_species_motif_deviation_112723/'
opt_dir <- 'opt0_motif_based_correlation_results_112823/'
targetClust <- 'reduce_resolution_annotation'

all_spe_meta_fl_list <- list.files(path= ipt_meta_dir)

for (i in (1:length(all_spe_meta_fl_list))){
  
  target_meta_fl_path_fl <- paste0(ipt_meta_dir,'/',all_spe_meta_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_meta_fl_path_fl)
  target_fl_path_nm <- gsub('_meta.txt','',target_fl_path_nm)

  target_motifdev_fl_path_fl <- paste0(ipt_motifdev_dir,'/',target_fl_path_nm,'.motif.deviation.txt')
  
  a <- read.table(target_motifdev_fl_path_fl)
  a <- t(a)
  #a <- read.table(target_motifdev_fl_path_fl,stringsAsFactors = T)
  
  b <- read.delim(target_meta_fl_path_fl,header=T,row.names = 1)

  #rownames(b) <- b$cellID
  
  
  # format
  #a <- sparseMatrix(i=as.numeric(a$V1),j=as.numeric(a$V2),x=as.numeric(a$V3),dimnames=list(levels(a$V1),levels(a$V2)))
  
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
  
  ##remove rows with NA
  dim(mat)
  mat <- mat[complete.cases(mat),]
  
  write.table(mat, file=paste0(opt_dir,'/opt_',target_fl_path_nm,'_motif_average_clusters.txt'), quote=F, row.names=T, col.names=T, sep="\t")


}


##conduct the correlation between two with the rice
ipt_rice_motif_avg_matrix_fl <- paste0(paste0(opt_dir,'/opt_','rice','_motif_average_clusters.txt'))
ipt_rice_motif_avg_matrix_dt <- read.delim(ipt_rice_motif_avg_matrix_fl)

all_motif_avg_fl_list = list.files(path = opt_dir,pattern = 'motif_average_clusters.txt')

outs <- lapply(all_motif_avg_fl_list, function(x){
  
  target_spe_fl <- paste0(opt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_spe_fl)
  target_fl_path_nm <- gsub('_motif_average_clusters.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  target_spe_dt <- read.delim(target_spe_fl)
  
  ##check shared motifs
  shared_motifs <- intersect(rownames(ipt_rice_motif_avg_matrix_dt),rownames(target_spe_dt))
  ipt_rice_motif_avg_matrix_reduce_dt <- ipt_rice_motif_avg_matrix_dt[shared_motifs,]
  target_spe_reduce_dt <- target_spe_dt[shared_motifs,]
  
  cor_res <- cor(ipt_rice_motif_avg_matrix_reduce_dt,target_spe_reduce_dt)
  head(cor_res)
  
  
  ##transfer to z score
  ##version 1
  dt_zscore_v1 <- apply(cor_res,1,function(x){
    (x - mean(x))/sd(x)
  })
  cor_res_zscore <- data.frame(t(dt_zscore_v1),check.names = F)
  head(cor_res_zscore)
  
  ##plot the cor_res
  pdf(paste0(opt_dir,'/','opt_',target_fl_path_nm,'_corr_heatmap.pdf'),width = 14,height = 14)
  pheatmap(cor_res_zscore,
           #scale="row",
           #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
           #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
           cluster_cols = T,
           cluster_rows = T,
           border_color='white',
           fontsize_col = 10,
           fontsize_row = 10
           #annotation_row = annotdf,
           #annotation_colors = mycolors,
           #annotation_col = annotdf_col,
           #color=colorRampPalette(c("navy", "white", "firebrick3"))(50)
           #color = c(hcl.colors(16, "BluYl"),'yellow2')
           #color=colorRampPalette(c("navy", "firebrick3"))(50)
           
           #annotation_row = annotdf,
           #annotation_colors = mycolors,
           
           #annotation_col = annotdf_col,
           #annotation_colors = mycolors_col,
           #color=colorRampPalette(c("white", "firebrick3"))(50)
  )
  dev.off()
  
  ##Here we will allow the cor_res to be three colum
  cor_res_sparse <- Matrix(cor_res, sparse = TRUE)
  ia <- as.data.frame(summary(cor_res_sparse))
  ia$i <- rownames(cor_res_sparse)[as.numeric(ia$i)]
  ia$j <- colnames(cor_res_sparse)[as.numeric(ia$j)]
  ia$spe <- target_fl_path_nm

  return(ia)
  
})

combine_dt <- do.call(rbind,outs)
write.csv(combine_dt,paste0(opt_dir,'/','opt_compare_samecelltype_crosspec_corr.csv'),quote = F)
  
  
 


##############################################################
##check the motifs and TF correlation for all the four species
##updating 112823
ipt_TFcorr_dir <- 'opt0_motif_TF_corr_112823//'
opt_dir <- 'opt0_motif_TF_corr_results_112823/'

all_spe_TFcorr_fl_list <- list.files(path= ipt_TFcorr_dir)

outs <- lapply(all_spe_TFcorr_fl_list, function(x){
  
  
  target_TFcorr_fl_path_fl <- paste0(ipt_TFcorr_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_TFcorr_fl_path_fl)
  target_fl_path_nm <- gsub('_corresponding_motif_TF.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_TFcorr_fl_path_fl,header = F)
  ipt_dt[1:3,1:3]
  #colnames(ipt_dt) <- c('motif','TFgene','corr')
  
  ipt_dt$spe <- target_fl_path_nm
  
  return(ipt_dt)
  
})

combine_dt <- do.call(rbind, outs)

dim(combine_dt)
combine_dt[1:3,1:3]

combine_dt$PosOrNeg <- ifelse(combine_dt$V3 > 0, 'Posi','Neg')
dim(combine_dt)                                       

colnames(combine_dt) <- c('motif','TF','final_corr','other1','other2','species','PosOrNeg')

combine_dt$species <- factor(combine_dt$species,levels = c('rice','maize','sorghum','Pm'))

##violin plot
p <- ggplot(combine_dt, aes(x=PosOrNeg,y=final_corr,fill=PosOrNeg)) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.5,fill='white') +
  
  geom_jitter(shape=6, position=position_jitter(0.2)) +
  labs(x="\nClass", y = 'Correlation\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35),
    axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(~species)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)

p

combine_dt$count <- 1

aggregate_dt <- aggregate(count ~ species + PosOrNeg , data = combine_dt, sum)

pdf(paste0(opt_dir,'/opt_corr_TF_motif.pdf'),width = 10,height = 10)
print(p)
dev.off()















################################################################
##compare the conserved broad and restricted evolutionary scores
##for broad to broad  looks like not work so many 0 and the updated data > 0.05
ipt_dt <- read.delim('temp_intersect_rice_acr_physcore_celltypecate_sorghum.txt',header = F)
ipt_dt$cate <- paste0(ipt_dt$V4,'_',ipt_dt$V5) 



ipt_broadbroad_ACR <- ipt_dt[ipt_dt$V4 == 'broadACR'& ipt_dt$V5 == 'broadACR',]

ipt_broadbroad_ACR$riceACR <- paste0(ipt_broadbroad_ACR$V1,'_',ipt_broadbroad_ACR$V2,'_',ipt_broadbroad_ACR$V3)
ipt_broadbroad_ACR_mean <- aggregate(V10 ~ riceACR,data = ipt_broadbroad_ACR,mean)
ipt_broadbroad_ACR_mean$cate <- 'BroadToBroad'

##for restricted to restricted per cell types
ipt_dt <- read.delim('./temp_intersect_rice_acr_physcore_percelltype_sorghum.txt',header = F)
head(ipt_dt)
ipt_dt$riceACR <- paste0(ipt_dt$V1,'_',ipt_dt$V2,'_',ipt_dt$V3)

ipt_acr_celltype <- unique(ipt_dt[c('V1','V2','V3','V4','V5')])
ipt_acr_celltype$riceACR <-  paste0(ipt_acr_celltype$V1,'_',ipt_acr_celltype$V2,'_',ipt_acr_celltype$V3)
head(ipt_acr_celltype)

merged_dt <- merge(ipt_dt_mean, ipt_acr_celltype,by.x = 'riceACR',by.y = 'riceACR')
head(merged_dt)

##check real restricted and same cell type
outs <- lapply( unique(merged_dt$V4), function(x){
  
  if (x != c('companion_cell')){
    ipt_target_dt <- merged_dt[merged_dt$V4 == x,]
    ipt_target_dt <- ipt_target_dt[ipt_target_dt$V5 == x,]
  }else{
    ipt_target_dt <- merged_dt[merged_dt$V4 == x,]
    ipt_target_dt <- ipt_target_dt[ipt_target_dt$V5 == 'companion_cells_sieve_elements',]
  }

  ipt_final_target_dt <- ipt_target_dt[c('riceACR','V10')]
  
  return(ipt_final_target_dt)
  
})

combine_dt <- do.call(rbind,outs)
combine_dt$cate <- 'CTtoCT'

merged_broad_ct <- rbind(ipt_broadbroad_ACR_mean,combine_dt)
head(ipt_broadbroad_ACR_mean)
head(combine_dt)
combine_dt_CTtoCT <- merged_broad_ct[merged_broad_ct$cate == 'CTtoCT',]
combine_dt_broadtobroad <- merged_broad_ct[merged_broad_ct$cate == 'BroadToBroad',]
t.test(combine_dt_CTtoCT$V10,combine_dt_broadtobroad$V10,alternative = 'greater')
mean(combine_dt_CTtoCT$V10)
mean(combine_dt_broadtobroad$V10)
length(combine_dt_CTtoCT$V10)
length(combine_dt_broadtobroad$V10)

##plot the broad and CT using violin
p <- ggplot(merged_broad_ct, aes(x=cate, y=V10,fill=cate)) + 
  geom_violin(position=position_dodge(1))+
  geom_boxplot(width=0.05,fatten = 2,position=position_dodge(0.75),outlier.shape = NA)+
  #geom_boxplot() +
  #geom_jitter(shape=16, position=position_jitter(0.2)) + 
  labs(x="\nCate", y = 'Conserved score\n')+
  #ggtitle(paste(title_nm,'\n',sep = ''))+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
#  scale_fill_manual(values=c('aliceblue','indianred1'))
p <- p + coord_cartesian(ylim = c(0, 1.5))
p <- p + scale_y_continuous(breaks=seq(0, 1.5, 0.25))

pdf('opt1_final_CS_compare_conserved_Broad_to_CT_violin.pdf',width = 10,height = 6 )
p
dev.off()


##plot the broad and CT using density
p<-ggplot(merged_broad_ct, aes(x=V10, color=cate)) +
  geom_density(alpha=0.4)
p <- p + coord_cartesian(xlim = c(0, 2))
p <- p + scale_x_continuous(breaks=seq(0,2, 0.25))
p <- p + theme(
  plot.title = element_text(size=20,hjust = 0.5),
  axis.title = element_text(size =20),
  #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
  axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
  axis.text.y = element_text(size=20,colour = "black"),
  axis.ticks = element_line(size = rel(2.5)),
  axis.ticks.length = unit(0.5, "cm"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  text = element_text(size = 15),
  strip.text = element_text(size=20),
  #legend.title = element_blank(),
  #legend.position = "none",
  plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
)
p
pdf('opt_CS_all_cates_density_addcate.pdf',width = 10,height = 6 )
p
dev.off()



################################
##plot the motif enrichment test
ipt_TFcorr_dir <- 'opt0_motif_enrichment_120323///'
opt_dir <- 'opt0_motif_enrichment_results_120323//'

fc_cutoff <- 1
FDR <- 0.01
min_beta_cutoff <- 0.05
min_beta_cutoff <- 0

all_spe_TFcorr_fl_list <- list.files(path= ipt_TFcorr_dir,pattern = 'opt_motif_enrichment')

outs <- lapply(all_spe_TFcorr_fl_list, function(x){
  
  target_fl_path_nm <- gsub('opt_motif_enrichment_cluster_','',x)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  meta_name_dt <- read.delim(paste0(ipt_TFcorr_dir,'/ipt_required_raw_data/opt_motif_common_name.txt'),header = F)
  head(meta_name_dt)
  


  TFmotif_dt <- read.delim(paste0(ipt_TFcorr_dir,'/',x),row.names = 1)
  head(TFmotif_dt)
  TFmotif_dt <- TFmotif_dt[TFmotif_dt$fc > fc_cutoff,]
  TFmotif_dt <- TFmotif_dt[TFmotif_dt$FDR < FDR,]
  TFmotif_dt <- TFmotif_dt[TFmotif_dt$beta > min_beta_cutoff,]
  head(TFmotif_dt)
  dim(TFmotif_dt)
  
  ##add the common name
  TFmotif_dt <- merge(TFmotif_dt,meta_name_dt,by.x='motif',by.y= 'V1')
  head(TFmotif_dt)
  dim(TFmotif_dt)

  ##check the ori name
  TFmotif_dt_check_ori_name_ver <- TFmotif_dt
  head(TFmotif_dt_check_ori_name_ver)
  TFmotif_dt_check_ori_name_ver <- TFmotif_dt_check_ori_name_ver[c('motif','organct','beta')]
  
  write.table(TFmotif_dt_check_ori_name_ver,paste0(opt_dir,'/','temp_',target_fl_path_nm,'_target_col_motif_enrichment_oriMotifName.txt'),quote = F,sep = '\t')
  head(TFmotif_dt_check_ori_name_ver)
  TFmotif_dt_check_ori_name_ver <- read.table(paste0(opt_dir,'/','temp_',target_fl_path_nm,'_target_col_motif_enrichment_oriMotifName.txt'),stringsAsFactors = T)
  dim(TFmotif_dt_check_ori_name_ver)
  
  TFmotif_mtx <- sparseMatrix(i=as.numeric(TFmotif_dt_check_ori_name_ver$motif),
                              j=as.numeric(TFmotif_dt_check_ori_name_ver$organct),
                              x=as.numeric(TFmotif_dt_check_ori_name_ver$beta),
                              dimnames=list(levels(TFmotif_dt_check_ori_name_ver$motif), levels(TFmotif_dt_check_ori_name_ver$organct)))
  mat <- as(TFmotif_mtx, "dgCMatrix")
  dim(mat)
  head(mat)
  
  ##filter the min_beta_cutoff
  mat <- as.data.frame(as.matrix(mat))
  mat$max <- apply(mat, 1, max, na.rm=TRUE)
  mat <- mat[mat$max > min_beta_cutoff,]
  dim(mat)
  head(mat)
  mat <- mat[,c(-1)]
  mat <- as.matrix(mat)
  mat <- mat[, !colnames(mat) %in% 'max']
  
  as.data.frame(rownames(mat))
  
  ##save the motif name
  write.table(as.data.frame(rownames(mat)),paste0(opt_dir,'/opt_',target_fl_path_nm,'_target_motif_names.txt'),quote = F, sep = '\t')
  

  ##use the new motif name to do the plotting
  dim(TFmotif_dt)
  TFmotif_dt <- TFmotif_dt[c('V2','organct','beta')]
  head(TFmotif_dt)
  write.table(TFmotif_dt,paste0(opt_dir,'/','opt_',target_fl_path_nm,'_target_col_motif_enrichment.txt'),quote = F,sep = '\t')
  head(TFmotif_dt)
  
  TFmotif_dt <- read.table(paste0(opt_dir,'/','opt_',target_fl_path_nm,'_target_col_motif_enrichment.txt'),stringsAsFactors = T)
  dim(TFmotif_dt)
  head(TFmotif_dt)

  TFmotif_mtx <- sparseMatrix(i=as.numeric(TFmotif_dt$V2),
                              j=as.numeric(TFmotif_dt$organct),
                              x=as.numeric(TFmotif_dt$beta),
                              dimnames=list(levels(TFmotif_dt$V2), levels(TFmotif_dt$organct)))
  mat <- as(TFmotif_mtx, "dgCMatrix")
  dim(mat)
  
  
  pdf(paste0(opt_dir,'/','opt_',target_fl_path_nm,'_FDR',FDR,'_fc',fc_cutoff,'_meta',min_beta_cutoff,'_heat_map_fc_motif_enrich_target_black.pdf'),width = 25, height = 40)
  dim(mat)
  head(mat)

  #dim(signal_enrichment_matrix_reorder)
  dim(mat)
  p <- pheatmap(mat,
                #scale="row",
                #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
                #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
                cluster_cols = F,
                cluster_rows = F,
                border_color = "black",
                #border_color='black',
                fontsize_col = 8,
                fontsize_row = 8,
                #annotation_row = annotdf,
                #annotation_colors = mycolors,
                #annotation_col = annotdf_col,
                #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
                #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
                na_col = "white",
                #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
                #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
                #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
                #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
                color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
                
                #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
                #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
                #color=colorRampPalette(viridis)(100),
                
                ##second candidate
                #color = viridis(n = 256, alpha = 1, 
                #                     begin = 0, end = 1, option = "viridis"),
                
                #annotation_row = annotdf,
                #annotation_colors = mycolors,
                
                #annotation_col = annotdf_col,
                #annotation_colors = mycolors_col,
                #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
                #color = inferno(length(mat_breaks) - 1),
                #breaks = mat_breaks,
                
                ##open it when necessary
                #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
                fontsize_number = 25
  ) 
  p
  dev.off()
  
  
  ##Here we will check family beta 
  motif_name_dt <- read.delim(paste0(ipt_TFcorr_dir,'/ipt_required_raw_data/ipt_target_motifs_list_all.txt'),header = F)
  head(motif_name_dt)
  table(motif_name_dt$V3)
  
  merged_dt <- merge(TFmotif_dt,motif_name_dt, by.x = 'V2',by.y = 'V3')
  head(merged_dt)
  colnames(merged_dt) <- c('motifnm','organct','beta','motifID','motifam','select')
  
  table(TFmotif_dt$V2)
  
  final_dt <- merged_dt[c('motifam','motifnm','organct','beta')]
  final_dt$spe <- target_fl_path_nm  
  
  ##we will modify different families of motif
  table(final_dt$motifam)
  
  final_dt$motifam <- gsub('Group D','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group A','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group B','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group C','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group G','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group H','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group I','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group K','bZIP',final_dt$motifam)
  final_dt$motifam <- gsub('Group S','bZIP',final_dt$motifam)
  
  final_dt$motifam <- gsub('Myb-related','Myb',final_dt$motifam)
  final_dt$motifam <- gsub('WRKY-like_FRS/FRF','WRKY',final_dt$motifam)
  
  final_dt <- final_dt[final_dt$motifam != 'unknown',]
  final_dt <- final_dt[final_dt$motifam != 'Unknown',]
  
  table(final_dt$motifam)
  
  
  
  
  aggregate_dt <- aggregate(beta ~ motifam + organct, data = final_dt, mean)
  aggregate_dt$spe <- target_fl_path_nm
  
  return(final_dt)
  
  
})

combine_dt <- unique(do.call(rbind,outs))

##updating 122923
##check how many motifs were above cutoff for all the five species
combine_dt_check <- combine_dt
combine_dt_check$count <- 1
combine_dt_check$motifnm_organct <- paste0(combine_dt_check$motifnm,'_',combine_dt_check$organct)
combine_dt_check$motifnm_organct  <- gsub('_ncell.+','',combine_dt_check$motifnm_organct)
combine_dt_check$motifnm_organct  <- gsub('companion_cells_sieve_elements','companion_cell',combine_dt_check$motifnm_organct)
aggregate_combine_dt_check <- aggregate(count ~ motifnm_organct, data = combine_dt_check, sum)
table(aggregate_combine_dt_check$count)
table(combine_dt$spe)
write.table(aggregate_combine_dt_check,paste0(opt_dir,'/opt_all_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_motif_overlap_count.txt'),quote = F,sep = '\t')
write.table(combine_dt_check,paste0(opt_dir,'/opt_all_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_motif_overlap.txt'),quote = F,sep = '\t')

##updating 123023
##plot the target motifs across five species
target_modifs <- c('ANL2','HDG1','HDG11','CCA1','JUB1','LHY','CDF3','CDF5','DOF1.8','KAN2')
combine_dt_check <- combine_dt_check[combine_dt_check$motifnm %in% target_modifs,]
combine_dt_check <- combine_dt_check[c('spe','motifnm','organct','beta')]

combine_dt_check$organct_spe <- paste0(combine_dt_check$organct,'__',combine_dt_check$spe)

combine_select_dt <- combine_dt_check[c('organct_spe','motifnm','beta')]
write.table(combine_select_dt,paste0(opt_dir,'/temp_avg_beta_all_celltype_spe.txt'),quote = F, sep = '\t',row.names = F)

combine_select_dt <- read.delim(paste0(opt_dir,'/temp_avg_beta_all_celltype_spe.txt'),stringsAsFactors = T)

#combine_select_dt[combine_select_dt$motifam == 'Click to view details',]

combine_select_mtx <- sparseMatrix(i=as.numeric(combine_select_dt$organct_spe),
                                   j=as.numeric(combine_select_dt$motifnm),
                                   x=as.numeric(combine_select_dt$beta),
                                   dimnames=list(levels(combine_select_dt$organct_spe), levels(combine_select_dt$motifnm)))



head(combine_select_mtx)
colnames(combine_select_mtx)


##reorder
rowname_dt <- as.data.frame(rownames(combine_select_mtx))
rowname_dt$celltype <- gsub('__.+','',rowname_dt$`rownames(combine_select_mtx)`)
rowname_dt$celltype <- gsub('_ncell_.+','',rowname_dt$celltype)
rowname_dt$spe <- gsub('.+__','',rowname_dt$`rownames(combine_select_mtx)`)
row.names(rowname_dt) <- rowname_dt$`rownames(combine_select_mtx)`
rowname_dt$celltype <- gsub('companion_cells_sieve_elements','companion_cell',rowname_dt$celltype)
rowname_dt$celltype <- gsub('companion_cell','companion_cells_sieve_elements',rowname_dt$celltype)

rowname_dt$spe <- factor(rowname_dt$spe,levels = c('rice','maize','sorghum','Pm','Uf'))
rowname_dt <- rowname_dt[order(rowname_dt$spe),]
rowname_dt$celltype <- factor(rowname_dt$celltype,levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements','epidermis'))
rowname_dt <- rowname_dt[order(rowname_dt$celltype),]

##decide the motif showing
mat <- as.matrix(combine_select_mtx)
z <- t(as.matrix(scale(t(as.matrix(mat)))))
dim(z)

o.order <- colnames(z)
col.clust <- hclust(as.dist(1-cor(z)))
col.o <- col.clust$order
col.dendro <- as.dendrogram(col.clust)


z <- z[,col.o]
z <- z[order(apply(z, 1, which.max), decreasing=F),]
z <- z[,o.order]
head(z)
n.range <- c(-4, -2, 0, 2, 4)

colnames(combine_select_mtx)
motif_order <- c('CCA1','JUB1','LHY','CDF3','CDF5','DOF1.8','KAN2','ANL2','HDG1','HDG11')
##MA0972.1  MA2017.2  MA1185.2  MA0974.3  MA1268.2   MA0981.2   MA1383.2  MA1375.2  MA1369.2   MA0990.2
 

combine_select_mtx <- combine_select_mtx[row.names(rowname_dt),]


combine_select_mtx <- as.data.frame(as.matrix(combine_select_mtx))
combine_select_mtx <- combine_select_mtx[row.names(rowname_dt),motif_order]


pdf(paste0(opt_dir,'/','opt_target_motif_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_allctspe_heat_map.pdf'),width = 10, height = 10)
p <- pheatmap(t(as.matrix(combine_select_mtx)),
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
p
dev.off()


combine_dt_check_target <- combine_dt_check[combine_dt_check$motifnm %in% target_modifs,]
combine_dt_check_target$spe_organct <- 
combine_dt_check_target$new_organct <- gsub('_ncell.+','',combine_dt_check_target$organct)
combine_dt_check_target$new_organct  <- gsub('companion_cells_sieve_elements','companion_cell',combine_dt_check_target$new_organct)

combine_dt_check_target <- combine_dt_check_target[c('new_organct','motifnm','beta')]
write.table(combine_dt_check_target,paste0(opt_dir,'/temp_target_motif_beta_all_celltype_spe.txt'),quote = F, sep = '\t',row.names = F)

combine_dt_check_target <- read.table(paste0(opt_dir,'/temp_target_motif_beta_all_celltype_spe.txt'),stringsAsFactors = T,header = T)

combine_select_target_mtx <- sparseMatrix(i=as.numeric(combine_dt_check_target$new_organct),
                                   j=as.numeric(combine_dt_check_target$motifnm),
                                   x=as.numeric(combine_dt_check_target$beta),
                                   dimnames=list(levels(combine_dt_check_target$new_organct), levels(combine_dt_check_target$motifnm)))





#################################
##check the average motif heatmap
aggregate_dt <- aggregate(beta ~ spe + motifam + organct, data = combine_dt, mean)


aggregate_dt$organct_spe <- paste0(aggregate_dt$organct,'__',aggregate_dt$spe)

combine_select_dt <- aggregate_dt[c('organct_spe','motifam','beta')]
write.table(combine_select_dt,paste0(opt_dir,'/temp_avg_beta_all_celltype_spe.txt'),quote = F, sep = '\t',row.names = F)

combine_select_dt <- read.delim(paste0(opt_dir,'/temp_avg_beta_all_celltype_spe.txt'),stringsAsFactors = T)

#combine_select_dt[combine_select_dt$motifam == 'Click to view details',]

combine_select_mtx <- sparseMatrix(i=as.numeric(combine_select_dt$organct_spe),
                            j=as.numeric(combine_select_dt$motifam),
                            x=as.numeric(combine_select_dt$beta),
                            dimnames=list(levels(combine_select_dt$organct_spe), levels(combine_select_dt$motifam)))



head(combine_select_mtx)
colnames(combine_select_mtx)


##reorder
rowname_dt <- as.data.frame(rownames(combine_select_mtx))
rowname_dt$celltype <- gsub('__.+','',rowname_dt$`rownames(combine_select_mtx)`)
rowname_dt$celltype <- gsub('_ncell_.+','',rowname_dt$celltype)
rowname_dt$spe <- gsub('.+__','',rowname_dt$`rownames(combine_select_mtx)`)
row.names(rowname_dt) <- rowname_dt$`rownames(combine_select_mtx)`
rowname_dt$celltype <- gsub('companion_cells_sieve_elements','companion_cell',rowname_dt$celltype)
rowname_dt$celltype <- gsub('companion_cell','companion_cells_sieve_elements',rowname_dt$celltype)

rowname_dt$spe <- factor(rowname_dt$spe,levels = c('rice','maize','sorghum','Pm','Uf'))
rowname_dt <- rowname_dt[order(rowname_dt$spe),]
rowname_dt$celltype <- factor(rowname_dt$celltype,levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements','epidermis'))
rowname_dt <- rowname_dt[order(rowname_dt$celltype),]

##decide the motif showing
mat <- as.matrix(combine_select_mtx)
z <- t(as.matrix(scale(t(as.matrix(mat)))))
dim(z)

o.order <- colnames(z)
col.clust <- hclust(as.dist(1-cor(z)))
col.o <- col.clust$order
col.dendro <- as.dendrogram(col.clust)


z <- z[,col.o]
z <- z[order(apply(z, 1, which.max), decreasing=F),]
z <- z[,o.order]
head(z)
n.range <- c(-4, -2, 0, 2, 4)

colnames(combine_select_mtx)
motif_order <- c('NAC','Myb','WRKY','HD-ZIP','PLINC','bZIP','ERF/DREB','SBP',
                 'C4-GATA-related','Trihelix',
                 'TCP','DOF',
                 'MIKC','GARP_ARR-B','GARP_G2-like',
                 'GRF','RAV',
                 'ABI3','ARF','E2F','BES/BZR','AP2','LBD','Click to view details')

combine_select_mtx <- combine_select_mtx[row.names(rowname_dt),]


combine_select_mtx <- as.data.frame(as.matrix(combine_select_mtx))
combine_select_mtx <- combine_select_mtx[row.names(rowname_dt),motif_order]


pdf(paste0(opt_dir,'/','opt_all_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_allctspe_heat_map.pdf'),width = 10, height = 10)
p <- pheatmap(t(as.matrix(combine_select_mtx)),
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
p
dev.off()

##transfer the the matrix to be three col
final_matrix <- t(as.matrix(combine_select_mtx))
head(final_matrix)

final_matrix <- as(final_matrix, "sparseMatrix")

ia <- as.data.frame(summary(final_matrix))
ia$i <- rownames(final_matrix)[as.numeric(ia$i)]
ia$j <- colnames(final_matrix)[as.numeric(ia$j)]
write.table(ia, file=paste0(opt_dir,'/','opt_all_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_allctspe_heat_map.sparse'), quote=F, row.names=F, col.names=F, sep="\t")

head(ia)

##updating 121323
##build a boxplot to show the average beta of all TF families in companion cells compared to the others
final_dt <- as.data.frame(ia)
colnames(final_dt) <- c('motif','celltype','score')
#companion_cell_dt <- final_dt[grepl ('companion_cells_sieve_elements',final_dt$j)|grepl ('companion_cell',final_dt$j),]
#colnames(companion_cell_dt) <- c('motif','celltype','score')


final_dt$celltype <- gsub('companion_cells_sieve_elements','companion_cell',final_dt$celltype)
final_dt$celltype <- gsub('companion_cell','companion_cells_sieve_elements',final_dt$celltype)
final_dt$celltype <- gsub('companion_cells_sieve_elements_ncell_1193','companion_cells_sieve_elements',final_dt$celltype)
final_dt$celltype <- gsub('bundle_sheath_ncell_7219','bundle_sheath',final_dt$celltype)
final_dt$celltype <- gsub('epidermis_ncell_5440','epidermis',final_dt$celltype)
final_dt$celltype <- gsub('mesophyll_ncell_5257','mesophyll',final_dt$celltype)

final_dt$spe <- gsub('.+__','',final_dt$celltype)
final_dt$celltype_true <- gsub('__.+','',final_dt$celltype)

final_dt$celltype_true <- factor(final_dt$celltype_true,levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements','epidermis'))

final_dt$spe <- factor(final_dt$spe,levels = c('rice','maize','sorghum','Pm','Uf'))

final_dt <- final_dt[order(final_dt$spe),]
final_dt <- final_dt[order(final_dt$celltype_true),]
unique(final_dt$celltype)
final_dt$celltype <- factor(final_dt$celltype,levels = unique(final_dt$celltype))


p<-ggplot(final_dt, aes(x=celltype, y=score, fill=celltype)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("bundle_sheath__rice"="#3084AF","bundle_sheath__maize"="#3084AF","bundle_sheath__sorghum"="#3084AF","bundle_sheath__Pm"="#3084AF","bundle_sheath__Uf"="#3084AF",
                             "companion_cells_sieve_elements__rice"="#eddc9f","companion_cells_sieve_elements__maize"="#eddc9f","companion_cells_sieve_elements__sorghum"="#eddc9f","companion_cells_sieve_elements__Pm"="#eddc9f","companion_cells_sieve_elements__Uf"="#eddc9f",
                            "epidermis__rice"="#C49B7B","epidermis__maize"="#C49B7B","epidermis__sorghum"="#C49B7B","epidermis__Pm"="#C49B7B","epidermis__Uf"="#C49B7B",
                            "mesophyll__rice"="#91e388", "mesophyll__maize"="#91e388", "mesophyll__sorghum"="#91e388", "mesophyll__Pm"="#91e388", "mesophyll__Uf"="#91e388",
                            "protoderm__rice"="#fcd2b1","protoderm__maize"="#fcd2b1","protoderm__sorghum"="#fcd2b1","protoderm__Pm"="#fcd2b1","protoderm__Uf"="#fcd2b1"))+
  theme(
    plot.title = element_text(size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
pdf(paste0(opt_dir,'/','opt_all_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_allctspe_boxplot.pdf'),width = 20,height = 12 )
print(p)
dev.off()

##calculate the number of motif per cell type per spe
final_dt$count <- 1

aggregate_dt <- aggregate(count ~ spe + celltype_true, data = final_dt, sum)
write.table(aggregate_dt,paste0(opt_dir,'/','opt_all_FDR',FDR,'_fc',fc_cutoff,'_beta',min_beta_cutoff,'_allctspe_boxplot_num.txt'),sep= '\t',quote = F)







##plot a bar plot to show the count of spe number for motif count per cell type 
ipt_dt <- read.delim(paste0(ipt_TFcorr_dir,'/opt_summary_TFfam_celltype_speCount_num.txt'),header = F)
head(ipt_dt)

colnames(ipt_dt) <- c('celltype','specount','number')


ipt_dt$celltype <- factor(ipt_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements','epidermis'))

p <- ggplot(data=ipt_dt, aes(x=celltype, y=number, fill = specount)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, vjust = 1,hjust =1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
pdf(paste0(opt_dir,'/opt_summary_TFfam_celltype_speCount_num.pdf'),width = 10,height = 7)
p
dev.off()





##check the correlation per epidermis not work
combine_dt_epidermis <- combine_dt[combine_dt$organct == 'companion_cells_sieve_elements'|combine_dt$organct == 'companion_cell',]
combine_dt_epidermis <- combine_dt[combine_dt$organct == 'bundle_sheath',]
combine_dt_epidermis <- combine_dt[combine_dt$organct == 'mesophyll',]
combine_dt_epidermis <- combine_dt_epidermis[c('spe','motifnm','beta')]
write.table(combine_dt_epidermis,paste0(opt_dir,'/temp_table.txt'),quote = F,sep = '\t')

combine_dt_epidermis <- read.table(paste0(opt_dir,'/temp_table.txt'),stringsAsFactors = T)
combine_dt_mtx <- sparseMatrix(i=as.numeric(combine_dt_epidermis$spe),
                                   j=as.numeric(combine_dt_epidermis$motifnm),
                                   x=as.numeric(combine_dt_epidermis$beta),
                                   dimnames=list(levels(combine_dt_epidermis$spe), levels(combine_dt_epidermis$motifnm)))
combine_dt_mtx <- t(combine_dt_mtx)
mydata.cor = cor(as.matrix(combine_dt_mtx))
mydata.cor

combine_dt_epidermis <- combine_dt_epidermis[c('spe','beta')]

combine_dt_maize_epidermis <- combine_dt_epidermis[combine_dt_epidermis$spe == 'maize',]
combine_dt_rice_epidermis <- combine_dt_epidermis[combine_dt_epidermis$spe == 'rice',]
dim(combine_dt_maize_epidermis)
dim(combine_dt_rice_epidermis)

cor(combine_dt_maize_epidermis$beta,combine_dt_rice_epidermis$beta)

library("Hmisc")
res2 <- rcorr(as.matrix(my_data))
res2


##########
##updating 121923
##build the ACR overlaping with the cns
ipt_dir <- 'opt0_check_cns_ACR_121923/Pablo_CNSs_version/'
opt_dir <- 'opt0_check_cns_ACR_results_121923/Pablo_CNSs_version/'
wide1 <- 15
wide2 <- 12

ipt_dir <- 'opt0_check_cns_ACR_121923/Database_CNSs_version/not_contain_UTR_ver_choose/'
opt_dir <- 'opt0_check_cns_ACR_results_121923/Database_CNSs_version/not_contain_UTR_choose/'
wide1 <- 10
wide2 <- 8

all_fl_list <- list.files(path= ipt_dir, pattern = '_cns_num.txt')


outs <- lapply(all_fl_list, function(x){
  
  target_fl_path_nm <- gsub('_celltype_cns_num.txt','',x)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(paste0(ipt_dir,'/',x),header = F)
  head(ipt_dt)
  colnames(ipt_dt) <- c('chr','st','ed','celltype_cate','CNS','geneCate')
  
  ipt_dt$acr <- paste0(ipt_dt$chr,'_',ipt_dt$st,'_',ipt_dt$ed)  
  ipt_dt$cate <- ifelse(ipt_dt$CNS == 0, 'notCoverCNS','CoverCNS')
  
  dim(ipt_dt)
  ipt_dt$count <- 1
  ipt_dt$celltype <- ifelse(ipt_dt$celltype_cate == 'broadly_accessible','Broad','CT')
  

  head(ipt_dt)
  aggregate_dt <- aggregate(count ~ celltype + cate,data = ipt_dt, sum)
  
  aggregate_celltype_dt <- aggregate(count ~ celltype, data = ipt_dt, sum)
  mergedt_dt <- merge(aggregate_dt,aggregate_celltype_dt,by.x = 'celltype',by.y = 'celltype')
  mergedt_dt$prop <- mergedt_dt$count.x/mergedt_dt$count.y
  
  mergedt_dt$spe <- target_fl_path_nm
  
  return(mergedt_dt)
    
})

combine_dt <- do.call(rbind,outs)



combine_dt$spe <- factor(combine_dt$spe,levels = c('rice','maize','sorghum','Pm','Uf'))
#combine_dt$spe <- factor(combine_dt$spe,levels = c('rice','maize'))

head(combine_dt)
p <- ggplot(combine_dt, aes(fill=cate, y=prop, x=celltype)) + 
  geom_bar(position="stack", stat="identity")
#p <- p + geom_text(aes(label = Proportion),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
p <- p+labs(x="\nCell type", y = "Proportion\n", cex.lab = 7) ## use \n to set the position
#p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 30), ##change all the text size
             axis.text.y = element_text(size=40,colour = "black"),
             axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
             axis.title.x = element_text(size=40,colour = "black"),
             axis.title.y = element_text(size=40,colour = "black"),
             plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
p <- p + facet_wrap(~spe,nrow = 1)
p

p <- p + coord_cartesian(ylim = c(0, 1)) + ##threshold
  scale_y_continuous(breaks=seq(0, 1, 0.1))
pdf(paste0(opt_dir,'/opt_all_spe_cns_prop.pdf'),width = wide1, height = 15)
p ##20 15
dev.off()


##updating 122223
##only plot the cns in CT proportion
combine_dt_CNS <- combine_dt[combine_dt$cate == 'CoverCNS',]
combine_dt_CNS$spe <- factor(combine_dt_CNS$spe, levels = c('rice','maize','sorghum','Pm','Uf'))
p <- ggplot(data=combine_dt_CNS, aes(x=spe, y=prop, fill=celltype)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 30), ##change all the text size
             axis.text.y = element_text(size=40,colour = "black"),
             axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
             axis.title.x = element_text(size=40,colour = "black"),
             axis.title.y = element_text(size=40,colour = "black"),
             plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
pdf(paste0(opt_dir,'/opt_all_spe_cns_prop_CNS.pdf'),width = wide2, height = 10)
p
dev.off()



# Use custom colors
#p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
#p + scale_fill_brewer(palette="Blues")





##updating 122223
##check if the CNS within CT more sig overlapping other than the broad
all_spe_list <- unique(combine_dt$spe)
outs <- lapply(all_spe_list, function(x){
  
  print(x)
  ipt_target_dt <- combine_dt[combine_dt$spe == x,]
  
  ipt_target_dt_broad <- ipt_target_dt[ipt_target_dt$celltype == 'Broad',]
  ipt_target_dt_broad_cns <- ipt_target_dt_broad[ipt_target_dt_broad$cate == 'CoverCNS',]
  ipt_target_dt_broad_not_cns <- ipt_target_dt_broad[ipt_target_dt_broad$cate == 'notCoverCNS',]
  
  ipt_target_dt_CT <- ipt_target_dt[ipt_target_dt$celltype == 'CT',]
  ipt_target_dt_CT_cns <- ipt_target_dt_CT[ipt_target_dt_CT$cate == 'CoverCNS',]
  ipt_target_dt_CT_not_cns <- ipt_target_dt_CT[ipt_target_dt_CT$cate == 'notCoverCNS',]
  
  ct_in_cns_num = ipt_target_dt_CT_cns$count.x
  ct_not_cns_num = ipt_target_dt_CT_not_cns$count.x
  not_ct_in_cns_num = ipt_target_dt_broad_cns$count.x
  not_ct_not_cns_num = ipt_target_dt_broad_not_cns$count.x
  
  my_data <- data.frame(
    Group1 = c(ct_in_cns_num, not_ct_in_cns_num),
    Group2 = c(ct_not_cns_num, not_ct_not_cns_num)
  )
  
  # Convert the data frame to a matrix
  data_matrix <- as.matrix(my_data)
  
  # Perform Fisher's exact test
  result <- fisher.test(data_matrix,alternative = 'greater')
  
  res <- c()
  res$pval = result$p.value
  res$spe = x
  
  res <- as.data.frame(res)
  return(res)

})

combine_enrich_dt <- do.call(rbind,outs)







##check number of rice atlas acr overlapping with the cns
ipt_dt <- read.delim(paste0(ipt_dir,'/opt_rice_atlas_acr_addPabloLeaf_addcns.txt'),header = F)
head(ipt_dt)
dim(ipt_dt)
length(unique(ipt_dt$V6))
ipt_dt$CNScate <- ifelse(ipt_dt$V7 != '0','CNS','notCNS')
ipt_dt$count <- 1
ipt_dt$leafcate <- ifelse(ipt_dt$V6 == 'none','OtherACR','LeafACR')
aggregate_dt <- aggregate(count ~ leafcate + CNScate, data = ipt_dt, sum)



##########
##updating 122323
##build the cns that could be captured by the ACR
ipt_dir <- 'opt0_check_cns_ACR_121923/Pablo_CNSs_version/cns_subject_122323/'
opt_dir <- 'opt0_check_cns_ACR_results_121923/Pablo_CNSs_version/'

ipt_dir <- 'opt0_check_cns_ACR_121923/Database_CNSs_version//not_contain_UTR_ver_choose/cns_subject_010524/'
opt_dir <- 'opt0_check_cns_ACR_results_121923/Database_CNSs_version/not_contain_UTR_choose/'





all_fl_list <- list.files(path= ipt_dir, pattern = 'celltype.txt')

generate_pie <- function(acr_num,total_acr_num) {
  #df <- dt[dt$celltype == target_celltype,]
  #pie(df$number)
  #library(RColorBrewer)
  other_acr_num = total_acr_num - acr_num
  pct <- round(acr_num/total_acr_num*100)
  lbls <- paste(c('ACRinSyntenic','OtherACR'),pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(c(acr_num,other_acr_num),labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path_nm <- gsub('_cns_subject_acr_celltype.txt','',x)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(paste0(ipt_dir,'/',x),header = F)
  head(ipt_dt)
  colnames(ipt_dt) <- c('chr','st','ed','acr','celltype_cate','celltype')
  
  ipt_dt$cns <- paste0(ipt_dt$chr,'_',ipt_dt$st,'_',ipt_dt$ed)  
  ipt_dt$cate <- ifelse(ipt_dt$acr == 'none', 'notCoverACR','CoverACR')
  
  
  dim(ipt_dt)
  ipt_dt$count <- 1
  
  
  table(ipt_dt$celltype)
 
  head(ipt_dt)
  aggregate_dt <- aggregate(count ~ celltype + cate,data = ipt_dt, sum)
  aggregate_dt <- aggregate_dt[aggregate_dt$cate == 'CoverACR',]

  aggregate_dt$spe <- target_fl_path_nm
  
  ##draw the pie chart  
  aggregate_cate <- aggregate(count ~ cate, data = ipt_dt,sum)
  
  total_cns <- aggregate_cate[aggregate_cate$cate == 'CoverACR',]$count + aggregate_cate[aggregate_cate$cate == 'notCoverACR',]$count
  coverACR_cns <- aggregate_cate[aggregate_cate$cate == 'CoverACR',]$count
  
  pdf(paste0(opt_dir,'/',target_fl_path_nm,'_cns_prop_coverACR_notcoverACR','.pdf'),width = 7,height = 7)
  print(generate_pie(coverACR_cns,total_cns))
  dev.off()
  
  print(paste0(target_fl_path_nm,'_total_cns_',total_cns))
  
  
  return(aggregate_dt)

})

combine_dt <- do.call(rbind,outs)


##updating 010524
##here we will manually generate a piechart figure to show the atlas overlapping
coverACR_cns <- 34667
total_cns <- 53253
target_fl_path_nm <- 'rice_atlas'
pdf(paste0(opt_dir,'/',target_fl_path_nm,'_cns_prop_coverACR_notcoverACR','.pdf'),width = 7,height = 7)
print(generate_pie(coverACR_cns,total_cns))
dev.off()






#################
##updating 123123
##we will plot the acr in syn and not in syn
ipt_dir <- 'opt0_rice_altas_rice_syn_nonsyn_123123//'
opt_dir <- 'opt0_rice_altas_rice_syn_nonsyn_results_123123//'

ipt_two_cate_list <- c('Syn','NoSyn')

outs <- lapply(ipt_two_cate_list, function(x){
  
  ipt_dt <- read.delim(paste0(ipt_dir,'/opt_', x,'_in_diffbins_plot_R.txt'),header = F)
  
  head(ipt_dt)
  colnames(ipt_dt) <- c('bin','cate','prop')
  ipt_dt <- ipt_dt[ipt_dt$bin != '[61,65]',]
  
  bin_list <- unique(ipt_dt$bin)
  outs2 <- lapply(bin_list, function(y){
    
    ipt_target_dt <- ipt_dt[ipt_dt$bin == y,]
    
    ipt_dt_simulate <- ipt_target_dt[ipt_target_dt$cate == 'control',]
    sd <- sd(ipt_dt_simulate$prop)
    avg <- mean(ipt_dt_simulate$prop)
    ipt_dt_real <- ipt_target_dt[ipt_target_dt$cate == 'real',]
    
    ipt_new_dt <- data.frame(
      cate = c('real', 'simulate'),  # Numeric values for Column1
      prop = c(ipt_dt_real$prop, avg),
      sd = c(0,sd)# Character values for Column2
    )
    
    ipt_new_dt$bin <- y
    return(ipt_new_dt)
    
  })
  
  combine_dt <- do.call(rbind,outs2)
  combine_dt$bin <- factor(combine_dt$bin,levels = bin_list)
  combine_real_dt <- combine_dt[combine_dt$cate == 'real',]
  
  combine_real_dt$SynOrNsynCate <- x
  
  return(combine_real_dt)
  
})

final_combine_dt <- do.call(rbind,outs)


p <- ggplot(data=final_combine_dt, aes(x=bin, y=prop, fill = SynOrNsynCate)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  #geom_bar(stat="identity", fill="steelblue")+
  #geom_text(aes(label=number), vjust=-0.3, size=3.5)+
  #geom_errorbar(aes(ymin=prop - sd, ymax=prop + sd), width=0.4, position=position_dodge(0.9)) +  # Add error bars
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    axis.text.x = element_text(colour = "black", size=15,angle = 45, vjust = 1,hjust =1),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  coord_cartesian(ylim = c(0, 0.25))


pdf(paste0(opt_dir,'/', 'opt_', prefix,'_ACR_diffbins_only_real.pdf'),width = 12,height = 5)
p
dev.off()

##plot the barplot to check the cell-type-specific and not in two cates
ipt_acr_syn_or_not_syn_dt <- read.delim('opt0_rice_altas_rice_syn_nonsyn_123123/opt_rice_atlas_acr_add_syn_or_not.txt',header = F)
dim(ipt_acr_syn_or_not_syn_dt)
table(ipt_acr_syn_or_not_syn_dt$V6)

not_syn_dt <- ipt_acr_syn_or_not_syn_dt[ipt_acr_syn_or_not_syn_dt$V6 == 'NonSyn',]
ct_not_syn_num <- nrow(not_syn_dt[not_syn_dt$V5 != 'none',])
ct_not_syn_num

syn_dt <- ipt_acr_syn_or_not_syn_dt[ipt_acr_syn_or_not_syn_dt$V6 != 'NonSyn',]
ct_syn_num <- nrow(syn_dt[syn_dt$V5 != 'none',])
ct_syn_num

broad_not_syn_num <- nrow(not_syn_dt[not_syn_dt$V5 == 'none',])
broad_syn_num <- nrow(syn_dt[syn_dt$V5 == 'none',])

data_sig <- matrix(c(ct_not_syn_num, broad_not_syn_num, ct_syn_num, broad_syn_num), nrow = 2)
result <- fisher.test(data_sig,alternative = 'greater')
result <- fisher.test(data_sig)
result$p.value
##4.817842e-215
fisher.test()

##plot the barplot
prop_data <- data.frame(
  cate <- c('CTinSynProp','CTinNsynProp'),
  prop <- c(ct_syn_num/nrow(syn_dt), ct_not_syn_num/nrow(not_syn_dt))
)

colnames(prop_data) <- c('cate','prop')

p <- ggplot(data=prop_data, aes(x=cate, y=prop, fill=cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_y_continuous(breaks=seq(0, 0.4, 0.1))
#facet_wrap(~Rice_CellType,nrow = 1)
pdf(paste0(opt_dir,'/opt_compare_CT_prop_insyn_to_noninsyn.pdf'),width = 6,height = 6 )
p
dev.off()


##updating 020124
##we will check the different categories of genic proximal and distal for the comaprision 
ipt_acr_syn_or_not_syn_dt <- read.delim('opt0_rice_altas_rice_syn_nonsyn_123123/opt_rice_atlas_acr_add_syn_or_not.txt',header = F)
head(ipt_acr_syn_or_not_syn_dt)

all_cate_list <- unique(ipt_acr_syn_or_not_syn_dt$V4)

##conduct the enrichment test
outs <- lapply(all_cate_list, function(x){
  
  ipt_target_dt <- ipt_acr_syn_or_not_syn_dt[ipt_acr_syn_or_not_syn_dt$V4 == x,]
  
  ##conduct the enrichment test
  not_syn_dt <- ipt_target_dt[ipt_target_dt$V6 == 'NonSyn',]
  ct_not_syn_num <- nrow(not_syn_dt[not_syn_dt$V5 != 'none',])
  ct_not_syn_num
  
  syn_dt <- ipt_target_dt[ipt_target_dt$V6 != 'NonSyn',]
  ct_syn_num <- nrow(syn_dt[syn_dt$V5 != 'none',])
  ct_syn_num
  
  broad_not_syn_num <- nrow(not_syn_dt[not_syn_dt$V5 == 'none',])
  broad_syn_num <- nrow(syn_dt[syn_dt$V5 == 'none',])
  
  data_sig <- matrix(c(ct_not_syn_num, broad_not_syn_num, ct_syn_num, broad_syn_num), nrow = 2)
  result <- fisher.test(data_sig,alternative = 'two.sided')
  result <- fisher.test(data_sig)
  
  res <- c()
  res$cate <- x
  res$pval <- result$p.value
  
  return(res)

})

combine_dt <- do.call(rbind,outs)


##plot the bar plot
outs <- lapply(all_cate_list, function(x){
  
  ipt_target_dt <- ipt_acr_syn_or_not_syn_dt[ipt_acr_syn_or_not_syn_dt$V4 == x,]
  
  ##conduct the enrichment test
  not_syn_dt <- ipt_target_dt[ipt_target_dt$V6 == 'NonSyn',]
  ct_not_syn_num <- nrow(not_syn_dt[not_syn_dt$V5 != 'none',])
  ct_not_syn_num
  
  syn_dt <- ipt_target_dt[ipt_target_dt$V6 != 'NonSyn',]
  ct_syn_num <- nrow(syn_dt[syn_dt$V5 != 'none',])
  ct_syn_num
  
  broad_not_syn_num <- nrow(not_syn_dt[not_syn_dt$V5 == 'none',])
  broad_syn_num <- nrow(syn_dt[syn_dt$V5 == 'none',])
  
  ##plot the barplot
  prop_data <- data.frame(
    prox_cate <- c(x,x),
    cate <- c('CTinSynProp','CTinNsynProp'),
    prop <- c(ct_syn_num/nrow(syn_dt), ct_not_syn_num/nrow(not_syn_dt))
  )

  colnames(prop_data) <- c('prox_cate','cate','prop')
  
  return(prop_data)
})

combine_dt <- do.call(rbind,outs)

combine_dt$prox_cate <- factor(combine_dt$prox_cate,levels = c('Genic','Proximal','Distal'))

p <- ggplot(data=combine_dt, aes(x=prox_cate, y=prop, fill=cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
#facet_wrap(~Rice_CellType,nrow = 1)
pdf(paste0(opt_dir,'/opt_compare_CT_prop_insyn_to_noninsyn_diff_proxCate.pdf'),width = 8,height = 6 )
p
dev.off()








########
##Step01
########
########################################
##build the pie chart to show the number 
ipt_dir <- 'opt1_acr_syntenic_summary_dir_111923/'
opt_dir <- 'opt1_acr_syntenic_summary_dir_results_111923/'

generate_pie <- function(acr_num,total_acr_num) {
  #df <- dt[dt$celltype == target_celltype,]
  #pie(df$number)
  #library(RColorBrewer)
  other_acr_num = total_acr_num - acr_num
  pct <- round(acr_num/total_acr_num*100)
  lbls <- paste(c('ACRinSyntenic','OtherACR'),pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(c(acr_num,other_acr_num),labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}

library(RColorBrewer)
myPalette <- brewer.pal(10, "Set3") 

all_fl_list <- list.files(path= ipt_dir)

rice_acr_num = 62669

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_H3K27me3_.+','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header=T)
  
  acr_in_syntenic_block_num <- length(unique(ipt_dt$rice_ACR))
  
  print(target_fl_path_nm)
  print(acr_in_syntenic_block_num)

  pdf(paste0(opt_dir,'/',target_fl_path_nm,'_prop_acr_in_syntenic','.pdf'),width = 7,height = 7)
  print(generate_pie(acr_in_syntenic_block_num,rice_acr_num))
  dev.off()
  
  
}

##updating 121823 
##build the composition bar plot to show the data
spe <- c("maize", "maize", "sorghum","sorghum","Pm","Pm",'Uf',"Uf")
syn <- c('Insyn', 'NotSyn','Insyn', 'NotSyn','Insyn', 'NotSyn','Insyn', 'NotSyn')
prop <- c(18, 82, 32, 68, 24, 76, 34, 66)

# Creating a dataframe
my_dataframe <- data.frame(Spe = spe, Syn = syn, Prop = prop)


my_dataframe$Spe <- factor(my_dataframe$Spe, levels = c('maize','sorghum','Pm','Uf'))

p <- ggplot(my_dataframe, aes(fill=Syn, y=Prop, x=Spe)) + 
  geom_bar(position="stack", stat="identity")
#p <- p + geom_text(aes(label = Prop),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
p <- p+labs(x="\nCell type", y = "Number of feature / Total number of feature\n", cex.lab = 7) ## use \n to set the position
#p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 30), ##change all the text size
             axis.text.y = element_text(size=40,colour = "black"),
             axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
             axis.title.x = element_text(size=40,colour = "black"),
             axis.title.y = element_text(size=40,colour = "black"),
             plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer

p <- p + coord_cartesian(ylim = c(0, 100)) + ##threshold
  scale_y_continuous(breaks=seq(0, 100, 20))

pdf(paste0(opt_dir,'/','opt_all_spe','_prop_acr_in_syntenic_barplot','.pdf'),width = 10,height = 12)
p
dev.off()


#################
##updating 123123
##Here we will check the non-syntenic ACRs for the cell type
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123//'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123//'

all_fl_list <- list.files(path= ipt_dir,pattern = '.txt')

spe1 <- 'rice'

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_H3K27me3_.+','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  if (spe1 == 'rice'){
    spe2_nm <- gsub('rice_','',target_fl_path_nm)
  }else{
    spe2_nm <- 'rice'
  }
  
  ipt_dt <- read.delim(target_fl_path,header=T)  
  
  ##forcus on the non-syntenic regions
  ipt_nonsyn_dt <- ipt_dt[ipt_dt[[paste0(spe2_nm,'_region')]] == 'none',]
  dim(ipt_nonsyn_dt)
  ipt_syn_dt <- ipt_dt[ipt_dt[[paste0(spe2_nm,'_region')]] != 'none',]
  dim(ipt_syn_dt)
  
  
  ipt_nonsyn_ACRnum <- length(unique(ipt_nonsyn_dt$rice_ACR))
  ipt_syn_ACRnum <- length(unique(ipt_syn_dt$rice_ACR))
  
  ##not syntenic broad the ct proportion
  ipt_nonsyn_broad_dt <- ipt_nonsyn_dt[ipt_nonsyn_dt$rice_CelltypeCate == 'broadly_accessible',]
  ipt_nonsyn_broad_ACRnum <- length(unique(ipt_nonsyn_broad_dt$rice_ACR))
  ipt_nonsyn_CT_dt <-  ipt_nonsyn_dt[ipt_nonsyn_dt$rice_CelltypeCate != 'broadly_accessible',]
  ipt_nonsyn_CT_ACRnum <- length(unique(ipt_nonsyn_CT_dt$rice_ACR))
  
  ipt_syn_broad_dt <- ipt_syn_dt[ipt_syn_dt$rice_CelltypeCate == 'broadly_accessible',]
  ipt_syn_broad_ACRnum <- length(unique(ipt_syn_broad_dt$rice_ACR))
  ipt_syn_CT_ct <-ipt_syn_dt[ipt_syn_dt$rice_CelltypeCate != 'broadly_accessible',]
  ipt_syn_CT_ACRnum <- length(unique(ipt_syn_CT_ct$rice_ACR))
  
  ##updating for check nonsyntenic and syntenic only for the CT
  CT_in_nonsyn_prop <- ipt_nonsyn_CT_ACRnum/ipt_nonsyn_ACRnum
  CT_in_syn_prop <- ipt_syn_CT_ACRnum/ipt_syn_ACRnum
  
  data_CT <- data.frame(
    cate <- c('CT_in_nonsyn_prop','CT_in_syn_prop'),
    prop <- c(CT_in_nonsyn_prop,CT_in_syn_prop)
  )
  colnames(data_CT) <- c('cate','prop')
  data_CT$spe <- spe2_nm
  
  ##print the sig test
  ##it is the same

  
  
  data <- data.frame(
    cate <- c('Broad_nonsyn_ACRprop','CT_nonsyn_ACRprop','Broad_syn_ACRprop','CT_syn_ACRprop'),
    prop <- c(ipt_nonsyn_broad_ACRnum/(ipt_nonsyn_broad_ACRnum+ipt_syn_broad_ACRnum),
            ipt_nonsyn_CT_ACRnum/(ipt_nonsyn_CT_ACRnum+ipt_syn_CT_ACRnum),
            ipt_syn_broad_ACRnum/(ipt_nonsyn_broad_ACRnum+ipt_syn_broad_ACRnum),
            ipt_syn_CT_ACRnum/(ipt_nonsyn_CT_ACRnum+ipt_syn_CT_ACRnum))
  )

  colnames(data) <- c('cate','prop')
  data$spe <- spe2_nm
  
  #print the sig test
  CT_nonsyn_ACRnum <- ipt_nonsyn_CT_ACRnum
  CT_syn_ACRnum <- ipt_syn_CT_ACRnum
  Broad_nonsyn_ACRnum <- ipt_nonsyn_broad_ACRnum
  Broad_syn_ACRnum <- ipt_syn_broad_ACRnum
  
  data_sig <- matrix(c(CT_nonsyn_ACRnum, Broad_nonsyn_ACRnum, CT_syn_ACRnum, Broad_syn_ACRnum), nrow = 2)
  result_for_nonsyn <- fisher.test(data_sig,alternative = 'greater')

  data$pval_nonsyn <- result_for_nonsyn$p.value
  
  Broad_syn_ACRnum <- ipt_syn_broad_ACRnum
  Broad_nonsyn_ACRnum <- ipt_nonsyn_broad_ACRnum
  CT_syn_ACRnum <- ipt_syn_CT_ACRnum
  CT_nonsyn_ACRnum <- ipt_nonsyn_CT_ACRnum
  
  data_sig <- matrix(c(Broad_syn_ACRnum, CT_syn_ACRnum, Broad_nonsyn_ACRnum, CT_nonsyn_ACRnum), nrow = 2)
  result_for_syn <- fisher.test(data_sig,alternative = 'greater')
  
  data$pval_syn <- result_for_syn$p.value
    
  data_CT$pval_syn <- result_for_syn$p.value
  ##Here we return the data_CT other than teh data
  return(data_CT)
})
  
combine_dt <- do.call(rbind,outs)

write.table(combine_dt,paste0(opt_dir,'/opt_syn_nonsyn_CTprop_ACR.txt'),quote = F,sep = '\t')

##we will plot the CT non syn prop version
#ipt_target_dt <- combine_dt[grepl('nonsyn',combine_dt$cate),]
combine_dt$spe <- factor(combine_dt$spe,levels = c('maize','sorghum','Pm','Uf'))
combine_dt$cate <- factor(combine_dt$cate, levels = c('CT_in_syn_prop','CT_in_nonsyn_prop'))


p <- ggplot(data=combine_dt, aes(x=spe, y=prop, fill=cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_y_continuous(breaks=seq(0, 0.2, 0.05))
#facet_wrap(~Rice_CellType,nrow = 1)
pdf(paste0(opt_dir,'/opt_check_non_syn_ct_prop.pdf'),width = 8,height = 6 )
p
dev.off()
  

#########
##we next check the cell type version
##updating 123123
##we will check the proprotion of each cell type under each category 
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123/////'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123///'
spe1 = 'rice'
cate_prefix = 'RiceToOthers'

#for others to rice
ipt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_others_to_rice_121023/'
opt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_others_to_rice_results_121023/'
spe1 = 'other'
cate_prefix = 'OthersToRice'

all_fl_list <- list.files(path= ipt_dir,pattern = '.txt')

outs <- lapply(all_fl_list, function(x){
  print(x)
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_H3K27me3.+','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  if (spe1 == 'rice'){
    spe2_nm <- gsub('rice_','',target_fl_path_nm)
  }else{
    spe2_nm <- 'rice'
  }

  ipt_dt <- read.delim(target_fl_path,header=T)
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[c(paste0(spe1,'_ACR'),paste0(spe2_nm,'_region'),paste0(spe1,'_CelltypeCate'))]
  dim(ipt_dt)
  ipt_dt <- unique(ipt_dt)
  ipt_dt_nonsyn <- unique(ipt_dt[ipt_dt[[paste0(spe2_nm,'_region')]] == 'none',][c(paste0(spe1,'_ACR'),paste0(spe1,'_CelltypeCate'))])
  dim(ipt_dt_nonsyn)
  ipt_dt_syn <- unique(ipt_dt[ipt_dt[[paste0(spe2_nm,'_region')]] != 'none',][c(paste0(spe1,'_ACR'),paste0(spe1,'_CelltypeCate'))])
  dim(ipt_dt_syn)
  
  ipt_dt_nonsyn$cate <- c('NonSyn')
  ipt_dt_syn$cate <- c('Syn')
  
  ipt_dt_combine <- rbind(ipt_dt_nonsyn,ipt_dt_syn)
  ipt_dt <- ipt_dt_combine
  
  colnames(ipt_dt) <- c('ACR','celltype','cate')
  
  ipt_dt <- as.data.frame(separate_rows(ipt_dt, celltype, sep = ","))
  
  ipt_dt$count <- 1  
  dim(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'broadly_accessible',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown_cells_1' & ipt_dt$celltype != 'unknown_cells_2',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown'&ipt_dt$celltype != 'procambial_meristem',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'procambium',]
  dim(ipt_dt)
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$cate != 'NotBlast',]
  
  aggregate_dt <- aggregate(count ~ cate + celltype, ipt_dt, sum)
  aggregate_celltype_dt <- aggregate(count ~ celltype, ipt_dt, sum)
  
  merged_dt <- merge(aggregate_dt,aggregate_celltype_dt, by.x = 'celltype',by.y = 'celltype')
  
  merged_dt$prop <- merged_dt$count.x/merged_dt$count.y
  
  merged_dt$spe <- target_fl_path_nm
  
  return(merged_dt)
  
})

combine_dt <- do.call(rbind,outs)

head(combine_dt)

combine_dt$cate <- factor(combine_dt$cate, levels = c('Syn','NonSyn'))

combine_dt <- combine_dt[combine_dt$cate != 'NotBlast',]

if (cate_prefix == 'RiceToOthers'){
  combine_dt$spe <- factor(combine_dt$spe, levels = c('rice_maize','rice_sorghum','rice_Pm','rice_Uf'))
  color_value <- c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                   "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1")
  combine_dt$celltype <- factor(combine_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cell',
                                                                'protoderm','epidermis'))
}
if (cate_prefix == 'OthersToRice'){
  combine_dt$spe <- factor(combine_dt$spe, levels = c('maize_rice','sorghum_rice','Pm_rice','Uf_rice'))
  color_value <- c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
                   "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1")
  combine_dt$celltype <- factor(combine_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements',
                                                                'protoderm','epidermis'))
}




p <- ggplot(data=combine_dt, aes(x=cate, y=prop, fill=celltype)) +
  scale_fill_manual(values=color_value) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_line(stat="identity", position=position_dodge())+
  #geom_point(stat="identity", position=position_dodge())+
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  facet_wrap(~spe,nrow = 1)
p

pdf(paste0(opt_dir,'/','opt_celltype_addTobe1_prop_for_syn_nonsyn_cate.pdf'),width = 14,height = 12 )
p
dev.off()



combine_dt_SS <- combine_dt[combine_dt$cate == 'NonSyn',]

#combine_dt_SS$celltype <- factor(combine_dt_SS$celltype,levels = rev(c('mesophyll','bundle_sheath','companion_cell',
#                                                                   'protoderm','epidermis')))

p <- ggplot(data=combine_dt_SS, aes(x=celltype, y=prop, fill=celltype)) +
  scale_fill_manual(values=color_value) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_line(stat="identity", position=position_dodge())+
  #geom_point(stat="identity", position=position_dodge())+
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  facet_wrap(~ spe, ncol= 1)+
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0,1))

p
pdf(paste0(opt_dir,'/','opt2_celltype_addTobe1_prop_for_threecate_only_NonSyn.pdf'),width = 9,height = 14 )
p
dev.off()













#########
##we next check the enrichment of cell type to other cell type using the fisher
ipt_dir <- 'opt1_acr_syntenic_summary_dir_add_nonsyntenic_123123/////'
opt_dir <- 'opt1_acr_syntenic_summary_dir_add_nonsyntenic_results_123123///'
spe1 = 'rice'
cate_prefix == 'RiceToOthers'

#for others to rice
ipt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_others_to_rice_121023/'
opt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_others_to_rice_results_121023/'
spe1 = 'other'
cate_prefix == 'OthersToRice'

all_fl_list <- list.files(path= ipt_dir)

outs <- lapply(all_fl_list, function(x){
  print(x)
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_H3K27me3.+','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  if (spe1 == 'rice'){
    spe2_nm <- gsub('rice_','',target_fl_path_nm)
  }else{
    spe2_nm <- 'rice'
  }
  
  ipt_dt <- read.delim(target_fl_path,header=T)
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[c(paste0(spe1,'_ACR'),paste0(spe2_nm,'_region'),paste0(spe1,'_CelltypeCate'))]
  dim(ipt_dt)
  ipt_dt <- unique(ipt_dt)
  ipt_dt_nonsyn <- unique(ipt_dt[ipt_dt[[paste0(spe2_nm,'_region')]] == 'none',][c(paste0(spe1,'_ACR'),paste0(spe1,'_CelltypeCate'))])
  dim(ipt_dt_nonsyn)
  ipt_dt_syn <- unique(ipt_dt[ipt_dt[[paste0(spe2_nm,'_region')]] != 'none',][c(paste0(spe1,'_ACR'),paste0(spe1,'_CelltypeCate'))])
  dim(ipt_dt_syn)
  
  ipt_dt_nonsyn$cate <- c('NonSyn')
  ipt_dt_syn$cate <- c('Syn')
  
  ipt_dt_combine <- rbind(ipt_dt_nonsyn,ipt_dt_syn)
  ipt_dt <- ipt_dt_combine
  
  colnames(ipt_dt) <- c('ACR','celltype','cate')
  
  ipt_dt <- as.data.frame(separate_rows(ipt_dt, celltype, sep = ","))
  
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'broadly_accessible',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown_cells_1' & ipt_dt$celltype != 'unknown_cells_2',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown'&ipt_dt$celltype != 'procambial_meristem',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'procambium',]
  
  ipt_dt$count <- 1
  
  all_celltype_list <- unique(ipt_dt$celltype)
  
  outs2 <- lapply(all_celltype_list, function(y){
    
    outs3 <- lapply(all_celltype_list, function(z){
      
      other_ct <- z
      
      target_ct_in_nonsyn_acrnum = nrow(unique(ipt_dt[ipt_dt$celltype == y & ipt_dt$cate == 'NonSyn',]))
      target_ct_in_syn_acrnum = nrow(unique(ipt_dt[ipt_dt$celltype == y & ipt_dt$cate == 'Syn',]))
      
      other_ct_in_nonsyn_acrnum = nrow(unique(ipt_dt[ipt_dt$celltype == other_ct & ipt_dt$cate == 'NonSyn',]))
      other_ct_in_syn_acrnum = nrow(unique(ipt_dt[ipt_dt$celltype == other_ct & ipt_dt$cate == 'Syn',]))
      
      data_sig <- matrix(c(target_ct_in_nonsyn_acrnum, other_ct_in_nonsyn_acrnum, target_ct_in_syn_acrnum, other_ct_in_syn_acrnum), nrow = 2)
      result <- fisher.test(data_sig,alternative = 'greater')
      
      res <- c()
      res$targetCT <- y
      res$otherCT <- other_ct
      res$pval <- result$p.value
      res <- as.data.frame(res)
      return(res)
      
    })
    
    combine_outs3_dt <- do.call(rbind,outs3)
      
  })
  
  combine_outs2_dt <- do.call(rbind,outs2)
  
  combine_outs2_dt$spe <- spe2_nm
  
  return(combine_outs2_dt)

})
  
combine_dt <- do.call(rbind,outs)  

write.table(combine_dt,paste0(opt_dir,'/opt_syn_nonsyn_enrichment_target_celltype_to_others.txt'),quote = F,sep = '\t')


###################
##make the heat map
spe_list <- unique(combine_dt$spe)

spe1 = 'rice'

outs <- lapply(spe_list, function(x){
  
  
  
  if (spe1 == 'rice'){
    spe2_nm <- x
    spe1_nm <- 'rice'
  }else{
    spe2_nm <- 'rice'
    spe1_nm <- x
  }
  
  opt_name <- paste0(spe1_nm,'_to_',spe2_nm)
  
  print(x)
  ipt_dt <- combine_dt[combine_dt$spe == x,]
  
  ipt_dt <- ipt_dt[c('targetCT','otherCT','pval')]

  #ipt_dt <- ipt_dt[ipt_dt$pval < 0.05,]
  ipt_dt$mlog10pval <- -log10(ipt_dt$pval)    
  ipt_dt <- ipt_dt[c('targetCT','otherCT','mlog10pval')]
  colnames(ipt_dt) <- c('target_celltype','other_celltype','mlog10pval')
  
  
  write.table(ipt_dt,paste0(opt_dir,'/temp_',opt_name,'.txt'),quote = F,sep = '\t')
  
  ipt_dt <- read.table(paste0(opt_dir,'/temp_',opt_name,'.txt'),stringsAsFactors = T)
  
  ipt_dt$target_celltype <- gsub('companion_cells_sieve_elements','companion_cell',ipt_dt$target_celltype)
  ipt_dt$other_celltype <- gsub('companion_cells_sieve_elements','companion_cell',ipt_dt$other_celltype)
  
  
  ipt_dt <- ipt_dt[ipt_dt$target_celltype != 'procambial_meristem',]
  ipt_dt <- ipt_dt[ipt_dt$other_celltype != 'procambial_meristem',]
  
  ipt_dt <- ipt_dt[ipt_dt$target_celltype != 'procambium',]
  ipt_dt <- ipt_dt[ipt_dt$other_celltype != 'procambium',]
  
  ipt_dt <- ipt_dt[ipt_dt$target_celltype != 'unknown',]
  ipt_dt <- ipt_dt[ipt_dt$other_celltype != 'unknown',]
  
  
  ipt_dt$target_celltype <- factor(ipt_dt$target_celltype,levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
  ipt_dt$other_celltype <- factor(ipt_dt$other_celltype,levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
  
  score_mtx <- sparseMatrix(i=as.numeric(ipt_dt$target_celltype),
                            j=as.numeric(ipt_dt$other_celltype),
                            x=as.numeric(ipt_dt$mlog10pval),
                            dimnames=list(levels(ipt_dt$target_celltype), levels(ipt_dt$other_celltype)))
  score_mtx <- as.matrix(score_mtx)
  
  pdf(paste0(opt_dir,'/opt_enrichment_celltype_on_nonsyntenic_species_specific_',opt_name,'.pdf'),width = 5,height = 3)
  p <- pheatmap(score_mtx,
                #scale="row",
                #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
                #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
                cluster_cols = F,
                cluster_rows = F,
                border_color='white',
                fontsize_col = 15,
                fontsize_row = 15,
                fontsize_number = 25
  ) 
  print(p)
  dev.off()
  
  
  
  
})



#################
##updating 010323
##we will check the gene cate composition between syn and nonsyn
ipt_dir <- 'opt1_acr_syntenic_summary_dir_add_nonsyntenic_123123/////'
opt_dir <- 'opt1_acr_syntenic_summary_dir_add_nonsyntenic_results_123123///opt_check_proximity_spe_prop_dir'
spe1 = 'rice'

ipt_acr_gene_cate_dt <- read.delim(paste0(ipt_dir,'/gene_cate_dir/opt_',spe1,'_acr_toGeneCate.txt'),header = F)
head(ipt_acr_gene_cate_dt)
ipt_acr_gene_cate_dt <- ipt_acr_gene_cate_dt[c('V1','V2')]
colnames(ipt_acr_gene_cate_dt) <- c('acrloc','genecate')

all_fl_list <- list.files(path= ipt_dir,pattern = '.txt')

outs <- lapply(all_fl_list, function(x){
  
  print(x)
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_H3K27me3.+','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  if (spe1 == 'rice'){
    spe2_nm <- gsub('rice_','',target_fl_path_nm)
  }else{
    spe2_nm <- 'rice'
  }

  ipt_dt <- read.delim(target_fl_path,header=T)
  head(ipt_dt)
  

  
  ipt_target_dt <- ipt_dt[c(paste0(spe1,'_ACRID'),paste0(spe1,'_ACR'),paste0(spe1,'_CelltypeCate'),'Syntenic_regionID_update')]
  
  ipt_target_dt <- merge(ipt_target_dt,ipt_acr_gene_cate_dt,by.x = paste0(spe1,'_ACR'),by.y = 'V1')
  
  ipt_target_dt$CTcate <- ifelse(ipt_target_dt[[paste0(spe1,'_CelltypeCate')]] == 'broadly_accessible','Broad','CT') 
  ipt_target_dt$SynOrNoSyn <- ifelse(ipt_target_dt$Syntenic_regionID_update == 'none','NSyn','Syn')
  head(ipt_target_dt)
  
  ipt_target_dt <- ipt_target_dt[c(paste0(spe1,'_ACR'),'CTcate','SynOrNoSyn','V2')]
  ipt_target_dt <- unique(ipt_target_dt)
  dim(ipt_target_dt)
  colnames(ipt_target_dt) <- c(paste0(spe1,'_ACR'),'CTcate','SynOrNoSyn','GeneCate')
  
  ipt_target_dt$count <- 1
  aggregate_dt <- aggregate(count ~ SynOrNoSyn + CTcate + GeneCate, data = ipt_target_dt,sum)  
  aggregate_broadCT_GeneCate_dt <- aggregate(count ~ CTcate + GeneCate, data = ipt_target_dt, sum)
  aggregate_broadCT_GeneCate_dt$CT_gene_Cate <- paste0(aggregate_broadCT_GeneCate_dt$CTcate,'_',aggregate_broadCT_GeneCate_dt$GeneCate)
  
  aggregate_dt$CT_gene_Cate <- paste0(aggregate_dt$CTcate,'_',aggregate_dt$GeneCate)
  
  mergedt_aggregate_broadCT_GeneCate_dt <- merge(aggregate_dt,aggregate_broadCT_GeneCate_dt,by.x = 'CT_gene_Cate',by.y = 'CT_gene_Cate')
  mergedt_aggregate_broadCT_GeneCate_dt$prop <- mergedt_aggregate_broadCT_GeneCate_dt$count.x/mergedt_aggregate_broadCT_GeneCate_dt$count.y
  
  ##sum the distal genic and proximal together per CTcate per SynOrNoSyn
  aggregate_CTcate_SynOrNoSyn <- aggregate(count ~ CTcate + SynOrNoSyn, data = ipt_target_dt, sum)
  aggregate_CTcate_SynOrNoSyn$CTcate_synnonsyn <- paste0(aggregate_CTcate_SynOrNoSyn$CTcate,'_',aggregate_CTcate_SynOrNoSyn$SynOrNoSyn)
  
  aggregate_dt$CTcate_synnonsyn <- paste0(aggregate_dt$CTcate,'_',aggregate_dt$SynOrNoSyn)
  
  mergedt_aggregate_CTcate_SynOrNoSyn <- merge(aggregate_dt,aggregate_CTcate_SynOrNoSyn,by.x = 'CTcate_synnonsyn',by.y = 'CTcate_synnonsyn')
  mergedt_aggregate_CTcate_SynOrNoSyn$prop <- mergedt_aggregate_CTcate_SynOrNoSyn$count.x/mergedt_aggregate_CTcate_SynOrNoSyn$count.y
  
  ##build the bar plot to show the proportion for the mergedt_aggregate_CTcate_SynOrNoSyn
  mergedt_aggregate_CTcate_SynOrNoSyn$spe <- spe2_nm
  
  return(mergedt_aggregate_CTcate_SynOrNoSyn)
})

final_combine_dt <- do.call(rbind,outs)
combine_dt <- final_combine_dt

spe_list <- unique(combine_dt$spe)

for (i in 1:length(spe_list)){
  
  combine_cate_dt <- combine_dt[combine_dt$spe == spe_list[i],]
  
  combine_cate_dt$GeneCate <- factor(combine_cate_dt$GeneCate,levels = c('Genic','Proximal','Distal'))
  
  combine_cate_dt$SynOrNoSyn.x <- factor(combine_cate_dt$SynOrNoSyn.x,levels = c('Syn','NSyn'))
  
  p <- ggplot(data=combine_cate_dt, aes(x=SynOrNoSyn.x , y=prop, fill=CTcate.x)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    facet_wrap(~GeneCate,nrow =  1) ## You can use face_wrap function only if you need it+
  
  pdf(paste0(opt_dir,'/','opt_',spe_list[i],'_gene_proximity_compare_prop','.pdf'),width = 12,height = 5 )
  print(p)
  dev.off()
  
}


####we will check the gene cate composition between syn and nonsyn
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123/////H3K27me3_summary_dir/'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123///opt_check_H3K27me3_prop_dir'
spe1 = 'rice'

#ipt_acr_gene_cate_dt <- read.delim(paste0(ipt_dir,'/gene_cate_dir/opt_',spe1,'_acr_toGeneCate.txt'),header = F)
#head(ipt_acr_gene_cate_dt)
#ipt_acr_gene_cate_dt <- ipt_acr_gene_cate_dt[c('V1','V2')]
#colnames(ipt_acr_gene_cate_dt) <- c('acrloc','genecate')

all_fl_list <- list.files(path= ipt_dir,pattern = '.txt')

outs <- lapply(all_fl_list, function(x){
  
  print(x)
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_summary_file_add_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  if (spe1 == 'rice'){
    spe2_nm <- gsub('rice_','',target_fl_path_nm)
  }else{
    spe2_nm <- 'rice'
  }
  
  ipt_dt <- read.delim(target_fl_path,header=T)
  head(ipt_dt)
  
  ipt_target_dt <- ipt_dt[c(paste0(spe1,'_ACR'),paste0(spe1,'_CelltypeCate'),paste0(spe1,'OverH3K27Cate'),paste0('Syntenic_regionID_update'))]
  
  ipt_target_dt$CTcate <- ifelse(ipt_target_dt[[paste0(spe1,'_CelltypeCate')]] == 'broadly_accessible','Broad','CT') 
  ipt_target_dt$SynOrNoSyn <- ifelse(ipt_target_dt$Syntenic_regionID_update == 'none','NSyn','Syn')
  head(ipt_target_dt)
  dim(ipt_target_dt)
  ipt_target_dt <- unique(ipt_target_dt[c(paste0(spe1,'_ACR'),paste0(spe1,'OverH3K27Cate'),'CTcate','SynOrNoSyn')])
  head(ipt_target_dt)
  dim(ipt_target_dt)
  
  ipt_target_dt$count <- 1
  
  aggregate_dt <- aggregate(count ~ riceOverH3K27Cate + CTcate + SynOrNoSyn, data = ipt_target_dt, sum)
  aggregate_dt$CTcate_synnonsyn <- paste0(aggregate_dt$CTcate,'_',aggregate_dt$SynOrNoSyn)
  
  aggregate_CTcate_SynOrNoSyn_dt <- aggregate(count ~ CTcate + SynOrNoSyn, data = aggregate_dt, sum)
  aggregate_CTcate_SynOrNoSyn_dt$CTcate_synnonsyn <- paste0(aggregate_CTcate_SynOrNoSyn_dt$CTcate,'_',aggregate_CTcate_SynOrNoSyn_dt$SynOrNoSyn)
  
  mergedt_CTcate_SynOrNoSyn_dt <- merge(aggregate_dt,aggregate_CTcate_SynOrNoSyn_dt,by.x= 'CTcate_synnonsyn', by.y = 'CTcate_synnonsyn')
  mergedt_CTcate_SynOrNoSyn_dt$prop <- mergedt_CTcate_SynOrNoSyn_dt$count.x/mergedt_CTcate_SynOrNoSyn_dt$count.y
  
  mergedt_CTcate_SynOrNoSyn_H3K27_dt <- mergedt_CTcate_SynOrNoSyn_dt[mergedt_CTcate_SynOrNoSyn_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  mergedt_CTcate_SynOrNoSyn_notH3K27_dt <- mergedt_CTcate_SynOrNoSyn_dt[mergedt_CTcate_SynOrNoSyn_dt$riceOverH3K27Cate != 'OverlapH3K27me3',]
  
  mergedt_CTcate_SynOrNoSyn_H3K27_dt$SynOrNoSyn.y <- factor(mergedt_CTcate_SynOrNoSyn_H3K27_dt$SynOrNoSyn.y,levels = c('Syn','NSyn'))
  
  p <- ggplot(data=mergedt_CTcate_SynOrNoSyn_H3K27_dt, aes(x=CTcate.x , y=prop, fill=SynOrNoSyn.x)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) 
  p
  
  pdf(paste0(opt_dir,'/','opt_',spe2_nm,'_H3K27me3_overlap_ACR','.pdf'),width = 6,height = 5 )
  print(p)
  dev.off()
  
  ##conduct a test to show if the CT or Broad enrich on the nonsyn
  ##broad acr in nsyn number
  ##broad acr in syn number
  ##ct acr in nonsyn number
  ##ct acr in syn number
  broad_acr_inH3K27me3_in_nsy_num <- mergedt_CTcate_SynOrNoSyn_H3K27_dt[mergedt_CTcate_SynOrNoSyn_H3K27_dt$CTcate_synnonsyn == 'Broad_NSyn',]$count.x
  broad_acr_inH3K27me3_in_syn_num <- mergedt_CTcate_SynOrNoSyn_H3K27_dt[mergedt_CTcate_SynOrNoSyn_H3K27_dt$CTcate_synnonsyn == 'Broad_Syn',]$count.x
  
  broad_acr_notInH3K27me3_in_nsyn_num <- mergedt_CTcate_SynOrNoSyn_notH3K27_dt[mergedt_CTcate_SynOrNoSyn_notH3K27_dt$CTcate_synnonsyn == 'Broad_NSyn',]$count.x
  broad_acr_notInH3K27me3_in_syn_num <- mergedt_CTcate_SynOrNoSyn_notH3K27_dt[mergedt_CTcate_SynOrNoSyn_notH3K27_dt$CTcate_synnonsyn == 'Broad_Syn',]$count.x
  
  data_broad <- matrix(c(broad_acr_inH3K27me3_in_nsy_num, broad_acr_notInH3K27me3_in_nsyn_num, broad_acr_inH3K27me3_in_syn_num, broad_acr_notInH3K27me3_in_syn_num), nrow = 2)
  result_broad <- fisher.test(data_broad,alternative = 'greater')
  
  CT_acr_in_inH3K27me3_nsyn_num <- mergedt_CTcate_SynOrNoSyn_H3K27_dt[mergedt_CTcate_SynOrNoSyn_H3K27_dt$CTcate_synnonsyn == 'CT_NSyn',]$count.x
  CT_acr_in_inH3K27me3_syn_num <- mergedt_CTcate_SynOrNoSyn_H3K27_dt[mergedt_CTcate_SynOrNoSyn_H3K27_dt$CTcate_synnonsyn == 'CT_Syn',]$count.x
  
  CT_acr_in_notInH3K27me3_nsyn_num <- mergedt_CTcate_SynOrNoSyn_notH3K27_dt[mergedt_CTcate_SynOrNoSyn_notH3K27_dt$CTcate_synnonsyn == 'CT_NSyn',]$count.x
  CT_acr_in_notInH3K27me3_syn_num <- mergedt_CTcate_SynOrNoSyn_notH3K27_dt[mergedt_CTcate_SynOrNoSyn_notH3K27_dt$CTcate_synnonsyn == 'CT_Syn',]$count.x
  
  data_CT <- matrix(c(broad_acr_inH3K27me3_in_nsy_num, broad_acr_notInH3K27me3_in_nsyn_num, broad_acr_inH3K27me3_in_syn_num, broad_acr_notInH3K27me3_in_syn_num), nrow = 2)
  result_CT <- fisher.test(data_CT,alternative = 'greater')
    
  res <- c()
  res$broadpval <- result_broad$p.value
  res$ctpval <- result_CT$p.value
  res$spe <- spe2_nm
  res <- as.data.frame(res)
  
  
  p <- ggplot(mergedt_CTcate_SynOrNoSyn_dt, aes(fill=riceOverH3K27Cate, y=prop, x=CTcate.x)) + 
    geom_bar(position="stack", stat="identity")
  #p <- p + geom_text(aes(label = Proportion),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
  p <- p+labs(x="\n", y = "Percentage\n", cex.lab = 7) ## use \n to set the position
  #p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
  p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size = 30), ##change all the text size
               axis.text.y = element_text(size=40,colour = "black"),
               axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
               axis.title.x = element_text(size=40,colour = "black"),
               axis.title.y = element_text(size=40,colour = "black"),
               plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
  
  p <- p + coord_cartesian(ylim = c(0, 1)) + ##threshold
    scale_y_continuous(breaks=seq(0, 1, 0.1))
  p <- p + facet_wrap(~ SynOrNoSyn.x)

  write.table(mergedt_CTcate_SynOrNoSyn_dt,paste0(opt_dir,'/','opt_',spe2_nm,'_H3K27me3_overlap_ACR_number','.txt'),quote = F,sep = '\t')
  
  pdf(paste0(opt_dir,'/','opt_',spe2_nm,'_H3K27me3_overlap_ACR_stack_Version','.pdf'),width = 12,height = 8 )
  print(p)
  dev.off()
  
  return(res)
  
})

combine_pval_dt <- do.call(rbind,outs)
combine_pval_dt
##     broadpval       ctpval     spe
#1 1.664156e-34 1.664156e-34   maize
#2 4.942617e-52 4.942617e-52 sorghum


##updating 010423
##check the enrichment of motifs 
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123/////collect_binomial_enrich_dir//'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123///opt_check_motif_enrich_dir/'
spe1 = 'rice'

target_cate <- 'nonsyn'
target_cate <- 'syn'


all_fl_list <- list.files(path= ipt_dir,pattern = '.txt')


##for the d
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_binomial_motif_enrich.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  ipt_dt$spe <- target_fl_path_nm 
  
  colnames(ipt_dt) <- c('cate','celltype','motifID','pval','num','total','prop','spe') 
  
  ipt_dt <- ipt_dt[c('cate','celltype','motifID','pval','spe')]
  return(ipt_dt)  
  
})


combine_dt <- do.call(rbind,outs)
combine_dt <- combine_dt[combine_dt$celltype != 'unknown_cells_1'&combine_dt$celltype != 'unknown_cells_2', ]

meta_name_dt <- read.delim(paste0(ipt_dir,'/ipt_required_raw_data/opt_motif_common_name.txt'),header = F)
head(meta_name_dt)

##add the common name
combine_dt <- merge(combine_dt,meta_name_dt,by.x='motifID',by.y= 'V1')
head(combine_dt)
combine_dt$qval <- p.adjust(combine_dt$pval, method="fdr")
combine_dt <- combine_dt[combine_dt$qval < 0.05,]

write.table(combine_dt,paste0(opt_dir,'/opt_enrichment_syn_nonsyn_dt.txt'),sep = '\t',quote = F)


##non syntenic regions 
combine_notshare_dt <- combine_dt[combine_dt$cate == target_cate,]
combine_notshare_dt$mlog10qval <- -log10(combine_notshare_dt$qval)
combine_notshare_dt$spe_pair_ct <- paste0(combine_notshare_dt$spe,'__',combine_notshare_dt$celltype)
combine_notshare_dt <- combine_notshare_dt[c('V2','spe_pair_ct','mlog10qval')]
colnames(combine_notshare_dt) <- c('motif','spe_pair_ct','mlog10qval')



write.table(combine_notshare_dt,paste0(opt_dir,'/temp_sparse.txt'),quote = F,sep = '\t')
combine_notshare_dt <- read.table(paste0(opt_dir,'/temp_sparse.txt'),stringsAsFactors = T)
score_mtx <- sparseMatrix(i=as.numeric(combine_notshare_dt$motif),
                          j=as.numeric(combine_notshare_dt$spe_pair),
                          x=as.numeric(combine_notshare_dt$mlog10qval),
                          dimnames=list(levels(combine_notshare_dt$motif), levels(combine_notshare_dt$spe_pair)))
score_mtx <- as.matrix(score_mtx)
dim(score_mtx)
head(score_mtx)
max(score_mtx)

##reorder the matrix
colnames(score_mtx)
target_order <- c('rice_maize__mesophyll','rice_sorghum__mesophyll','rice_Pm__mesophyll','rice_Uf__mesophyll',
                  'rice_maize__bundle_sheath','rice_sorghum__bundle_sheath','rice_Pm__bundle_sheath','rice_Uf__bundle_sheath',
                  'rice_maize__companion_cell','rice_sorghum__companion_cell','rice_Pm__companion_cell','rice_Uf__companion_cell',
                  'rice_maize__protoderm','rice_sorghum__protoderm','rice_Pm__protoderm','rice_Uf__protoderm',
                  'rice_maize__epidermis','rice_sorghum__epidermis','rice_Pm__epidermis','rice_Uf__epidermis'
)

intersect_colnm <- intersect(target_order,colnames(score_mtx))

score_mtx <- score_mtx[,intersect_colnm]

# Replace Inf with the maximum value of the matrix
max_non_inf <- max(score_mtx[!is.infinite(score_mtx)], na.rm = TRUE)
score_mtx[is.infinite(score_mtx)] <- max_non_inf
max(score_mtx)

p <- pheatmap(score_mtx,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
pdf(paste0(opt_dir,'/opt_',target_cate,'_celltype_enrich_motif.pdf'),width = 15,height = 15)
p
dev.off()


##we will plot a dot plot to show the family with avg pval and number of TF motif
##Here we will check family beta 
motif_name_dt <- read.delim(paste0(ipt_dir,'/ipt_required_raw_data/ipt_target_motifs_list_all.txt'),header = F)
head(motif_name_dt)
table(motif_name_dt$V3)

merged_dt <- merge(combine_notshare_dt,motif_name_dt, by.x = 'motif',by.y = 'V3')
head(merged_dt)
colnames(merged_dt) <- c('motifnm','spe_pair_ct','mlog10qval','motifID','motifam','select')



final_dt <- merged_dt[c('motifam','motifnm','spe_pair_ct','mlog10qval')]
#final_dt$spe <- target_fl_path_nm  

##we will modify different families of motif
table(final_dt$motifam)

final_dt$motifam <- gsub('Group D','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group A','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group B','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group C','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group G','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group H','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group I','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group K','bZIP',final_dt$motifam)
final_dt$motifam <- gsub('Group S','bZIP',final_dt$motifam)

final_dt$motifam <- gsub('Myb-related','Myb',final_dt$motifam)
final_dt$motifam <- gsub('WRKY-like_FRS/FRF','WRKY',final_dt$motifam)

final_dt <- final_dt[final_dt$motifam != 'unknown',]
final_dt <- final_dt[final_dt$motifam != 'Unknown',]

table(final_dt$motifam)

final_dt$spe <- gsub('__.+','',final_dt$spe_pair_ct)

all_spe <- unique(final_dt$spe)

table(final_dt$motifam)

final_dt <- final_dt[final_dt$motifam != 'Click to view details',]

outs <- lapply(all_spe, function(x){
  
  ipt_target_dt <- final_dt[final_dt$spe == x,]
  ipt_target_dt$celltype <- gsub('.+__','',ipt_target_dt$spe_pair_ct)
  
  celltype_list <- unique(ipt_target_dt$celltype)
  
  outs2 <- lapply(celltype_list, function(y){
    
    ipt_target_ct_dt <- unique(ipt_target_dt[ipt_target_dt$celltype == y,])
    
    ipt_target_ct_dt$count <- 1
    
    aggregate_avg_fam_dt <- aggregate(mlog10qval ~ motifam, data = ipt_target_ct_dt, mean)
    aggregate_count_fam_dt <- aggregate(count ~ motifam, data = ipt_target_ct_dt, sum)
    
    aggregate_merge <- merge(aggregate_avg_fam_dt,aggregate_count_fam_dt,by.x = 'motifam',by.y= 'motifam')
    
    rownames(aggregate_merge) <- aggregate_merge$motifam
    aggregate_merge$log10count <- log10(aggregate_merge$count)
    aggregate_merge <- aggregate_merge[c('count','mlog10qval')]
    aggregate_merge
    
    points_to_label <- rownames(aggregate_merge)
    #make the dot plot
    aggregate_merge <- as.data.frame(aggregate_merge)
    p <- ggplot(aggregate_merge, aes(count, mlog10qval)) +
      geom_point() +
      geom_point(data = aggregate_merge, color = 'black', size = 5) +
      geom_text(data = aggregate_merge, aes(label = rownames(aggregate_merge)), hjust = -0.2,color='red',cex= 7) + 
      theme(
        plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
        axis.title = element_text(size =20),
        #axis.text.x = element_text(colour = "black", size=30,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
        axis.text.x = element_text(size=20,colour = "black"),
        axis.text.y = element_text(size=20,colour = "black"),
        
        axis.ticks = element_line(size = rel(2.5)),
        axis.ticks.length = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
      )
    
    p <- p +
      scale_x_continuous(breaks=seq(0, 30, 5)) +
      #scale_y_continuous(breaks=seq(-2, 2, 0.5)) +
      #coord_cartesian(xlim = c(-2, 2), ylim = c(-2,2))
      coord_cartesian(xlim = c(0, 30))  
    
    target_name <- paste0(x,'_',y,'_',target_cate)
    pdf(paste0(opt_dir,'/opt_dotplot_',target_name,'_celltype_enrich_motif.pdf'),width = 8,height = 8)
    print(p)
    dev.off()
  
    
  })
  
})



#################
##updating 011124
##we will check the TEs and their enrichment
##these are for the rice to others
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123/TE_rice_to_others_dir/collect_all_TEs_cate_dir/'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123/opt_check_rice_to_others_TE_dir/'

##the followings are the others to rice
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123/TE_others_to_rice_dir/collect_all_TE_cate_prop_dir/'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123/opt_check_others_to_rice_TE_dir//'


all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_final_TE_prop_in_eachcate')

##for the d
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_TE_prop_in_eachcate_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)

  
  ipt_dt <- read.delim(target_fl_path, header = F)
  table(ipt_dt$V6)
  ipt_dt$V6 <- gsub('nLTRothers.+','nLTRothers',ipt_dt$V6)
  ipt_dt$V6 <- gsub('DNAothers.+','unclassifiedDNA',ipt_dt$V6)
  ipt_dt <- ipt_dt[ipt_dt$V6 != 'Others_centromeric_repeat_Centro.tandem'& ipt_dt$V6 != 'SateliteDNA' & ipt_dt$V6 != 'Others_low_complexity_Simple_repeat'
                   &ipt_dt$V6 != 'Others_knob_knob.knob180' & ipt_dt$V6 != 'Others_centromeric_repeat_Cent.CentC',]
  ipt_dt$V6 <- gsub('Others_repeat_region_Unknown','Unknown', ipt_dt$V6)
  ipt_dt$V6 <- gsub('DNAnona','DNAnona_others',ipt_dt$V6)
  ipt_dt$V6 <- gsub('MITE','DNAnona_MITE',ipt_dt$V6)
  
  ipt_target_dt <- ipt_dt[c('V1','V4','V6')]
  colnames(ipt_target_dt) <- c('Cate','Num','TEfam')
  
  aggregate_eachcate <- aggregate(Num ~ Cate, data = ipt_target_dt, sum)
  
  ipt_aggregate_target_dt <- aggregate(Num ~ TEfam + Cate , data= ipt_target_dt ,sum)

  mergedt_dt <- merge(ipt_aggregate_target_dt,aggregate_eachcate, by.x = 'Cate', by.y = 'Cate')
  mergedt_dt$prop <- mergedt_dt$Num.x/mergedt_dt$Num.y  
  
  ##compare number
  ##ggplot
  
  mergedt_dt$TEfam <- factor(mergedt_dt$TEfam,levels = c('LTR','LINE','SINE','nLTRothers','DNAauto','DNAnona_MITE','DNAnona_others','Helitron','unclassifiedDNA','Unknown'))

  
  mergedt_dt <- mergedt_dt[mergedt_dt$Cate != 'nonsyn_all' & mergedt_dt$Cate != 'syn_all',]
  mergedt_dt$Cate <- factor(mergedt_dt$Cate,levels = c('syn_broad','nonsyn_broad','syn_CT','nonsyn_CT'))
  
  
  p <- ggplot(mergedt_dt, aes(fill=TEfam, y=prop, x=Cate)) + 
    geom_bar(position="stack", stat="identity")
  #p <- p + geom_text(aes(label = Proportion),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
  p <- p+labs(x="\n", y = "\n", cex.lab = 7) ## use \n to set the position
  #p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
  p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size = 30), ##change all the text size
               axis.text.y = element_text(size=40,colour = "black"),
               axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
               axis.title.x = element_text(size=40,colour = "black"),
               axis.title.y = element_text(size=40,colour = "black"),
               plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
  
  p <- p + coord_cartesian(ylim = c(0, 1)) + ##threshold
    scale_y_continuous(breaks=seq(0, 1, 0.1))
  
  pdf(paste0(opt_dir,'/opt_TEfam_cate_', target_fl_path_nm,'.pdf'),width = 10, height = 10)
  print(p) ##20 15
  dev.off()
  
  ##updating 020124
  ##only plot the LTR for the comparison
  mergedt_LTR_dt <- mergedt_dt[mergedt_dt$TEfam == 'LTR',]
  
  p <- ggplot(mergedt_LTR_dt, aes(fill=TEfam, y=prop, x=Cate)) + 
    geom_bar(position="stack", stat="identity")
  #p <- p + geom_text(aes(label = Proportion),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
  p <- p+labs(x="\n", y = "\n", cex.lab = 7) ## use \n to set the position
  #p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
  p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size = 30), ##change all the text size
               axis.text.y = element_text(size=40,colour = "black"),
               axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
               axis.title.x = element_text(size=40,colour = "black"),
               axis.title.y = element_text(size=40,colour = "black"),
               plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
  
  p <- p + coord_cartesian(ylim = c(0, 0.15)) + ##threshold
    scale_y_continuous(breaks=seq(0, 0.15, 0.05))
  
  pdf(paste0(opt_dir,'/opt_TEfam_cate_', target_fl_path_nm,'_LTR.pdf'),width = 10, height = 10)
  print(p) ##20 15
  dev.off()
  
  
  
  ##updating 020124
  ##only consider the LTR
  ##we will conduct the enrichment test for the LTR
  ##LTR in nonsyn_broad
  LTR_dt <- ipt_aggregate_target_dt[ipt_aggregate_target_dt$TEfam == 'LTR',]
  LTR_nonsyn_broad <- LTR_dt[LTR_dt$Cate == 'nonsyn_broad',]$Num
  LTR_syn_broad <- LTR_dt[LTR_dt$Cate == 'syn_broad',]$Num
  
  not_LTR_dt <- ipt_aggregate_target_dt[ipt_aggregate_target_dt$TEfam != 'LTR',]
  aggregate_not_LTR_dt <- aggregate(Num ~ Cate, data = not_LTR_dt,sum)
  not_LTR_nonsyn_broad <- aggregate_not_LTR_dt[aggregate_not_LTR_dt$Cate == 'nonsyn_broad',]$Num
  
  not_LTR_syn_broad <- aggregate_not_LTR_dt[aggregate_not_LTR_dt$Cate == 'syn_broad', ]$Num
  
  ##do the enrichment test
  data <- matrix(c(LTR_nonsyn_broad, not_LTR_nonsyn_broad, LTR_syn_broad, not_LTR_syn_broad), nrow = 2)
  result <- fisher.test(data,alternative = 'two.sided')
  
  
  res <- c()
  res$spe <- target_fl_path_nm
  res$pval <- result$p.value
  res <- as.data.frame(res)
  return(res)
  
})

combine_dt <- do.call(rbind,outs)



########################
##conduct the enrichment test for the syn and nonsyn broad and CT
##but it looks like does not work
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123//TE_rice_to_others_dir/collect_broad_enrich_dir/'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123/opt_check_rice_to_others_TE_dir/'
target_cate <- 'broad'
spe1 = 'rice'

ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123//TE_rice_to_others_dir/collect_broad_enrich_dir/'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123/opt_check_rice_to_others_TE_dir/'
target_cate <- 'CT'
spe1 = 'rice'

all_fl_list <- list.files(path= ipt_dir,pattern = '*')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_nonsyn_vs_syn_TE_enrich_fam_broad_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  ipt_dt
  
  ipt_dt <- ipt_dt[c('V1','V2')]
  colnames(ipt_dt) <- c('TEfam','pval')
  
  ipt_dt$qval <- p.adjust(ipt_dt$pval, method="fdr")
  ipt_dt <- ipt_dt[ipt_dt$qval < 0.05,]
  
  ipt_dt$spe <- target_fl_path_nm
  ipt_dt$cate <- target_cate
  
  return(ipt_dt)
  
})

combine_dt_broad <- do.call(rbind,outs)
combine_dt_CT <- do.call(rbind,outs)

all_combine_dt <- rbind(combine_dt_broad,combine_dt_CT)

all_combine_dt$newqval <- ifelse(all_combine_dt$qval < 2e-16, 2e-16,all_combine_dt$qval )

all_combine_dt$spe_cate <- paste0(all_combine_dt$spe,'_',all_combine_dt$cate)
all_combine_dt$mlog10qval <- -log10(all_combine_dt$newqval)

all_combine_final_dt <- all_combine_dt[c('TEfam','mlog10qval','spe_cate')]


write.table(all_combine_final_dt, paste0(opt_dir,'/temp_TEfam_qval_enrich.txt'),quote = F,sep = '\t')

ipt_count_dt <- read.table(paste0(opt_dir,'/temp_TEfam_qval_enrich.txt'),stringsAsFactors = T)
head(ipt_count_dt)

##transfer to matrix
score_mtx <- sparseMatrix(i=as.numeric(ipt_count_dt$TEfam),
                          j=as.numeric(ipt_count_dt$spe_cate),
                          x=as.numeric(ipt_count_dt$mlog10qval),
                          dimnames=list(levels(ipt_count_dt$TEfam), levels(ipt_count_dt$spe_cate)))
score_mtx <- as.matrix(score_mtx)
score_mtx <- t(score_mtx)

dt_zscore <- apply(score_mtx,1,function(x){
  (x - mean(x))/sd(x)
})

dt_zscore <- data.frame(t(dt_zscore))

#cate_levels <- c('Distal_Distal','Proximal_Proximal','Genic_Genic','Distal_Proximal','Distal_Genic',
#                 'Proximal_Distal','Proximal_Genic','Genic_Distal','Genic_Proximal')
#colnames(dt_zscore) <- factor(colnames(dt_zscore),levels = cate_levels)

dt_zscore <- as.matrix(dt_zscore)
#dt_zscore <- dt_zscore[,cate_levels]

spe_cate <- c('maize_broad','maize_CT','sorghum_broad','sorghum_CT','Pm_broad','Pm_CT','Uf_broad','Uf_CT')
dt_zscore <- dt_zscore[spe_cate,]


pdf(paste0(opt_dir,'/opt_all_spe_TEfam_enrich.pdf'),width = 12,height = 12)
p <- pheatmap(dt_zscore,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color='white',
              fontsize_col = 15,
              fontsize_row = 15,
              fontsize_number = 25
) 
print(p)
dev.off()

#################
##updating 011424
##conduct enrichment for cell-type TE family enrichment
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123//TE_rice_to_others_dir//collect_allcelltype_enrich_dir/'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123/opt_check_rice_to_others_TE_dir/'


all_fl_list <- list.files(path= ipt_dir,pattern = '*')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_nonsyn_vs_syn_TE_enrich_fam_combineAllCT_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  ipt_dt
  
  ipt_dt <- ipt_dt[c('V1','V2','V3')]
  colnames(ipt_dt) <- c('Celltype','TE','pval')
  
  ipt_dt <- ipt_dt[!grepl('unknown',ipt_dt$Celltype),]
  
  #ipt_dt$qval <- p.adjust(ipt_dt$pval, method="fdr")
  ipt_dt <- ipt_dt[ipt_dt$pval < 0.05,]
  
  ipt_dt$spe <- target_fl_path_nm
  
  return(ipt_dt)
  
})

combine_dt <- do.call(rbind,outs)


combine_dt$spe_celltype <- paste0('rice_',combine_dt$spe,'__',combine_dt$Celltype)
combine_dt$mlog10pval <- -log10(combine_dt$pval)

all_combine_final_dt <- combine_dt[c('TE','mlog10pval','spe_celltype','pval')]


write.table(all_combine_final_dt, paste0(opt_dir,'/temp_TEcelltype_pval_enrich.txt'),quote = F,sep = '\t')

ipt_count_dt <- read.table(paste0(opt_dir,'/temp_TEcelltype_pval_enrich.txt'),stringsAsFactors = T)
head(ipt_count_dt)

##transfer to matrix
score_mtx <- sparseMatrix(i=as.numeric(ipt_count_dt$TE),
                          j=as.numeric(ipt_count_dt$spe_celltype),
                          x=as.numeric(ipt_count_dt$mlog10pval),
                          dimnames=list(levels(ipt_count_dt$TE), levels(ipt_count_dt$spe_celltype)))
score_mtx <- as.matrix(score_mtx)
score_mtx <- t(score_mtx)

#dt_zscore <- apply(score_mtx,1,function(x){
#  (x - mean(x))/sd(x)
#})

#dt_zscore <- data.frame(t(dt_zscore))

#cate_levels <- c('Distal_Distal','Proximal_Proximal','Genic_Genic','Distal_Proximal','Distal_Genic',
#                 'Proximal_Distal','Proximal_Genic','Genic_Distal','Genic_Proximal')
#colnames(dt_zscore) <- factor(colnames(dt_zscore),levels = cate_levels)

#dt_zscore <- as.matrix(dt_zscore)
#dt_zscore <- dt_zscore[,cate_levels]

#spe_cate <- c('maize_broad','maize_CT','sorghum_broad','sorghum_CT','Pm_broad','Pm_CT','Uf_broad','Uf_CT')
#dt_zscore <- dt_zscore[spe_cate,]

rownames(score_mtx)
target_order <- c('rice_maize__mesophyll','rice_sorghum__mesophyll','rice_Pm__mesophyll','rice_Uf__mesophyll',
                  'rice_maize__bundle_sheath','rice_sorghum__bundle_sheath','rice_Pm__bundle_sheath','rice_Uf__bundle_sheath',
                  'rice_maize__companion_cell','rice_sorghum__companion_cell','rice_Pm__companion_cell','rice_Uf__companion_cell',
                  'rice_maize__protoderm','rice_sorghum__protoderm','rice_Pm__protoderm','rice_Uf__protoderm',
                  'rice_maize__epidermis','rice_sorghum__epidermis','rice_Pm__epidermis','rice_Uf__epidermis'
)

intersect_order <- intersect(target_order,rownames(score_mtx))

score_mtx <- score_mtx[intersect_order,]

p <- pheatmap(score_mtx,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              color=colorRampPalette(c("white",rev(c("#440154","#3b528b","#21918c","#5ec962","#fde725","#fde725","#fde725","#fde725"))))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              
              ##ori true
              #color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
pdf(paste0(opt_dir,'/opt_all_spe_TE_celltype_enrich.pdf'),width = 10,height = 15)
p
dev.off()


#################
##updating 011724
##we will check the TE enrichment
ipt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123//TE_rice_to_others_dir//collect_bionomial_enrich_specific_TEfam_dir//'
opt_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_results_123123/opt_check_rice_to_others_TE_dir/'

ipt_motif_required_dir <- 'opt1_opt2_acr_syntenic_summary_dir_add_nonsyntenic_123123/ipt_required_raw_data/'

all_fl_list <- list.files(path= ipt_dir,pattern = '*')


outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_cate_motif_enrichment_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  
  meta_name_dt <- read.delim(paste0(ipt_motif_required_dir,'/opt_motif_common_name.txt'),header = F)
  head(meta_name_dt)
  
  ##add the common name
  combine_dt <- merge(ipt_dt,meta_name_dt,by.x='V3',by.y= 'V1')
  head(combine_dt)
  combine_dt$pval <- combine_dt$V4
  dim(combine_dt)
  
  #combine_dt <- combine_dt[combine_dt$pval < 0.3,]
  
  combine_dt <- combine_dt[order(combine_dt$V4,decreasing = F),]
  
  #combine_flt_dt <-  combine_dt[combine_dt$V7 > 0.1,]
  
  
  ##no result if we use the qval
  #combine_dt$qval <- p.adjust(combine_dt$V4, method="fdr")
  #combine_dt <- combine_dt[combine_dt$qval < 0.05,]
  
  ##modify the pval to double as we want the two sides
  combine_dt$pval_right <- combine_dt$pval*2
  
  
  combine_dt <- combine_dt[combine_dt$pval_right < 0.05,]
  
  combine_dt$mlog10pval <- -log10(combine_dt$pval_right)
  
  combine_dt <- combine_dt[c('V2.y','mlog10pval','V2.x','pval_right')]
  colnames(combine_dt) <- c('motif','mlog10pval','celltype','pval_right')
  
  #combine_dt <- combine_dt[order(combine_dt[[paste0('mlog10','pval')]],decreasing = T),]
  #combine_dt$motif <- factor(combine_dt$motif ,levels = combine_dt$motif)
  
  combine_dt$spe <- target_fl_path_nm
  

  combine_dt <- unique(combine_dt)
  return(combine_dt)
})

final_combine_dt <- do.call(rbind,outs)
#final_combine_dt$qval <- p.adjust(final_combine_dt$pval, method="fdr")
#final_combine_dt$mlog10pval <- -log10(final_combine_dt$pval)
#final_combine_dt <- final_combine_dt[c('spe','motif','mlog10pval')]




write.table(final_combine_dt, paste0(opt_dir,'/temp_spe_motif_target_TEfam_enrich.txt'),quote = F,sep = '\t')

ipt_count_dt <- read.table(paste0(opt_dir,'/temp_spe_motif_target_TEfam_enrich.txt'),stringsAsFactors = T)
head(ipt_count_dt)

##transfer to matrix
score_mtx <- sparseMatrix(i=as.numeric(ipt_count_dt$spe),
                          j=as.numeric(ipt_count_dt$motif),
                          x=as.numeric(ipt_count_dt$mlog10pval),
                          dimnames=list(levels(ipt_count_dt$spe), levels(ipt_count_dt$motif)))
score_mtx <- as.matrix(score_mtx)
score_mtx <- t(score_mtx)

score_mtx <- score_mtx[,c('maize','sorghum','Pm','Uf')]
p <- pheatmap(score_mtx,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              color=colorRampPalette(c("white",rev(c("#440154","#3b528b","#21918c","#5ec962","#fde725","#fde725","#fde725","#fde725"))))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              
              ##ori true
              #color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
p
pdf(paste0(opt_dir,'/opt_allspe_target_TE_enrich_motif_heatmap.pdf'),width = 8,height = 8)
print(p)
dev.off()




########
##Step02
########
#################################################################
##build the barplot of shared-acc proportion between broad and ct
##the fisher test it is significant
ipt_dt <- data.frame(
  Species = c("Os:Zm", "Os:Zm","Os:Sb", "Os:Sb","Os:Pm", "Os:Pm","Os:Uf", "Os:Uf"),
  Prop = c(2520/(2520+6921), 342/(342+1572), 5340/(5340 + 11174), 778/(778 + 2649),
           3367/12496,370/2607,4699/17582,553/3672),
  Cate = c("Broad_sharedA", "CT_sharedA", "Broad_sharedA", "CT_sharedA",
           "Broad_sharedA", "CT_sharedA", "Broad_sharedA", "CT_sharedA")
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Os:Zm','Os:Sb','Os:Pm','Os:Uf'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
  #facet_wrap(~Rice_CellType,nrow = 1)
pdf(paste0('opt2_Broad_sharedA_prop_compare_CT_sharedA_prop.pdf'),width = 8,height = 6 )
p
dev.off()

##updating 120423
##we will check shared-inacc ACRs enriched for the CT or not
##for rice to maize
shared_InAcc_CT_num <- 549
shared_InAcc_Broad_num <- 2301
Other_ACR_CT_num <- 1023 + 342
Other_ACR_Broad_num <- 2520 + 4520

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
##p-value = 0.0002467

##for rice to sorghum
shared_InAcc_CT_num <- 1172
shared_InAcc_Broad_num <- 4523
Other_ACR_CT_num <- 2255
Other_ACR_Broad_num <- 11991

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 2.639e-15

##for rice to Pm
shared_InAcc_CT_num <- 1087
shared_InAcc_Broad_num <- 4189
Other_ACR_CT_num <- 1520
Other_ACR_Broad_num <- 8307

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 3.538e-15

##for rice to Uf
shared_InAcc_CT_num <- 1544
shared_InAcc_Broad_num <- 6089
Other_ACR_CT_num <- 2128
Other_ACR_Broad_num <- 11493

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value < 2.2e-16

##We will plot the proportion 
ipt_dt <- data.frame(
  Species = c("Os:Zm", "Os:Zm","Os:Sb", "Os:Sb","Os:Pm", "Os:Pm","Os:Uf", "Os:Uf"),
  Prop = c(2301/(9441), 549/(1914), 4523/(16514), 1172/(3427),
           4189/12496,1087/2607,6089/17582,1544/3672),
  Cate = c("Broad_sharedInA", "CT_sharedInA", "Broad_sharedInA", "CT_sharedInA",
           "Broad_sharedInA", "CT_sharedInA", "Broad_sharedInA", "CT_sharedInA")
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Os:Zm','Os:Sb','Os:Pm','Os:Uf'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
#facet_wrap(~Rice_CellType,nrow = 1)
pdf(paste0('opt2_Broad_sharedInA_prop_compare_CT_sharedInA_prop.pdf'),width = 8,height = 6 )
p
dev.off()



##updating 120423
##we will check not-shared ACRs enriched for the CT or not
not_shared_CT_num <- 1023
not_shared_Broad_num <- 4620
Other_ACR_CT_num <- 891
Other_ACR_Broad_num <- 4821

data <- matrix(c(not_shared_CT_num, Other_ACR_CT_num, not_shared_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 0.000336


##for rice to sorghum
not_shared_CT_num <- 1477
not_shared_Broad_num <- 6651
Other_ACR_CT_num <- 1950
Other_ACR_Broad_num <- 9863

data <- matrix(c(not_shared_CT_num, Other_ACR_CT_num, not_shared_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 0.002243


##for rice to Pm
not_shared_CT_num <- 1150
not_shared_Broad_num <- 4940
Other_ACR_CT_num <- 1457
Other_ACR_Broad_num <- 7556

data <- matrix(c(not_shared_CT_num, Other_ACR_CT_num, not_shared_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 1.683e-05

##for rice to Uf
not_shared_CT_num <- 1575
not_shared_Broad_num <- 6794
Other_ACR_CT_num <- 2097
Other_ACR_Broad_num <- 10788

data <- matrix(c(not_shared_CT_num, Other_ACR_CT_num, not_shared_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 1.812e-06

##We will plot the proportion 
ipt_dt <- data.frame(
  Species = c("Os:Zm", "Os:Zm","Os:Sb", "Os:Sb","Os:Pm", "Os:Pm","Os:Uf", "Os:Uf"),
  Prop = c(4620/(9441), 1023/(1914), 6651/(16514), 1477/(3427),
           4940/12496,1150/2607,6794/17582,1575/3672),
  Cate = c("Broad_notshared", "CT_notshared", "Broad_notshared", "CT_notshared",
           "Broad_notshared", "CT_notshared", "Broad_notshared", "CT_notshared")
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Os:Zm','Os:Sb','Os:Pm','Os:Uf'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
#facet_wrap(~Rice_CellType,nrow = 1)
pdf(paste0('opt2_Broad_notshared_prop_compare_CT_notshared_prop.pdf'),width = 8,height = 6 )
p
dev.off()


##updating 120823
#################
##we get the flip and to see if the results are the same
ipt_dt <- data.frame(
  Species = c("Zm:Os", "Zm:Os","Sb:Os", "Sb:Os","Pm:Os", "Pm:Os","Uf:Os", "Uf:Os"),
  Prop = c(4002/(10446), 705/(2571), 5472/(18708), 653/(3113),
           3170/7970,687/2426,3695/8953,1280/4449),
  Cate = c("Broad_sharedA", "CT_sharedA", "Broad_sharedA", "CT_sharedA",
           "Broad_sharedA", "CT_sharedA", "Broad_sharedA", "CT_sharedA")
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Zm:Os','Sb:Os','Pm:Os','Uf:Os'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
  
#facet_wrap(~Rice_CellType,nrow = 1)
p
pdf(paste0('opt2_Broad_sharedA_prop_compare_CT_sharedA_prop_others_to_rice.pdf'),width = 8,height = 8 )
p
dev.off()

##for the shared A
shared_InAcc_CT_num <- 549
shared_InAcc_Broad_num <- 2301
Other_ACR_CT_num <- 1023 + 342
Other_ACR_Broad_num <- 2520 + 4520

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
##p-value = 0.0002467

##for rice to sorghum
shared_InAcc_CT_num <- 1172
shared_InAcc_Broad_num <- 4523
Other_ACR_CT_num <- 2255
Other_ACR_Broad_num <- 11991

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 2.639e-15

##for rice to Pm
shared_InAcc_CT_num <- 1087
shared_InAcc_Broad_num <- 4189
Other_ACR_CT_num <- 1520
Other_ACR_Broad_num <- 8307

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value = 3.538e-15

##for rice to Uf
shared_InAcc_CT_num <- 1544
shared_InAcc_Broad_num <- 6089
Other_ACR_CT_num <- 2128
Other_ACR_Broad_num <- 11493

data <- matrix(c(shared_InAcc_CT_num, Other_ACR_CT_num, shared_InAcc_Broad_num, Other_ACR_Broad_num), nrow = 2)
result <- fisher.test(data)
result
##p-value < 2.2e-16





#################
##for the sharedI
ipt_dt <- data.frame(
  Species = c("Zm:Os", "Zm:Os","Sb:Os", "Sb:Os","Pm:Os", "Pm:Os","Uf:Os", "Uf:Os"),
  Prop = c(2480/(10446), 752/(2571), 4224/(18708), 819/(3113),
           1283/7970,520/2426,1541/8953,968/4449),
  Cate = c("Broad_sharedInA", "CT_sharedInA", "Broad_sharedInA", "CT_sharedInA",
           "Broad_sharedInA", "CT_sharedInA", "Broad_sharedInA", "CT_sharedInA")
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Zm:Os','Sb:Os','Pm:Os','Uf:Os'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
#facet_wrap(~Rice_CellType,nrow = 1)
p
pdf(paste0('opt2_Broad_sharedA_prop_compare_CT_sharedI_prop_others_to_rice.pdf'),width = 8,height = 8 )
p
dev.off()




####################
##for the not shared
ipt_dt <- data.frame(
  Species = c("Zm:Os", "Zm:Os","Sb:Os", "Sb:Os","Pm:Os", "Pm:Os","Uf:Os", "Uf:Os"),
  Prop = c(3964/(10446), 1114/(2571), 9012/(18708), 1641/(3113),
           3517/7970,1219/2426,3717/8953,2201/4449),
  Cate = c("Broad_notshared", "CT_notshared", "Broad_notshared", "CT_notshared",
           "Broad_notshared", "CT_notshared", "Broad_notshared", "CT_notshared")
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Zm:Os','Sb:Os','Pm:Os','Uf:Os'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
#facet_wrap(~Rice_CellType,nrow = 1)
p
pdf(paste0('opt2_Broad_sharedA_prop_compare_CT_not_shared_prop_others_to_rice.pdf'),width = 8,height = 8 )
p
dev.off()






###################################################################################
##check the cell type specific ACR per cell type corresponds to the maize cell type
ipt_maize_dt <- read.delim('./opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_112023/opt2_syntenic_region_rice_in_maize_ACRnum_EachCelltype_adjustPair.txt',header = F)
ipt_maize_dt$species <- 'maize'
ipt_maize_dt <- ipt_maize_dt[ipt_maize_dt$V2 != 'broadly_accessible',]
head(ipt_maize_dt)
ipt_maize_dt <- ipt_maize_dt[ipt_maize_dt$V2 != 'procambial_meristem',]
aggregate_sum_dt <- aggregate(V3 ~ V1, data= ipt_maize_dt, sum)
merged_maize_dt <- merge(ipt_maize_dt,aggregate_sum_dt,by.x = 'V1',by.y = 'V1')
merged_maize_dt$prop <- merged_maize_dt$V3.x/merged_maize_dt$V3.y


ipt_sorghum_dt <- read.delim('./opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_112023/opt2_syntenic_region_rice_in_sorghum_ACRnum_EachCelltype_adjustPair.txt',header = F)
ipt_sorghum_dt$species <- 'sorghum'
ipt_sorghum_dt <- ipt_sorghum_dt[ipt_sorghum_dt$V2 != 'broadly_accessible',]
ipt_sorghum_dt <- ipt_sorghum_dt[ipt_sorghum_dt$V2 != 'unknown',]
ipt_sorghum_dt <- ipt_sorghum_dt[ipt_sorghum_dt$V2 != 'procambial_meristem',]
aggregate_sum_dt <- aggregate(V3 ~ V1, data= ipt_sorghum_dt, sum)
merged_sorghum_dt <- merge(ipt_sorghum_dt,aggregate_sum_dt,by.x = 'V1',by.y = 'V1')
merged_sorghum_dt$prop <- merged_sorghum_dt$V3.x/merged_sorghum_dt$V3.y


ipt_Pm_dt <- read.delim('./opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_112023/opt2_syntenic_region_rice_in_Pm_ACRnum_EachCelltype_adjustPair.txt',header = F)
ipt_Pm_dt$species <- 'Pm'
ipt_Pm_dt <- ipt_Pm_dt[ipt_Pm_dt$V2 != 'broadly_accessible',]
ipt_Pm_dt <- ipt_Pm_dt[ipt_Pm_dt$V2 != 'unknown',]
ipt_Pm_dt <- ipt_Pm_dt[ipt_Pm_dt$V2 != 'procambium',]
aggregate_sum_dt <- aggregate(V3 ~ V1, data= ipt_Pm_dt, sum)
merged_Pm_dt <- merge(ipt_Pm_dt,aggregate_sum_dt,by.x = 'V1',by.y = 'V1')
merged_Pm_dt$prop <- merged_Pm_dt$V3.x/merged_Pm_dt$V3.y

ipt_Uf_dt <- read.delim('./opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_112023/opt2_syntenic_region_rice_in_Uf_ACRnum_EachCelltype_adjustPair.txt',header = F)
ipt_Uf_dt$species <- 'Uf'
ipt_Uf_dt <- ipt_Uf_dt[ipt_Uf_dt$V2 != 'broadly_accessible',]
ipt_Uf_dt <- ipt_Uf_dt[ipt_Uf_dt$V2 != 'unknown',]
ipt_Uf_dt <- ipt_Uf_dt[ipt_Uf_dt$V2 != 'procambium',]
aggregate_sum_dt <- aggregate(V3 ~ V1, data= ipt_Uf_dt, sum)
merged_Uf_dt <- merge(ipt_Uf_dt,aggregate_sum_dt,by.x = 'V1',by.y = 'V1')
merged_Uf_dt$prop <- merged_Uf_dt$V3.x/merged_Uf_dt$V3.y


##merge the two ipt files
merged_two_dt <- rbind(merged_maize_dt,merged_sorghum_dt,merged_Pm_dt,merged_Uf_dt)

head(merged_two_dt)
dim(merged_two_dt)
colnames(merged_two_dt) <- c('Rice_CellType','other_CellType','Number','Species','Total','Prop')
merged_two_dt <- merged_two_dt[merged_two_dt$Rice_CellType != 'unknown_cells_1' & merged_two_dt$Rice_CellType != 'unknown_cells_2',]

##Have broad accessible version
##do not have broad accessible version
merged_two_dt <- merged_two_dt[merged_two_dt$other_CellType != 'broadly_accessible',]



merged_two_dt$Species <- factor(merged_two_dt$Species ,levels = c('maize','sorghum','Pm','Uf'))
merged_two_dt$Rice_CellType <- factor(merged_two_dt$Rice_CellType ,levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
merged_two_dt$other_CellType <- factor(merged_two_dt$other_CellType ,levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements','protoderm','epidermis'))



##Here we will use the composition bar plot to show
p <- ggplot(data=merged_two_dt, aes(x=Species, y=Prop, fill=other_CellType)) +
  scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
                             "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  #geom_bar(stat="identity", position=position_dodge()) +
  geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  facet_wrap(~Rice_CellType,nrow = 1)
            
pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_all_ratio.pdf'),width = 12,height = 6 )
p
dev.off()


##updating 011624
##we will check the bar plot of number other than the prop
p <- ggplot(data=merged_two_dt, aes(x=Species, y=Number, fill=other_CellType)) +
  scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
                             "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  #geom_bar(stat="identity", position=position_dodge()) +
  geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_ACR_num.pdf'),width = 12,height = 6 )
p
dev.off()


##updating 011624
##we will check the number of corresponding and not corresponding for each species pair
##epidermis how many of them are corresponding to the different species
merged_two_modi_CT_dt <-  merged_two_dt
merged_two_modi_CT_dt$other_CellType <- gsub('companion_cells_sieve_elements','companion_cell',merged_two_modi_CT_dt$other_CellType)
merged_two_modi_CT_dt$Rice_CellType <- gsub('epidermis','EpiCell',merged_two_modi_CT_dt$Rice_CellType)
merged_two_modi_CT_dt$Rice_CellType <- gsub('protoderm','EpiCell',merged_two_modi_CT_dt$Rice_CellType)
merged_two_modi_CT_dt$other_CellType <- gsub('epidermis','EpiCell',merged_two_modi_CT_dt$other_CellType)
merged_two_modi_CT_dt$other_CellType <- gsub('protoderm','EpiCell',merged_two_modi_CT_dt$other_CellType)

merged_two_modi_CT_dt$check_corres <- ifelse(merged_two_modi_CT_dt$Rice_CellType == merged_two_modi_CT_dt$other_CellType, 'yes','no')

aggregate_dt <- aggregate(Number ~ Species + check_corres,data = merged_two_modi_CT_dt,sum)

aggregate_spe_dt <- aggregate(Number ~ Species, data = aggregate_dt,sum)

aggregate_dt$check_corres <- factor(aggregate_dt$check_corres,levels = c('yes','no'))
p <- ggplot(data=aggregate_dt, aes(x=Species, y=Number, fill=check_corres)) +
  #scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
  #                           "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
p
#  facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_all_num_barplot_summary_all.pdf'),width = 8,height = 8 )
p
dev.off()

##updating 011824
##we will build the pie charts for showing the compositions 
generate_pie <- function(df) {
  #df <- dt[dt$celltype == target_celltype,]
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df$cate,pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}

spe_list <- unique(aggregate_dt$Species)
for (i in (1:length(spe_list))){
  
  target_spe = spe_list[i]
  aggregate_spe_dt <- aggregate_dt[aggregate_dt$Species == target_spe,]
  aggregate_spe_dt <- aggregate_spe_dt[c('check_corres','Number')]
  colnames(aggregate_spe_dt) <- c('cate','count')
  
  pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_all_num_piechar_',target_spe,'.pdf'),width = 8,height = 8 )
  print(generate_pie(aggregate_spe_dt))
  dev.off()
  
  
}





##check shifting of meosphyll and bundle sheath
merged_two_mesophyll_dt <- merged_two_dt[merged_two_dt$Rice_CellType == 'mesophyll',]
merged_two_mesophyll_dt$corr <- ifelse(merged_two_mesophyll_dt$other_CellType == 'mesophyll', 'yes','no')
merged_two_mesophyll_dt$corr_bundle <- ifelse(merged_two_mesophyll_dt$other_CellType == 'bundle_sheath','yes','no')
merged_two_mesophyll_nomeso_dt <- merged_two_mesophyll_dt[merged_two_mesophyll_dt$corr == 'no',]
aggregate_dt <- aggregate(Number ~ Species + corr_bundle, data = merged_two_mesophyll_nomeso_dt, sum)


aggregate_dt$corr_bundle <- factor(aggregate_dt$corr_bundle,levels = c('yes','no'))
p <- ggplot(data=aggregate_dt, aes(x=Species, y=Number, fill=corr_bundle)) +
  #scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
  #                           "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 

#  facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_all_num_barplot_mesop_shift_to_bundle_andothers.pdf'),width = 8,height = 8 )
p
dev.off()

##build the piecharts
merged_two_mesophyll_dt <- merged_two_dt[merged_two_dt$Rice_CellType == 'mesophyll',]
spe_list <- unique(merged_two_mesophyll_dt$Species)
for (i in (1:length(spe_list))){
  
  target_spe = spe_list[i]
  merged_two_mesophyll_spe_dt <- merged_two_mesophyll_dt[merged_two_mesophyll_dt$Species == target_spe,]
  merged_two_mesophyll_spe_dt <- merged_two_mesophyll_spe_dt[c('other_CellType','Number')]
  colnames(merged_two_mesophyll_spe_dt) <- c('cate','count')
  
  pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/store_mesophyll_to_bundlesheath_piechart/opt2_final_celltype_in_rice_piechar_',target_spe,'.pdf'),width = 8,height = 8 )
  print(generate_pie(merged_two_mesophyll_spe_dt))
  dev.off()
  
  
}






##check shift of bundle sheath to mesophyll
merged_two_bundle_dt <- merged_two_dt[merged_two_dt$Rice_CellType == 'bundle_sheath',]
merged_two_bundle_dt$corr <- ifelse(merged_two_bundle_dt$other_CellType == 'bundle_sheath', 'yes','no')
merged_two_bundle_dt$corr_meso <- ifelse(merged_two_bundle_dt$other_CellType == 'mesophyll','yes','no')
merged_two_bundel_nobundle_dt <- merged_two_bundle_dt[merged_two_bundle_dt$corr == 'no',]
aggregate_dt <- aggregate(Number ~ Species + corr_meso, data = merged_two_bundel_nobundle_dt, sum)


aggregate_dt$corr_meso <- factor(aggregate_dt$corr_meso,levels = c('yes','no'))
p <- ggplot(data=aggregate_dt, aes(x=Species, y=Number, fill=corr_meso)) +
  #scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
  #                           "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 

#  facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_all_num_barplot_bundle_shift_to_meso_andothers.pdf'),width = 8,height = 8 )
p
dev.off()


##build the piecharts
merged_two_mesophyll_dt <- merged_two_dt[merged_two_dt$Rice_CellType == 'bundle_sheath',]
spe_list <- unique(merged_two_mesophyll_dt$Species)
for (i in (1:length(spe_list))){
  
  target_spe = spe_list[i]
  merged_two_mesophyll_spe_dt <- merged_two_mesophyll_dt[merged_two_mesophyll_dt$Species == target_spe,]
  merged_two_mesophyll_spe_dt <- merged_two_mesophyll_spe_dt[c('other_CellType','Number')]
  colnames(merged_two_mesophyll_spe_dt) <- c('cate','count')
  
  pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/store_bundlesheat_to_mesophyll_piechart/opt2_final_celltype_in_rice_piechar_',target_spe,'.pdf'),width = 8,height = 8 )
  print(generate_pie(merged_two_mesophyll_spe_dt))
  dev.off()
  
  
}




##check companion cell shared
merged_two_cc_dt <- merged_two_dt[merged_two_dt$Rice_CellType == 'companion_cell',]
merged_two_cc_dt$other_CellType <- gsub('companion_cells_sieve_elements','companion_cell',merged_two_cc_dt$other_CellType)
merged_two_cc_dt$corrCC <- ifelse(merged_two_cc_dt$Rice_CellType == merged_two_cc_dt$other_CellType, 'yes','no')

merged_two_cc_dt$corrCC <- factor(aggregate_dt$corrCC,levels = c('yes','no'))

p <- ggplot(data=merged_two_cc_dt, aes(x=Species, y=Number, fill=corrCC)) +
  #scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
  #                           "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  
  #facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_CC_to_CC_andothers_barplot.pdf'),width = 8,height = 8 )
p
dev.off()








##updating 011624
##we will build a bar plot to show the number
mergedt_two_target_dt <- unique(merged_two_dt[c('Rice_CellType','Species','Total')])
mergedt_two_target_dt$Species <- factor(mergedt_two_target_dt$Species ,levels = c('maize','sorghum','Pm','Uf'))
mergedt_two_target_dt$Rice_CellType <- factor(mergedt_two_target_dt$Rice_CellType ,levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))

p <- ggplot(data=mergedt_two_target_dt, aes(x=Species, y=Total, fill=Rice_CellType)) +
  scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                             "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0('opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023/','/opt2_final_celltype_in_rice_to_all_num_barplot.pdf'),width = 12,height = 4 )
p
dev.off()



##updating 020124
##we will calculate the how many ct in rice corresponding to all species 
##in the shared ACRs
ipt_dir <- 'opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_112023/ipt_ACR_rice_to_other_correspondings_020124//////'
opt_dir <- 'opt2_syntenic_CT_ACR_correspond_prop_EachCelltype_results_112023//'


all_fl_list <- list.files(path = ipt_dir)

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('opt1_syntenic_region_','',x)
  target_fl_path_nm <- gsub('_SharedNotSharedCate_celltypes.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  colnames(ipt_dt) <- c('ACR','CelltypeCate','Category')
  
  head(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$CelltypeCate != 'broadly_accessible',]
  ipt_dt <- ipt_dt[!grepl('unknown',ipt_dt$CelltypeCate),]
  ipt_dt <- ipt_dt[ipt_dt$Category == 'SharedAcc',]
  dim(ipt_dt)
  
  ipt_dt$spe <- target_fl_path_nm
  
  return(ipt_dt)

})

combine_dt <- do.call(rbind,outs)
combine_dt$count <- 1

aggregate_dt <- aggregate(count ~ ACR + CelltypeCate, data = combine_dt, sum)

aggregate_dt$countCate <- paste0('count', aggregate_dt$count)
aggregate_dt$count2 <- 1
aggregate_dt2 <- aggregate(count2 ~ countCate + CelltypeCate,  data = aggregate_dt,sum)

write.csv(aggregate_dt2, paste0(opt_dir,'/opt_shared_ACR_count_eachcelltype_count.csv'),quote = F)

#################
##updating 111523
##we will check the distal proximal and genic for these cases
##updating 112023 for all the species
ipt_dir <- 'opt2_different_cate_num_zscore_riceTomazie_dir_112023////'
opt_heatmap_dir <- 'opt2_different_cate_num_zscore_riceTomazie_heatmap_dir_results_112023//'
opt_barplot_dir <- 'opt2_different_cate_num_zscore_riceTomazie_barplot_dir_results_112023/'

all_fl_list <- list.files(path= ipt_dir)

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt2_forBothACR_syntenic_region_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_ACRName_CelltypeCate.txt','',target_fl_path_nm)
  sp_nm <- gsub('rice_in_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  
  head(ipt_dt)
  ipt_dt$GeneCatePair <- paste0(ipt_dt$riceGeneCate,'_',ipt_dt[[paste0(sp_nm,'GeneCate')]])
  head(ipt_dt)
  ipt_dt$CTcatePair <- paste0(ipt_dt$riceCTcate,'_',ipt_dt[[paste0(sp_nm,'CTcate')]])
  head(ipt_dt)
  
  ipt_dt$count <- 1
  
  ##first we will build a bar plot to show the number
  aggregate_CTcatePair_dt <- aggregate(count ~ CTcatePair, ipt_dt, sum)
  
  p <- ggplot(data=aggregate_CTcatePair_dt, aes(x=CTcatePair, y=count)) +
    #geom_bar(stat="identity", position=position_dodge()) +
    geom_bar(stat="identity") +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=10),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  
  pdf(paste0(opt_barplot_dir,'/opt2_different_cate_num_totalCTcatePair_',target_fl_path_nm,'.pdf'),width = 5,height = 8 )
  print(p)
  dev.off()
  
  
  aggregate_dt <- aggregate(count ~ CTcatePair + GeneCatePair, data = ipt_dt, sum)
  
  write.table(aggregate_dt, paste0(opt_heatmap_dir,'/temp_CTcate_Genecate_pair_count.txt'),quote = F,sep = '\t')
  
  ipt_count_dt <- read.table(paste0(opt_heatmap_dir,'/temp_CTcate_Genecate_pair_count.txt'),stringsAsFactors = T)
  head(ipt_count_dt)
  
  ##transfer to matrix
  score_mtx <- sparseMatrix(i=as.numeric(ipt_count_dt$CTcatePair),
                            j=as.numeric(ipt_count_dt$GeneCatePair),
                            x=as.numeric(ipt_count_dt$count),
                            dimnames=list(levels(ipt_count_dt$CTcatePair), levels(ipt_count_dt$GeneCatePair)))
  score_mtx <- as.matrix(score_mtx)
  
  
  dt_zscore <- apply(score_mtx,1,function(x){
    (x - mean(x))/sd(x)
  })
  
  dt_zscore <- data.frame(t(dt_zscore))
  
  cate_levels <- c('Distal_Distal','Proximal_Proximal','Genic_Genic','Distal_Proximal','Distal_Genic',
                   'Proximal_Distal','Proximal_Genic','Genic_Distal','Genic_Proximal')
  colnames(dt_zscore) <- factor(colnames(dt_zscore),levels = cate_levels)
  
  dt_zscore <- as.matrix(dt_zscore)
  dt_zscore <- dt_zscore[,cate_levels]
  
  
  pdf(paste0(opt_heatmap_dir,'/opt2_different_cate_num_zscore_',target_fl_path_nm,'.pdf'),width = 5,height = 3)
  p <- pheatmap(dt_zscore,
           #scale="row",
           #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
           #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
           cluster_cols = F,
           cluster_rows = F,
           border_color='white',
           fontsize_col = 15,
           fontsize_row = 15,
           fontsize_number = 25
  ) 
  print(p)
  dev.off()
  

   
}


#################
##updating 120423
##we will check the sharedInAcc ACR cate
ipt_dir <- 'opt2_different_cate_num_zscore_riceTomazie_dir_sharedInAcc_proximity_120423//'
opt_dir <- 'opt2_different_cate_num_zscore_riceTomazie_dir_sharedInAcc_proximity_results_120423//'

ipt_dir <- 'opt2_different_cate_num_zscore_riceTomazie_dir_Notshared_proximity_120423/////'
opt_dir <- 'opt2_different_cate_num_zscore_riceTomazie_dir_Notshared_proximity_results_120423//'


all_fl_list <- list.files(path = ipt_dir)
outs <- lapply(all_fl_list, function(x){
  
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('opt2_for_ACR_sharedInAcc_syntenic_region_rice_in_','',x)
  target_fl_path_nm <- gsub('_ACRName_CelltypeCate.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  ipt_dt$count <- 1
  colnames(ipt_dt) <- c('Cate','ACRloc','Proximity','count')
  aggregate_dt <- aggregate(count ~ Cate + Proximity, data = ipt_dt, sum)
  aggregate_dt$spe <- target_fl_path_nm
  
  return(aggregate_dt)

})

combine_dt <- do.call(rbind,outs)

write.table(combine_dt,paste0(opt_dir,'/opt_proximity_num_all_spe.txt'),quote = F,sep= '\t')















ipt_dt <- read.delim('./opt2_forBothACR_syntenic_region_rice_in_maize_ACRName_CelltypeCate.txt')
prefix <- 'riceTomazie'
head(ipt_dt)
ipt_dt$GeneCatePair <- paste0(ipt_dt$riceGeneCate,'_',ipt_dt$maizeGeneCate)
head(ipt_dt)
ipt_dt$CTcatePair <- paste0(ipt_dt$riceCTcate,'_',ipt_dt$maizeCTcate)
head(ipt_dt)

ipt_dt <- read.delim('./opt2_forBothACR_syntenic_region_rice_in_sorghum_ACRName_CelltypeCate.txt')
prefix <- 'riceTosorghum'
head(ipt_dt)
ipt_dt$GeneCatePair <- paste0(ipt_dt$riceGeneCate,'_',ipt_dt$sorghumGeneCate)
head(ipt_dt)
ipt_dt$CTcatePair <- paste0(ipt_dt$riceCTcate,'_',ipt_dt$sorghumCTcate)
head(ipt_dt)


ipt_dt$count <- 1


##first we will build a bar plot to show the number
aggregate_CTcatePair_dt <- aggregate(count ~ CTcatePair, ipt_dt, sum)

p <- ggplot(data=aggregate_CTcatePair_dt, aes(x=CTcatePair, y=count)) +
  #geom_bar(stat="identity", position=position_dodge()) +
  geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )

pdf(paste0('opt2_different_cate_num_totalCTcatePair_',prefix,'.pdf'),width = 5,height = 8 )
p
dev.off()



aggregate_dt <- aggregate(count ~ CTcatePair + GeneCatePair, data = ipt_dt, sum)

write.table(aggregate_dt, 'temp_CTcate_Genecate_pair_count.txt',quote = F,sep = '\t')

ipt_count_dt <- read.table('temp_CTcate_Genecate_pair_count.txt',stringsAsFactors = T)
head(ipt_count_dt)

##transfer to matrix
score_mtx <- sparseMatrix(i=as.numeric(ipt_count_dt$CTcatePair),
                          j=as.numeric(ipt_count_dt$GeneCatePair),
                          x=as.numeric(ipt_count_dt$count),
                          dimnames=list(levels(ipt_count_dt$CTcatePair), levels(ipt_count_dt$GeneCatePair)))
score_mtx <- as.matrix(score_mtx)


dt_zscore <- apply(score_mtx,1,function(x){
  (x - mean(x))/sd(x)
})

dt_zscore <- data.frame(t(dt_zscore))

cate_levels <- c('Distal_Distal','Proximal_Proximal','Genic_Genic','Distal_Proximal','Distal_Genic',
                 'Proximal_Distal','Proximal_Genic','Genic_Distal','Genic_Proximal')
colnames(dt_zscore) <- factor(colnames(dt_zscore),levels = cate_levels)

dt_zscore <- as.matrix(dt_zscore)
dt_zscore <- dt_zscore[,cate_levels]


pdf(paste0('opt2_different_cate_num_zscore_',prefix,'.pdf'),width = 5,height = 3)
pheatmap(dt_zscore,
         #scale="row",
         #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
         #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
         cluster_cols = F,
         cluster_rows = F,
         border_color='white',
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize_number = 25
) 
dev.off()

















##updating 111723
########################################################
##check the syntenic ACR overlapping with the atlas ACRs
ipt_dir <- 'opt2_syntenic_ACR_add_atlasACR_111723///'
opt_dir <- 'opt2_syntenic_ACR_add_atlasACR_results_111723/'

all_fl_list <- list.files(path= ipt_dir)

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_syntenic_ACR_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_addAtlasACRcelltype.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  ipt_dt <- ipt_dt[c('V4','V10')]
  colnames(ipt_dt) <- c('Cate','Num')
  
  ipt_dt <- ipt_dt[ipt_dt$Cate != 'All',]
  
  ipt_dt$Cate <- factor(ipt_dt$Cate,levels =c('NotShared','SharedAcc','SharedInAcc'))
  
  ##we first build the density plot
  p<-ggplot(ipt_dt, aes(x=Num, color = Cate,fill = Cate)) +
    geom_density(alpha=0.6) + 
    scale_fill_manual(values=c("SharedAcc"="#ED1C24", "SharedInAcc"="#FFDE17", "NotShared"="#1CBDC2")) +
    scale_color_manual(values=c("SharedAcc"="#ED1C24", "SharedInAcc"="#FFDE17", "NotShared"="#1CBDC2")) 
  #p <- p + coord_cartesian(xlim = c(0, 2))
  #p <- p + scale_x_continuous(breaks=seq(0,2, 0.25))
  p <- p + theme(
    plot.title = element_text(size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
  
  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'.pdf'),width = 10,height = 6 )
  print(p)
  dev.off()
  
  ##we will also make a bar plot to compare the significant
  ipt_dt$Cate <- factor(ipt_dt$Cate,levels =c('SharedAcc','SharedInAcc','NotShared'))
  p<-ggplot(ipt_dt, aes(x=Cate, y=Num, fill=Cate)) +
    geom_boxplot() + 
    theme(
      plot.title = element_text(size=20,hjust = 0.5),
      axis.title = element_text(size =20),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black"),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'_boxplot.pdf'),width = 6,height = 6 )
  print(p)
  dev.off()
  
  
  ipt_sharedAcc_dt <- ipt_dt[ipt_dt$Cate == 'SharedAcc',]
  ipt_sharedInAcc_dt <- ipt_dt[ipt_dt$Cate == 'SharedInAcc',]
  ipt_noshared_dt <- ipt_dt[ipt_dt$Cate == 'NotShared',]
  test_sharedAcc_sharedInAcc = wilcox.test(ipt_sharedAcc_dt$Num,ipt_sharedInAcc_dt$Num)
  test_sharedAcc_notshared = wilcox.test(ipt_sharedAcc_dt$Num,ipt_noshared_dt$Num)
  test_shareInAcc_notshared = wilcox.test(ipt_sharedInAcc_dt$Num,ipt_noshared_dt$Num)
  
  df <- data.frame(
    compare = c('sharedAcc_sharedInAcc','sharedAcc_notshared','sharedInAcc_notshared'),
    pval = c(test_sharedAcc_sharedInAcc$p.value,test_sharedAcc_notshared$p.value,test_shareInAcc_notshared$p.value)
  )
  
  write.table(df,paste0(opt_dir,'/opt_',target_fl_path_nm,'_pval.txt'),quote = F,sep = '\t')
  

}

##updating 112523
#############################################################################
##heck the syntenic ACR overlapping with the atlas ACRs with different CTcate
ipt_dir <- 'opt2_syntenic_ACR_add_atlasACR_addCTcate_112523////'
opt_dir <- 'opt2_syntenic_ACR_add_atlasACR_addCTcate_results_112523///'

all_fl_list <- list.files(path= ipt_dir)

##two options
##1) Broad
##2) CT
prefix <- 'Broad'
prefix <- 'CT'

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_syntenic_ACR_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_addAtlasACRcelltype_addCTcate.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  
  if (prefix == 'Broad'){
    ipt_dt <- ipt_dt[ipt_dt$V5 == 'broadly_accessible',]
  }
  if (prefix == 'CT'){
    ipt_dt <- ipt_dt[ipt_dt$V5 != 'broadly_accessible',]
  } 
  
  ipt_dt <- ipt_dt[c('V4','V11')]
  colnames(ipt_dt) <- c('Cate','Num')
  
  ipt_dt <- ipt_dt[ipt_dt$Cate != 'All',]
  
  ipt_dt$Cate <- factor(ipt_dt$Cate,levels =c('NotShared','SharedAcc','SharedInAcc'))
  
  ##we first build the density plot
  p<-ggplot(ipt_dt, aes(x=Num, color = Cate,fill = Cate)) +
    geom_density(alpha=0.6) + 
    scale_fill_manual(values=c("SharedAcc"="#ED1C24", "SharedInAcc"="#FFDE17", "NotShared"="#1CBDC2")) +
    scale_color_manual(values=c("SharedAcc"="#ED1C24", "SharedInAcc"="#FFDE17", "NotShared"="#1CBDC2")) 
  #p <- p + coord_cartesian(xlim = c(0, 2))
  #p <- p + scale_x_continuous(breaks=seq(0,2, 0.25))
  p <- p + theme(
    plot.title = element_text(size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
  
  pdf(paste0(opt_dir,'/opt_',prefix,'_',target_fl_path_nm,'.pdf'),width = 10,height = 6 )
  print(p)
  dev.off()
  
  ##we will also make a bar plot to compare the significant
  ipt_dt$Cate <- factor(ipt_dt$Cate,levels =c('SharedAcc','SharedInAcc','NotShared'))
  p<-ggplot(ipt_dt, aes(x=Cate, y=Num, fill=Cate)) +
    geom_boxplot() + 
    theme(
      plot.title = element_text(size=20,hjust = 0.5),
      axis.title = element_text(size =20),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black"),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  pdf(paste0(opt_dir,'/opt_',prefix,'_',target_fl_path_nm,'_boxplot.pdf'),width = 6,height = 6 )
  print(p)
  dev.off()
  
  
  ipt_sharedAcc_dt <- ipt_dt[ipt_dt$Cate == 'SharedAcc',]
  ipt_sharedInAcc_dt <- ipt_dt[ipt_dt$Cate == 'SharedInAcc',]
  ipt_noshared_dt <- ipt_dt[ipt_dt$Cate == 'NotShared',]
  test_sharedAcc_sharedInAcc = wilcox.test(ipt_sharedAcc_dt$Num,ipt_sharedInAcc_dt$Num)
  test_sharedAcc_notshared = wilcox.test(ipt_sharedAcc_dt$Num,ipt_noshared_dt$Num)
  test_shareInAcc_notshared = wilcox.test(ipt_sharedInAcc_dt$Num,ipt_noshared_dt$Num)
  
  df <- data.frame(
    compare = c('sharedAcc_sharedInAcc','sharedAcc_notshared','sharedInAcc_notshared'),
    pval = c(test_sharedAcc_sharedInAcc$p.value,test_sharedAcc_notshared$p.value,test_shareInAcc_notshared$p.value)
  )
  
  write.table(df,paste0(opt_dir,'/opt_',prefix,'_',target_fl_path_nm,'_pval.txt'),quote = F,sep = '\t')
  
  
}

##udpating 120823
##############################################################################
##check the syntenic ACR overlapping with the atlas ACRs with different CTcate
ipt_dir <- 'opt2_syntenic_ACR_add_atlasACR_addCTcate_112523////'
opt_dir <- 'opt2_syntenic_ACR_add_atlasACR_addCTcate_results_112523///'

all_fl_list <- list.files(path= ipt_dir)

##two options
##1) Broad
##2) CT
#prefix <- 'Broad'
prefix <- 'CT'

target_cate <- 'NotShared'

for (i in (1:length(all_fl_list))){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_syntenic_ACR_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_addAtlasACRcelltype_addCTcate.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  
  if (prefix == 'Broad'){
    ipt_dt <- ipt_dt[ipt_dt$V5 == 'broadly_accessible',]
  }
  if (prefix == 'CT'){
    ipt_dt <- ipt_dt[ipt_dt$V5 != 'broadly_accessible',]
  } 
  
  if (target_cate == 'NotShared'){
    ipt_dt <- ipt_dt[ipt_dt$V4 == target_cate,]
  }
  
  
  ipt_dt <- ipt_dt[c('V5','V11')]
  colnames(ipt_dt) <- c('Cate','Num')
  head(ipt_dt)
  
  ipt_dt <- as.data.frame(separate_rows(ipt_dt, Cate, sep = ","))
  
  
  ipt_dt <- ipt_dt[ipt_dt$Cate != 'unknown_cells_1' & ipt_dt$Cate != 'unknown_cells_2',]
  

  
  ##we first build the density plot
  p<-ggplot(ipt_dt, aes(x=Num, color = Cate,fill = Cate)) +
    geom_density(alpha=0.6) + 
    scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                               "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
    scale_color_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                               "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))
  #p <- p + coord_cartesian(xlim = c(0, 2))
  #p <- p + scale_x_continuous(breaks=seq(0,2, 0.25))
  p <- p + theme(
    plot.title = element_text(size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
  p
  
  pdf(paste0(opt_dir,'/opt_',prefix,'_',target_fl_path_nm,'_',target_cate,'_celltype.pdf'),width = 10,height = 6 )
  print(p)
  dev.off()
  
  ##we will also make a bar plot to compare the significant
  ipt_dt$Cate <- factor(ipt_dt$Cate,levels =c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
  p<-ggplot(ipt_dt, aes(x=Cate, y=Num, fill=Cate)) +
    scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                               "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
    geom_boxplot() + 
    theme(
      plot.title = element_text(size=20,hjust = 0.5),
      axis.title = element_text(size =20),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black"),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  pdf(paste0(opt_dir,'/opt_',prefix,'_',target_fl_path_nm,'_',target_cate,'_celltype_boxplot.pdf'),width = 6,height = 6 )
  print(p)
  dev.off()
  
  
  
  
  
}
















##updating 111723 
###########################################
##check the syntenic ACR overlap three cate
ipt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_111723////'
opt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_results_111723//'

all_fl_list <- list.files(path= ipt_dir)

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt1_syntenic_region_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_SharedNotSharedCate_celltypes.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header=F)
  head(ipt_dt)
  colnames(ipt_dt) <- c('ACR','celltype','cate')
  ipt_dt$count <- 1  
  dim(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'broadly_accessible',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown_cells_1' & ipt_dt$celltype != 'unknown_cells_2',]
  dim(ipt_dt)
  
  aggregate_dt <- aggregate(count ~ cate + celltype, ipt_dt, sum)
  aggregate_cate_dt <- aggregate(count ~ cate, ipt_dt, sum)
  
  merged_dt <- merge(aggregate_dt,aggregate_cate_dt, by.x = 'cate',by.y = 'cate')
  
  merged_dt$prop <- merged_dt$count.x/merged_dt$count.y
  
  merged_dt$spe <- target_fl_path_nm
  
  return(merged_dt)
  
})

combine_dt <- do.call(rbind,outs)

head(combine_dt)

combine_dt$cate <- factor(combine_dt$cate, levels = c('SharedAcc','SharedInAcc','NotShared','NotBlast'))
combine_dt$celltype <- factor(combine_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cell',
                                                              'epidermis','protoderm'))

combine_dt <- combine_dt[combine_dt$cate != 'NotBlast',]

combine_dt$spe <- factor(combine_dt$spe, levels = c('rice_maize','rice_sorghum','rice_Pm','rice_Uf'))

p <- ggplot(data=combine_dt, aes(x=celltype, y=prop, fill=cate)) +
  scale_fill_manual(values=c("SharedAcc"="#ED1C24", "SharedInAcc"="#FFDE17", "NotShared"="#72CAC4","NotBlast" = "#A7A9AC")) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  facet_wrap(~spe,nrow = 1)

pdf(paste0(opt_dir,'/','opt2_celltype_cate_prop.pdf'),width = 12,height = 6 )
p
dev.off()

##updating 120523
##we will plot the figures
all_spe_list <- unique(combine_dt$spe)

lapply(all_spe_list, function(x){
  
  ipt_target_spe_dt <- combine_dt[combine_dt$spe == x,]
  
  p <- ggplot(data=ipt_target_spe_dt, aes(x=celltype, y=prop,fill=celltype)) +
    scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                               "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1")) + 
    geom_bar(stat="identity", position=position_dodge()) +
    #geom_bar(stat="identity") +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=10),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) + 
    facet_wrap(~cate,nrow = 1)
  
  pdf(paste0(opt_dir,'/','opt2_',x,'_celltype_cate_prop.pdf'),width = 12,height = 6 )
  print(p)
  dev.off()
  
})


##updating 120623
##we will check the proprotion of each cell type under each category 
ipt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_111723////'
opt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_results_111723//'
cate_prefix <- 'RiceToOthers'

#for others to rice
ipt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_others_to_rice_121023/'
opt_dir <- 'opt2_syntenic_ACR_Shared_IA_NotShared_celltypes_others_to_rice_results_121023/'
cate_prefix <- 'OthersToRice'


all_fl_list <- list.files(path= ipt_dir)

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt1_syntenic_region_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_SharedNotSharedCate_celltypes.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header=F)
  head(ipt_dt)
  colnames(ipt_dt) <- c('ACR','celltype','cate')
  ipt_dt$count <- 1  
  dim(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'broadly_accessible',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown_cells_1' & ipt_dt$celltype != 'unknown_cells_2',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'unknown'&ipt_dt$celltype != 'procambial_meristem',]
  ipt_dt <- ipt_dt[ipt_dt$celltype != 'procambium',]
  dim(ipt_dt)
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$cate != 'NotBlast',]
  

  
  aggregate_dt <- aggregate(count ~ cate + celltype, ipt_dt, sum)
  aggregate_celltype_dt <- aggregate(count ~ celltype, ipt_dt, sum)
  
  merged_dt <- merge(aggregate_dt,aggregate_celltype_dt, by.x = 'celltype',by.y = 'celltype')
  
  merged_dt$prop <- merged_dt$count.x/merged_dt$count.y
  
  merged_dt$spe <- target_fl_path_nm
  
  return(merged_dt)
  
})

combine_dt <- do.call(rbind,outs)

head(combine_dt)

combine_dt$cate <- factor(combine_dt$cate, levels = c('SharedAcc','SharedInAcc','NotShared','NotBlast'))





combine_dt <- combine_dt[combine_dt$cate != 'NotBlast',]

if (cate_prefix == 'RiceToOthers'){
  combine_dt$spe <- factor(combine_dt$spe, levels = c('rice_maize','rice_sorghum','rice_Pm','rice_Uf'))
  color_value <- c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                   "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1")
  combine_dt$celltype <- factor(combine_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cell',
                                                                'protoderm','epidermis'))
}
if (cate_prefix == 'OthersToRice'){
  combine_dt$spe <- factor(combine_dt$spe, levels = c('maize_rice','sorghum_rice','Pm_rice','Uf_rice'))
  color_value <- c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
                   "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1")
  combine_dt$celltype <- factor(combine_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cells_sieve_elements',
                                                                'protoderm','epidermis'))
}




p <- ggplot(data=combine_dt, aes(x=cate, y=prop, fill=celltype)) +
  scale_fill_manual(values=color_value) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_line(stat="identity", position=position_dodge())+
  #geom_point(stat="identity", position=position_dodge())+
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + 
  facet_wrap(~spe,nrow = 1)
p

pdf(paste0(opt_dir,'/','opt2_celltype_addTobe1_prop_for_threecate.pdf'),width = 14,height = 6 )
p
dev.off()


##updating 010624
##we will plot the epidermis to each cell types
combine_dt_SS <- combine_dt[combine_dt$cate == 'NotShared',]



  
#combine_dt_SS$celltype <- factor(combine_dt_SS$celltype,levels = rev(c('mesophyll','bundle_sheath','companion_cell',
#                                                                   'protoderm','epidermis')))
  
p <- ggplot(data=combine_dt_SS, aes(x=celltype, y=prop, fill=celltype)) +
  scale_fill_manual(values=color_value) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_line(stat="identity", position=position_dodge())+
  #geom_point(stat="identity", position=position_dodge())+
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  facet_wrap(~ spe, ncol= 1)+
  scale_y_continuous(breaks=seq(0, 0.8, 0.2)) +
  coord_cartesian(ylim = c(0,0.8))
  
p
pdf(paste0(opt_dir,'/','opt2_celltype_addTobe1_prop_for_threecate_only_Species_specific.pdf'),width = 9,height = 14 )
p
dev.off()

  



combine_dt$new_cate <- ifelse(combine_dt$cate == 'NotShared', 'SpeS','notSS')
aggregate_add_prop_dt <- aggregate(prop ~ spe + celltype + new_cate, data = combine_dt, sum)

aggregate_add_prop_dt <- aggregate_add_prop_dt[aggregate_add_prop_dt$spe == 'rice_maize',]







##updating 121123
####################################################################
##we will use fisher exact test to show whether the cell type specific ACRs
##are enriched on the not shared ACR which cell types are more likely enriched on the not shared ACRs
ipt_dir <- 'opt2_enrichment_celltype_for_not_shared_group_121123/'
opt_dir <- 'opt2_enrichment_celltype_for_not_shared_group_results_121123/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt2_')


##for the d
outs <- lapply(all_fl_list, function(x){
  
  print(x)
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt2_enrichment_celltype_diff_cate_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$V1 != 'broadly_accessible',]
  ipt_dt <- ipt_dt[ipt_dt$V2 != 'broadly_accessible',]
  
  ipt_dt <- ipt_dt[c('V1','V2','V3')]
  colnames(ipt_dt) <- c('target_celltype','other_celltype','pval')
  #ipt_dt <- ipt_dt[ipt_dt$pval < 0.05,]
  ipt_dt$mlog10pval <- -log10(ipt_dt$pval)    
  ipt_dt <- ipt_dt[c('target_celltype','other_celltype','mlog10pval')]
  
  write.table(ipt_dt,paste0(opt_dir,'/temp_',target_fl_path_nm,'.txt'),quote = F,sep = '\t')
  
  ipt_dt <- read.table(paste0(opt_dir,'/temp_',target_fl_path_nm,'.txt'),stringsAsFactors = T)
  
  ipt_dt$target_celltype <- gsub('companion_cells_sieve_elements','companion_cell',ipt_dt$target_celltype)
  ipt_dt$other_celltype <- gsub('companion_cells_sieve_elements','companion_cell',ipt_dt$other_celltype)
  
  
  ipt_dt <- ipt_dt[ipt_dt$target_celltype != 'procambial_meristem',]
  ipt_dt <- ipt_dt[ipt_dt$other_celltype != 'procambial_meristem',]
  
  ipt_dt <- ipt_dt[ipt_dt$target_celltype != 'procambium',]
  ipt_dt <- ipt_dt[ipt_dt$other_celltype != 'procambium',]
  
  ipt_dt <- ipt_dt[ipt_dt$target_celltype != 'unknown',]
  ipt_dt <- ipt_dt[ipt_dt$other_celltype != 'unknown',]
  
  
  ipt_dt$target_celltype <- factor(ipt_dt$target_celltype,levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
  ipt_dt$other_celltype <- factor(ipt_dt$other_celltype,levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
  
  score_mtx <- sparseMatrix(i=as.numeric(ipt_dt$target_celltype),
                            j=as.numeric(ipt_dt$other_celltype),
                            x=as.numeric(ipt_dt$mlog10pval),
                            dimnames=list(levels(ipt_dt$target_celltype), levels(ipt_dt$other_celltype)))
  score_mtx <- as.matrix(score_mtx)
  
  pdf(paste0(opt_dir,'/opt_enrichment_celltype_on_not_shared_',target_fl_path_nm,'.pdf'),width = 5,height = 3)
  p <- pheatmap(score_mtx,
                #scale="row",
                #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
                #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
                cluster_cols = F,
                cluster_rows = F,
                border_color='white',
                fontsize_col = 15,
                fontsize_row = 15,
                fontsize_number = 25
  ) 
  print(p)
  dev.off()
  
  
  
    
})










##updating 120723
####################################################################
##we will conduct the motif enrichment for the lineage specific ACRs 
ipt_dir <- 'opt2_motif_enrichment_cate_per_celltype_120723//////'
opt_dir <- 'opt2_motif_enrichment_cate_per_celltype_results_120723////'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_')

target_cate <- 'notshare'
target_cate <- 'shareI'


##for the d
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_celltype_enrich.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  #ipt_dt$spe <- target_fl_path_nm 
  
  colnames(ipt_dt) <- c('cate','celltype','motifID','pval','num','total','prop','spe') 
  
  ipt_dt <- ipt_dt[c('cate','celltype','motifID','pval','spe')]
  return(ipt_dt)  
  
})


combine_dt <- do.call(rbind,outs)
combine_dt <- combine_dt[combine_dt$celltype != 'unknown_cells_1'&combine_dt$celltype != 'unknown_cells_2', ]

meta_name_dt <- read.delim(paste0(ipt_dir,'/ipt_required_raw_data/opt_motif_common_name.txt'),header = F)
head(meta_name_dt)

##add the common name
combine_dt <- merge(combine_dt,meta_name_dt,by.x='motifID',by.y= 'V1')
head(combine_dt)
combine_dt$qval <- p.adjust(combine_dt$pval, method="fdr")
combine_dt <- combine_dt[combine_dt$qval < 0.05,]

##lineage specific
combine_notshare_dt <- combine_dt[combine_dt$cate == target_cate,]
combine_notshare_dt$mlog10qval <- -log10(combine_notshare_dt$qval)
combine_notshare_dt$spe_pair_ct <- paste0(combine_notshare_dt$spe,'__',combine_notshare_dt$celltype)
combine_notshare_dt <- combine_notshare_dt[c('V2','spe_pair_ct','mlog10qval')]
colnames(combine_notshare_dt) <- c('motif','spe_pair_ct','mlog10qval')



write.table(combine_notshare_dt,paste0(opt_dir,'/temp_sparse.txt'),quote = F,sep = '\t')
combine_notshare_dt <- read.table(paste0(opt_dir,'/temp_sparse.txt'),stringsAsFactors = T)
score_mtx <- sparseMatrix(i=as.numeric(combine_notshare_dt$motif),
                          j=as.numeric(combine_notshare_dt$spe_pair),
                          x=as.numeric(combine_notshare_dt$mlog10qval),
                          dimnames=list(levels(combine_notshare_dt$motif), levels(combine_notshare_dt$spe_pair)))
score_mtx <- as.matrix(score_mtx)
dim(score_mtx)
head(score_mtx)
max(score_mtx)

##reorder the matrix
colnames(score_mtx)
target_order <- c('rice_maize__mesophyll','rice_sorghum__mesophyll','rice_Pm__mesophyll','rice_Uf__mesophyll',
                  'rice_maize__bundle_sheath','rice_sorghum__bundle_sheath','rice_Pm__bundle_sheath','rice_Uf__bundle_sheath',
                  'rice_maize__companion_cell','rice_sorghum__companion_cell','rice_Pm__companion_cell','rice_Uf__companion_cell',
                  'rice_maize__protoderm','rice_sorghum__protoderm','rice_Pm__protoderm','rice_Uf__protoderm',
                  'rice_maize__epidermis','rice_sorghum__epidermis','rice_Pm__epidermis','rice_Uf__epidermis'
                  )

intersect_colnm <- intersect(target_order,colnames(score_mtx))

score_mtx <- score_mtx[,intersect_colnm]

# Replace Inf with the maximum value of the matrix
max_non_inf <- max(score_mtx[!is.infinite(score_mtx)], na.rm = TRUE)
score_mtx[is.infinite(score_mtx)] <- max_non_inf
max(score_mtx)

p <- pheatmap(score_mtx,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
pdf(paste0(opt_dir,'/opt_',target_cate,'_celltype_enrich_motif.pdf'),width = 15,height = 15)
p
dev.off()

##Here we will build a plot to show 
write.table(score_mtx,paste0(opt_dir,'/opt_', target_cate,'_celltype_enrich_mtx.txt'), quote = F, sep = '\t')
score_mtx <- as(score_mtx, "sparseMatrix")
ia <- as.data.frame(summary(score_mtx))
ia$i <- rownames(score_mtx)[as.numeric(ia$i)]
ia$j <- colnames(score_mtx)[as.numeric(ia$j)]
write.table(ia, file=paste0(opt_dir,'/opt_', target_cate,'_celltype_enrich.sparse'), quote=F, row.names=F, col.names=F, sep="\t")

##read the new analysed data
ipt_motif_Sepcount_qval_fl <- 'opt2_motif_enrichment_cate_per_celltype_120723/add_downstream_analysis_for_s2_s4_120823/opt_not_shared_motif_specount_avglog10qval.txt'
ipt_motif_Sepcount_qval_dt <- read.delim(ipt_motif_Sepcount_qval_fl,header = F)
colnames(ipt_motif_Sepcount_qval_dt) <- c('motif','celltype','SpeCount','Avglog10qval')

p <- ggplot(ipt_motif_Sepcount_qval_dt, aes(Avglog10qval, SpeCount, color = celltype, fill = celltype,shape= celltype)) +
  geom_point(size = 7) +
  scale_color_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                             "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cell"="#eddc9f",
                              "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  scale_shape_manual(values = c("bundle_sheath" = 16, "companion_cell" = 25,
                                "epidermis" = 18, "mesophyll" = 15, "protoderm" = 17))+
  #geom_point(data = points_to_label, color = 'red', size = 3) +
  #geom_text(data = points_to_label, aes(label = rownames(points_to_label)), hjust = -0.2,color='red',cex= 7) + 
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text.x = element_text(colour = "black", size=30,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
    axis.text.x = element_text(size=20,colour = "black"),
    axis.text.y = element_text(size=20,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
  #scale_x_continuous(breaks=seq(-2, 2, 0.5)) +
  #scale_y_continuous(breaks=seq(-2, 2, 0.5)) +
  #coord_cartesian(xlim = c(-2, 2), ylim = c(-2,2))

pdf(paste0(opt_dir,'/opt_',target_cate,'_motif_enrich_speCount_avglog10.pdf'),width = 12,height = 6)
print(p)
dev.off()

write.table(ipt_motif_Sepcount_qval_dt,paste0(opt_dir,'/opt_',target_cate,'_motif_enrich_speCount_avglog10.txt'),sep = '\t',quote = F)




##updating 120823
##we will check other species to rice
ipt_ori_dir <- 'opt2_motif_enrichment_cate_per_celltype_120723//////'
ipt_dir <- 'opt2_motif_enrichment_cate_per_celltype_120723/for_other_species_to_rice/'
opt_dir <- 'opt2_motif_enrichment_cate_per_celltype_results_120723////'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_')

target_cate <- 'notshare'
target_cate <- 'shareI'

##for the d
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_celltype_enrich.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  #ipt_dt$spe <- target_fl_path_nm 
  
  colnames(ipt_dt) <- c('cate','celltype','motifID','pval','num','total','prop','spe') 
  
  ipt_dt <- ipt_dt[c('cate','celltype','motifID','pval','spe')]
  return(ipt_dt)  
  
})

combine_dt <- do.call(rbind,outs)
combine_dt <- combine_dt[combine_dt$celltype != 'unknown'&combine_dt$celltype != 'procambium'&combine_dt$celltype != 'procambial_meristem', ]

meta_name_dt <- read.delim(paste0(ipt_ori_dir,'/ipt_required_raw_data/opt_motif_common_name.txt'),header = F)
head(meta_name_dt)

##add the common name
combine_dt <- merge(combine_dt,meta_name_dt,by.x='motifID',by.y= 'V1')
head(combine_dt)
#combine_dt <- combine_dt[combine_dt$pval < 0.0001,]
combine_dt$qval <- p.adjust(combine_dt$pval, method="fdr")
combine_dt <- combine_dt[combine_dt$qval < 0.05,]

table(combine_dt$celltype)

##lineage specific
combine_notshare_dt <- combine_dt[combine_dt$cate == target_cate,]
combine_notshare_dt$mlog10qval <- -log10(combine_notshare_dt$qval)
combine_notshare_dt$spe_pair_ct <- paste0(combine_notshare_dt$spe,'__',combine_notshare_dt$celltype)
combine_notshare_dt <- combine_notshare_dt[c('V2','spe_pair_ct','mlog10qval')]
colnames(combine_notshare_dt) <- c('motif','spe_pair_ct','mlog10qval')

write.table(combine_notshare_dt,paste0(opt_dir,'/temp_sparse.txt'),quote = F,sep = '\t')
combine_notshare_dt <- read.table(paste0(opt_dir,'/temp_sparse.txt'),stringsAsFactors = T)
score_mtx <- sparseMatrix(i=as.numeric(combine_notshare_dt$motif),
                          j=as.numeric(combine_notshare_dt$spe_pair),
                          x=as.numeric(combine_notshare_dt$mlog10qval),
                          dimnames=list(levels(combine_notshare_dt$motif), levels(combine_notshare_dt$spe_pair)))
score_mtx <- as.matrix(score_mtx)
dim(score_mtx)
head(score_mtx)
max(score_mtx)
colnames(score_mtx)

##reorder the matrix
colnames(score_mtx)
target_order <- c('maize_rice__mesophyll','sorghum_rice__mesophyll','Pm_rice__mesophyll','Uf_rice__mesophyll',
                  'maize_rice__bundle_sheath','sorghum_rice__bundle_sheath','Pm_rice__bundle_sheath','Uf_rice__bundle_sheath',
                  'maize_rice__companion_cells_sieve_elements','sorghum_rice__companion_cells_sieve_elements','Pm_rice__companion_cells_sieve_elements','Uf_rice__companion_cells_sieve_elements',
                  'maize_rice__protoderm','sorghum_rice__protoderm','Pm_rice__protoderm','Uf_rice__protoderm',
                  'maize_rice__epidermis','sorghum_rice__epidermis','Pm_rice__epidermis','Uf_rice__epidermis'
)

intersect_colnm <- intersect(target_order,colnames(score_mtx))

score_mtx <- score_mtx[,intersect_colnm]

# Replace Inf with the maximum value of the matrix
max_non_inf <- max(score_mtx[!is.infinite(score_mtx)], na.rm = TRUE)
score_mtx[is.infinite(score_mtx)] <- max_non_inf
max(score_mtx)

p <- pheatmap(score_mtx,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
pdf(paste0(opt_dir,'/opt_',target_cate,'_celltype_enrich_motif_others_to_rice.pdf'),width = 15,height = 35)
p
dev.off()

##Here we will build a plot to show 
write.table(score_mtx,paste0(opt_dir,'/opt_', target_cate,'_celltype_enrich_mtx_otherspe_to_rice.txt'), quote = F, sep = '\t')
score_mtx <- as(score_mtx, "sparseMatrix")
ia <- as.data.frame(summary(score_mtx))
ia$i <- rownames(score_mtx)[as.numeric(ia$i)]
ia$j <- colnames(score_mtx)[as.numeric(ia$j)]
write.table(ia, file=paste0(opt_dir,'/opt_', target_cate,'_celltype_enrich_otherspe_to_rice.sparse'), quote=F, row.names=F, col.names=F, sep="\t")


##read the new analysed data
ipt_motif_Sepcount_qval_fl <- 'opt2_motif_enrichment_cate_per_celltype_120723/add_downstream_analysis_for_s2_s4_120823/opt_not_shared_others_to_rice_motif_specount_avglog10qval.txt'
ipt_motif_Sepcount_qval_dt <- read.delim(ipt_motif_Sepcount_qval_fl,header = F)
colnames(ipt_motif_Sepcount_qval_dt) <- c('motif','celltype','SpeCount','Avglog10qval')

p <- ggplot(ipt_motif_Sepcount_qval_dt, aes(Avglog10qval, SpeCount, color = celltype, fill = celltype,shape= celltype)) +
  geom_point(size = 7) +
  scale_color_manual(values=c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
                              "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
                             "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
  scale_shape_manual(values = c("bundle_sheath" = 16, "companion_cells_sieve_elements" = 25,
                                "epidermis" = 18, "mesophyll" = 15, "protoderm" = 17))+
  #geom_point(data = points_to_label, color = 'red', size = 3) +
  #geom_text(data = points_to_label, aes(label = rownames(points_to_label)), hjust = -0.2,color='red',cex= 7) + 
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text.x = element_text(colour = "black", size=30,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
    axis.text.x = element_text(size=20,colour = "black"),
    axis.text.y = element_text(size=20,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
#scale_x_continuous(breaks=seq(-2, 2, 0.5)) +
#scale_y_continuous(breaks=seq(-2, 2, 0.5)) +
#coord_cartesian(xlim = c(-2, 2), ylim = c(-2,2))

pdf(paste0(opt_dir,'/opt_',target_cate,'_motif_enrich_speCount_avglog10_othersTorice.pdf'),width = 12,height = 6)
print(p)
dev.off()

write.table(ipt_motif_Sepcount_qval_dt,paste0(opt_dir,'/opt_',target_cate,'_motif_enrich_speCount_avglog10_othersTorice.txt'),sep = '\t',quote = F)






##updating 120723
######################################
##we will check the cns number per ACR
ipt_dir <- 'opt2_check_cns_num_per_acr_120723/Pablo_version/'
opt_dir <- 'opt2_check_cns_num_per_acr_results_120723/Pablo_version/'

##choose the database version
ipt_dir <- 'opt2_check_cns_num_per_acr_120723/data_base_version/not_contain_UTR_version_choose/'
opt_dir <- 'opt2_check_cns_num_per_acr_results_120723/database_version/not_contain_UTR_version_choose/'


##check the others to rice which would be the Pablo's ACRs 
ipt_dir <- 'opt2_check_cns_num_per_acr_others_to_rice_121823/Pablo_CNSs_version/'
opt_dir <- 'opt2_check_cns_num_per_acr_results_others_to_rice_121823/Pablo_CNSs_version/'

ipt_dir <- 'opt2_check_cns_num_per_acr_others_to_rice_121823/data_base_version/contain_UTR_version//'
opt_dir <- 'opt2_check_cns_num_per_acr_results_others_to_rice_121823/data_base_version//'




all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_final_cate')

##for the d
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_cate_acr_cns_num_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)  
  table(ipt_dt$V1)
  ipt_dt_three_cate_dt <- ipt_dt[ipt_dt$V1 == 'shared_A_Broad' | ipt_dt$V1 == 'shared_A_CT' |
                                 ipt_dt$V1 == 'shared_I_Broad' | ipt_dt$V1 == 'shared_I_CT' | 
                                 ipt_dt$V1 == 'not_shared_Broad' | ipt_dt$V1 == 'not_shared_CT',]
  
  ipt_dt_three_cate_dt$cate <- ifelse(ipt_dt_three_cate_dt$V3 == 0,'notCoverCNS','CoverCNS')
  ipt_dt_three_cate_dt$count <- 1
  
  table(ipt_dt_three_cate_dt$V1)
  
  aggregate_dt_three_cate <- aggregate(count ~ V1 + cate, data = ipt_dt_three_cate_dt, sum)
  aggregate_dt_sum_cate <- aggregate(count ~ V1 ,data = ipt_dt_three_cate_dt,sum)
  merge_three_cate = merge(aggregate_dt_three_cate,aggregate_dt_sum_cate,by.x = 'V1',by.y = 'V1')
  merge_three_cate$prop <- merge_three_cate$count.x/merge_three_cate$count.y
  merge_three_cate$V1 <- factor(merge_three_cate$V1,levels = c('shared_A_Broad','shared_A_CT',
                                                               'shared_I_Broad','shared_I_CT',
                                                               'not_shared_Broad','not_shared_CT'))
  

  
  p <- ggplot(data=merge_three_cate, aes(x=V1, y=prop, fill=cate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )

  pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_cns_prop','.pdf'),width = 8,height = 8 )
  print(p)
  dev.off()
  
  merge_three_cate$spe <- target_fl_path_nm
  
  return(merge_three_cate)
  
  
})

combine_dt <- do.call(rbind,outs)

combine_coverCNS_dt <- combine_dt[combine_dt$cate == 'CoverCNS',]
head(combine_coverCNS_dt)
colnames(combine_coverCNS_dt) <- c('syn_cate','cate','countx','county','prop','spe')


write.csv(combine_coverCNS_dt,paste0(opt_dir,'/opt_threeclass_acr_coverCNS_num.csv'),quote = F)

p <- ggplot(combine_coverCNS_dt, aes(x=syn_cate, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=syn_cate, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=1, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
p <- p + coord_cartesian(ylim = c(0, 0.5)) + ##threshold
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
p

pdf(paste0(opt_dir,'/opt2_rice_to_allspe_syntenic_cate_coverCNS_prop.pdf'),width = 10,height = 10)
p
dev.off()
##other option for allspe to rice
pdf(paste0(opt_dir,'/opt2_allspe_to_rice_syntenic_cate_coverCNS_prop.pdf'),width = 10,height = 10)
p
dev.off()


##check the significant test
shared_acr_broad = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'shared_A_Broad',]$prop
shared_acr_CT = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'shared_A_CT',]$prop
variable_acr_broad = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'shared_I_Broad',]$prop
variable_acr_CT = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'shared_I_CT',]$prop
not_shared_acr_broad = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'not_shared_Broad',]$prop
not_shared_acr_CT = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'not_shared_CT',]$prop

t.test(shared_acr_CT,shared_acr_broad,alternative = 'greater')
t.test(variable_acr_CT,variable_acr_broad,alternative = 'greater')
t.test(not_shared_acr_CT,not_shared_acr_broad,alternative = 'greater')


##check the significant test for the others to rice
not_shared_broad_merged_all_dt <- combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'not_shared_Broad',]
not_shared_CT_merged_all_dt <- combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'not_shared_CT',]
t.test(not_shared_CT_merged_all_dt$prop,not_shared_broad_merged_all_dt$prop, paired = T)


##updating 121923 check the proprotion of CT across all three categories compared to broad across all three categories
generate_pie <- function(dt,target_celltype) {
  df <- dt[dt$celltype == target_celltype,]
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df$cate,pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}


myPalette <- brewer.pal(10, "Set3") 

spe_list <- unique(combine_coverCNS_dt$spe)
for (i in (1:length(spe_list))){
  
  combine_coverCNS_dt_spe <- combine_coverCNS_dt[combine_coverCNS_dt$spe == spe_list[i],]
  combine_coverCNS_dt_spe$BroadCTcate <- ifelse(grepl('_Broad',combine_coverCNS_dt_spe$syn_cate),'Broad','CT')
  
  aggregate_syn_cate <- aggregate(countx ~ BroadCTcate,data = combine_coverCNS_dt_spe, sum)
  merged_syn_cate <- merge(combine_coverCNS_dt_spe,aggregate_syn_cate,by.x = 'BroadCTcate',by.y = 'BroadCTcate')
  
  ##for the broad
  merged_syn_cate_broad <- merged_syn_cate[merged_syn_cate$BroadCTcate == 'Broad',]
  df <- merged_syn_cate_broad
  pct <- round(df$countx.x/sum(df$countx.x)*100)
  lbls <- paste(df$syn_cate,pct)
  lbls <- paste(lbls,"%",sep="") 
  pdf(paste0(opt_dir,'/opt_piechart_',spe_list[i],'_broad.pdf'), width = 10,height = 10)
  print(pie(df$countx.x,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1))
  dev.off()
  
  ##for the cell-type-specific
  merged_syn_cate_CT <- merged_syn_cate[merged_syn_cate$BroadCTcate == 'CT',]
  df <- merged_syn_cate_CT
  pct <- round(df$countx.x/sum(df$countx.x)*100)
  lbls <- paste(df$syn_cate,pct)
  lbls <- paste(lbls,"%",sep="") 
  pdf(paste0(opt_dir,'/opt_piechart_',spe_list[i],'_CT.pdf'), width = 10,height = 10)
  print(pie(df$countx.x,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1))
  dev.off()
  
}

##updating 020224
##we will plot the total shared and not shared
spe_list <- unique(combine_coverCNS_dt$spe)

generate_pie <- function(dt) {
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(dt$count/sum(dt$count)*100)
  lbls <- paste(dt$cate,pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(dt$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}


for (i in (1:length(spe_list))){
 
  combine_coverCNS_dt_spe <- combine_coverCNS_dt[combine_coverCNS_dt$spe == spe_list[i],]
  
  not_shared_Broad <- combine_coverCNS_dt_spe[combine_coverCNS_dt_spe$syn_cate == 'not_shared_Broad',]$countx
  not_shared_CT <- combine_coverCNS_dt_spe[combine_coverCNS_dt_spe$syn_cate == 'not_shared_CT',]$countx
  not_shared_coverCNS_num <- not_shared_Broad + not_shared_CT
  
  shared_A_Broad <- combine_coverCNS_dt_spe[combine_coverCNS_dt_spe$syn_cate == 'shared_A_Broad',]$countx
  shared_A_CT <- combine_coverCNS_dt_spe[combine_coverCNS_dt_spe$syn_cate == 'shared_A_CT',]$countx
  shared_A_coverCNS_num <- shared_A_Broad + shared_A_CT
  
  shared_I_Broad <- combine_coverCNS_dt_spe[combine_coverCNS_dt_spe$syn_cate == 'shared_I_Broad',]$countx
  shared_I_CT <- combine_coverCNS_dt_spe[combine_coverCNS_dt_spe$syn_cate == 'shared_I_CT',]$countx
  shared_I_coverCNS_num <- shared_I_Broad + shared_I_CT
  
  dt <- data.frame(
    cate <- c('not_shared_coverCNS_num','shared_A_coverCNS_num','shared_I_coverCNS_num'),
    number <- c(not_shared_coverCNS_num,shared_A_coverCNS_num,shared_I_coverCNS_num)
  )
  
  colnames(dt) <- c('cate','count')
  
  pdf(paste0(opt_dir,'/opt_piechart_',spe_list[i],'_total.pdf'), width = 10,height = 10)
  print(generate_pie(dt))
  dev.off()
  
}







##updating 121823 plot the different cate per cell type
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_cate_acr_cns_num_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  spe_nm <- gsub('rice_to_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice', '',spe_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)  
  table(ipt_dt$V1)
  ipt_dt_three_cate_dt <- ipt_dt[ipt_dt$V1 != 'shared_A_Broad' & ipt_dt$V1 != 'shared_A_CT' &
                                   ipt_dt$V1 != 'shared_I_Broad' & ipt_dt$V1 != 'shared_I_CT' & 
                                   ipt_dt$V1 != 'not_shared_Broad' & ipt_dt$V1 != 'not_shared_CT',]
  table(ipt_dt_three_cate_dt$V1)
  
  ##Here we will modify a little bit for the name of ipt_dt_three_cate_dt
  ipt_dt_three_cate_dt$V1 <- gsub('companion_cells_sieve_elements','companion_cell',ipt_dt_three_cate_dt$V1)
  ipt_dt_three_cate_dt <- ipt_dt_three_cate_dt[!grepl('procambial_meristem', ipt_dt_three_cate_dt$V1),]
  ipt_dt_three_cate_dt <- ipt_dt_three_cate_dt[!grepl('procambium', ipt_dt_three_cate_dt$V1),]
  
  ipt_dt_three_cate_dt$cate <- ifelse(ipt_dt_three_cate_dt$V3 == 0,'notCoverCNS','CoverCNS')
  ipt_dt_three_cate_dt$count <- 1
  aggregate_dt_three_cate <- aggregate(count ~ V1 + cate, data = ipt_dt_three_cate_dt, sum)
  aggregate_dt_sum_cate <- aggregate(count ~ V1 ,data = ipt_dt_three_cate_dt,sum)
  merge_three_cate = merge(aggregate_dt_three_cate,aggregate_dt_sum_cate,by.x = 'V1',by.y = 'V1')
  merge_three_cate$prop <- merge_three_cate$count.x/merge_three_cate$count.y
  
  three_cate <- c('shared_A','shared_I','not_shared')
  outs2 <- lapply(three_cate, function(x){
    
    print(x)
    merge_three_cate_target <- merge_three_cate[grepl(x,merge_three_cate$V1),]
    merge_three_cate_target$syn_cate <- x
    
    merge_three_cate_target <- merge_three_cate_target[!grepl('unknown',merge_three_cate_target$V1),]
    merge_three_cate_target$celltype <- gsub('shared_A_','',merge_three_cate_target$V1)
    merge_three_cate_target$celltype <- gsub('shared_I_','',merge_three_cate_target$celltype)
    merge_three_cate_target$celltype <- gsub('not_shared_','',merge_three_cate_target$celltype)
    celltype_levels <- c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis')
    merge_three_cate_target$celltype <- factor(merge_three_cate_target$celltype, levels = celltype_levels)
    
    p <- ggplot(data=merge_three_cate_target, aes(x=celltype, y=prop, fill=cate)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme(
        plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
        axis.title = element_text(size =20, face="bold"),
        #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
        axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
        axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
        axis.ticks = element_line(size = rel(2.5)),
        axis.ticks.length = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 15),
        strip.text = element_text(size=20),
        #legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
      )
    #pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_cns_prop_celltype_',x,'.pdf'),width = 15,height = 8 )
    #print(p)
    #dev.off()
    
    
    return(merge_three_cate_target) 
    
  })
  
  merged_dt <- do.call(rbind,outs2)
  
  #merge_three_cate$V1 <- factor(merge_three_cate$V1,levels = c('shared_A_Broad','shared_A_CT',
  #                                                             'shared_I_Broad','shared_I_CT',
  #                                                             'not_shared_Broad','not_shared_CT'))
  
  
  merged_dt$syn_cate <- factor(merged_dt$syn_cate,levels =  c('shared_A','shared_I','not_shared'))
  
  class(merged_dt)
  p <- ggplot(data=merged_dt, aes(x=celltype, y=prop, fill=cate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) + 
    facet_wrap(~syn_cate)
  
  pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_cns_prop_celltype_all_syn_cate','.pdf'),width = 15,height = 8 )
  print(p)
  dev.off()
  
  merged_dt$spe <- spe_nm
  
  return(merged_dt)  
  
})


merged_all_dt <- do.call(rbind,outs)
##conduct the t test to compare not shared overlapping with the CNS for all species
not_shared_epidermis_merged_all_dt <- merged_all_dt[merged_all_dt$V1 == 'not_shared_epidermis'&merged_all_dt$cate =='notCoverCNS',]
not_shared_mesophyll_merged_all_dt <- merged_all_dt[merged_all_dt$V1 == 'not_shared_mesophyll'&merged_all_dt$cate =='notCoverCNS',]
not_shared_bundle_sheath_merged_all_dt <- merged_all_dt[merged_all_dt$V1 == 'not_shared_bundle_sheath'&merged_all_dt$cate =='notCoverCNS',]
not_shared_companion_cell_all_dt <- merged_all_dt[merged_all_dt$V1 == 'not_shared_companion_cell'&merged_all_dt$cate =='notCoverCNS',]
not_shared_protoderm_all_dt <- merged_all_dt[merged_all_dt$V1 == 'not_shared_protoderm'&merged_all_dt$cate =='notCoverCNS',]
t.test(not_shared_epidermis_merged_all_dt$prop,not_shared_mesophyll_merged_all_dt$prop, paired = T, alternative = 'greater')
t.test(not_shared_epidermis_merged_all_dt$prop,not_shared_bundle_sheath_merged_all_dt$prop, paired = T, alternative = 'greater')
t.test(not_shared_epidermis_merged_all_dt$prop,not_shared_companion_cell_all_dt$prop, paired = T, alternative = 'greater')
t.test(not_shared_epidermis_merged_all_dt$prop,not_shared_protoderm_all_dt$prop, paired = T, alternative = 'greater')


##use the box plot to plot all proportion for all the species 
merged_all_not_shared_dt <- merged_all_dt[grepl('not_shared_',merged_all_dt$V1),]
merged_all_not_shared_dt <- merged_all_not_shared_dt[merged_all_not_shared_dt$cate == 'CoverCNS',]
p <- ggplot(merged_all_not_shared_dt, aes(x=celltype, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=celltype, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=0.75, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 

pdf(paste0(opt_dir,'/opt2_rice_to_allspe_cover_cns_prop_for_not_shared.pdf'),width = 10,height = 14)
p
dev.off()


merged_all_shared_A_dt <- merged_all_dt[grepl('shared_A',merged_all_dt$V1),]
merged_all_shared_A_dt <- merged_all_shared_A_dt[merged_all_shared_A_dt$cate == 'CoverCNS',]
p <- ggplot(merged_all_shared_A_dt, aes(x=celltype, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=celltype, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=0.75, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 

pdf(paste0(opt_dir,'/opt2_rice_to_allspe_cover_cns_prop_for_shared_A.pdf'),width = 10,height = 14)
p
dev.off()


merged_all_shared_I_dt <- merged_all_dt[grepl('shared_I',merged_all_dt$V1),]
merged_all_shared_I_dt <- merged_all_shared_I_dt[merged_all_shared_I_dt$cate == 'CoverCNS',]
p <- ggplot(merged_all_shared_I_dt, aes(x=celltype, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=celltype, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=0.75, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 

pdf(paste0(opt_dir,'/opt2_rice_to_allspe_cover_cns_prop_for_shared_I.pdf'),width = 10,height = 14)
p
dev.off()






##check how many species-specific ACRs are shared across all four species pairs
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_cate_acr_cns_num_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  spe_nm <- gsub('rice_to_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  
  #ipt_dt_epider_mis <- ipt_dt[ipt_dt$V1 == 'not_shared_epidermis',]
  ipt_dt_epider_mis <- ipt_dt[ipt_dt$V1 == 'not_shared_protoderm',]
  
  ipt_dt_epider_mis$spe <- spe_nm
  
  
  ipt_dt_epider_mis$cate <- ifelse(ipt_dt_epider_mis$V3 == 0,'notCoverCNS','CoverCNS')
  ipt_dt_epider_mis$count <- 1
  aggregate_dt_three_cate <- aggregate(count ~ V1 + cate, data = ipt_dt_epider_mis, sum)
  aggregate_dt_sum_cate <- aggregate(count ~ V1 ,data = ipt_dt_epider_mis,sum)
  
  merge_three_cate = merge(aggregate_dt_three_cate,aggregate_dt_sum_cate,by.x = 'V1',by.y = 'V1')
  merge_three_cate$prop <- merge_three_cate$count.x/merge_three_cate$count.y
  
  
  return(ipt_dt_epider_mis)
  
})

combine_dt <- do.call(rbind,outs)
combine_dt_maize <- combine_dt[combine_dt$spe == 'maize',]
combine_dt_sorghum <- combine_dt[combine_dt$spe == 'sorghum',]
combine_dt_Pm <- combine_dt[combine_dt$spe == 'Pm',]
combine_dt_Uf <- combine_dt[combine_dt$spe == 'Uf',]

maize_sorghum_inter <- intersect(combine_dt_maize$V2,combine_dt_sorghum$V2)
maize_sorghum_Uf_inter <- intersect(maize_sorghum_inter,combine_dt_Uf$V2)
maize_sorghum_Pm_Uf_inter <- intersect(combine_dt_Pm$V2,maize_sorghum_Uf_inter)

combine_unique <- unique(combine_dt[c('V1','V2','V3')])
combine_dt_all_have_epidermis <- combine_unique[combine_unique$V2 %in% maize_sorghum_Pm_Uf_inter,]
table(combine_dt_all_have_epidermis$V3)







##check the cell type enrichment compare to each other
all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_plot_R')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_plot_R_prop_allcelltype_enrich_on_CNS_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)

  ipt_dt <- read.delim(target_fl_path,header = F)
  colnames(ipt_dt) <- c('cate','celltype','groupCate','prop')
  ipt_dt$celltype <- gsub('shared_I_','',ipt_dt$celltype)
  ipt_dt$celltype <- gsub('not_shared_','',ipt_dt$celltype)
  ipt_dt$celltype <- gsub('shared_A_','',ipt_dt$celltype)
  
  ipt_dt$celltype <- factor(ipt_dt$celltype, levels = c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis'))
  ipt_dt$groupCate <- factor(ipt_dt$groupCate, levels = c('TargetCT','NotTargetCT'))
  ipt_dt$cate <- factor(ipt_dt$cate, levels = c('shared_A','shared_I','not_shared'))
  
  p <- ggplot(data=ipt_dt, aes(x=celltype , y=prop, fill=groupCate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    facet_wrap(~cate) ## You can use face_wrap function only if you need it+
  pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_cns_prop_TargetCellTP_versus_otherCellTP','.pdf'),width = 12,height = 8 )
  print(p)
  dev.off()

})



##upadting 010524
##check syn and nonsyn for the CNS
ipt_dir <- 'opt2_check_cns_num_per_acr_120723/data_base_version/not_contain_UTR_version_choose/'
opt_dir <- 'opt2_check_cns_num_per_acr_results_120723/database_version/not_contain_UTR_version_choose/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_final_cate')

##for the d
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_cate_acr_cns_num_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)  
  table(ipt_dt$V1)
  ipt_dt_three_cate_dt <- ipt_dt[ipt_dt$V1 == 'nonsyntenic_Broad' | ipt_dt$V1 == 'nonsyntenic_CT' |
                                   ipt_dt$V1 == 'syntenic_Broad' | ipt_dt$V1 == 'syntenic_CT',]
  
  ipt_dt_three_cate_dt$cate <- ifelse(ipt_dt_three_cate_dt$V3 == 0,'notCoverCNS','CoverCNS')
  ipt_dt_three_cate_dt$count <- 1
  aggregate_dt_three_cate <- aggregate(count ~ V1 + cate, data = ipt_dt_three_cate_dt, sum)
  aggregate_dt_sum_cate <- aggregate(count ~ V1 ,data = ipt_dt_three_cate_dt,sum)
  merge_three_cate = merge(aggregate_dt_three_cate,aggregate_dt_sum_cate,by.x = 'V1',by.y = 'V1')
  merge_three_cate$prop <- merge_three_cate$count.x/merge_three_cate$count.y
  merge_three_cate$V1 <- factor(merge_three_cate$V1,levels = c('syntenic_Broad','syntenic_CT',
                                                               'nonsyntenic_Broad','nonsyntenic_CT'))
  
  
  
  p <- ggplot(data=merge_three_cate, aes(x=V1, y=prop, fill=cate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  
  pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_synOrnonsyn_cns_prop','.pdf'),width = 8,height = 8 )
  print(p)
  dev.off()
  
  merge_three_cate$spe <- target_fl_path_nm
  
  return(merge_three_cate)
  
  
})

combine_dt <- do.call(rbind,outs)

combine_coverCNS_dt <- combine_dt[combine_dt$cate == 'CoverCNS',]
head(combine_coverCNS_dt)
colnames(combine_coverCNS_dt) <- c('syn_cate','cate','countx','county','prop','spe')

p <- ggplot(combine_coverCNS_dt, aes(x=syn_cate, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=syn_cate, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=1.5, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
  p <- p + coord_cartesian(ylim = c(0, 0.5)) + ##threshold
  scale_y_continuous(breaks=seq(0, 0.5, 0.1))
p

pdf(paste0(opt_dir,'/opt2_rice_to_allspe_nonsynOrsyn_cate_coverCNS_prop.pdf'),width = 8,height = 10)
p
dev.off()
##other option for allspe to rice
pdf(paste0(opt_dir,'/opt2_allspe_to_rice_nonsynOrsyn_cate_coverCNS_prop.pdf'),width = 8,height = 10)
p
dev.off()


##check the significant test
syntenic_broad = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'syntenic_Broad',]$prop
syntenic_CT = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'syntenic_CT',]$prop
nonsyntenic_broad = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'nonsyntenic_Broad',]$prop
nonsyntenic_CT = combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'nonsyntenic_CT',]$prop

t.test(syntenic_CT,syntenic_broad,alternative = 'greater') ##0.0007467
t.test(nonsyntenic_CT,nonsyntenic_broad,alternative = 'greater') ##0.5334



##check the significant test for the others to rice
not_shared_broad_merged_all_dt <- combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'not_shared_Broad',]
not_shared_CT_merged_all_dt <- combine_coverCNS_dt[combine_coverCNS_dt$syn_cate == 'not_shared_CT',]
t.test(not_shared_CT_merged_all_dt$prop,not_shared_broad_merged_all_dt$prop, paired = T)





##updating 121823 plot the different cate per cell type
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_cate_acr_cns_num_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  spe_nm <- gsub('rice_to_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice', '',spe_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)  
  table(ipt_dt$V1)
  ipt_dt_three_cate_dt <- ipt_dt[ipt_dt$V1 != 'nonsyntenic_Broad' & ipt_dt$V1 != 'nonsyntenic_CT' &
                                   ipt_dt$V1 != 'syntenic_Broad' & ipt_dt$V1 != 'syntenic_CT' &
                                   !grepl('shared',ipt_dt$V1),]
  

  table(ipt_dt_three_cate_dt$V1)
  
  ##Here we will modify a little bit for the name of ipt_dt_three_cate_dt
  ipt_dt_three_cate_dt$V1 <- gsub('companion_cells_sieve_elements','companion_cell',ipt_dt_three_cate_dt$V1)
  ipt_dt_three_cate_dt <- ipt_dt_three_cate_dt[!grepl('procambial_meristem', ipt_dt_three_cate_dt$V1),]
  ipt_dt_three_cate_dt <- ipt_dt_three_cate_dt[!grepl('procambium', ipt_dt_three_cate_dt$V1),]
  
  ipt_dt_three_cate_dt$cate <- ifelse(ipt_dt_three_cate_dt$V3 == 0,'notCoverCNS','CoverCNS')
  ipt_dt_three_cate_dt$count <- 1
  aggregate_dt_three_cate <- aggregate(count ~ V1 + cate, data = ipt_dt_three_cate_dt, sum)
  aggregate_dt_sum_cate <- aggregate(count ~ V1 ,data = ipt_dt_three_cate_dt,sum)
  merge_three_cate = merge(aggregate_dt_three_cate,aggregate_dt_sum_cate,by.x = 'V1',by.y = 'V1')
  merge_three_cate$prop <- merge_three_cate$count.x/merge_three_cate$count.y
  
  
  merge_three_cate$newcate <- gsub('nonsyntenic','nonyntenic',merge_three_cate$V1)
  
  three_cate <- c('nonyntenic','syntenic')
  outs2 <- lapply(three_cate, function(y){
    
    print(y)
    merge_three_cate_target <- merge_three_cate[grepl(y,merge_three_cate$newcate),]
    merge_three_cate_target$syn_cate <- y
    
    merge_three_cate_target <- merge_three_cate_target[!grepl('unknown',merge_three_cate_target$newcate),]
    merge_three_cate_target$celltype <- gsub('nonyntenic_','',merge_three_cate_target$newcate)
    merge_three_cate_target$celltype <- gsub('syntenic_','',merge_three_cate_target$celltype)
    celltype_levels <- c('mesophyll','bundle_sheath','companion_cell','protoderm','epidermis')
    merge_three_cate_target$celltype <- factor(merge_three_cate_target$celltype, levels = celltype_levels)
    
    p <- ggplot(data=merge_three_cate_target, aes(x=celltype, y=prop, fill=cate)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme(
        plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
        axis.title = element_text(size =20, face="bold"),
        #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
        axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
        axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
        axis.ticks = element_line(size = rel(2.5)),
        axis.ticks.length = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 15),
        strip.text = element_text(size=20),
        #legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
      )
    #pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_cns_prop_celltype_',x,'.pdf'),width = 15,height = 8 )
    #print(p)
    #dev.off()
    
    
    return(merge_three_cate_target) 
    
  })
  
  merged_dt <- do.call(rbind,outs2)
  
  #merge_three_cate$V1 <- factor(merge_three_cate$V1,levels = c('shared_A_Broad','shared_A_CT',
  #                                                             'shared_I_Broad','shared_I_CT',
  #                                                             'not_shared_Broad','not_shared_CT'))
  
  
  merged_dt$syn_cate <- factor(merged_dt$syn_cate,levels =  c('syntenic','nonyntenic'))
  
  class(merged_dt)
  p <- ggplot(data=merged_dt, aes(x=celltype, y=prop, fill=cate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) + 
    facet_wrap(~syn_cate)
  
  pdf(paste0(opt_dir,'/','opt2_',target_fl_path_nm,'_cns_prop_celltype_all_synOrnonsyn_cate','.pdf'),width = 15,height = 8 )
  print(p)
  dev.off()
  
  merged_dt$spe <- spe_nm
  
  return(merged_dt)  
  
})


merged_all_dt <- do.call(rbind,outs)
##conduct the t test to compare not shared overlapping with the CNS for all species
nonsyntenic_epidermis_merged_all_dt <- merged_all_dt[merged_all_dt$newcate == 'nonyntenic_epidermis'&merged_all_dt$cate =='CoverCNS',]
nonsyntenic_mesophyll_merged_all_dt <- merged_all_dt[merged_all_dt$newcate == 'nonyntenic_mesophyll'&merged_all_dt$cate =='CoverCNS',]
nonsyntenic_bundle_sheath_merged_all_dt <- merged_all_dt[merged_all_dt$newcate == 'nonyntenic_bundle_sheath'&merged_all_dt$cate =='CoverCNS',]
nonsyntenic_companion_cell_all_dt <- merged_all_dt[merged_all_dt$newcate == 'nonyntenic_companion_cell'&merged_all_dt$cate =='CoverCNS',]
nonsyntenic_protoderm_all_dt <- merged_all_dt[merged_all_dt$newcate == 'nonyntenic_protoderm'&merged_all_dt$cate =='CoverCNS',]
t.test(nonsyntenic_epidermis_merged_all_dt$prop,nonsyntenic_mesophyll_merged_all_dt$prop, paired = T, alternative = 'less')
t.test(nonsyntenic_epidermis_merged_all_dt$prop,nonsyntenic_bundle_sheath_merged_all_dt$prop, paired = T, alternative = 'less')
t.test(nonsyntenic_epidermis_merged_all_dt$prop,nonsyntenic_companion_cell_all_dt$prop, paired = T, alternative = 'less')
t.test(nonsyntenic_epidermis_merged_all_dt$prop,nonsyntenic_protoderm_all_dt$prop, paired = T, alternative = 'less')

t.test(nonsyntenic_protoderm_all_dt$prop,nonsyntenic_mesophyll_merged_all_dt$prop, paired = T, alternative = 'less')
t.test(nonsyntenic_protoderm_all_dt$prop,nonsyntenic_bundle_sheath_merged_all_dt$prop, paired = T, alternative = 'less')
t.test(nonsyntenic_protoderm_all_dt$prop,nonsyntenic_companion_cell_all_dt$prop, paired = T, alternative = 'less')
t.test(nonsyntenic_protoderm_all_dt$prop,nonsyntenic_epidermis_merged_all_dt$prop, paired = T, alternative = 'less')




##use the box plot to plot all proportion for all the species 
merged_all_not_shared_dt <- merged_all_dt[grepl('nonyntenic_',merged_all_dt$newcate),]
merged_all_not_shared_dt <- merged_all_not_shared_dt[merged_all_not_shared_dt$cate == 'CoverCNS',]
p <- ggplot(merged_all_not_shared_dt, aes(x=celltype, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=celltype, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=0.75, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 

pdf(paste0(opt_dir,'/opt2_rice_to_allspe_cover_cns_prop_for_synOrnonsyn.pdf'),width = 10,height = 14)
p
dev.off()









##updating 121923
##we will check how many variable ACRs covering CNSs in other species syntenic to rice altas ACRs
##not useful
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923///'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_')

build_dt <- function(ipt_dt,spe_nm){
  
  ##select the variable ACR containing the CNS
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  
  #ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'no',]
  
  
  #ipt_dt <- ipt_dt[ipt_dt[[paste0('rice','_region_capture_CNS')]] == 'yes',]
  
  maize_acr_num <- length(unique(ipt_dt[[paste0(spe_nm,'_ACR')]]))
  
  ##how many of them could be detect in having rice ACR
  head(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$Atlas_ACR != 'none',]
  maize_acr_capture_atlas_acr_num <-  length(unique(ipt_dt[[paste0(spe_nm,'_ACR')]]))
  
  ##make a dataframe
  dt <- data.frame(
    spe = c(spe_nm,spe_nm),
    cate = c('VarACRCNS_atlasACR','VarACRCNS_notatlasACR'),
    number = c(maize_acr_capture_atlas_acr_num,maize_acr_num-maize_acr_capture_atlas_acr_num),
    prop = c(maize_acr_capture_atlas_acr_num/maize_acr_num, (maize_acr_num-maize_acr_capture_atlas_acr_num)/maize_acr_num)
  )
  
  return(dt)
  
}

build_organ_num_dt <- function(ipt_dt,spe_nm){
  
  ##select the variable ACR containing the CNS
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  
  ##how many of them could be detect in having rice ACR
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$Atlas_ACR != 'none',]
  #ipt_dt <- 
  
}


outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_summary_add_rice_atlas_acr.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path, header = T)
  head(ipt_dt)
  
  ipt_CT_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_CelltypeCate')]] != 'broadly_accessible',]
  CT_dt <- build_dt(ipt_CT_dt,spe_nm)
  CT_dt$celltype <- 'CT'
  
  ipt_Broad_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_CelltypeCate')]] == 'broadly_accessible',]
  Broad_dt <- build_dt(ipt_Broad_dt,spe_nm)
  Broad_dt$celltype <- 'Broad'
  
  combine_dt <- rbind(CT_dt,Broad_dt)
  
  return(combine_dt)

})

combine_dt <- do.call(rbind,outs)

##do the composition plot
combine_dt$spe <- factor(combine_dt$spe,levels = c('maize','sorghum','Pm','Uf'))


p <- ggplot(combine_dt, aes(fill=cate, y=prop, x=spe)) + 
  geom_bar(position="stack", stat="identity")
#p <- p + geom_text(aes(label = prop),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
p <- p+labs(x="\nCell type", y = "Proportion\n", cex.lab = 7) ## use \n to set the position
#p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 30), ##change all the text size
             axis.text.y = element_text(size=40,colour = "black"),
             axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
             axis.title.x = element_text(size=40,colour = "black"),
             axis.title.y = element_text(size=40,colour = "black"),
             plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
p <- p + facet_wrap(~celltype)

p <- p + coord_cartesian(ylim = c(0, 1)) + ##threshold
  scale_y_continuous(breaks=seq(0, 1, 0.1))
pdf(paste0(opt_dir,'/opt_variable_acr_cns_capture_atlas_rice_acr.pdf'),width = 16,height = 10)
p
dev.off()

##build the boxplot barplot
combine_capture_atlasACR_dt <- combine_dt[combine_dt$cate == 'VarACRCNS_atlasACR',]
p <- ggplot(combine_capture_atlasACR_dt, aes(x=celltype, y=prop)) +
  #geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width=0.5,fatten = 2,position=position_dodge(0.75)) +
  geom_dotplot(aes(x=celltype, y=prop, fill = spe, color = spe), binaxis='y', stackdir='center', dotsize=3, position=position_dodge(0.75),alpha=1,
               color = 'white')+
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) 
p <- p + coord_cartesian(ylim = c(0, 1)) + ##threshold
  scale_y_continuous(breaks=seq(0, 1, 0.1))
pdf(paste0(opt_dir,'/opt_variable_acr_cns_capture_atlas_rice_acr_boxplot_update.pdf'),width = 5,height = 6)
p
dev.off()

##not overlapversion
pdf(paste0(opt_dir,'/opt_variable_acr_cns_capture_atlas_rice_acr_boxplot_not_overlapCNS_update.pdf'),width = 5,height = 6)
p
dev.off()


##check the statistic
overlap_CNS_broad = c(0.4258475,0.4594595,0.5305720,0.4869110)
overlap_CNS_CT = c(0.5607940,0.6744186,0.6917808,0.5723684)

not_overlap_CNS_broad = c(0.2603383,0.3978495,0.3833737,0.3770370)
not_overlap_CNS_CT = c(0.3610315,0.4758910,0.5319465,0.4852941)

t.test(overlap_CNS_broad,not_overlap_CNS_broad)

t.test(overlap_CNS_CT,not_overlap_CNS_CT)


##updating 010524
##check the not UTR regions database only maize to rice
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//database_notcontain_UTR/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_')

build_dt <- function(ipt_dt,spe_nm,captureCNSorNot){
  
  ##select the variable ACR containing the CNS
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  
  if (captureCNSorNot == 'yes'){
    ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  }else{
    ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'no',]
  }
  
  #ipt_dt <- ipt_dt[ipt_dt[[paste0('rice','_region_capture_CNS')]] == 'yes',]
  
  maize_acr_num <- length(unique(ipt_dt[[paste0(spe_nm,'_ACR')]]))
  
  ##how many of them could be detect in having rice ACR
  head(ipt_dt)
  ipt_dt <- ipt_dt[ipt_dt$Atlas_ACR != 'none',]
  maize_acr_capture_atlas_acr_num <-  length(unique(ipt_dt[[paste0(spe_nm,'_ACR')]]))
  
  ##make a dataframe
  dt <- data.frame(
    spe = c(spe_nm,spe_nm),
    cate = c('VarACRCNS_atlasACR','VarACRCNS_notatlasACR'),
    number = c(maize_acr_capture_atlas_acr_num,maize_acr_num-maize_acr_capture_atlas_acr_num),
    prop = c(maize_acr_capture_atlas_acr_num/maize_acr_num, (maize_acr_num-maize_acr_capture_atlas_acr_num)/maize_acr_num)
  )
  
  return(dt)
  
}

build_organ_num_dt <- function(ipt_dt,spe_nm){
  
  ##select the variable ACR containing the CNS
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  
  ##how many of them could be detect in having rice ACR
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$Atlas_ACR != 'none',]
  #ipt_dt <- 
  
}


outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_summary_add_rice_atlas_acr.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path, header = T)
  head(ipt_dt)
  
  ipt_CT_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_CelltypeCate')]] != 'broadly_accessible',]
  CT_dt_captureCNS <- build_dt(ipt_CT_dt,spe_nm,'yes')
  CT_dt_notcaptureCNS <- build_dt(ipt_CT_dt,spe_nm,'no')
  CT_dt_captureCNS$captureCNS <- 'overlapCNS'
  CT_dt_notcaptureCNS$captureCNS <- 'notCNS'
  CT_dt <- rbind(CT_dt_captureCNS,CT_dt_notcaptureCNS)
  CT_dt$celltype <- 'CT'
  
  ipt_Broad_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_CelltypeCate')]] == 'broadly_accessible',]
  Broad_dt_captureCNS <- build_dt(ipt_Broad_dt,spe_nm,'yes')
  Broad_dt_notcaptureCNS <- build_dt(ipt_Broad_dt,spe_nm,'no')
  Broad_dt_captureCNS$captureCNS <- 'overlapCNS'
  Broad_dt_notcaptureCNS$captureCNS <- 'notCNS'
  Broad_dt <- rbind(Broad_dt_captureCNS,Broad_dt_notcaptureCNS)
  Broad_dt$celltype <- 'Broad'
  
  combine_dt <- rbind(CT_dt,Broad_dt)
  
  return(combine_dt)
  
})

combine_dt <- do.call(rbind,outs)

##do the composition plot
#combine_dt$spe <- factor(combine_dt$spe,levels = c('maize','sorghum','Pm','Uf'))

combine_varACRCNS_prop_dt <- combine_dt[combine_dt$cate == 'VarACRCNS_atlasACR',]

p <- ggplot(combine_varACRCNS_prop_dt, aes(fill=celltype , y=prop, x=captureCNS)) + 
  geom_bar(stat="identity", position=position_dodge())
#p <- p + geom_text(aes(label = prop),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
p <- p+labs(x="\nCell type", y = "Proportion\n", cex.lab = 7) ## use \n to set the position
#p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 30), ##change all the text size
             axis.text.y = element_text(size=40,colour = "black"),
             axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
             axis.title.x = element_text(size=40,colour = "black"),
             axis.title.y = element_text(size=40,colour = "black"),
             plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
#p <- p + facet_wrap(~celltype)

p <- p + coord_cartesian(ylim = c(0, 1)) + ##threshold
  scale_y_continuous(breaks=seq(0, 1, 0.1))
pdf(paste0(opt_dir,'/opt_variable_acr_cns_capture_atlas_rice_acr.pdf'),width = 16,height = 10)
p
dev.off()

##Here we will consider the sig test
combine_broad_dt <- combine_dt[combine_dt$celltype == 'Broad',]
combine_CT_dt <- combine_dt[combine_dt$celltype == 'CT',]

combine_broad_atlas_dt <- combine_broad_dt[combine_broad_dt$cate == 'VarACRCNS_atlasACR',]
combine_broad_atlas_CNS_dt <- combine_broad_atlas_dt[combine_broad_atlas_dt$captureCNS == 'overlapCNS',]
combine_broad_atlas_not_CNS_dt <- combine_broad_atlas_dt[combine_broad_atlas_dt$captureCNS == 'notCNS',]

combine_broad_not_atlas_dt <- combine_broad_dt[combine_broad_dt$cate == 'VarACRCNS_notatlasACR',]
combine_broad_not_atlas_CNS_dt <- combine_broad_not_atlas_dt[combine_broad_not_atlas_dt$captureCNS == 'overlapCNS',]
combine_broad_not_atlas_not_CNS_dt <- combine_broad_not_atlas_dt[combine_broad_not_atlas_dt$captureCNS == 'notCNS',]

overlapCNS_in_capture_atlasACR_num <- combine_broad_atlas_CNS_dt$number
overlapCNS_not_capture_atlasACR_num <- combine_broad_not_atlas_CNS_dt$number
notoverlapCNS_in_capture_atlasACR_num <- combine_broad_atlas_not_CNS_dt$number
notoverlapCNS_not_capture_atlasACR_num <- combine_broad_not_atlas_not_CNS_dt$number

data <- matrix(c(overlapCNS_in_capture_atlasACR_num, notoverlapCNS_in_capture_atlasACR_num, overlapCNS_not_capture_atlasACR_num, notoverlapCNS_not_capture_atlasACR_num), nrow = 2)

fisher_result <- fisher.test(data,alternative = 'greater')
print(fisher_result$p.value)


combine_CT_atlas_dt <- combine_CT_dt[combine_CT_dt$cate == 'VarACRCNS_atlasACR',]
combine_CT_atlas_CNS_dt <- combine_CT_atlas_dt[combine_CT_atlas_dt$captureCNS == 'overlapCNS',]
combine_CT_atlas_not_CNS_dt <- combine_CT_atlas_dt[combine_CT_atlas_dt$captureCNS == 'notCNS',]

combine_CT_not_atlas_dt <- combine_CT_dt[combine_CT_dt$cate == 'VarACRCNS_notatlasACR',]
combine_CT_not_atlas_CNS_dt <- combine_CT_not_atlas_dt[combine_CT_not_atlas_dt$captureCNS == 'overlapCNS',]
combine_CT_not_atlas_not_CNS_dt <- combine_CT_not_atlas_dt[combine_CT_not_atlas_dt$captureCNS == 'notCNS',]

overlapCNS_in_capture_atlasACR_num <- combine_CT_atlas_CNS_dt$number
overlapCNS_not_capture_atlasACR_num <- combine_CT_not_atlas_CNS_dt$number
notoverlapCNS_in_capture_atlasACR_num <- combine_CT_atlas_not_CNS_dt$number
notoverlapCNS_not_capture_atlasACR_num <- combine_CT_not_atlas_not_CNS_dt$number

data <- matrix(c(overlapCNS_in_capture_atlasACR_num, notoverlapCNS_in_capture_atlasACR_num, overlapCNS_not_capture_atlasACR_num, notoverlapCNS_not_capture_atlasACR_num), nrow = 2)

fisher_result <- fisher.test(data,alternative = 'greater')
print(fisher_result$p.value)


#################
##updating 011124
##we will check how many organs and ct would be captured for each ACR in between CNS and not CNS
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//database_notcontain_UTR/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_maize_to_rice_summary_add_rice_atlas_acr_rmleaf')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_summary_add_rice_atlas_acr_rmleaf.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = T)
  head(ipt_dt)
  dim(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$capture_altas_celltypenum_noleaf != 'none',]  

  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  ipt_coverCNS_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  ipt_coverCNS_dt <- ipt_coverCNS_dt[ipt_coverCNS_dt$rice_region_capture_CNS == 'yes',]
  
  ipt_not_coverCNS_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] != 'yes',]
  ipt_not_coverCNS_dt <- ipt_not_coverCNS_dt[ipt_not_coverCNS_dt$rice_region_capture_CNS == 'yes',]
  
  ##updating 020224
  ##we will allow the rice region to contain the CNS
  
  
  
  
  
  ##check the cell states
  head(ipt_coverCNS_dt)
  ipt_coverCNS_dt <- ipt_coverCNS_dt[c('capture_altas_celltypenum_noleaf','capture_atlas_organ_num')]
  ipt_coverCNS_dt$cate <- c('CoverCNS')
  
  ipt_not_coverCNS_dt <- ipt_not_coverCNS_dt[c('capture_altas_celltypenum_noleaf','capture_atlas_organ_num')]
  ipt_not_coverCNS_dt$cate <- c('NotCoverCNS')
  
  combine_dt <- rbind(ipt_coverCNS_dt,ipt_not_coverCNS_dt)
  
  combine_dt$capture_altas_celltypenum_noleaf <- as.numeric(combine_dt$capture_altas_celltypenum_noleaf)
  combine_dt$capture_atlas_organ_num <- as.numeric(combine_dt$capture_atlas_organ_num)
  
  
  
  p<-ggplot(combine_dt, aes(x=capture_altas_celltypenum_noleaf, fill=cate, color = cate)) +
    geom_histogram(position="identity", alpha=0.5) + 
    #scale_color_manual(values=c("Broad"="#ED1C24", "CT"="#FFDE17")) +
    #scale_fill_manual(values=c("Broad"="#ED1C24", "CT"="#FFDE17")) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
   # coord_cartesian(ylim = c(0,7000)) +
  #  scale_y_continuous(breaks=seq(0,7000,2000))
  pdf(paste0(opt_dir,'/opt_density_variable_acr_cnsnotcns_altas_cellstate_num.pdf'),width = 8,height = 8)
  print(p)
  dev.off()
  
  combine_dt$capture_altas_celltypenum_noleaf <- as.numeric(combine_dt$capture_altas_celltypenum_noleaf)
  #max(combine_dt$capture_altas_celltypenum_noleaf)
  
  ##use the density plot
  p<-ggplot(combine_dt, aes(x=capture_altas_celltypenum_noleaf, color=cate,fill=cate)) +
    geom_density(alpha=0.4)
  #p <- p + coord_cartesian(xlim = c(0, 2))
  #p <- p + scale_x_continuous(breaks=seq(0,2, 0.25))
  p <- p + theme(
    plot.title = element_text(size=20,hjust = 0.5),
    axis.title = element_text(size =20),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black"),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) + coord_cartesian(xlim = c(0, 120)) + scale_x_continuous(breaks=seq(0,120, 10))

  pdf(paste0(opt_dir,'/opt_density_real_plot_variable_acr_cnsnotcns_altas_cellstate_num_rmleaf.pdf'),width = 8,height = 8)
  print(p)
  dev.off()
  
  
  ##plot the difference or organ capture number
  p<-ggplot(combine_dt, aes(x=capture_atlas_organ_num, fill=cate, color = cate)) +
    geom_histogram(position="identity", alpha=0.5) + 
    #scale_color_manual(values=c("Broad"="#ED1C24", "CT"="#FFDE17")) +
    #scale_fill_manual(values=c("Broad"="#ED1C24", "CT"="#FFDE17")) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  # coord_cartesian(ylim = c(0,7000)) +
  #  scale_y_continuous(breaks=seq(0,7000,2000))
  pdf(paste0(opt_dir,'/opt_density_variable_acr_cnsnotcns_altas_organ_num.pdf'),width = 8,height = 8)
  print(p)
  dev.off()
  
  
  
  
})









#################
##updating 011124
##check number of acr for cell type 
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//database_notcontain_UTR/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'relativeDiffProp')


##define a function first
plot_barplot <- function(ipt_flt_dt,prefix){
  
  ##order based on the avg value
  aggregate_median <- aggregate(OverlapCNS_ACRprop ~ FinalCelltype, data = ipt_flt_dt, mean)
  rownames(aggregate_median) <- aggregate_median$FinalCelltype
  aggregate_median <- aggregate_median[order(aggregate_median$OverlapCNS_ACRprop,decreasing = T),]
  
  ipt_flt_dt$FinalCelltype <- factor(ipt_flt_dt$FinalCelltype,levels = rownames(aggregate_median))
  
  p <- ggplot(ipt_flt_dt, aes(x=FinalCelltype,y=OverlapCNS_ACRprop,fill=FinalCelltype)) + 
    #geom_violin(trim=FALSE) + 
    geom_boxplot(position=position_dodge(1)) +
    #geom_boxplot(width=0.5,fill='white') +
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,angle = 45, vjust = 1,hjust =1),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks=seq(0, 1, 0.2))
  pdf(paste0(opt_dir,'/opt_all_organ_overlapCNS_prop_',prefix,'.pdf'),width = 10,height = 10)
  print(p)
  dev.off()
  
  
  
}


outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_celltype_ACRcount_CNSaltasACR_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_relativeDiffProp.txt','',target_fl_path_nm)
  spe_nm <- gsub('_rice','',target_fl_path_nm)

  ipt_dt <- read.delim(target_fl_path,header = T)
  head(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$Organ != 'leaf',]

  ##plot organ
  ipt_dt$Organ <- factor(ipt_dt$Organ,levels = c('bud','Eseedling','seedling','panicle','crownroot','semroot',
                                                 'Eseed','Lseed'))
  
  p <- ggplot(ipt_dt, aes(x=Organ,y=OverlapCNS_ACRprop,fill=Organ)) + 
    #geom_violin(trim=FALSE) + 
    geom_boxplot(position=position_dodge(1)) +
    #geom_boxplot(width=0.5,fill='white') +
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,angle = 45, vjust = 1,hjust =1),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks=seq(0, 1, 0.2))
  #scale_fill_brewer(palette="Blues") +
  #facet_wrap(~spe,ncol = 1)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)
  
    

  pdf(paste0(opt_dir,'/opt_all_organ_overlapCNS_prop.pdf'),width = 10,height = 10)
  print(p)
  dev.off()
  
  ##check the statistical sig 
  ipt_dt$leafRelated <- ifelse(ipt_dt$Organ == 'Eseedling'|ipt_dt$Organ == 'seedling' | ipt_dt$Organ == 'bud','LeafRelated','Other')
  ipt_leafrelated_dt <- ipt_dt[ipt_dt$leafRelated == 'LeafRelated',]
  ipt_other_dt <- ipt_dt[ipt_dt$leafRelated != 'LeafRelated',]
  compare_leaf_not_leaf_test <- t.test(ipt_leafrelated_dt$OverlapCNS_ACRprop,ipt_other_dt$OverlapCNS_ACRprop,alternative = 'greater')
  
  
  ##we will plot for the organ for the specific cell types
  ipt_flt_dt <- ipt_dt[ipt_dt$FinalCelltype != 'notInclude',]

  ipt_flt_leafrelated_dt <- ipt_flt_dt[ipt_flt_dt$leafRelated == 'LeafRelated',]
  plot_barplot(ipt_flt_leafrelated_dt,'LeafRelated')
  
  ##for the root
  ipt_flt_rootrelated_dt <- ipt_flt_dt[ipt_flt_dt$Organ == 'semroot'|ipt_flt_dt$Organ == 'crownroot',]
  plot_barplot(ipt_flt_rootrelated_dt,'RootRelated')
  
  ##for the panicle
  ipt_flt_panicle_dt <- ipt_flt_dt[ipt_flt_dt$Organ == 'panicle',]
  plot_barplot(ipt_flt_panicle_dt,'Panicle')
  
  ##for the seed
  ipt_flt_seed_dt <- ipt_flt_dt[ipt_flt_dt$Organ == 'Lseed'| ipt_flt_dt$Organ == 'Eseed',]
  plot_barplot(ipt_flt_seed_dt,'SeedRelated')
  

  
})


#################
##updating 012424
##we will plot the exact number of cell type specific ACRs
##also try to plot the total number of cell type specific ACRs or proportion of the cell type specific ACRs
##Here we will try to plot a target cell types specific ACR number
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//database_notcontain_UTR/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'relativeDiffProp')

outs <- lapply(all_fl_list, function(x){
  
  ##read the file
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_celltype_ACRcount_CNSaltasACR_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('_relativeDiffProp.txt','',target_fl_path_nm)
  spe_nm <- gsub('_rice','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = T)
  head(ipt_dt)
  
  ipt_Eseed_dt <- ipt_dt[ipt_dt$Organ == 'Eseed',]
  
  ##for the meta
  ipt_meta_dt <- read.delim('/Users/haidong/Desktop/Repeat_analysis/Post_doc_projects/rice_altlas_scATAC_seq_042121/collect_final_annotation_012222/opt_final_comb_annot_ver8.4.txt',header = T)
  head(ipt_meta_dt)
  
  ##for the total DA peak
  target_fl_path <- paste0(ipt_dir,'/opt_final_summary_DApeak_tn5_cellnum.txt')
  ipt_final_DA_peak_dt <- read.delim(target_fl_path,header = F)
  head(ipt_final_DA_peak_dt)
  ipt_final_DA_peak_dt$organ <- gsub('\\..+','',ipt_final_DA_peak_dt$V2)
  
  #################################################
  ##we will check the total number of specific acrs
  ipt_summary_dt <- read.delim(paste0(ipt_dir,'/opt_maize_to_rice_summary_add_rice_atlas_acr_rmleaf.txt'))
  head(ipt_summary_dt)
  ipt_summary_target_dt <- ipt_summary_dt[grepl('.',ipt_summary_dt$Atlas_ACR_celltype),]
  ipt_summary_broad_dt <- ipt_summary_dt[!grepl('.',ipt_summary_dt$Atlas_ACR_celltype),]
  head(ipt_summary_broad_dt)
  
  ipt_summary_target_dt <- ipt_summary_target_dt[ipt_summary_target_dt$Atlas_ACR_celltype != 'none',]
  ipt_summary_target_dt <- ipt_summary_target_dt[ipt_summary_target_dt$maize_ACRID != 'none',]
  ipt_summary_target_dt <- ipt_summary_target_dt[ipt_summary_target_dt$rice_ACR == 'none',]
  ipt_summary_target_dt <- ipt_summary_target_dt[ipt_summary_target_dt$maize_ACR_capture_CNS == 'yes',]
  ipt_summary_target_dt <- ipt_summary_target_dt[ipt_summary_target_dt$rice_region_capture_CNS == 'yes',]
  dim(ipt_summary_target_dt)
  
  ipt_summary_broad_dt <- ipt_summary_broad_dt[ipt_summary_broad_dt$maize_ACRID != 'none',]
  ipt_summary_broad_dt <- ipt_summary_broad_dt[ipt_summary_broad_dt$rice_ACR == 'none',]
  ipt_summary_broad_dt <- ipt_summary_broad_dt[ipt_summary_broad_dt$maize_ACR_capture_CNS == 'yes',]
  ipt_summary_broad_dt <- ipt_summary_broad_dt[ipt_summary_broad_dt$rice_region_capture_CNS == 'yes',]
  dim(ipt_summary_broad_dt)
  length(unique(ipt_summary_broad_dt$Atlas_ACR))
  
  ipt_summary_target_split_dt <- as.data.frame(separate_rows(ipt_summary_target_dt, Atlas_ACR_celltype, sep = ","))
  ipt_summary_target_split_dt <- ipt_summary_target_split_dt[grepl('.',ipt_summary_target_split_dt$Atlas_ACR_celltype),]
  ipt_summary_target_split_dt$organ <- gsub('\\..+','',ipt_summary_target_split_dt$Atlas_ACR_celltype)
  
  all_organ_list <- c('bud','Eseedling','seedling','semroot','crownroot','panicle','Eseed','Lseed')
  
  outs <- lapply(all_organ_list, function(x){
    
    ipt_summary_target_split_torgan_dt <- ipt_summary_target_split_dt[ipt_summary_target_split_dt$organ == x,]
    ACRnum <- length(unique(ipt_summary_target_split_torgan_dt$Atlas_ACR))
    
    res <- c()
    res$organ <- x
    res$num <- ACRnum
    res <- as.data.frame(res)
    return(res)
  
  })
  
  combine_dt <- do.call(rbind,outs)
  
  broad_dt <- data.frame(
    organ <- 'broad',
    num <- length(unique(ipt_summary_broad_dt$Atlas_ACR))
  )
  colnames(broad_dt) <- c('organ','num')
  combine_dt <- rbind(combine_dt, broad_dt)
  
  ##we do not fill that is too messy
  combine_dt$organ <- factor(combine_dt$organ,levels = c('broad','semroot','panicle','Eseed','bud','Eseedling','seedling','crownroot','Lseed'))
  p <- ggplot(data=combine_dt, aes(x=organ, y=num)) +
    #scale_fill_manual(values=color_list)+
    #geom_bar(stat="identity", position=position_dodge()) +
    geom_bar(stat="identity") +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=10),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) 
  #  facet_wrap(~Rice_CellType,nrow = 1)

  pdf(paste0(opt_dir,'/opt_all_organ_CNSatlasACR_broad_CT_num_barplot.pdf'),width = 12,height = 6 )
  print(p)
  dev.off()
  
  ########################################
  ##check CT acr number for each cell type
  
  target_organ_list <- c('broad','bud','Eseedling','seedling','semroot','crownroot','panicle','Eseed','Lseed')
  
  outs2 <- lapply(target_organ_list, function(y){
    
    ipt_target_organ_dt <- ipt_dt[ipt_dt$Organ == y,]
    
    ipt_target_organ_DA_peak_dt <- ipt_final_DA_peak_dt[ipt_final_DA_peak_dt$organ == y,]
    
    ##order the ACR count
    ipt_target_organ_dt <- ipt_target_organ_dt[order(ipt_target_organ_dt$OverlapCNS_ACRcount,decreasing = T),]
    celltype_list = ipt_target_organ_dt$Celltype
    ipt_target_organ_dt$Celltype <- factor(ipt_target_organ_dt$Celltype,levels = ipt_target_organ_dt$Celltype)
    
    shared_celltypes <- intersect(ipt_target_organ_dt$Celltype ,ipt_target_organ_DA_peak_dt$V2)
    ipt_target_organ_dt <- ipt_target_organ_dt[ipt_target_organ_dt$Celltype %in% shared_celltypes,]
    
    ##find the color of cell types
    #ipt_meta_organ_dt <- ipt_meta_dt[ipt_meta_dt$tissue == y,]
    #ipt_meta_organ_dt <- unique(ipt_meta_organ_dt[c('color_TCP','Final_annotation_TCP_up')])

    #color_list <- ipt_meta_organ_dt$color_TCP
    #names(color_list) <- ipt_meta_organ_dt$Final_annotation_TCP_up 


    ##we do not fill that is too messy
    p <- ggplot(data=ipt_target_organ_dt, aes(x=Celltype, y=OverlapCNS_ACRcount)) +
      #scale_fill_manual(values=color_list)+
      #geom_bar(stat="identity", position=position_dodge()) +
      geom_bar(stat="identity") +
      theme(
        plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
        axis.title = element_text(size =20, face="bold"),
        #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
        axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
        axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
        axis.ticks = element_line(size = rel(2.5)),
        axis.ticks.length = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 15),
        strip.text = element_text(size=10),
        #legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
      ) 
    #  facet_wrap(~Rice_CellType,nrow = 1)

    
    pdf(paste0(opt_dir,'/opt_',y,'_all_celltype_num_orderby_ACRnum_CNSoverlapVariableACR.pdf'),width = 12,height = 6 )
    print(p)
    dev.off()
    
    
    ##build a bar plot for the DA peaks
    ipt_target_organ_DA_peak_dt <- ipt_target_organ_DA_peak_dt[ipt_target_organ_DA_peak_dt$V2 %in% shared_celltypes,]
    ipt_target_organ_DA_peak_dt$V2 <- factor(ipt_target_organ_DA_peak_dt$V2,levels =shared_celltypes)
  
    p <- ggplot(data=ipt_target_organ_DA_peak_dt, aes(x=V2, y=V3)) +
      #scale_fill_manual(values=c("bundle_sheath"="#3084AF", "companion_cells_sieve_elements"="#eddc9f",
      #                           "epidermis"="#C49B7B","mesophyll"="#91e388","protoderm"="#fcd2b1"))+
      #geom_bar(stat="identity", position=position_dodge()) +
      geom_bar(stat="identity") +
      theme(
        plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
        axis.title = element_text(size =20, face="bold"),
        #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
        axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
        axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
        axis.ticks = element_line(size = rel(2.5)),
        axis.ticks.length = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 15),
        strip.text = element_text(size=10),
        #legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
      ) 
    #  facet_wrap(~Rice_CellType,nrow = 1)
   
    pdf(paste0(opt_dir,'/opt_',y,'_all_celltype_num_orderby_ACRnum_totalACR.pdf'),width = 12,height = 6 )
    print(p)
    dev.off()
    
    
  })
  

})


















#################
##updating 011224
##check the cell type composition in the variable ACRs shown the maize
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//database_notcontain_UTR/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_maize_to_rice_summary_add_rice_atlas_acr_rmleaf')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_summary_add_rice_atlas_acr_rmleaf.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = T)
  ipt_dt <- ipt_dt[ipt_dt$capture_altas_celltypenum != 'none',]  
  
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  ipt_coverCNS_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  length(unique(ipt_coverCNS_dt$maize_ACR))
  
  ipt_coverCNS_dt <- ipt_coverCNS_dt[ipt_coverCNS_dt[[paste0('rice_region_capture_CNS')]] == 'yes',]
  length(unique(ipt_coverCNS_dt$maize_ACR))
  
  ##updating 020224 do not considering the not cover the CNS
  #ipt_not_coverCNS_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] != 'yes',]

  ipt_coverCNS_maize_celltype_dt <- unique(ipt_coverCNS_dt[c('maize_ACR','maize_CelltypeCate')])
  ipt_coverCNS_maize_celltype_dt <- as.data.frame(separate_rows(ipt_coverCNS_maize_celltype_dt, maize_CelltypeCate, sep = ","))
  coverCNS_dt <- as.data.frame(table(ipt_coverCNS_maize_celltype_dt$maize_CelltypeCate))
  coverCNS_dt$cate <- 'overlap_CNS'
  
  #ipt_not_coverCNS_maize_celltype_dt <- unique(ipt_not_coverCNS_dt[c('maize_ACR','maize_CelltypeCate')])
  #ipt_not_coverCNS_maize_celltype_dt <- as.data.frame(separate_rows(ipt_not_coverCNS_maize_celltype_dt, maize_CelltypeCate, sep = ","))
  #notcoverCNS_dt <- as.data.frame(table(ipt_not_coverCNS_maize_celltype_dt$maize_CelltypeCate))
  #notcoverCNS_dt$cate <- 'notoverlap_CNS'
  
  #merged_dt <- rbind(coverCNS_dt,notcoverCNS_dt)  
  
  p <- ggplot(data=coverCNS_dt, aes(x=Var1, y=Freq, fill=cate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    #geom_bar(stat="identity") +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=10),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )  
  #  coord_cartesian(ylim = c(0, 0.4)) +
  #  scale_y_continuous(breaks=seq(0, 0.4, 0.1))

  
  #facet_wrap(~Rice_CellType,nrow = 1)
  pdf(paste0(opt_dir,'/opt_check_',spe_nm,'_variable_ACR_CT_barplot_makesureRiceRegionCNS.pdf'),width = 8,height = 8 )
  print(p)
  dev.off()
  
  
  

})


##updating 020924
##check how many CNS ACR were found to be 
##ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
ipt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_121923//database_version_notContainUTR//'
opt_dir <- 'opt2_check_cns_num_per_acr_variableACR_overlap_atlas_ACR_results_121923//database_notcontain_UTR/'

all_fl_list <- list.files(path= ipt_dir,pattern = 'opt_maize_to_rice_summary_add_rice_atlas_acr_rmleaf')

generate_pie <- function(df,cate) {
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df[[cate]],pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}


outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_summary_add_rice_atlas_acr_rmleaf.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  spe_nm <- gsub('_to_rice','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = T)

  ##should have blast  
  ipt_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACRID')]] != 'none',]
  ipt_dt <- ipt_dt[ipt_dt$rice_ACR == 'none',]
  ##shoud capture cns for ACR
  ipt_coverCNS_dt <- ipt_dt[ipt_dt[[paste0(spe_nm,'_ACR_capture_CNS')]] == 'yes',]
  length(unique(ipt_coverCNS_dt$maize_ACR))
  
  ##should capture cns for region
  ipt_coverCNS_dt <- ipt_coverCNS_dt[ipt_coverCNS_dt[[paste0('rice_region_capture_CNS')]] == 'yes',]
  length(unique(ipt_coverCNS_dt$maize_ACR))
  
  ipt_coverCNS_capture_atlas_dt <- ipt_coverCNS_dt[ipt_coverCNS_dt$capture_altas_celltypenum != 'none',]  
  length(unique(ipt_coverCNS_capture_atlas_dt$maize_ACR))
  
  ##generate the pie chart
  target_dt <- data.frame(
    cate <- c('capture','not_capture'),
    number <- c(248,84)
  )
  
  colnames(target_dt) <- c('cate','count')
  
  pdf(paste0(opt_dir,'/opt_variable_acr_cns_capture_atlas_rice_acr_piechart.pdf'),width = 6,height = 6)
  generate_pie(target_dt,'cate')
  dev.off()
  
  
})






##updating 120823
#########################################################################
##check the genes and expression under the stress condition per cell type
##this part looks like the mesophyll and bundle are more likely reponsive to stress related genes
##so we will check the gene expression for snRNAseq
##we will not use this part
ipt_dir <- 'opt2_nearby_genes_120623/////'
opt_dir <- 'opt2_nearby_genes_results_120623///'

all_fl_list <- list.files(path = ipt_dir,pattern = 'celltypeCpm_addGeneCate')

#check total of species-specific ACRs
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene_add_celltypeCpm_addGeneCate.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  
  #ipt_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none',]
  #head(ipt_dt)
  
  ipt_dt <- ipt_dt[ipt_dt$rice_CelltypeCate != 'broadly_accessible',]
  
  ipt_dt_target <- as.data.frame(ipt_dt[c('rice_ACR')])
  ipt_dt_target$spe <- target_fl_path_nm
  head(ipt_dt_target)
  
  return(ipt_dt_target)

})

combine_dt <- do.call(rbind,outs)
length(unique(combine_dt$rice_ACR)) ##2452

combine_dt_maize <- combine_dt[combine_dt$spe == 'maize',]
combine_dt_sorghum <- combine_dt[combine_dt$spe == 'sorghum',]
combine_dt_Pm <- combine_dt[combine_dt$spe == 'Pm',]
combine_dt_Uf <- combine_dt[combine_dt$spe == 'Uf',]
maize_sorghum <- intersect(combine_dt_maize$rice_ACR,combine_dt_sorghum$rice_ACR)
maize_sorghum_Pm <- intersect(maize_sorghum,combine_dt_Pm$rice_ACR)
maize_sorghum_Pm_Uf <- intersect(maize_sorghum_Pm,combine_dt_Uf$rice_ACR)



outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene_add_celltypeCpm_addGeneCate.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)

  ipt_dt <- read.delim(target_fl_path)
  
  #ipt_dt <- ipt_dt[c('rice_CelltypeCate','geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')]
  ipt_dt <- as.data.frame(separate_rows(ipt_dt, rice_CelltypeCate, sep = ","))
  head(ipt_dt)
  
  ##for the not shared
  ipt_target_noshared_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none',]
  ipt_target_noshared_dt <- unique(ipt_target_noshared_dt[c('rice_ACR','rice_CelltypeCate','geneCate','Dist','geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')])
  ipt_target_noshared_dt <- ipt_target_noshared_dt[ipt_target_noshared_dt$seedling.Epidermis != 'none',]
  ipt_target_noshared_dt <- ipt_target_noshared_dt[ipt_target_noshared_dt$rice_CelltypeCate == 'mesophyll' | 
                                                     ipt_target_noshared_dt$rice_CelltypeCate == 'epidermis'|
                                                     ipt_target_noshared_dt$rice_CelltypeCate == 'companion_cell',]
  ipt_target_noshared_dt$cate <- 'not_shared'
  
  ##for the shared I
  ipt_target_sharedI_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_fl_path_nm,'_ACR')]] == 'none',]
  ipt_target_sharedI_dt <- unique(ipt_target_sharedI_dt[c('rice_ACR','rice_CelltypeCate','geneCate','Dist','geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')])
  ipt_target_sharedI_dt <- ipt_target_sharedI_dt[ipt_target_sharedI_dt$seedling.Epidermis != 'none',]
  ipt_target_sharedI_dt <- ipt_target_sharedI_dt[ipt_target_sharedI_dt$rice_CelltypeCate == 'mesophyll' | 
                                                   ipt_target_sharedI_dt$rice_CelltypeCate == 'epidermis'|
                                                   ipt_target_sharedI_dt$rice_CelltypeCate == 'companion_cell',]
  ipt_target_sharedI_dt$cate <- 'sharedI'
  
  ##for the shared A
  ipt_target_sharedA_dt <- ipt_dt[ipt_dt$rice_ACR != 'none' & ipt_dt[[paste0(target_fl_path_nm,'_ACR')]] != 'none',]
  ipt_target_sharedA_dt <- unique(ipt_target_sharedA_dt[c('rice_ACR','rice_CelltypeCate','geneCate','Dist','geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')])
  ipt_target_sharedA_dt <- ipt_target_sharedA_dt[ipt_target_sharedA_dt$seedling.Epidermis != 'none',]
  ipt_target_sharedA_dt <- ipt_target_sharedA_dt[ipt_target_sharedA_dt$rice_CelltypeCate == 'mesophyll' | 
                                                   ipt_target_sharedA_dt$rice_CelltypeCate == 'epidermis'|
                                                   ipt_target_sharedA_dt$rice_CelltypeCate == 'companion_cell',]
  
  ipt_target_sharedA_dt$cate <- 'sharedA'
  
  merged_dt <- rbind(ipt_target_noshared_dt,ipt_target_sharedI_dt,ipt_target_sharedA_dt)
  merged_dt$spe <- target_fl_path_nm
  return(merged_dt)
  
})

combine_dt <- do.call(rbind,outs)
dim(combine_dt)
length(unique(combine_dt$rice_ACR))

##we will split them into different cate and species
combine_dt_epidermis <- combine_dt[combine_dt$rice_CelltypeCate == 'epidermis',]
combine_dt_epidermis <- combine_dt_epidermis[c('rice_ACR','geneCate','Dist','seedling.Epidermis','cate','spe')]

##build the 
spe_list <- unique(combine_dt_epidermis$spe)
cate_list <- unique(combine_dt_epidermis$cate)
genecate_list <- unique(combine_dt_epidermis$geneCate)

##do the plotting
combine_dt_epidermis$geneCate <- factor(combine_dt_epidermis$geneCate,levels = c('Genic','Proximal','Distal'))
combine_dt_epidermis$cate <- factor(combine_dt_epidermis$cate,levels = c('sharedA','sharedI','not_shared'))

combine_dt_epidermis$seedling.Epidermis <- as.numeric(combine_dt_epidermis$seedling.Epidermis)

combine_dt_epidermis_proximal <- combine_dt_epidermis[combine_dt_epidermis$geneCate == 'Proximal',]
combine_dt_epidermis_proximal <- combine_dt_epidermis_proximal[combine_dt_epidermis_proximal$Dist < 500,]
combine_dt_epidermis_proximal_not <- combine_dt_epidermis[combine_dt_epidermis$geneCate != 'Proximal',]
combine_dt_epidermis <- rbind(combine_dt_epidermis_proximal,combine_dt_epidermis_proximal_not)

p <- ggplot(combine_dt_epidermis, aes(x=geneCate,y=seedling.Epidermis,fill=cate)) + 
  #geom_violin(trim=FALSE) + 
  geom_boxplot(position=position_dodge(1)) +
  #geom_boxplot(width=0.5,fill='white') +
  #geom_jitter(shape=6, position=position_jitter(0.2)) +
  labs(x="\nClass", y = 'Correlation\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35),
    axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  #scale_fill_brewer(palette="Blues") +
  facet_wrap(~spe,ncol = 1)## You can use face_wrap function only if you need it+
#geom_text(data =tukey_letters, 
#          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
#          size = 8,position=position_dodge(.5)  ##change the letter
#)
pdf(paste0(opt_dir,'/opt_epidermis_genecate_exp.pdf'),width = 10,height = 15)
p
dev.off()

##check the test
outs <- lapply(spe_list, function(x){
  
  outs2 <- lapply(genecate_list, function(y){
    
 
    
    target_dt <- combine_dt_epidermis[combine_dt_epidermis$spe == x,]
    
    if (y == 'Proximal'){
      target_dt <- target_dt[target_dt$geneCate == y,]
      target_dt <- target_dt[target_dt$Dist > 500,]
    }else{
      target_dt <- target_dt[target_dt$geneCate == y,]
    }
        
        
    target_dt$seedling.Epidermis <- as.numeric(target_dt$seedling.Epidermis)
    
    ##for sharedA to sharedI
    target_sharedA_dt <- target_dt[target_dt$cate == 'sharedA',]
    target_sharedI_dt <- target_dt[target_dt$cate == 'sharedI',]
    target_not_shared_dt <- target_dt[target_dt$cate == 'not_shared',]
    
    sharedA_list = target_sharedA_dt$seedling.Epidermis
    sharedI_list = target_sharedI_dt$seedling.Epidermis
    not_shared_list = target_not_shared_dt$seedling.Epidermis
    
    if (mean(sharedA_list) > mean(sharedI_list)){
      res <- t.test( target_sharedA_dt$seedling.Epidermis,target_sharedI_dt$seedling.Epidermis, alternative = 'greater',paired = F)
      pval_sharedA_sharedI <- res$p.value
    }else{
      res <- t.test( target_sharedA_dt$seedling.Epidermis,target_sharedI_dt$seedling.Epidermis, alternative = 'less',paired = F)
      pval_sharedA_sharedI <- res$p.value
    }
    
    if (mean(sharedA_list) > mean(not_shared_list)){
      res <- t.test( target_sharedA_dt$seedling.Epidermis,target_not_shared_dt$seedling.Epidermis, alternative = 'greater',paired = F)
      pval_sharedA_not_shared <- res$p.value
    }else{
      res <- t.test( target_sharedA_dt$seedling.Epidermis,target_not_shared_dt$seedling.Epidermis, alternative = 'less',paired = F)
      pval_sharedA_not_shared <- res$p.value
    }
      
    if (mean(sharedI_list) > mean(not_shared_list)){
      res <- t.test( target_sharedI_dt$seedling.Epidermis,target_not_shared_dt$seedling.Epidermis, alternative = 'greater',paired = F)
      pval_sharedI_not_shared <- res$p.value
    }else{
      res <- t.test( target_sharedI_dt$seedling.Epidermis,target_not_shared_dt$seedling.Epidermis, alternative = 'less',paired = F)
      pval_sharedI_not_shared <- res$p.value
    }

    data <- data.frame(
      Name = c("sharedA_sharedI", "sharedA_not_shared", "sharedI_not_shared"),
      Age = c(pval_sharedA_sharedI, pval_sharedA_not_shared, pval_sharedI_not_shared)
    )
    data$genecate <- y
  
    return(data)
    
  })
  
  combine_dt <- do.call(rbind,outs2)
  combine_dt$spe <- x
  return(combine_dt)
})

final_combine_dt <- do.call(rbind,outs)
write.table(final_combine_dt,paste0(opt_dir,'/opt_epidermis_genecate_exp_ttest.txt'),quote = F,sep = '\t')








#combine_dt_maize <- combine_dt[combine_dt$spe == 'Pm',]

combine_dt_maize <- combine_dt
##for the companion_cell
combine_dt_maize <- combine_dt_maize[combine_dt_maize$rice_CelltypeCate == 'companion_cell',]

combine_dt_maize_notshared <- combine_dt_maize[combine_dt_maize$cate == 'not_shared',]
table(combine_dt_maize_notshared$geneCate)
combine_dt_maize_sharedA <- combine_dt_maize[combine_dt_maize$cate == 'sharedA',]
table(combine_dt_maize_sharedA$geneCate)

combine_dt_maize_notshared_distal <- combine_dt_maize_notshared[combine_dt_maize_notshared$geneCate == 'Genic',]
combine_dt_maize_notshared_distal <- combine_dt_maize_notshared_distal[c('seedling.CompanionCell','cate')]
combine_dt_maize_notshared_distal$seedling.CompanionCell <- as.numeric(combine_dt_maize_notshared_distal$seedling.CompanionCell)

combine_dt_maize_sharedA_distal <- combine_dt_maize_sharedA[combine_dt_maize_sharedA$geneCate == 'Genic',]
combine_dt_maize_sharedA_distal <- combine_dt_maize_sharedA_distal[c('seedling.CompanionCell','cate')]
combine_dt_maize_sharedA_distal$seedling.CompanionCell <- as.numeric(combine_dt_maize_sharedA_distal$seedling.CompanionCell)

wilcox.test(combine_dt_maize_notshared_distal$seedling.CompanionCell,combine_dt_maize_sharedA_distal$seedling.CompanionCell,alternative = 'less',paired = F)

##genic less

combine_dt_maize <- combine_dt_maize[c('seedling.CompanionCell','cate')]
combine_dt_maize$seedling.CompanionCell <- as.numeric(combine_dt_maize$seedling.CompanionCell)
wilcox.test(combine_dt_maize[combine_dt_maize$cate == 'not_shared',]$seedling.CompanionCell,combine_dt_maize[combine_dt_maize$cate == 'sharedA',]$seedling.CompanionCell,alternative = 'less',paired = F)

mean(combine_dt_maize[combine_dt_maize$cate == 'not_shared',]$seedling.CompanionCell)
mean(combine_dt_maize[combine_dt_maize$cate == 'sharedA',]$seedling.CompanionCell)



##for the mesophyll
combine_dt_maize <- combine_dt
combine_dt_maize <- combine_dt_maize[combine_dt_maize$rice_CelltypeCate == 'mesophyll',]

combine_dt_maize_notshared <- combine_dt_maize[combine_dt_maize$cate == 'not_shared',]
table(combine_dt_maize_notshared$geneCate)
combine_dt_maize_sharedA <- combine_dt_maize[combine_dt_maize$cate == 'sharedA',]
table(combine_dt_maize_sharedA$geneCate)


combine_dt_maize_notshared_distal <- combine_dt_maize_notshared[combine_dt_maize_notshared$geneCate == 'Genic',]
combine_dt_maize_notshared_distal <- combine_dt_maize_notshared_distal[c('seedling.Mesophyll','cate')]
combine_dt_maize_notshared_distal$seedling.Mesophyll <- as.numeric(combine_dt_maize_notshared_distal$seedling.Mesophyll)

combine_dt_maize_sharedA_distal <- combine_dt_maize_sharedA[combine_dt_maize_sharedA$geneCate == 'Genic',]
combine_dt_maize_sharedA_distal <- combine_dt_maize_sharedA_distal[c('seedling.Mesophyll','cate')]
combine_dt_maize_sharedA_distal$seedling.Mesophyll <- as.numeric(combine_dt_maize_sharedA_distal$seedling.Mesophyll)

wilcox.test(combine_dt_maize_notshared_distal$seedling.Mesophyll,combine_dt_maize_sharedA_distal$seedling.Mesophyll,alternative = 'less',paired = F)

##distal greater
##proximal less




table(combine_dt_maize$geneCate)
combine_dt_maize <- combine_dt_maize[c('seedling.Mesophyll','cate')]
combine_dt_maize$seedling.Mesophyll <- as.numeric(combine_dt_maize$seedling.Mesophyll)
wilcox.test(combine_dt_maize[combine_dt_maize$cate == 'not_shared',]$seedling.Mesophyll,combine_dt_maize[combine_dt_maize$cate == 'sharedA',]$seedling.Mesophyll,alternative = 'less',paired = F)

mean(combine_dt_maize[combine_dt_maize$cate == 'not_shared',]$seedling.Mesophyll)
mean(combine_dt_maize[combine_dt_maize$cate == 'sharedA',]$seedling.Mesophyll)



##for the epidermis
combine_dt_maize <- combine_dt
#combine_dt_maize <- combine_dt_maize[combine_dt_maize$Dist > 250,]
combine_dt_maize <- combine_dt_maize[combine_dt_maize$rice_CelltypeCate == 'epidermis',]


combine_dt_maize_notshared <- combine_dt_maize[combine_dt_maize$cate == 'not_shared',]
table(combine_dt_maize_notshared$geneCate)
combine_dt_maize_sharedA <- combine_dt_maize[combine_dt_maize$cate == 'sharedA',]
table(combine_dt_maize_sharedA$geneCate)
combine_dt_maize_sharedI <- combine_dt_maize[combine_dt_maize$cate == 'sharedI',]
table(combine_dt_maize_sharedI$geneCate)


combine_dt_maize_notshared_distal <- combine_dt_maize_notshared[combine_dt_maize_notshared$geneCate == 'Proximal',]
combine_dt_maize_notshared_distal <- combine_dt_maize_notshared_distal[c('seedling.Epidermis','cate')]
combine_dt_maize_notshared_distal$seedling.Epidermis <- as.numeric(combine_dt_maize_notshared_distal$seedling.Epidermis)

combine_dt_maize_sharedA_distal <- combine_dt_maize_sharedA[combine_dt_maize_sharedA$geneCate == 'Proximal',]
combine_dt_maize_sharedA_distal <- combine_dt_maize_sharedA_distal[c('seedling.Epidermis','cate')]
combine_dt_maize_sharedA_distal$seedling.Epidermis <- as.numeric(combine_dt_maize_sharedA_distal$seedling.Epidermis)

combine_dt_maize_sharedI <- combine_dt_maize_sharedI[combine_dt_maize_sharedI$geneCate == 'Proximal',]
combine_dt_maize_sharedI <- combine_dt_maize_sharedI[c('seedling.Epidermis','cate')]
combine_dt_maize_sharedI$seedling.Epidermis <- as.numeric(combine_dt_maize_sharedI$seedling.Epidermis)


wilcox.test(combine_dt_maize_notshared_distal$seedling.Epidermis,combine_dt_maize_sharedA_distal$seedling.Epidermis,alternative = 'greater',paired = F)
wilcox.test(combine_dt_maize_notshared_distal$seedling.Epidermis,combine_dt_maize_sharedI$seedling.Epidermis,alternative = 'greater',paired = F)

##Proximal has greater 

##we will build the plot for the epidermis


#################
##updating 122823
##intersect expression to the H3K27me3
ipt_dir <- 'opt2_nearby_genes_check_exp_122723//////'
opt_dir <- 'opt2_nearby_genes_check_exp_results_122723////'

ipt_H3K27me3_dir <- 'opt3_check_H3K27me3_in_syntenic_121223/'


all_fl_list <- list.files(path = ipt_dir,pattern = '')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene_add_celltypeCpm_two_spe_species_specificVer.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  head(ipt_dt)
  
  ipt_summary_H3K27me3_dt <- read.delim(paste0(ipt_H3K27me3_dir,'/opt_final_summary_file_add_rice_maize.txt'))
  head(ipt_summary_H3K27me3_dt)

  merged_dt <- merge(ipt_dt,ipt_summary_H3K27me3_dt,by.x = 'rice_ACR',by.y = 'rice_ACR')
  head(merged_dt)
  
  write.table(merged_dt,paste0(opt_dir,'/opt_merged_H3K27me3_rice_maize_exp.txt'),sep = '\t',quote = F,row.names = F)

  ##check the different cell type and see nearby gene expression
  ##check epidermal rice
  #merged_dt <- merged_dt[merged_dt$riceOverH3K27Cate == 'none',]
  
  merged_dt <- merge(ipt_dt,ipt_summary_H3K27me3_dt,by.x = 'rice_ACR',by.y = 'rice_ACR')
  
  #merged_dt <- merged_dt[merged_dt$riceGeneDist < 5000,]
  
  ##best
  merged_dt_epidermis <- merged_dt[merged_dt$rice_celltype == 'epidermis' |merged_dt$rice_celltype == 'protoderm' | merged_dt$rice_celltype == 'epidermis,protoderm',]
  #merged_dt_epidermis <- merged_dt[merged_dt$rice_celltype == 'epidermis' | merged_dt$rice_celltype == 'epidermis,protoderm',]
  
  ##bad
  #merged_dt_epidermis <- merged_dt[grepl('epidermis',merged_dt$rice_celltype) | grepl('protoderm',merged_dt$rice_celltype),]
  
  ##not best
  #merged_dt_epidermis <- merged_dt[grepl('epidermis',merged_dt$rice_celltype),]
  
  
  merged_dt_epidermis_exp <- merged_dt_epidermis[c('rice_geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')]
  merged_dt_epidermis_exp <- merged_dt_epidermis_exp[merged_dt_epidermis_exp$seedling.Epidermis != 'none',]
  dim(merged_dt_epidermis_exp)
  merged_dt_epidermis_exp$targetct <- c('epidermis')
  
  merged_dt_epi_epi <- merged_dt_epidermis_exp[c('rice_geneID','seedling.Epidermis','targetct')]
  merged_dt_epi_epi$cate <- c('epidermis_epidermis')
  colnames(merged_dt_epi_epi) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_epi_epi <- unique(merged_dt_epi_epi)
  
  merged_dt_epi_meso <- merged_dt_epidermis_exp[c('rice_geneID','seedling.Mesophyll','targetct')]
  merged_dt_epi_meso$cate <- c('epidermis_meso')
  colnames(merged_dt_epi_meso) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_epi_meso <- unique(merged_dt_epi_meso)
  
  merged_dt_epi_cc <- merged_dt_epidermis_exp[c('rice_geneID','seedling.CompanionCell','targetct')]
  merged_dt_epi_cc$cate <- c('epidermis_cc')
  colnames(merged_dt_epi_cc) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_epi_cc <- unique(merged_dt_epi_cc)
  
  ttest_epidermis_meso <- t.test(as.numeric(merged_dt_epi_epi$exp),as.numeric(merged_dt_epi_meso$exp),paired = T,alternative = 'greater')
  ttest_epidermis_companion <- t.test(as.numeric(merged_dt_epi_epi$exp),as.numeric(merged_dt_epi_cc$exp),paired = T,alternative = 'greater')
  mean(as.numeric(merged_dt_epi_epi$exp))
  mean(as.numeric(merged_dt_epi_meso$exp))
  
  ##check the mesophyll
  #merged_dt_mesophyll <- merged_dt[merged_dt$rice_celltype == 'mesophyll',]
  merged_dt_mesophyll <- merged_dt[grepl('mesophyll',merged_dt$rice_celltype),]
  merged_dt_mesophyll_exp <- merged_dt_mesophyll[c('rice_geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')]
  merged_dt_mesophyll_exp <- merged_dt_mesophyll_exp[merged_dt_mesophyll_exp$seedling.Mesophyll != 'none',]
  merged_dt_mesophyll_exp$targetct <- c('mesophyll')
  dim(merged_dt_mesophyll_exp)
  
  merged_dt_meso_epi <- merged_dt_mesophyll_exp[c('rice_geneID','seedling.Epidermis','targetct')]
  merged_dt_meso_epi$cate <- c('meso_epidermis')
  colnames(merged_dt_meso_epi) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_meso_epi <- unique(merged_dt_meso_epi)
  
  merged_dt_meso_meso <- merged_dt_mesophyll_exp[c('rice_geneID','seedling.Mesophyll','targetct')]
  merged_dt_meso_meso$cate <- c('meso_meso')
  colnames(merged_dt_meso_meso) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_meso_meso <- unique(merged_dt_meso_meso)
  
  merged_dt_meso_cc <- merged_dt_mesophyll_exp[c('rice_geneID','seedling.CompanionCell','targetct')]
  merged_dt_meso_cc$cate <- c('meso_cc')
  colnames(merged_dt_meso_cc) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_meso_cc <- unique(merged_dt_meso_cc)
  
  ttest_meso_epidermis <- t.test(as.numeric(merged_dt_meso_meso$exp),as.numeric(merged_dt_meso_epi$exp),paired = T,alternative = 'greater')
  ttest_meso_companion <- t.test(as.numeric(merged_dt_meso_meso$exp),as.numeric(merged_dt_meso_cc$exp),paired = T,alternative = 'greater')


  ##for the companion cell
  #merged_dt_companioncell <- merged_dt[merged_dt$rice_celltype == 'companion_cell',]
  
  merged_dt_companioncell <- merged_dt[grepl('companion_cell',merged_dt$rice_celltype),]
  
  merged_dt_companioncell_exp <- merged_dt_companioncell[c('rice_geneID','seedling.Epidermis','seedling.Mesophyll','seedling.CompanionCell')]
  merged_dt_companioncell_exp <- merged_dt_companioncell_exp[merged_dt_companioncell_exp$seedling.CompanionCell != 'none',]
  merged_dt_companioncell_exp$targetct <- c('cc')
  dim(merged_dt_companioncell_exp)
  
  
  merged_dt_cc_epi <- merged_dt_companioncell_exp[c('rice_geneID','seedling.Epidermis','targetct')]
  merged_dt_cc_epi$cate <- c('cc_epidermis')
  colnames(merged_dt_cc_epi) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_cc_epi <- unique(merged_dt_cc_epi)
  
  merged_dt_cc_meso <- merged_dt_companioncell_exp[c('rice_geneID','seedling.Mesophyll','targetct')]
  merged_dt_cc_meso$cate <- c('cc_meso')
  colnames(merged_dt_cc_meso) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_cc_meso <- unique(merged_dt_cc_meso)
  
  merged_dt_cc_cc <- merged_dt_companioncell_exp[c('rice_geneID','seedling.CompanionCell','targetct')]
  merged_dt_cc_cc$cate <- c('cc_cc')
  colnames(merged_dt_cc_cc) <- c('rice_geneID','exp','targetct','cate')
  merged_dt_cc_cc <- unique(merged_dt_cc_cc)
  
  ttest_companion_epidermis <- t.test(as.numeric(merged_dt_cc_cc$exp),as.numeric(merged_dt_cc_epi$exp),paired = T,alternative = 'greater')
  ttest_companion_meso <- t.test(as.numeric(merged_dt_cc_cc$exp),as.numeric(merged_dt_cc_meso$exp),paired = T,alternative = 'greater')


  #make a dataframe
  data_sig_ttest <- data.frame(
    cate <- c('epidermis_meso','epidermis_cc','meso_epidermis','meso_cc',
              'cc_epidermis','cc_meso'),
    pval <- c(ttest_epidermis_meso$p.value,ttest_epidermis_companion$p.value,
              ttest_meso_epidermis$p.value,ttest_meso_companion$p.value,
              ttest_companion_epidermis$p.value,ttest_companion_meso$p.value)
  )
  
  colnames(data_sig_ttest) <- c('cate','pval')
  data_sig_ttest
  
  write.table(data_sig_test,paste0(opt_dir,'/opt_sig_ttest_species_celltype_specific_ACRs_geneexp.txt'),quote = F,sep = '\t')
  
  
  ##merge all the exp data
  all_exp_dt <- rbind(merged_dt_epi_epi,merged_dt_epi_meso,merged_dt_epi_cc,
                      merged_dt_meso_meso,merged_dt_meso_epi,merged_dt_meso_cc,
                      merged_dt_cc_cc,merged_dt_cc_epi,merged_dt_cc_meso)
  
  all_exp_dt$targetct <- gsub('_.+','',all_exp_dt$cate)
  all_exp_dt$otherct <- gsub('.+_','',all_exp_dt$cate)

  
  #all_exp_dt$otherct <- factor(all_exp_dt$otherct,levels = rev(c('meso','cc','epidermis')))
  all_exp_dt$targetct <- factor(all_exp_dt$targetct,levels = rev(c('meso','cc','epidermis')))
  all_exp_dt$exp <- as.numeric(all_exp_dt$exp)
  
  all_exp_dt_cate <- unique(all_exp_dt[c('cate','targetct','otherct')])
  rownames(all_exp_dt_cate) <- all_exp_dt_cate$cate
  all_exp_dt_cate$targetct <- factor(all_exp_dt_cate$targetct,levels = rev(c('meso','cc','epidermis')))
  all_exp_dt_cate$otherct <- factor(all_exp_dt_cate$otherct,levels = rev(c('meso','cc','epidermis')))
  all_exp_dt_cate <- all_exp_dt_cate[order(all_exp_dt_cate$targetct),]
  
  right_level <- c('meso_meso','meso_cc','meso_epidermis','cc_cc','cc_meso','','cc_epidermis','epidermis_epidermis','epidermis_meso','epidermis_cc')
  
  all_exp_dt$cate <- factor(all_exp_dt$cate,levels =right_level)
  
  
  p <- ggplot(all_exp_dt, aes(x=cate,y=exp,fill=targetct)) + 
    #geom_violin(trim=FALSE) + 
    #geom_boxplot(position=position_dodge(1)) +
    geom_boxplot()+
    #geom_boxplot(width=0.5,fill='white') +
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Expression\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,angle = 90,vjust = 0.5,hjust =0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) 
    #scale_fill_brewer(palette="Blues") +
    #facet_wrap(~targetct,nrow = 1)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)
  
  pdf(paste0(opt_dir,'/opt_three_target_ct_exp_comp.pdf'),width = 10,height = 15)
  p
  dev.off()

  ##we will next plot exp with and without H3K27me3
  ##for the epidermis
  merged_dt <- merge(ipt_dt,ipt_summary_H3K27me3_dt,by.x = 'rice_ACR',by.y = 'rice_ACR')
  merged_dt_overlapH3K27me3 <- merged_dt[merged_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  merged_dt_not_overlapH3K27me3 <- merged_dt[merged_dt$riceOverH3K27Cate == 'none',]

  ##epidermis
  merged_dt_overlapH3K27me3_epi <- merged_dt_overlapH3K27me3[merged_dt_overlapH3K27me3$rice_celltype == 'epidermis' |merged_dt_overlapH3K27me3$rice_celltype == 'protoderm' | merged_dt_overlapH3K27me3$rice_celltype == 'epidermis,protoderm',]
  merged_dt_overlapH3K27me3_epi <- merged_dt_overlapH3K27me3_epi[c('rice_geneID','seedling.Epidermis')]
  merged_dt_overlapH3K27me3_epi <- unique(merged_dt_overlapH3K27me3_epi[merged_dt_overlapH3K27me3_epi$seedling.Epidermis != 'none',])
  merged_dt_overlapH3K27me3_epi$ct <- c('epidermis')
  merged_dt_overlapH3K27me3_epi$h3cate <- c('H3K27')
  merged_dt_overlapH3K27me3_epi$cate <- c('CT')
  colnames(merged_dt_overlapH3K27me3_epi) <- c('gene','exp','ct','h3cate','cate')
  ##mesophyll
  merged_dt_overlapH3K27me3_meso <- merged_dt_overlapH3K27me3[grepl('mesophyll',merged_dt_overlapH3K27me3$rice_celltype),]
  merged_dt_overlapH3K27me3_meso <- merged_dt_overlapH3K27me3_meso[c('rice_geneID','seedling.Mesophyll')]
  merged_dt_overlapH3K27me3_meso <- unique(merged_dt_overlapH3K27me3_meso[merged_dt_overlapH3K27me3_meso$seedling.Mesophyll != 'none',])
  merged_dt_overlapH3K27me3_meso$ct <- c('meso')
  merged_dt_overlapH3K27me3_meso$h3cate <- c('H3K27')
  merged_dt_overlapH3K27me3_meso$cate <- c('CT')
  colnames(merged_dt_overlapH3K27me3_meso) <- c('gene','exp','ct','h3cate','cate')
  ##cc
  merged_dt_overlapH3K27me3_cc <- merged_dt_overlapH3K27me3[grepl('companion_cell',merged_dt_overlapH3K27me3$rice_celltype),]
  merged_dt_overlapH3K27me3_cc <- merged_dt_overlapH3K27me3_cc[c('rice_geneID','seedling.CompanionCell')]
  merged_dt_overlapH3K27me3_cc <- unique(merged_dt_overlapH3K27me3_cc[merged_dt_overlapH3K27me3_cc$seedling.CompanionCell != 'none',])
  merged_dt_overlapH3K27me3_cc$ct <- c('cc')
  merged_dt_overlapH3K27me3_cc$h3cate <- c('H3K27')
  merged_dt_overlapH3K27me3_cc$cate <- c('CT')
  colnames(merged_dt_overlapH3K27me3_cc) <- c('gene','exp','ct','h3cate','cate')
  
  ##epidermis
  merged_dt_not_overlapH3K27me3_epi <- merged_dt_not_overlapH3K27me3[merged_dt_not_overlapH3K27me3$rice_celltype == 'epidermis' |merged_dt_not_overlapH3K27me3$rice_celltype == 'protoderm' | merged_dt_not_overlapH3K27me3$rice_celltype == 'epidermis,protoderm',]
  merged_dt_not_overlapH3K27me3_epi <- merged_dt_not_overlapH3K27me3_epi[c('rice_geneID','seedling.Epidermis')]
  merged_dt_not_overlapH3K27me3_epi <- unique(merged_dt_not_overlapH3K27me3_epi[merged_dt_not_overlapH3K27me3_epi$seedling.Epidermis != 'none',])
  merged_dt_not_overlapH3K27me3_epi$ct <- c('epidermis')
  merged_dt_not_overlapH3K27me3_epi$h3cate <- c('noH3K27')
  merged_dt_not_overlapH3K27me3_epi$cate <- c('CT')
  colnames(merged_dt_not_overlapH3K27me3_epi) <- c('gene','exp','ct','h3cate','cate')
  ##mesophyll
  merged_dt_not_overlapH3K27me3_meso <- merged_dt_not_overlapH3K27me3[grepl('mesophyll',merged_dt_not_overlapH3K27me3$rice_celltype),]
  merged_dt_not_overlapH3K27me3_meso <- merged_dt_not_overlapH3K27me3_meso[c('rice_geneID','seedling.Mesophyll')]
  merged_dt_not_overlapH3K27me3_meso <- unique(merged_dt_not_overlapH3K27me3_meso[merged_dt_not_overlapH3K27me3_meso$seedling.Mesophyll != 'none',])
  merged_dt_not_overlapH3K27me3_meso$ct <- c('meso')
  merged_dt_not_overlapH3K27me3_meso$h3cate <- c('noH3K27')
  merged_dt_not_overlapH3K27me3_meso$cate <- c('CT')
  colnames(merged_dt_not_overlapH3K27me3_meso) <- c('gene','exp','ct','h3cate','cate')
  ##cc
  merged_dt_not_overlapH3K27me3_cc <- merged_dt_not_overlapH3K27me3[grepl('companion_cell',merged_dt_not_overlapH3K27me3$rice_celltype),]
  merged_dt_not_overlapH3K27me3_cc <- merged_dt_not_overlapH3K27me3_cc[c('rice_geneID','seedling.CompanionCell')]
  merged_dt_not_overlapH3K27me3_cc <- unique(merged_dt_not_overlapH3K27me3_cc[merged_dt_not_overlapH3K27me3_cc$seedling.CompanionCell != 'none',])
  merged_dt_not_overlapH3K27me3_cc$ct <- c('cc')
  merged_dt_not_overlapH3K27me3_cc$h3cate <- c('noH3K27')
  merged_dt_not_overlapH3K27me3_cc$cate <- c('CT')
  colnames(merged_dt_not_overlapH3K27me3_cc) <- c('gene','exp','ct','h3cate','cate')
  
  #sig test
  epdiermis_H3K27me3 <- t.test(as.numeric(merged_dt_overlapH3K27me3_epi$exp),as.numeric(merged_dt_not_overlapH3K27me3_epi$exp),alternative = 'less')
  meso_H3K27me3 <- t.test(as.numeric(merged_dt_overlapH3K27me3_meso$exp),as.numeric(merged_dt_not_overlapH3K27me3_meso$exp),alternative = 'less')
  cc_H3K27me3 <- t.test(as.numeric(merged_dt_overlapH3K27me3_cc$exp),as.numeric(merged_dt_not_overlapH3K27me3_cc$exp),alternative = 'less')
  
  broad_H3K27me3 <- merged_dt[merged_dt$rice_celltype == 'broadly_accessible'&merged_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  broad_not_H3K27me3 <- merged_dt[merged_dt$rice_celltype == 'broadly_accessible'&merged_dt$riceOverH3K27Cate == 'none',]
  
  ##for the broad ACR
  #epidermis
  broad_H3K27me3_epi <- broad_H3K27me3[c('rice_geneID','seedling.Epidermis')]
  broad_H3K27me3_epi <- unique(broad_H3K27me3_epi[broad_H3K27me3_epi$seedling.Epidermis != 'none',])
  broad_H3K27me3_epi$ct <- c('epidermis')
  broad_H3K27me3_epi$h3cate <- c('H3K27')
  broad_H3K27me3_epi$cate <- c('Broad')
  colnames(broad_H3K27me3_epi) <- c('gene','exp','ct','h3cate','cate')
  #meso
  broad_H3K27me3_meso <- broad_H3K27me3[c('rice_geneID','seedling.Mesophyll')]
  broad_H3K27me3_meso <- unique(broad_H3K27me3_meso[broad_H3K27me3_meso$seedling.Mesophyll != 'none',])
  broad_H3K27me3_meso$ct <- c('meso')
  broad_H3K27me3_meso$h3cate <- c('H3K27')
  broad_H3K27me3_meso$cate <- c('Broad')
  colnames(broad_H3K27me3_meso) <- c('gene','exp','ct','h3cate','cate')
  #cc
  broad_H3K27me3_cc <- broad_H3K27me3[c('rice_geneID','seedling.CompanionCell')]
  broad_H3K27me3_cc <- unique(broad_H3K27me3_cc[broad_H3K27me3_cc$seedling.CompanionCell != 'none',])
  broad_H3K27me3_cc$ct <- c('cc')
  broad_H3K27me3_cc$h3cate <- c('H3K27')
  broad_H3K27me3_cc$cate <- c('Broad')
  colnames(broad_H3K27me3_cc) <- c('gene','exp','ct','h3cate','cate')
  
  #epidermis
  broad_not_H3K27me3_epi <- broad_not_H3K27me3[c('rice_geneID','seedling.Epidermis')]
  broad_not_H3K27me3_epi <- unique(broad_not_H3K27me3_epi[broad_not_H3K27me3_epi$seedling.Epidermis != 'none',])
  broad_not_H3K27me3_epi$ct <- c('epidermis')
  broad_not_H3K27me3_epi$h3cate <- c('noH3K27')
  broad_not_H3K27me3_epi$cate <- c('Broad')
  colnames(broad_not_H3K27me3_epi) <- c('gene','exp','ct','h3cate','cate')
  #meso
  broad_not_H3K27me3_meso <- broad_not_H3K27me3[c('rice_geneID','seedling.Mesophyll')]
  broad_not_H3K27me3_meso <- unique(broad_not_H3K27me3_meso[broad_not_H3K27me3_meso$seedling.Mesophyll != 'none',])
  broad_not_H3K27me3_meso$ct <- c('meso')
  broad_not_H3K27me3_meso$h3cate <- c('noH3K27')
  broad_not_H3K27me3_meso$cate <- c('Broad')
  colnames(broad_not_H3K27me3_meso) <- c('gene','exp','ct','h3cate','cate')
  #cc
  broad_not_H3K27me3_cc <- broad_not_H3K27me3[c('rice_geneID','seedling.CompanionCell')]
  broad_not_H3K27me3_cc <- unique(broad_not_H3K27me3_cc[broad_not_H3K27me3_cc$seedling.CompanionCell != 'none',])
  broad_not_H3K27me3_cc$ct <- c('cc')
  broad_not_H3K27me3_cc$h3cate <- c('noH3K27')
  broad_not_H3K27me3_cc$cate <- c('Broad')
  colnames(broad_not_H3K27me3_cc) <- c('gene','exp','ct','h3cate','cate')
  
  broad_epdiermis_H3K27me3 <- t.test(as.numeric(broad_H3K27me3_epi$exp),as.numeric(broad_not_H3K27me3_epi$exp),alternative = 'less')
  broad_meso_H3K27me3 <- t.test(as.numeric(broad_H3K27me3_meso$exp),as.numeric(broad_not_H3K27me3_meso$exp),alternative = 'less')
  broad_cc_H3K27me3 <- t.test(as.numeric(broad_H3K27me3_cc$exp),as.numeric(broad_not_H3K27me3_cc$exp),alternative = 'less')
  
  ##merge
  broad_mergedt_dt <- rbind(merged_dt_overlapH3K27me3_epi,merged_dt_overlapH3K27me3_meso,merged_dt_overlapH3K27me3_cc,
                            merged_dt_not_overlapH3K27me3_epi,merged_dt_not_overlapH3K27me3_meso,merged_dt_not_overlapH3K27me3_cc,
                            broad_H3K27me3_epi,broad_H3K27me3_meso,broad_H3K27me3_cc,
                            broad_not_H3K27me3_epi,broad_not_H3K27me3_meso,broad_not_H3K27me3_cc)
  
  broad_mergedt_dt_test <- broad_mergedt_dt[broad_mergedt_dt$cate == 'Broad',] 
  
  head(broad_mergedt_dt)
  

  broad_mergedt_dt$h3cate <- as.factor(broad_mergedt_dt$h3cate)
  broad_mergedt_dt$exp <- as.numeric(broad_mergedt_dt$exp)
  broad_mergedt_dt$cate <- factor(broad_mergedt_dt$cate,levels = c('CT','Broad'))
  broad_mergedt_dt$ct <- factor(broad_mergedt_dt$ct,levels = c('meso','cc','epidermis'))
  p <- ggplot(broad_mergedt_dt, aes(x=ct,y=exp,fill=h3cate)) + 
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.2,fatten = 2,position=position_dodge(0.9)) +
    #geom_boxplot(width=0.5,fill='white') +
    
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    scale_fill_brewer(palette="Blues") +
    facet_wrap(~cate)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)
  
  pdf(paste0(opt_dir,'/opt_three_target_ct_exp_comp_H3K27me3.pdf'),width = 12,height = 6)
  p
  dev.off()
  
  
  
  

})

##updating 010924
##we will check the exp patterns diff for the species specific ACRs under the syntenic and non syntenic regions
ipt_dir <- 'opt2_nearby_genes_check_exp_122723/'
opt_dir <- 'opt2_nearby_genes_check_exp_results_122723/'

all_fl_list <- list.files(path = ipt_dir,pattern = 'opt_rice_maize.+collapes.txt')


#check total of species-specific ACRs
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr_collapes.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = T)
  head(ipt_dt)
  
  ipt_syn_dt <- ipt_dt[ipt_dt$Syntenic_regionID_update != 'none',]
  ipt_syn_SS_dt <- ipt_syn_dt[ipt_syn_dt$maize_region == 'none',]
  
  ipt_non_syn_dt <- ipt_dt[ipt_dt$Syntenic_regionID_update == 'none',]
  
  ipt_synSS_nonsyn_dt <- rbind(ipt_syn_SS_dt,ipt_non_syn_dt)
  
  ##check high exp epidermis gene
  ipt_synSS_nonsyn_dt_highEpi <- ipt_synSS_nonsyn_dt[ipt_synSS_nonsyn_dt$rice_highestCT == 'seedling.Epidermis' & 
                                                       ipt_synSS_nonsyn_dt$maize_highestCT != 'protoderm',]
  dim(ipt_synSS_nonsyn_dt_highEpi)
  
  ipt_synSS_nonsyn_dt_highEpi_acrEpiProto <- ipt_synSS_nonsyn_dt_highEpi[grepl('epidermis',ipt_synSS_nonsyn_dt_highEpi$rice_CelltypeCate)|
                                                                                 grepl('protoderm',ipt_synSS_nonsyn_dt_highEpi$rice_CelltypeCate),]
  
  ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_ACR','geneID','rice_seedling.Epidermis','rice_seedling.Mesophyll',
                                                                                              'rice_seedling.CompanionCell','maize_ortho','maize_companion_cells_sieve_elements','maize_mesophyll','maize_protoderm')]
  write.csv(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target,paste0(opt_dir,'/opt_final_table_rice_to_maize_geneExp.csv'),quote = F)
  
  length(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$geneID))
  write.csv(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$geneID,paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto.csv'),quote = F)
  write.csv(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$maize_ortho),paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto_maizeOthor.csv'),quote = F)
  ##87
  length(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$geneID))
  
  length(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$rice_ACR))
  ##102
  write.csv(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$rice_ACR),paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto_ACRlist.csv'),quote = F)
  
  
  ##upating 013124
  ##we will check the H3K27me3 
  ipt_H3K27me3_dt <- read.delim(paste0(ipt_dir,'/opt_acr_in_and_not_in_H3K27peak_add_ct_addnum_addflank150bp_addcombineCate_addGeneRatio_addGeneID.txt'),header = F)
  head(ipt_H3K27me3_dt)
  ipt_H3K27me3_dt$geneID_ACR <- paste0(ipt_H3K27me3_dt$V11,'__',ipt_H3K27me3_dt$V1,'_',ipt_H3K27me3_dt$V2,'_',ipt_H3K27me3_dt$V3)

  ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target_col <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('geneID','rice_InorNotH3K27m3','rice_ACR')]
  ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target_col$geneID_ACR <- paste0(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target_col$geneID,'__',ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target_col$rice_ACR)
  
  
  ipt_H3K27me3_merge_all_dt <- merge(ipt_H3K27me3_dt,ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target_col,by.x = 'geneID_ACR',by.y = 'geneID_ACR',all.x = TRUE)
  dim(ipt_H3K27me3_merge_all_dt)
  
  write.table(ipt_H3K27me3_merge_all_dt,paste0(opt_dir,'/opt_H3K27me3_target_87_rice_to_maize_genes_ACRs.txt'),quote = F,sep = '\t')
  
  ##we will extract the whole gene information
  ipt_H3K27me3_merge_dt <- merge(ipt_H3K27me3_dt,ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target_col,by.x = 'geneID_ACR',by.y = 'geneID_ACR')
  dim(ipt_H3K27me3_merge_dt)
  alltarget_genes <- unique(ipt_H3K27me3_merge_dt$geneID)
  
  ipt_H3K27me3_merge_all_target_genes_dt <- ipt_H3K27me3_merge_all_dt[ipt_H3K27me3_merge_all_dt$V11 %in% alltarget_genes,]
  write.table(ipt_H3K27me3_merge_all_target_genes_dt,paste0(opt_dir,'/opt_H3K27me3_target_87_rice_to_maize_genes_ACRs_filter_genes.txt'),quote = F,sep = '\t')
  
  ##we would add the gene annotaiton of these 87 genes
  ipt_rice_symbol_dt <- read.delim(paste0(ipt_dir,'/opt_rice_symbol.txt'),header = F)
  head(ipt_rice_symbol_dt)
  
  ipt_rice_annot_dt <- read.delim(paste0(ipt_dir,'/Osativa_323_v7.0.annotation_info_concise.txt'),header = T)
  head(ipt_rice_annot_dt)
  
  ipt_H3K27me3_merge_all_target_genes_symbol_dt <- merge(ipt_H3K27me3_merge_all_target_genes_dt,ipt_rice_symbol_dt,by.x = 'V11',by.y = 'V1', all.x = TRUE)
  
  ipt_H3K27me3_merge_all_target_genes_symbol_annot_dt <- merge(ipt_H3K27me3_merge_all_target_genes_symbol_dt,ipt_rice_annot_dt,by.x = 'V11',by.y = 'locusName', all.x = TRUE)
  
  write.table(ipt_H3K27me3_merge_all_target_genes_symbol_annot_dt,paste0(opt_dir,'/opt_H3K27me3_target_87_rice_to_maize_genes_ACRs_filter_genes_add_Annot.txt'),quote = F,sep = '\t')
  
  
  
  ##check the expression boxplot
  ##for the rice
  rice_epi <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_seedling.Epidermis')]
  colnames(rice_epi) <- c('exp')
  rice_epi$ct <- 'epidermis'
  rice_meso <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_seedling.Mesophyll')]
  colnames(rice_meso) <- c('exp')
  rice_meso$ct <- 'meso'
  rice_cc <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_seedling.CompanionCell')]
  colnames(rice_cc) <- c('exp')
  rice_cc$ct <- 'cc'
  rice_all <- rbind(rice_epi,rice_meso,rice_cc)
  rice_all$exp <- as.numeric(rice_all$exp)
  rice_all$ct <- factor(rice_all$ct,levels = c('epidermis','meso','cc'))
  p <- ggplot(rice_all, aes(x=ct,y=exp,fill=ct)) + 
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.2,fatten = 2,position=position_dodge(0.9)) +
    #geom_boxplot(width=0.5,fill='white') +
    
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    scale_fill_brewer(palette="Blues") 
    #facet_wrap(~cate)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)
  pdf(paste0(opt_dir,'/opt_rice_',target_fl_path_nm,'_synSS_Nsyn_highEpi_acrEpiProto_geneExp_rice.pdf'),width = 8,height = 8)
  p
  dev.off()
  
  ##check the sig test
  rice_epi <- rice_all[rice_all$ct == 'epidermis',]$exp
  rice_meso <- rice_all[rice_all$ct == 'meso',]$exp
  rice_cc <- rice_all[rice_all$ct == 'cc',]$exp
  t.test(rice_epi,rice_meso,alternative = "greater")
  t.test(rice_epi,rice_cc,alternative = 'greater')
  
  ##for the maize
  maize_epi <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_protoderm')]
  colnames(maize_epi) <- c('exp')
  maize_epi$ct <- 'epidermis'
  maize_meso <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_mesophyll')]
  colnames(maize_meso) <- c('exp')
  maize_meso$ct <- 'meso'
  maize_cc <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_companion_cells_sieve_elements')]
  colnames(maize_cc) <- c('exp')
  maize_cc$ct <- 'cc'
  maize_all <- rbind(maize_epi,maize_meso,maize_cc)
  maize_all$exp <- as.numeric(maize_all$exp)
  maize_all$ct <- factor(maize_all$ct,levels = c('epidermis','meso','cc'))
  
  maize_epi <- maize_all[maize_all$ct == 'epidermis',]$exp
  maize_meso <- maize_all[maize_all$ct == 'meso',]$exp
  maize_cc <- maize_all[maize_all$ct == 'cc',]$exp
  t.test(maize_meso,maize_epi,alternative = "greater")
  t.test(maize_cc,maize_epi,alternative = "greater")
  
  
  p <- ggplot(maize_all, aes(x=ct,y=exp,fill=ct)) + 
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.2,fatten = 2,position=position_dodge(0.9)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    scale_fill_brewer(palette="Blues") 
  pdf(paste0(opt_dir,'/opt_rice_',target_fl_path_nm,'_synSS_Nsyn_highEpi_acrEpiProto_geneExp_maize.pdf'),width = 8,height = 8)
  p
  dev.off()
  
})


##updating 011724
##for the maize to rice

all_fl_list <- list.files(path = ipt_dir,pattern = 'opt_maize_rice.+collapes.txt')

#check total of species-specific ACRs
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr_collapes.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_maize_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = T)
  head(ipt_dt)
  
  ipt_syn_dt <- ipt_dt[ipt_dt$Syntenic_regionID_update != 'none',]
  ipt_syn_SS_dt <- ipt_syn_dt[ipt_syn_dt$rice_region == 'none',]
  
  ipt_non_syn_dt <- ipt_dt[ipt_dt$Syntenic_regionID_update == 'none',]
  
  ipt_synSS_nonsyn_dt <- rbind(ipt_syn_SS_dt,ipt_non_syn_dt)
  table( ipt_synSS_nonsyn_dt$rice_highestCT)
  ##check high exp epidermis gene
  
  ipt_synSS_nonsyn_dt_highEpi <- ipt_synSS_nonsyn_dt[ipt_synSS_nonsyn_dt$maize_highestCT == 'protoderm' & 
                                                       ipt_synSS_nonsyn_dt$rice_highestCT != 'seedling.Epidermis',]
  dim(ipt_synSS_nonsyn_dt_highEpi)
  
  table(ipt_synSS_nonsyn_dt_highEpi$maize_CelltypeCat)
  
  ipt_synSS_nonsyn_dt_highEpi_acrEpiProto <- ipt_synSS_nonsyn_dt_highEpi[grepl('epidermis',ipt_synSS_nonsyn_dt_highEpi$maize_CelltypeCate)|
                                                                           grepl('protoderm',ipt_synSS_nonsyn_dt_highEpi$maize_CelltypeCate),]
  
  ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_ACR','geneID','rice_seedling.Epidermis','rice_seedling.Mesophyll',
                                                                                              'rice_seedling.CompanionCell','rice_ortho','maize_companion_cells_sieve_elements','maize_mesophyll','maize_protoderm')]
  write.csv(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto_target,paste0(opt_dir,'/opt_final_table_maize_to_rice_geneExp.csv'),quote = F)
  
  write.csv(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$geneID,paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto_maizeTorice_maizeGenes.csv'),quote = F)
  ##102
  write.csv(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$rice_ortho),paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto_maizeTorice_maizeGenes_riceOthor.csv'),quote = F)
  
  length(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$geneID))
  
  length(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$maize_ACR))
  write.csv(unique(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$maize_ACR),paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto_maizeTorice_maizeGenes_ACRlist.csv'),quote = F)
  
  geneID_dt <- read.delim(paste0(ipt_dir,'/geneIDs_v3_v4_v5.txt'),header = F)
  target_maize_ID_dt <- as.data.frame(ipt_synSS_nonsyn_dt_highEpi_acrEpiProto$geneID)
  colnames(target_maize_ID_dt) <- c('geneID')
  
  mergedt_gene_dt <- merge(geneID_dt,unique(target_maize_ID_dt),by.x = 'V3',by.y = 'geneID')
  write.csv(mergedt_gene_dt$V2,paste0(opt_dir,'/temp_synSS_nonsyn_dt_highEpi_acrEpiProto_maizeTorice_maizeGenes_V4version.csv'),quote = F)
  
  ##check the expression boxplot
  ##for the maize
  rice_epi <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_protoderm')]
  colnames(rice_epi) <- c('exp')
  rice_epi$ct <- 'epidermis'
  rice_meso <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_mesophyll')]
  colnames(rice_meso) <- c('exp')
  rice_meso$ct <- 'meso'
  rice_cc <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('maize_companion_cells_sieve_elements')]
  colnames(rice_cc) <- c('exp')
  rice_cc$ct <- 'cc'
  rice_all <- rbind(rice_epi,rice_meso,rice_cc)
  rice_all$exp <- as.numeric(rice_all$exp)
  rice_all$ct <- factor(rice_all$ct,levels = c('epidermis','meso','cc'))
  p <- ggplot(rice_all, aes(x=ct,y=exp,fill=ct)) + 
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.2,fatten = 2,position=position_dodge(0.9)) +
    #geom_boxplot(width=0.5,fill='white') +
    
    #geom_jitter(shape=6, position=position_jitter(0.2)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    scale_fill_brewer(palette="Blues") 
  #facet_wrap(~cate)## You can use face_wrap function only if you need it+
  #geom_text(data =tukey_letters, 
  #          aes(x=xpos, y=ymax+offset_asterisk,label=groups), 
  #          size = 8,position=position_dodge(.5)  ##change the letter
  #)
  pdf(paste0(opt_dir,'/opt_maize_',target_fl_path_nm,'_synSS_Nsyn_highEpi_acrEpiProto_geneExp_maize.pdf'),width = 8,height = 8)
  p
  dev.off()
  
  ##check the sig test
  rice_epi <- rice_all[rice_all$ct == 'epidermis',]$exp
  rice_meso <- rice_all[rice_all$ct == 'meso',]$exp
  rice_cc <- rice_all[rice_all$ct == 'cc',]$exp
  t.test(rice_epi,rice_meso,alternative = "greater") ##0.0002543
  t.test(rice_epi,rice_cc,alternative = 'greater') ##0.0001733
  
  ##for the rice
  maize_epi <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_seedling.Epidermis')]
  colnames(maize_epi) <- c('exp')
  maize_epi$ct <- 'epidermis'
  maize_meso <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_seedling.Mesophyll')]
  colnames(maize_meso) <- c('exp')
  maize_meso$ct <- 'meso'
  maize_cc <- ipt_synSS_nonsyn_dt_highEpi_acrEpiProto[c('rice_seedling.CompanionCell')]
  colnames(maize_cc) <- c('exp')
  maize_cc$ct <- 'cc'
  maize_all <- rbind(maize_epi,maize_meso,maize_cc)
  maize_all$exp <- as.numeric(maize_all$exp)
  maize_all$ct <- factor(maize_all$ct,levels = c('epidermis','meso','cc'))
  
  maize_epi <- maize_all[maize_all$ct == 'epidermis',]$exp
  maize_meso <- maize_all[maize_all$ct == 'meso',]$exp
  maize_cc <- maize_all[maize_all$ct == 'cc',]$exp
  t.test(maize_meso,maize_epi,alternative = "greater") ##0.2315
  t.test(maize_cc,maize_epi,alternative = "greater") ##0.3807
  
  
  p <- ggplot(maize_all, aes(x=ct,y=exp,fill=ct)) + 
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.2,fatten = 2,position=position_dodge(0.9)) +
    labs(x="\nClass", y = 'Correlation\n')+
    theme(
      plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
      axis.title = element_text(size =35),
      axis.text.x = element_text(colour = "black", size=30,hjust = 0.5),  ##change the text to italic
      axis.text.y = element_text(size=30,colour = "black"),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    scale_fill_brewer(palette="Blues") 
  pdf(paste0(opt_dir,'/opt_maize_',target_fl_path_nm,'_synSS_Nsyn_highEpi_acrEpiProto_geneExp_rice.pdf'),width = 8,height = 8)
  p
  dev.off()
  
})












##updating 011024
##we will check the enrichment of motifs
ipt_dir <- 'opt2_nearby_genes_check_exp_122723/'
opt_dir <- 'opt2_nearby_genes_check_exp_results_122723/'

##this is for the motif enrichment from rice to maize
combine_dt <- read.delim(paste0(ipt_dir,'/opt_final_cate_motif_enrichment.txt'), header = F)
head(combine_dt)
opt_prefix <- 'rice_to_maize'

##this is for the motif enrichment from maize to rice
combine_dt <- read.delim(paste0(ipt_dir,'/opt_final_cate_motif_enrichment_maize_to_rice.txt'), header = F)
head(combine_dt)
opt_prefix <- 'maize_to_rice'


colnames(combine_dt) <- c('TACR','celltype','motifID','pval','num','totalnum','avgprop')
combine_dt <- combine_dt[c('celltype','motifID','pval')]

meta_name_dt <- read.delim(paste0(ipt_dir,'/ipt_required_raw_data/opt_motif_common_name.txt'),header = F)
head(meta_name_dt)

##add the common name
combine_dt <- merge(combine_dt,meta_name_dt,by.x='motifID',by.y= 'V1')
head(combine_dt)
#combine_dt <- combine_dt[order(combine_dt$qval),]

##only for the maize to rice
combine_less005_dt <- combine_dt[combine_dt$pval < 0.35,]
combine_less005_dt$qval <- p.adjust(combine_less005_dt$pval, method="fdr")
combine_dt <- combine_less005_dt

##for rice to maize
combine_dt$qval <- p.adjust(combine_dt$pval, method="fdr")

combine_dt$mlog10qval <- -log10(combine_dt$qval)
combine_dt$mlog10pval <- -log10(combine_dt$pval)




ipt_two_kinds_of_values <- c('pval','qval')

outs2 <- lapply(ipt_two_kinds_of_values, function(sig){
  
  print(sig)
  if (sig == 'pval'){
    
    combine_dt <- combine_dt[combine_dt$pval < 0.05,]
    
    combine_dt <- combine_dt[c('V2','celltype','mlog10pval')]
    colnames(combine_dt) <- c('motif','celltype','mlog10pval')
    
    write.table(combine_dt,paste0(opt_dir,'/temp_sparse.txt'),quote = F,sep = '\t')
    combine_dt <- read.table(paste0(opt_dir,'/temp_sparse.txt'),stringsAsFactors = T)
    score_mtx <- sparseMatrix(i=as.numeric(combine_dt$motif),
                              j=as.numeric(combine_dt$celltype),
                              x=as.numeric(combine_dt$mlog10pval),
                              dimnames=list(levels(combine_dt$motif), levels(combine_dt$celltype)))
    score_mtx <- as.matrix(score_mtx)
    dim(score_mtx)
    head(score_mtx)
    max(score_mtx)
    
  }else{
    
    combine_dt <- combine_dt[combine_dt$qval < 0.05,]
    
    write.csv(combine_dt,paste0(opt_dir,'/temp_',opt_prefix,'_motif_enrichment_qval.csv'),quote = F)
    
    combine_dt <- combine_dt[c('V2','celltype','mlog10qval')]
    colnames(combine_dt) <- c('motif','celltype','mlog10qval')
    
    
    write.table(combine_dt,paste0(opt_dir,'/temp_sparse.txt'),quote = F,sep = '\t')
    combine_dt <- read.table(paste0(opt_dir,'/temp_sparse.txt'),stringsAsFactors = T)
    score_mtx <- sparseMatrix(i=as.numeric(combine_dt$motif),
                              j=as.numeric(combine_dt$celltype),
                              x=as.numeric(combine_dt$mlog10qval),
                              dimnames=list(levels(combine_dt$motif), levels(combine_dt$celltype)))
    score_mtx <- as.matrix(score_mtx)
    dim(score_mtx)
    head(score_mtx)
    max(score_mtx)
  }

  
  
  
  ##reorder the matrix
  colnames(score_mtx)
  #target_order <- c('rice_maize__mesophyll','rice_sorghum__mesophyll','rice_Pm__mesophyll','rice_Uf__mesophyll',
  #                  'rice_maize__bundle_sheath','rice_sorghum__bundle_sheath','rice_Pm__bundle_sheath','rice_Uf__bundle_sheath',
  #                  'rice_maize__companion_cell','rice_sorghum__companion_cell','rice_Pm__companion_cell','rice_Uf__companion_cell',
  #                  'rice_maize__protoderm','rice_sorghum__protoderm','rice_Pm__protoderm','rice_Uf__protoderm',
  #                  'rice_maize__epidermis','rice_sorghum__epidermis','rice_Pm__epidermis','rice_Uf__epidermis'
  #)
  
  #intersect_colnm <- intersect(target_order,colnames(score_mtx))
  
  #score_mtx <- score_mtx[,intersect_colnm]
  
  # Replace Inf with the maximum value of the matrix
  max_non_inf <- max(score_mtx[!is.infinite(score_mtx)], na.rm = TRUE)
  score_mtx[is.infinite(score_mtx)] <- max_non_inf
  max(score_mtx)
  
  p <- pheatmap(score_mtx,
                #scale="row",
                #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
                #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
                cluster_cols = F,
                cluster_rows = F,
                border_color = "black",
                #border_color='black',
                fontsize_col = 8,
                fontsize_row = 8,
                #annotation_row = annotdf,
                #annotation_colors = mycolors,
                #annotation_col = annotdf_col,
                #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
                #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
                na_col = "white",
                #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
                #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
                #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
                #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
                color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
                
                #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
                #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
                #color=colorRampPalette(viridis)(100),
                
                ##second candidate
                #color = viridis(n = 256, alpha = 1, 
                #                     begin = 0, end = 1, option = "viridis"),
                
                #annotation_row = annotdf,
                #annotation_colors = mycolors,
                
                #annotation_col = annotdf_col,
                #annotation_colors = mycolors_col,
                #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
                #color = inferno(length(mat_breaks) - 1),
                #breaks = mat_breaks,
                
                ##open it when necessary
                #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
                fontsize_number = 25
  ) 
  pdf(paste0(opt_dir,'/opt_synSSorNsyn_higestCT_acrEpiProto_geneExp_motifEnrich_',sig,'_',opt_prefix,'.pdf'),width = 6,height = 8)
  print(p)
  dev.off()
  
  ##build the dot plot for each cell type 
  combine_dt
  celltype_list <- unique(combine_dt$celltype)
  
  outs <- lapply(celltype_list, function(x){
    
    target_dt <- combine_dt[combine_dt$celltype == x,]
    
    target_dt <- target_dt[order(target_dt[[paste0('mlog10',sig)]],decreasing = T),]
    target_dt$motif <- factor(target_dt$motif ,levels = target_dt$motif)
    
    if (sig == 'pval'){
      p <- ggplot(target_dt, aes(motif,mlog10pval, color = celltype, fill = celltype,shape= celltype)) 
    }else{
      p <- ggplot(target_dt, aes(motif,mlog10qval, color = celltype, fill = celltype,shape= celltype)) 
    }
      
    p <-  p + geom_point(size = 7) +
      scale_color_manual(values=c("bundle_sheath"="#3084AF", "seedling.CompanionCell"="#eddc9f",
                                  "seedling.Epidermis"="#C49B7B","seedling.Mesophyll"="#91e388","protoderm"="#fcd2b1"))+
      scale_fill_manual(values=c("bundle_sheath"="#3084AF", "seedling.CompanionCell"="#eddc9f",
                                 "seedling.Epidermis"="#C49B7B","seedling.Mesophyll"="#91e388","protoderm"="#fcd2b1"))+
      scale_shape_manual(values = c("bundle_sheath" = 16, "seedling.CompanionCell" = 25,
                                    "seedling.Epidermis" = 18, "seedling.Mesophyll" = 15, "protoderm" = 17))+
      #geom_point(data = points_to_label, color = 'red', size = 3) +
      #geom_text(data = points_to_label, aes(label = rownames(points_to_label)), hjust = -0.2,color='red',cex= 7) + 
      theme(
        plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
        axis.title = element_text(size =20),
        axis.text.x = element_text(colour = "black", size=30,angle = 45, vjust = 1,hjust = 1,face = 'bold'),  ##change the text to italic
        #axis.text.x = element_text(size=20,colour = "black"),
        axis.text.y = element_text(size=20,colour = "black"),
        
        axis.ticks = element_line(size = rel(2.5)),
        axis.ticks.length = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20),
        #legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
      ) +
    #scale_x_continuous(breaks=seq(-2, 2, 0.5)) +
    #scale_y_continuous(breaks=seq(0, 2.5, 0.5)) +
    scale_y_continuous(breaks=seq(0, 7, 1)) +
    coord_cartesian(ylim = c(0,7))
 
   
    pdf(paste0(opt_dir,'/opt_',x,'_synSSorNsyn_higestCT_acrEpiProto_geneExp_motifEnrich_dotplot_',sig,'_',opt_prefix,'.pdf'),width = 8,height = 6)
    print(p)
    dev.off()
    
    
    
    
    
    
    
  })

})
















##########
##previous for the stress related genes does not work
##load other
ipt_target_gene_dt <- read.delim('ipt_nearby_gene_090323.txt',header = F)
ipt_target_genes <- ipt_target_gene_dt$V1

ipt_target_libID_groups_dt <- read.delim(paste0(ipt_dir,'/','ipt_target_libID_annot_add_groups_090323.txt'),header = T)
rownames(ipt_target_libID_groups_dt) <- ipt_target_libID_groups_dt$libID

ipt_target_lib_fpkm_dt <- read.delim(paste0(ipt_dir,'/','./opt_target_lib_fpkm_nochangeName.txt'),row.names = 1)
head(ipt_target_lib_fpkm_dt)


##
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  head(ipt_dt)

  ##select the target
  if (target_cate == 'not_shared'){
    
    ipt_target_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none',]
    unique(ipt_target_dt$geneID)
    colnames(ipt_target_dt)
    
    ipt_target_dt <- unique(ipt_target_dt[c('geneID','rice_CelltypeCate')])
  
  }
  
  ipt_target_genes <- ipt_target_dt$geneID
  
  ipt_target_genes <- gsub('Os','OS',ipt_target_genes)
  ipt_target_genes <- gsub('g','G',ipt_target_genes)
  
  ipt_target_dt$modiGene <- ipt_target_genes
  
  ipt_target_lib_fpkm_targetGene_dt <- ipt_target_lib_fpkm_dt[ipt_target_genes,]
  head(ipt_target_lib_fpkm_targetGene_dt)
  
  all_treats_name <- colnames(ipt_target_libID_groups_dt)
  all_treats_name <- all_treats_name[3:length(all_treats_name)]
  all_treats_name
  
  all_geneID <- ipt_target_genes
  
  outs2 <- lapply(all_geneID, function(x){
    
    message(x)
    
    outs <- lapply(all_treats_name, function(y){
      
      ##find two groups of name
      treat_IDs <- ipt_target_libID_groups_dt[ipt_target_libID_groups_dt[[y]] == paste0(y,'_treat'),]$libID
      control_IDs <- ipt_target_libID_groups_dt[ipt_target_libID_groups_dt[[y]] == paste0(y,'_control'),]$libID
      
      treat_exp <- as.numeric(ipt_target_lib_fpkm_targetGene_dt[x,treat_IDs])
      control_exp <- as.numeric(ipt_target_lib_fpkm_targetGene_dt[x,control_IDs])
      
      res = wilcox.test(treat_exp,control_exp,paired = FALSE)
      pval <- res$p.value
      
      res <- c()
      res$treat <- y
      res$pvalue <- pval
      
      return(res)    
      
    })
    
    combine_treat_dt <- as.data.frame(do.call(rbind,outs))
    
    combine_treat_dt$gene <- x
    
    return(combine_treat_dt)
    
  })
  
  combine_all_dt <- as.data.frame(do.call(rbind,outs2))
  combine_all_dt$response <- ifelse(combine_all_dt$pvalue < 0.05, 'yes','no')
  combine_all_dt[is.na(combine_all_dt$response), "response"] <- 'no'
  combine_all_dt$count <- 1
  head(combine_all_dt)
  
  ##add the cell type information
  combine_ct_gene <- merge(ipt_target_dt,combine_all_dt,by.x = 'modiGene',by.y = 'gene')
  
  unique(combine_ct_gene$rice_CelltypeCate)
  
  combine_ct_gene_dt <- as.data.frame(separate_rows(combine_ct_gene, rice_CelltypeCate, sep = ","))
  
  unique(combine_ct_gene_dt$rice_CelltypeCate)
  
  combine_ct_gene_dt <- combine_ct_gene_dt[combine_ct_gene_dt$pvalue != 'NaN',]
  combine_ct_gene_dt$pvalue <- as.numeric(combine_ct_gene_dt$pvalue)
  
  combine_ct_gene_dt <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate != 'unknown_cells_2',]
  combine_ct_gene_dt <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate != 'unknown_cells_1',]
  
  combine_ct_gene_dt_epidermis <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate == 'epidermis',]
  aggregate(count ~ response, data = combine_ct_gene_dt_epidermis, sum)
  
  combine_ct_gene_dt_mesophyll <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate == 'mesophyll',]
  aggregate(count ~ response, data = combine_ct_gene_dt_mesophyll, sum)
  
  combine_ct_gene_dt_protoderm <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate == 'protoderm',]
  aggregate(count ~ response, data = combine_ct_gene_dt_protoderm, sum)
  
  combine_ct_gene_dt_companion_cell <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate == 'companion_cell',]
  aggregate(count ~ response, data = combine_ct_gene_dt_companion_cell, sum)
  
  combine_ct_gene_dt_bundle_sheath <- combine_ct_gene_dt[combine_ct_gene_dt$rice_CelltypeCate == 'bundle_sheath',]
  aggregate(count ~ response, data = combine_ct_gene_dt_bundle_sheath, sum)
  
  aggregate(pvalue ~ response, data = combine_ct_gene_dt_epidermis, mean)
  aggregate(pvalue ~ response, data = combine_ct_gene_dt_mesophyll, mean)
  
  combine_ct_gene_dt <- combine_ct_gene_dt[complete.cases(combine_ct_gene_dt), ]
  
  
  
  
  
  
  
  class(combine_all_dt)
  combine_all_dt$treat <- as.factor(combine_all_dt$treat)
  table(combine_all_dt$treat)
  tail(combine_all_dt)
  combine_all_dt <- as.matrix(combine_all_dt)
  write.table(combine_all_dt,'opt_s2_gene_check_sig_diff_all_treat.txt',quote = F,sep = '\t')
  
  combine_all_dt <- read.table('opt_s2_gene_check_sig_diff_all_treat.txt')
  aggregate_dt <- aggregate(count ~ treat + response, data = combine_all_dt, sum)
  
  ##change to the matrix 
  combine_all_dt <- read.table('opt_s2_gene_check_sig_diff_all_treat.txt',stringsAsFactors = T)
  head(combine_all_dt)
  
  
})






##previous for chekcing the gene overlapping
combine_dt_sharedA <- do.call(rbind,outs)
maize_gene_list <- combine_dt_sharedA[combine_dt_sharedA$spe == 'maize',]$geneID
sorghum_gene_list <- combine_dt_sharedA[combine_dt_sharedA$spe == 'sorghum',]$geneID
Pm_gene_list <- combine_dt_sharedA[combine_dt_sharedA$spe == 'Pm',]$geneID
Uf_gene_list <- combine_dt_sharedA[combine_dt_sharedA$spe == 'Uf',]$geneID
maize_sorghum <- intersect(maize_gene_list,sorghum_gene_list)
maize_sorghum_Pm <- intersect(maize_sorghum,Pm_gene_list)
maize_sorghum_Pm_Uf <- intersect(maize_sorghum_Pm,Uf_gene_list)


##for the share I 
outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  head(ipt_dt)
  
  ##check the shared_acc genes
  ipt_dt_flt <- ipt_dt[ipt_dt$rice_ACRID != 'none',]
  ipt_dt_flt <- ipt_dt_flt[ipt_dt_flt[[paste0(target_fl_path_nm,'_ACR')]] == 'none',]
  
  shared_acc_genes <- unique(ipt_dt_flt$geneID)
  res <- c()
  res$geneID <- shared_acc_genes
  res$spe <- target_fl_path_nm
  res_dt <- as.data.frame(res)
  return(res_dt)
  
})

combine_dt_sharedI <- do.call(rbind,outs)
maize_gene_list <- combine_dt_sharedI[combine_dt_sharedI$spe == 'maize',]$geneID
sorghum_gene_list <- combine_dt_sharedI[combine_dt_sharedI$spe == 'sorghum',]$geneID
Pm_gene_list <- combine_dt_sharedI[combine_dt_sharedI$spe == 'Pm',]$geneID
Uf_gene_list <- combine_dt_sharedI[combine_dt_sharedI$spe == 'Uf',]$geneID
maize_sorghum <- intersect(maize_gene_list,sorghum_gene_list)
maize_sorghum_Pm <- intersect(maize_sorghum,Pm_gene_list)
maize_sorghum_Pm_Uf <- intersect(maize_sorghum_Pm,Uf_gene_list)




##updating 121323
##we will check if the species-specific ACRs were more likely in the distal compared to the others
ipt_dir <- 'opt2_nearby_genes_120623/////'
opt_dir <- 'opt2_nearby_genes_results_120623///'

all_fl_list <- list.files(path = ipt_dir,pattern = 'celltypeCpm_addGeneCate')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_final_blast_summary_add_gene_add_celltypeCpm_addGeneCate.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)

  ##for the not shared
  ipt_target_noshared_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none',]
  ipt_target_noshared_dt <- unique(ipt_target_noshared_dt[c('rice_ACR','geneCate','rice_CelltypeCate')])
  ipt_target_noshared_dt$celltypeCate <- ifelse(ipt_target_noshared_dt$rice_CelltypeCate == 'broadly_accessible', 'Broad','CT')
  ipt_target_noshared_dt$count <- 1
  aggregate_dt_notshared <- aggregate(count ~ geneCate + celltypeCate, data = ipt_target_noshared_dt, sum)
  aggregate_dt_notshared$cate <- 'not_shared'
  
  aggregate_celltypecate <- aggregate(count ~ celltypeCate, data = aggregate_dt_notshared, sum)
  combine_aggregate_notshared <- merge(aggregate_dt_notshared, aggregate_celltypecate, by.x = 'celltypeCate',by.y = 'celltypeCate')
  combine_aggregate_notshared$prop <- combine_aggregate_notshared$count.x/combine_aggregate_notshared$count.y
  
  ##for the shared I
  ipt_target_sharedI_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_fl_path_nm,'_ACR')]] == 'none',]
  ipt_target_sharedI_dt <- unique(ipt_target_sharedI_dt[c('rice_ACR','geneCate','rice_CelltypeCate')])
  ipt_target_sharedI_dt$celltypeCate <- ifelse(ipt_target_sharedI_dt$rice_CelltypeCate == 'broadly_accessible', 'Broad','CT')
  ipt_target_sharedI_dt$count <- 1
  aggregate_dt_sharedI <- aggregate(count ~ geneCate + celltypeCate, data = ipt_target_sharedI_dt, sum)
  aggregate_dt_sharedI$cate <- 'shared_I'
  
  aggregate_celltypecate <- aggregate(count ~ celltypeCate, data = aggregate_dt_sharedI, sum)
  combine_aggregate_sharedI <- merge(aggregate_dt_sharedI, aggregate_celltypecate, by.x = 'celltypeCate',by.y = 'celltypeCate')
  combine_aggregate_sharedI$prop <- combine_aggregate_sharedI$count.x/combine_aggregate_sharedI$count.y
  
  ##for the shared A
  ipt_target_sharedA_dt <- ipt_dt[ipt_dt$rice_ACR != 'none' & ipt_dt[[paste0(target_fl_path_nm,'_ACR')]] != 'none',]
  ipt_target_sharedA_dt <- unique(ipt_target_sharedA_dt[c('rice_ACR','geneCate','rice_CelltypeCate')])
  ipt_target_sharedA_dt$celltypeCate <- ifelse(ipt_target_sharedA_dt$rice_CelltypeCate == 'broadly_accessible', 'Broad','CT')
  ipt_target_sharedA_dt$count <- 1
  aggregate_dt_sharedA <- aggregate(count ~ geneCate + celltypeCate, data = ipt_target_sharedA_dt, sum)
  aggregate_dt_sharedA$cate <- 'shared_A'
  
  aggregate_celltypecate <- aggregate(count ~ celltypeCate, data = aggregate_dt_sharedA, sum)
  combine_aggregate_sharedA <- merge(aggregate_dt_sharedA, aggregate_celltypecate, by.x = 'celltypeCate',by.y = 'celltypeCate')
  combine_aggregate_sharedA$prop <- combine_aggregate_sharedA$count.x/combine_aggregate_sharedA$count.y
  
  
  merged_dt <- rbind(combine_aggregate_notshared,combine_aggregate_sharedI,combine_aggregate_sharedA)
  merged_dt$spe <- target_fl_path_nm
  return(merged_dt)
  
})

combine_dt <- do.call(rbind,outs)
dim(combine_dt)
head(combine_dt)

##do the plotting
cate_list <- unique(combine_dt$cate)

for (i in 1:length(cate_list)){
  
  combine_cate_dt <- combine_dt[combine_dt$cate == cate_list[i],]
  
  combine_cate_dt$geneCate <- factor(combine_cate_dt$geneCate,levels = c('Genic','Proximal','Distal'))

  p <- ggplot(data=combine_cate_dt, aes(x=geneCate , y=prop, fill=celltypeCate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    facet_wrap(~spe,nrow =  1) ## You can use face_wrap function only if you need it+
  
  pdf(paste0(opt_dir,'/','opt_',cate_list[i],'_',target_fl_path_nm,'_gene_proximity_compare_prop','.pdf'),width = 12,height = 5 )
  print(p)
  dev.off()

}


spe_list <- unique(combine_dt$spe)

for (i in 1:length(spe_list)){
  
  combine_cate_dt <- combine_dt[combine_dt$spe == spe_list[i],]
  
  combine_cate_dt$geneCate <- factor(combine_cate_dt$geneCate,levels = c('Genic','Proximal','Distal'))
  
  combine_cate_dt$cate <- factor(combine_cate_dt$cate,levels = c('shared_A','shared_I','not_shared'))
  
  p <- ggplot(data=combine_cate_dt, aes(x=cate , y=prop, fill=celltypeCate)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    ) +
    facet_wrap(~geneCate,nrow =  1) ## You can use face_wrap function only if you need it+
  
  pdf(paste0(opt_dir,'/','opt_',spe_list[i],'_gene_proximity_compare_prop','.pdf'),width = 12,height = 5 )
  print(p)
  dev.off()
  
}








########
##step03 check the H3K27me3 overlapping with ACR across different species
########
ipt_dt <- data.frame(
  Species = c("Rice:Maize", "Rice:Maize", "Rice:Sorghum", "Rice:Sorghum"),
  Prop = c(0.319587629, 0.213781418, 0.401197605, 0.296886142),
  Cate = c("H3K27me3", "H3K27me3-absent", "H3K27me3", "H3K27me3-absent")
)

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )
pdf('opt3_final_compare_BtoB_underH3K27me3_H3K27absent_in_sorghum_maize.pdf',width = 8,height = 8 )
p
dev.off()


#################################
##check the distance to the genes
ipt_dt <- read.delim('opt_final_summary_file_add_rice_maize.txt',header = T)
head(ipt_dt)
prefix <- 'maize'

ipt_dt <- read.delim('opt_final_summary_file_add_rice_sorghum.txt',header = T)
head(ipt_dt)
prefix <- 'sorghum'


##For the maize
##Under the H3K27me3
##check broad to broad
ipt_broadTobroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' & ipt_dt$maize_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]
head(ipt_broadTobroad_dt)
ipt_broadTobroad_target_dt <- ipt_broadTobroad_dt[c('rice_ACR','maize_ACR','riceGeneCate','maizeGeneCate')]
ipt_broadTobroad_target_dt$final_cate <- paste0(ipt_broadTobroad_target_dt$riceGeneCate,'_',ipt_broadTobroad_target_dt$maizeGeneCate)
head(ipt_broadTobroad_target_dt)
ipt_broadTobroad_target_dt <- unique(ipt_broadTobroad_target_dt)
opt_broadTobroad_summary_num_dt <- as.data.frame(table(ipt_broadTobroad_target_dt$final_cate))
opt_broadTobroad_summary_num_dt
colnames(opt_broadTobroad_summary_num_dt) <- c('Cate','BtoBH3k27')

##broad to not broad
ipt_broadToNbroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$maize_InorNotH3K27m3 == 'none'| 
                                             ipt_broadToNbroad_dt$maize_InorNotH3K27m3 == 'Others' | 
                                             ipt_broadToNbroad_dt$maize_InorNotH3K27m3 == 'RestrictedInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$maizeOVerH3K27Cate == 'OverlapH3K27me3' |
                                               ipt_broadToNbroad_dt$maizeRegionOverH3K27Cate == 'OverlapH3K27me3',]


head(ipt_broadToNbroad_dt)
ipt_broadToNbroad_target_dt <- ipt_broadToNbroad_dt[c('rice_ACR','maize_ACR','riceGeneCate','maizeGeneCate')]
ipt_broadToNbroad_target_dt$final_cate <- paste0(ipt_broadToNbroad_target_dt$riceGeneCate,'_',ipt_broadToNbroad_target_dt$maizeGeneCate)
head(ipt_broadToNbroad_target_dt)
ipt_broadToNbroad_target_dt <- unique(ipt_broadToNbroad_target_dt)
opt_broadToNbroad_summary_num_dt <- as.data.frame(table(ipt_broadToNbroad_target_dt$final_cate))
opt_broadToNbroad_summary_num_dt
colnames(opt_broadToNbroad_summary_num_dt) <- c('Cate','BtoNBH3K27')

##not under the H3K27me3
##check broad to broad
ipt_broadTobroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                       ipt_dt$maizeOVerH3K27Cate == 'none' &
                                       ipt_dt$rice_BroOrRes == 'broad' &
                                       ipt_dt$maize_BroOrRes == 'broad',]

head(ipt_broadTobroad_nH3K27_dt)
ipt_broadTobroad_target_nH3K27_dt <- ipt_broadTobroad_nH3K27_dt[c('rice_ACR','maize_ACR','riceGeneCate','maizeGeneCate')]
ipt_broadTobroad_target_nH3K27_dt$final_cate <- paste0(ipt_broadTobroad_target_nH3K27_dt$riceGeneCate,'_',ipt_broadTobroad_target_nH3K27_dt$maizeGeneCate)
head(ipt_broadTobroad_target_nH3K27_dt)
ipt_broadTobroad_target_nH3K27_dt <- unique(ipt_broadTobroad_target_nH3K27_dt)
opt_broadTobroad_nH3K27_summary_num_dt <- as.data.frame(table(ipt_broadTobroad_target_nH3K27_dt$final_cate))
opt_broadTobroad_nH3K27_summary_num_dt
colnames(opt_broadTobroad_nH3K27_summary_num_dt) <- c('Cate','BtoBnoH3K27')


##broad to not broad
ipt_broadToNbroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                       ipt_dt$maizeOVerH3K27Cate == 'none' &
                                       ipt_dt$rice_BroOrRes == 'broad' &
                                       ipt_dt$maize_BroOrRes != 'broad'&
                                        ipt_dt$maizeRegionOverH3K27Cate == 'none',]

head(ipt_broadToNbroad_nH3K27_dt)
ipt_broadToNbroad_nH3K27_target_dt <- ipt_broadToNbroad_nH3K27_dt[c('rice_ACR','maize_ACR','riceGeneCate','maizeGeneCate')]
ipt_broadToNbroad_nH3K27_target_dt$final_cate <- paste0(ipt_broadToNbroad_nH3K27_target_dt$riceGeneCate,'_',ipt_broadToNbroad_nH3K27_target_dt$maizeGeneCate)
head(ipt_broadToNbroad_nH3K27_target_dt)
ipt_broadToNbroad_nH3K27_target_dt <- unique(ipt_broadToNbroad_nH3K27_target_dt)
opt_broadToNbroad_nH3K27_summary_num_dt <- as.data.frame(table(ipt_broadToNbroad_nH3K27_target_dt$final_cate))
opt_broadToNbroad_nH3K27_summary_num_dt
colnames(opt_broadToNbroad_nH3K27_summary_num_dt) <- c('Cate','BtoNBnoH3K27')

##Here we will build a heatmap matrix to show the categories
merged_dt_H3K27 <- merge(x = opt_broadTobroad_summary_num_dt, y = opt_broadToNbroad_summary_num_dt,by.x = 'Cate',by.y ='Cate',all=T)
merged_dt_noH3K27 <- merge(x = opt_broadTobroad_nH3K27_summary_num_dt, y = opt_broadToNbroad_nH3K27_summary_num_dt,by.x = 'Cate',by.y = 'Cate', all=T)

merged_dt <- merge(x = merged_dt_H3K27, y = merged_dt_noH3K27, by.x = 'Cate', by.y = 'Cate',all=T)
merged_dt[is.na(merged_dt)] <- 0

rownames(merged_dt) <- merged_dt$Cate
merged_dt <- merged_dt[,2:5]
merged_t_dt <- t(merged_dt)


dt_zscore <- apply(merged_t_dt,1,function(x){
  (x - mean(x))/sd(x)
})

dt_zscore <- data.frame(t(dt_zscore))

cate_levels <- c('Distal_Distal','Proximal_Proximal','Genic_Genic','Distal_Proximal','Distal_Genic',
                 'Proximal_Distal','Proximal_Genic','Genic_Distal','Genic_Proximal',
                 'Distal_none','Proximal_none','Genic_none')
colnames(dt_zscore) <- factor(colnames(dt_zscore),levels = cate_levels)

dt_zscore <- as.matrix(dt_zscore)
dt_zscore <- dt_zscore[,cate_levels]


pdf(paste0('opt3_different_cate_num_zscore_maize.pdf'),width = 5,height = 3)
pheatmap(dt_zscore,
         #scale="row",
         #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
         #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
         cluster_cols = F,
         cluster_rows = F,
         border_color='white',
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize_number = 25
) 
dev.off()



##For the sorghum
##Under the H3K27me3
##check broad to broad
ipt_broadTobroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' & ipt_dt$sorghum_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]
head(ipt_broadTobroad_dt)
ipt_broadTobroad_target_dt <- ipt_broadTobroad_dt[c('rice_ACR','sorghum_ACR','riceGeneCate','sorghumGeneCate')]
ipt_broadTobroad_target_dt$final_cate <- paste0(ipt_broadTobroad_target_dt$riceGeneCate,'_',ipt_broadTobroad_target_dt$sorghumGeneCate)
head(ipt_broadTobroad_target_dt)
ipt_broadTobroad_target_dt <- unique(ipt_broadTobroad_target_dt)
opt_broadTobroad_summary_num_dt <- as.data.frame(table(ipt_broadTobroad_target_dt$final_cate))
opt_broadTobroad_summary_num_dt
colnames(opt_broadTobroad_summary_num_dt) <- c('Cate','BtoBH3k27')

##broad to not broad
ipt_broadToNbroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$sorghum_InorNotH3K27m3 == 'none'| 
                                               ipt_broadToNbroad_dt$sorghum_InorNotH3K27m3 == 'Others' | 
                                               ipt_broadToNbroad_dt$sorghum_InorNotH3K27m3 == 'RestrictedInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$sorghumOVerH3K27Cate == 'OverlapH3K27me3' |
                                               ipt_broadToNbroad_dt$sorghumRegionOverH3K27Cate == 'OverlapH3K27me3',]


head(ipt_broadToNbroad_dt)
ipt_broadToNbroad_target_dt <- ipt_broadToNbroad_dt[c('rice_ACR','sorghum_ACR','riceGeneCate','sorghumGeneCate')]
ipt_broadToNbroad_target_dt$final_cate <- paste0(ipt_broadToNbroad_target_dt$riceGeneCate,'_',ipt_broadToNbroad_target_dt$sorghumGeneCate)
head(ipt_broadToNbroad_target_dt)
ipt_broadToNbroad_target_dt <- unique(ipt_broadToNbroad_target_dt)
opt_broadToNbroad_summary_num_dt <- as.data.frame(table(ipt_broadToNbroad_target_dt$final_cate))
opt_broadToNbroad_summary_num_dt
colnames(opt_broadToNbroad_summary_num_dt) <- c('Cate','BtoNBH3K27')

##not under the H3K27me3
##check broad to broad
ipt_broadTobroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                       ipt_dt$sorghumOVerH3K27Cate == 'none' &
                                       ipt_dt$rice_BroOrRes == 'broad' &
                                       ipt_dt$sorghum_BroOrRes == 'broad',]

head(ipt_broadTobroad_nH3K27_dt)
ipt_broadTobroad_target_nH3K27_dt <- ipt_broadTobroad_nH3K27_dt[c('rice_ACR','sorghum_ACR','riceGeneCate','sorghumGeneCate')]
ipt_broadTobroad_target_nH3K27_dt$final_cate <- paste0(ipt_broadTobroad_target_nH3K27_dt$riceGeneCate,'_',ipt_broadTobroad_target_nH3K27_dt$sorghumGeneCate)
head(ipt_broadTobroad_target_nH3K27_dt)
ipt_broadTobroad_target_nH3K27_dt <- unique(ipt_broadTobroad_target_nH3K27_dt)
opt_broadTobroad_nH3K27_summary_num_dt <- as.data.frame(table(ipt_broadTobroad_target_nH3K27_dt$final_cate))
opt_broadTobroad_nH3K27_summary_num_dt
colnames(opt_broadTobroad_nH3K27_summary_num_dt) <- c('Cate','BtoBnoH3K27')


##broad to not broad
ipt_broadToNbroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                        ipt_dt$sorghumOVerH3K27Cate == 'none' &
                                        ipt_dt$rice_BroOrRes == 'broad' &
                                        ipt_dt$sorghum_BroOrRes != 'broad'&
                                        ipt_dt$sorghumRegionOverH3K27Cate == 'none',]

head(ipt_broadToNbroad_nH3K27_dt)
ipt_broadToNbroad_nH3K27_target_dt <- ipt_broadToNbroad_nH3K27_dt[c('rice_ACR','sorghum_ACR','riceGeneCate','sorghumGeneCate')]
ipt_broadToNbroad_nH3K27_target_dt$final_cate <- paste0(ipt_broadToNbroad_nH3K27_target_dt$riceGeneCate,'_',ipt_broadToNbroad_nH3K27_target_dt$sorghumGeneCate)
head(ipt_broadToNbroad_nH3K27_target_dt)
ipt_broadToNbroad_nH3K27_target_dt <- unique(ipt_broadToNbroad_nH3K27_target_dt)
opt_broadToNbroad_nH3K27_summary_num_dt <- as.data.frame(table(ipt_broadToNbroad_nH3K27_target_dt$final_cate))
opt_broadToNbroad_nH3K27_summary_num_dt
colnames(opt_broadToNbroad_nH3K27_summary_num_dt) <- c('Cate','BtoNBnoH3K27')

##Here we will build a heatmap matrix to show the categories
merged_dt_H3K27 <- merge(x = opt_broadTobroad_summary_num_dt, y = opt_broadToNbroad_summary_num_dt,by.x = 'Cate',by.y ='Cate',all=T)
merged_dt_noH3K27 <- merge(x = opt_broadTobroad_nH3K27_summary_num_dt, y = opt_broadToNbroad_nH3K27_summary_num_dt,by.x = 'Cate',by.y = 'Cate', all=T)

merged_dt <- merge(x = merged_dt_H3K27, y = merged_dt_noH3K27, by.x = 'Cate', by.y = 'Cate',all=T)
merged_dt[is.na(merged_dt)] <- 0

rownames(merged_dt) <- merged_dt$Cate
merged_dt <- merged_dt[,2:5]
merged_t_dt <- t(merged_dt)


dt_zscore <- apply(merged_t_dt,1,function(x){
  (x - mean(x))/sd(x)
})

dt_zscore <- data.frame(t(dt_zscore))

cate_levels <- c('Distal_Distal','Proximal_Proximal','Genic_Genic','Distal_Proximal','Distal_Genic',
                 'Proximal_Distal','Proximal_Genic','Genic_Distal','Genic_Proximal',
                 'Distal_none','Proximal_none','Genic_none')
colnames(dt_zscore) <- factor(colnames(dt_zscore),levels = cate_levels)

dt_zscore <- as.matrix(dt_zscore)
dt_zscore <- dt_zscore[,cate_levels]


pdf(paste0('opt3_different_cate_num_zscore_sorghum.pdf'),width = 5,height = 3)
pheatmap(dt_zscore,
         #scale="row",
         #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
         #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
         cluster_cols = F,
         cluster_rows = F,
         border_color='white',
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize_number = 25
) 
dev.off()




###############################
##calculate the rice ACR number
ipt_dt <- read.delim('opt_final_summary_file_add_rice_maize.txt',header = T)
head(ipt_dt)
prefix <- 'maize'

ipt_dt <- read.delim('opt_final_summary_file_add_rice_sorghum.txt',header = T)
head(ipt_dt)
prefix <- 'sorghum'
###############
##for the maize
##under H3K27me3
##broadtobroad 
ipt_broadTobroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' & ipt_dt$maize_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]
head(ipt_broadTobroad_dt)
length(unique(ipt_broadTobroad_dt$rice_ACR))
#62

##broad to nbroad
ipt_broadToNbroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$maize_InorNotH3K27m3 == 'none'| 
                                               ipt_broadToNbroad_dt$maize_InorNotH3K27m3 == 'Others' | 
                                               ipt_broadToNbroad_dt$maize_InorNotH3K27m3 == 'RestrictedInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$maizeOVerH3K27Cate == 'OverlapH3K27me3' |
                                               ipt_broadToNbroad_dt$maizeRegionOverH3K27Cate == 'OverlapH3K27me3',]
length(unique(ipt_broadToNbroad_dt$rice_ACR))
##132


##not under the H3K27me3
##broad to broad
ipt_broadTobroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                       ipt_dt$maizeOVerH3K27Cate == 'none' &
                                       ipt_dt$rice_BroOrRes == 'broad' &
                                       ipt_dt$maize_BroOrRes == 'broad',]

length(unique(ipt_broadTobroad_nH3K27_dt$rice_ACR))
##1514

##broad to Nbroad
ipt_broadToNbroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                        ipt_dt$maizeOVerH3K27Cate == 'none' &
                                        ipt_dt$rice_BroOrRes == 'broad' &
                                        ipt_dt$maize_BroOrRes != 'broad'&
                                        ipt_dt$maizeRegionOverH3K27Cate == 'none',]

length(unique(ipt_broadToNbroad_nH3K27_dt$rice_ACR))
##5568

##conduct the fisher exact test
BtoB_underH3K27me3 <- length(unique(ipt_broadTobroad_dt$rice_ACR))
BtoB_notH3K27me3 <- length(unique(ipt_broadTobroad_nH3K27_dt$rice_ACR))
BtoNB_underH3K27me3 <- length(unique(ipt_broadToNbroad_dt$rice_ACR))
BtoNB_notH3K27me3 <- length(unique(ipt_broadToNbroad_nH3K27_dt$rice_ACR))
  
data <- matrix(c(BtoB_underH3K27me3, BtoB_notH3K27me3, BtoNB_underH3K27me3, BtoNB_notH3K27me3), nrow = 2)

# Perform Fisher's exact test
fisher_result <- fisher.test(data,alternative = 'greater')
fisher_result
##0.0004789


#################
##for the sorghum
ipt_broadTobroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' & ipt_dt$sorghum_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]
head(ipt_broadTobroad_dt)
length(unique(ipt_broadTobroad_dt$rice_ACR))
#134

##broad to nbroad
ipt_broadToNbroad_dt <- ipt_dt[ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$sorghum_InorNotH3K27m3 == 'none'| 
                                               ipt_broadToNbroad_dt$sorghum_InorNotH3K27m3 == 'Others' | 
                                               ipt_broadToNbroad_dt$sorghum_InorNotH3K27m3 == 'RestrictedInflankH3K27me3peak',]

ipt_broadToNbroad_dt <- ipt_broadToNbroad_dt[ipt_broadToNbroad_dt$sorghumOVerH3K27Cate == 'OverlapH3K27me3' |
                                               ipt_broadToNbroad_dt$sorghumRegionOverH3K27Cate == 'OverlapH3K27me3',]
length(unique(ipt_broadToNbroad_dt$rice_ACR))
##200


##not under the H3K27me3
##broad to broad
ipt_broadTobroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                       ipt_dt$sorghumOVerH3K27Cate == 'none' &
                                       ipt_dt$rice_BroOrRes == 'broad' &
                                       ipt_dt$sorghum_BroOrRes == 'broad',]

length(unique(ipt_broadTobroad_nH3K27_dt$rice_ACR))
##3747

##broad to Nbroad
ipt_broadToNbroad_nH3K27_dt <- ipt_dt[ipt_dt$riceOverH3K27Cate == 'none' & 
                                        ipt_dt$sorghumOVerH3K27Cate == 'none' &
                                        ipt_dt$rice_BroOrRes == 'broad' &
                                        ipt_dt$sorghum_BroOrRes != 'broad'&
                                        ipt_dt$sorghumRegionOverH3K27Cate == 'none',]

length(unique(ipt_broadToNbroad_nH3K27_dt$rice_ACR))
##8874



##conduct the fisher exact test
BtoB_underH3K27me3 <- length(unique(ipt_broadTobroad_dt$rice_ACR))
BtoB_notH3K27me3 <- length(unique(ipt_broadTobroad_nH3K27_dt$rice_ACR))
BtoNB_underH3K27me3 <- length(unique(ipt_broadToNbroad_dt$rice_ACR))
BtoNB_notH3K27me3 <- length(unique(ipt_broadToNbroad_nH3K27_dt$rice_ACR))

data <- matrix(c(BtoB_underH3K27me3, BtoB_notH3K27me3, BtoNB_underH3K27me3, BtoNB_notH3K27me3), nrow = 2)

# Perform Fisher's exact test
fisher_result <- fisher.test(data,alternative = 'greater')
##3.899e-05




#################
##updating 112823
##check the organ overlapping
ipt_dir <- 'opt3_different_organs_H3K27me3_112823/////'
opt_dir <- 'opt3_different_organs_H3K27me3_results_112823///'

all_fl_list <- list.files(path= ipt_dir)

for (i in 1:length(all_fl_list)){
  
  target_fl_path <- paste0(ipt_dir,'/',all_fl_list[i])
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('H3K27me3_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)

  ipt_dt <- read.delim(target_fl_path,header = F)
  ipt_dt$riceACR <- paste0(ipt_dt$V1,'_',ipt_dt$V2,'_',ipt_dt$V3)  
  head(ipt_dt)
  
  ##total broadflankIn in rice
  ipt_dt_BFI_rice <- ipt_dt[ipt_dt$V4 == 'BroadInflankH3K27me3peak',]
  length(unique(ipt_dt_BFI_rice$riceACR))
  ##716
  
  ##895
  
  ##broadFin in rice leaf and in rice other organ
  ipt_dt_BFI_rice_otherorgan <- ipt_dt_BFI_rice[ipt_dt_BFI_rice$V13 == 'BroadInflankH3K27me3peak',]
  length(unique(ipt_dt_BFI_rice_otherorgan$riceACR))
  ##225
  
  ##390
  
  ##broadFin in rice leaf and in rice other organ in maize 
  ipt_dt_BFI_rice_otherorgan_maize <- ipt_dt_BFI_rice_otherorgan[ipt_dt_BFI_rice_otherorgan$V8 == 'BroadInflankH3K27me3peak',]
  length(unique(ipt_dt_BFI_rice_otherorgan_maize$riceACR))
  ##17
  
  ##26
  
  ipt_dt_BFI_rice_NotInotherorgan <- ipt_dt_BFI_rice[ipt_dt_BFI_rice$V13 != 'BroadInflankH3K27me3peak',]
  length(unique(ipt_dt_BFI_rice_NotInotherorgan$riceACR))
  ##496
  
  ##513
  
  
  ipt_dt_BFI_rice_NotInotherorgan_maize <- ipt_dt_BFI_rice_NotInotherorgan[ipt_dt_BFI_rice_NotInotherorgan$V8 == 'BroadInflankH3K27me3peak',]
  length(unique(ipt_dt_BFI_rice_NotInotherorgan_maize$riceACR))
  ##25
  
  ##32
  
  ipt_dt_BFI_rice_maize <- ipt_dt_BFI_rice[ipt_dt_BFI_rice$V8 == 'BroadInflankH3K27me3peak',]
  length(unique(ipt_dt_BFI_rice_maize$riceACR))

  ##use the fisher 
  data <- matrix(c(17, 25, 225-17,496-25), nrow = 2)
  
  # Assign row and column names
  rownames(data) <- c("Group A", "Group B")
  colnames(data) <- c("Outcome 1", "Outcome 2")
  data
  result <- fisher.test(data)
  ##p = 0.2
  
  data <- matrix(c(26, 32, 390-26,513-32), nrow = 2)
  
  # Assign row and column names
  rownames(data) <- c("Group A", "Group B")
  colnames(data) <- c("Outcome 1", "Outcome 2")
  data
  result <- fisher.test(data)
  result
  ##0.7862
}



#################
##updating 121223
##this function is to check whether we find something for the lineage specific ACRs under the H3K27me3
ipt_dir <- 'opt3_check_H3K27me3_in_syntenic_121223/'
opt_dir <- 'opt3_check_H3K27me3_in_syntenic_results_121223/'

all_fl_list <- list.files(path= ipt_dir)

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_summary_file_add_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  target_species <- gsub('rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  
  colnames(ipt_dt)
  
  ##we will compare the broadACR under H3K27me3 is corresponding to region without H3K27me3
  #ipt_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none',]
  
  #shared_A_broad_bothspe_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt$maize_ACR != 'none' & ipt_dt$rice_BroOrRes == 'broad' & ipt_dt$maize_BroOrRes == 'broad',]
  #length(unique(shared_A_broad_bothspe_dt$rice_ACR))
  ##1732
  
  #shared_I_broad_rice_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt$maize_ACR == 'none' & ipt_dt$rice_BroOrRes == 'broad',]
  #length(unique(shared_I_broad_rice_dt$rice_ACR))
  ##1706
  
  #shared_not_shared_broad_rice_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none' & ipt_dt$rice_BroOrRes == 'broad',]
  #length(unique(shared_not_shared_broad_rice_dt$rice_ACR))
  ##4309
  
  
  
  ##check the proportion under the H3K27me3
  shared_A_broad_bothspe_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] != 'none' & ipt_dt$rice_BroOrRes == 'broad' & ipt_dt[[paste0(target_species,'_BroOrRes')]] == 'broad' &
                                                 ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' ,]
  length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  
  
  shared_I_broad_rice_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                              ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' ,]
  length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  
  
  shared_not_shared_broad_rice_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                                       ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]
  length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  
  
  ##check the proportion not under the H3K27me3
  shared_A_broad_bothspe_not_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] != 'none' & ipt_dt$rice_BroOrRes == 'broad' & ipt_dt[[paste0(target_species,'_BroOrRes')]] == 'broad' &
                                                     ipt_dt$riceOverH3K27Cate == 'none' ,]
  length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  
  
  shared_I_broad_rice_not_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                                  ipt_dt$riceOverH3K27Cate == 'none',]
  length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  
  shared_not_shared_broad_rice_not_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                                           ipt_dt$riceOverH3K27Cate == 'none',]
  length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  
  
  #3733/(1578+1499+3733)
  #576/(154 + 207 + 576)
  
  ##for the not shared
  In_H3K27me3_not_shared_broad_num <- length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  In_H3K27me3_Not_notshared_broad_num <- length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  
  NotIn_H3K27me3_not_shared_broad_num <- length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  NotIn_H3K27me3_Not_notshared_broad_num <-  length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  data <- matrix(c(In_H3K27me3_not_shared_broad_num, NotIn_H3K27me3_not_shared_broad_num, In_H3K27me3_Not_notshared_broad_num, NotIn_H3K27me3_Not_notshared_broad_num), nrow = 2)
  
  # Perform Fisher's exact test
  result_not_shared <- fisher.test(data, alternative = 'two.sided')
  
  ##for the shared A
  In_H3K27me3_shared_A_broad_num <- length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  In_H3K27me3_Not_shared_A_broad_num <- length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  
  NotIn_H3K27me3_shared_A_broad_num <- length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  NotIn_H3K27me3_Not_shared_A_broad_num <-  length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  data <- matrix(c(In_H3K27me3_shared_A_broad_num, NotIn_H3K27me3_shared_A_broad_num, In_H3K27me3_Not_shared_A_broad_num, NotIn_H3K27me3_Not_shared_A_broad_num), nrow = 2)
  
  # Perform Fisher's exact test
  result_shared_A <- fisher.test(data, alternative = 'two.sided')
  
  ##for the shared I
  In_H3K27me3_shared_I_broad_num <- length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  In_H3K27me3_Not_shared_I_broad_num <- length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR)) + length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  
  NotIn_H3K27me3_shared_I_broad_num <- length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  NotIn_H3K27me3_Not_shared_I_broad_num <-  length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR)) + length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  
  data <- matrix(c(In_H3K27me3_shared_I_broad_num, NotIn_H3K27me3_shared_I_broad_num, In_H3K27me3_Not_shared_I_broad_num, NotIn_H3K27me3_Not_shared_I_broad_num), nrow = 2)
  
  # Perform Fisher's exact test
  result_shared_I <- fisher.test(data, alternative = 'two.sided')
  

  
  #make the barplot for the proportion compared to the not under H3K27me3
  SS_ACR_H3K27_rice_maize = length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  total_ACR_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR)) + 
                               length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR)) + 
                               length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  
  SS_ACR_not_H3K27_rice_maize = length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  total_ACR_not_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR)) + 
                                   length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR)) + 
                                   length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  SH_ACR_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  SH_ACR_not_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  
  SI_ACR_H3K27_rice_maize = length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  SI_ACR_not_H3K27_rice_maize = length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
    
  
  
  ##make a data frame
  my_data <- data.frame(
    SS_ACR_H3K27 = SS_ACR_H3K27_rice_maize,
    total_ACR_H3K27 = total_ACR_H3K27_rice_maize,
    SS_ACR_not_H3K27 = SS_ACR_not_H3K27_rice_maize,
    total_ACR_not_H3K27 = total_ACR_not_H3K27_rice_maize,
    SH_ACR_H3K27 = SH_ACR_H3K27_rice_maize,
    SH_ACR_not_H3K27 = SH_ACR_not_H3K27_rice_maize,
    SI_ACR_H3K27 = SI_ACR_H3K27_rice_maize,
    SI_ACR_not_H3K27 = SI_ACR_not_H3K27_rice_maize
  )
  
  my_data <- t(my_data)
  colnames(my_data) <- c('Count')
  my_data <- as.data.frame(my_data)
  my_data$spe <- target_fl_path_nm
  my_data$not_shared_pval <- result_not_shared$p.value
  my_data$shared_A_pval <- result_shared_A$p.value
  my_data$shared_I_pval <- result_shared_I$p.value
  
  
  
  return(my_data)
  
})

combine_dt <- do.call(rbind,outs)

write.csv(combine_dt,paste0(opt_dir,'/opt_proportion_under_H3K27_no_H3K27_pval.csv'),quote = F)

combine_maize_dt <- combine_dt[combine_dt$spe == 'rice_maize',]
combine_sorghum_dt <- combine_dt[combine_dt$spe == 'rice_sorghum',]

##for the species-specific-ACR
ipt_dt <- data.frame(
  Species = c("Os:Zm", "Os:Zm","Os:Sb", "Os:Sb"),
  Prop = c(combine_maize_dt['SS_ACR_H3K27','Count']/combine_maize_dt['total_ACR_H3K27','Count'],
           combine_maize_dt['SS_ACR_not_H3K27','Count']/combine_maize_dt['total_ACR_not_H3K27','Count'],
           combine_sorghum_dt['SS_ACR_H3K27','Count']/combine_sorghum_dt['total_ACR_H3K27','Count'],
           combine_sorghum_dt['SS_ACR_not_H3K27','Count']/combine_sorghum_dt['total_ACR_not_H3K27','Count']),
  Cate = c("SS-ACR-H3K27", "SS-ACR-no-H3K27", "SS-ACR-H3K27", "SS-ACR-no-H3K27"),
  Target = c(combine_maize_dt['SS_ACR_H3K27','Count'], combine_maize_dt['SS_ACR_not_H3K27','Count'],
             combine_sorghum_dt['SS_ACR_H3K27','Count'], combine_sorghum_dt['SS_ACR_not_H3K27','Count'])
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Os:Zm','Os:Sb'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  
#facet_wrap(~Rice_CellType,nrow = 1)
p
pdf(paste0(opt_dir,'/','opt_not_shared_proportion_under_H3K27_no_H3K27.pdf'),width = 8,height = 6 )
p
dev.off()


##for the shared-ACR
ipt_dt <- data.frame(
  Species = c("Os:Zm", "Os:Zm","Os:Sb", "Os:Sb"),
  Prop = c(combine_maize_dt['SH_ACR_H3K27','Count']/combine_maize_dt['total_ACR_H3K27','Count'],
           combine_maize_dt['SH_ACR_not_H3K27','Count']/combine_maize_dt['total_ACR_not_H3K27','Count'],
           combine_sorghum_dt['SH_ACR_H3K27','Count']/combine_sorghum_dt['total_ACR_H3K27','Count'],
           combine_sorghum_dt['SH_ACR_not_H3K27','Count']/combine_sorghum_dt['total_ACR_not_H3K27','Count']),
    
  Cate = c("SH-ACR-H3K27", "SH-ACR-no-H3K27", "SH-ACR-H3K27", "SH-ACR-no-H3K27"),
  
  Target = c(combine_maize_dt['SH_ACR_H3K27','Count'],combine_maize_dt['SH_ACR_not_H3K27','Count'],
             combine_sorghum_dt['SH_ACR_H3K27','Count'],combine_sorghum_dt['SH_ACR_not_H3K27','Count'])         
           
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Os:Zm','Os:Sb'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  
#facet_wrap(~Rice_CellType,nrow = 1)
p
pdf(paste0(opt_dir,'/','opt_shared_proportion_under_H3K27_no_H3K27.pdf'),width = 8,height = 6 )
p
dev.off()


##for the shared-I-ACR
ipt_dt <- data.frame(
  Species = c("Os:Zm", "Os:Zm","Os:Sb", "Os:Sb"),
  Prop = c(combine_maize_dt['SI_ACR_H3K27','Count']/combine_maize_dt['total_ACR_H3K27','Count'],
           combine_maize_dt['SI_ACR_not_H3K27','Count']/combine_maize_dt['total_ACR_not_H3K27','Count'],
           combine_sorghum_dt['SI_ACR_H3K27','Count']/combine_sorghum_dt['total_ACR_H3K27','Count'],
           combine_sorghum_dt['SI_ACR_not_H3K27','Count']/combine_sorghum_dt['total_ACR_not_H3K27','Count']),
  
  Cate = c("SI-ACR-H3K27", "SI-ACR-no-H3K27", "SI-ACR-H3K27", "SI-ACR-no-H3K27"),
  
  Target = c(combine_maize_dt['SI_ACR_H3K27','Count'],combine_maize_dt['SI_ACR_not_H3K27','Count'],
             combine_sorghum_dt['SI_ACR_H3K27','Count'],combine_sorghum_dt['SI_ACR_not_H3K27','Count'])
  
)

ipt_dt$Species <- factor(ipt_dt$Species, levels = c('Os:Zm','Os:Sb'))

p <- ggplot(data=ipt_dt, aes(x=Species, y=Prop, fill=Cate)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  
#facet_wrap(~Rice_CellType,nrow = 1)
p
pdf(paste0(opt_dir,'/','opt_shared_I_proportion_under_H3K27_no_H3K27.pdf'),width = 8,height = 6 )
p
dev.off()

##updating 020124
##we will plot the pie chart
combine_maize_dt <- combine_dt[combine_dt$spe == 'rice_maize',]
combine_sorghum_dt <- combine_dt[combine_dt$spe == 'rice_sorghum',]

##for the rice to maize
generate_pie <- function(df,cate) {
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df[[cate]],pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}

combine_maize_target_dt <- combine_maize_dt[c('SS_ACR_H3K27','SH_ACR_H3K27','SI_ACR_H3K27'),]
combine_maize_target_dt <- combine_maize_target_dt[c('Count')]
combine_maize_target_dt$cate <- rownames(combine_maize_target_dt)
combine_maize_target_dt$count <- combine_maize_target_dt$Count

pdf(paste0(opt_dir,'/opt_all_H3K27me3_syntenic_ACR_composition_piechart_','maize','.pdf'),width = 8,height = 8 )
print(generate_pie(combine_maize_target_dt,'cate'))
dev.off()

combine_sorghum_target_dt <- combine_sorghum_dt[c('SS_ACR_H3K27','SH_ACR_H3K27','SI_ACR_H3K27'),]
combine_sorghum_target_dt <- combine_sorghum_target_dt[c('Count')]
combine_sorghum_target_dt$cate <- rownames(combine_sorghum_target_dt)
combine_sorghum_target_dt$count <- combine_sorghum_target_dt$Count

pdf(paste0(opt_dir,'/opt_all_H3K27me3_syntenic_ACR_composition_piechart_','sorghum','.pdf'),width = 8,height = 8 )
print(generate_pie(combine_sorghum_target_dt,'cate'))
dev.off()


#################
##updating 122023
##we will check the overlapping with the CNS data
ipt_dir <- 'opt3_check_H3K27me3_in_syntenic_CNS_122023///'
opt_dir <- 'opt3_check_H3K27me3_in_syntenic_CNS_results_122023//'

all_fl_list <- list.files(path= ipt_dir)

analyze_ACR_syntenic_H3K27me3 <- function(ipt_dt,target_species,alternative_str){
  
  ##check the proportion under the H3K27me3
  shared_A_broad_bothspe_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] != 'none' & ipt_dt$rice_BroOrRes == 'broad' & ipt_dt[[paste0(target_species,'_BroOrRes')]] == 'broad' &
                                                 ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' ,]
  length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  
  
  shared_I_broad_rice_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                              ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak' ,]
  length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  
  
  shared_not_shared_broad_rice_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                                       ipt_dt$rice_InorNotH3K27m3 == 'BroadInflankH3K27me3peak',]
  length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  
  
  ##check the proportion not under the H3K27me3
  shared_A_broad_bothspe_not_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] != 'none' & ipt_dt$rice_BroOrRes == 'broad' & ipt_dt[[paste0(target_species,'_BroOrRes')]] == 'broad' &
                                                     ipt_dt$riceOverH3K27Cate == 'none' ,]
  length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  
  
  shared_I_broad_rice_not_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID != 'none' & ipt_dt[[paste0(target_species,'_ACR')]] == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                                  ipt_dt$riceOverH3K27Cate == 'none',]
  length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  
  shared_not_shared_broad_rice_not_H3K27me3_dt <- ipt_dt[ipt_dt$rice_ACRID == 'none' & ipt_dt$rice_BroOrRes == 'broad' &
                                                           ipt_dt$riceOverH3K27Cate == 'none',]
  length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  
  
  #3733/(1578+1499+3733)
  #576/(154 + 207 + 576)
  
  ##for the not shared
  In_H3K27me3_not_shared_broad_num <- length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  In_H3K27me3_Not_notshared_broad_num <- length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  
  NotIn_H3K27me3_not_shared_broad_num <- length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  NotIn_H3K27me3_Not_notshared_broad_num <-  length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  data <- matrix(c(In_H3K27me3_not_shared_broad_num, NotIn_H3K27me3_not_shared_broad_num, In_H3K27me3_Not_notshared_broad_num, NotIn_H3K27me3_Not_notshared_broad_num), nrow = 2)
  
  # Perform Fisher's exact test
  result_not_shared <- fisher.test(data, alternative = alternative_str)
  
  ##for the shared A
  In_H3K27me3_shared_A_broad_num <- length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  In_H3K27me3_Not_shared_A_broad_num <- length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  
  NotIn_H3K27me3_shared_A_broad_num <- length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  NotIn_H3K27me3_Not_shared_A_broad_num <-  length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR)) + length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  data <- matrix(c(In_H3K27me3_shared_A_broad_num, NotIn_H3K27me3_shared_A_broad_num, In_H3K27me3_Not_shared_A_broad_num, NotIn_H3K27me3_Not_shared_A_broad_num), nrow = 2)
  
  # Perform Fisher's exact test
  result_shared_A <- fisher.test(data, alternative = alternative_str)
  
  ##for the shared I
  In_H3K27me3_shared_I_broad_num <- length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  In_H3K27me3_Not_shared_I_broad_num <- length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR)) + length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  
  NotIn_H3K27me3_shared_I_broad_num <- length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  NotIn_H3K27me3_Not_shared_I_broad_num <-  length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR)) + length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  
  data <- matrix(c(In_H3K27me3_shared_I_broad_num, NotIn_H3K27me3_shared_I_broad_num, In_H3K27me3_Not_shared_I_broad_num, NotIn_H3K27me3_Not_shared_I_broad_num), nrow = 2)
  
  # Perform Fisher's exact test
  result_shared_I <- fisher.test(data, alternative = alternative_str)
  
  
  
  #make the barplot for the proportion compared to the not under H3K27me3
  SS_ACR_H3K27_rice_maize = length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  total_ACR_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR)) + 
    length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR)) + 
    length(unique(shared_not_shared_broad_rice_H3K27me3_dt$rice_ACR))
  
  SS_ACR_not_H3K27_rice_maize = length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  total_ACR_not_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR)) + 
    length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR)) + 
    length(unique(shared_not_shared_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  SH_ACR_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_H3K27me3_dt$rice_ACR))
  SH_ACR_not_H3K27_rice_maize = length(unique(shared_A_broad_bothspe_not_H3K27me3_dt$rice_ACR))
  
  SI_ACR_H3K27_rice_maize = length(unique(shared_I_broad_rice_H3K27me3_dt$rice_ACR))
  SI_ACR_not_H3K27_rice_maize = length(unique(shared_I_broad_rice_not_H3K27me3_dt$rice_ACR))
  
  
  
  ##make a data frame
  my_data <- data.frame(
    SS_ACR_H3K27 = SS_ACR_H3K27_rice_maize,
    total_ACR_H3K27 = total_ACR_H3K27_rice_maize,
    SS_ACR_not_H3K27 = SS_ACR_not_H3K27_rice_maize,
    total_ACR_not_H3K27 = total_ACR_not_H3K27_rice_maize,
    SH_ACR_H3K27 = SH_ACR_H3K27_rice_maize,
    SH_ACR_not_H3K27 = SH_ACR_not_H3K27_rice_maize,
    SI_ACR_H3K27 = SI_ACR_H3K27_rice_maize,
    SI_ACR_not_H3K27 = SI_ACR_not_H3K27_rice_maize
  )
  
  my_data <- t(my_data)
  colnames(my_data) <- c('Count')
  my_data <- as.data.frame(my_data)
  my_data$spe <- target_fl_path_nm
  my_data$not_shared_pval <- result_not_shared$p.value
  my_data$shared_A_pval <- result_shared_A$p.value
  my_data$shared_I_pval <- result_shared_I$p.value
  
  
  
  return(my_data)
  

}



outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_summary_file_add_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  target_species <- gsub('CNS_rice_','',target_fl_path_nm)

  
  ipt_dt <- read.delim(target_fl_path)
  
  colnames(ipt_dt)
  
  ##for the data overlapping with the CNS
  ipt_CNS_dt <- ipt_dt[ipt_dt$OverlapCNS == 'yes',]
  
  dt_CNS <- analyze_ACR_syntenic_H3K27me3(ipt_CNS_dt,target_species,'greater')
  
  ipt_not_CNS_dt <- ipt_dt[ipt_dt$OverlapCNS == 'no',]
  dt_CNS <- analyze_ACR_syntenic_H3K27me3(ipt_not_CNS_dt,target_species,'greater')
  
  
})



#################
##updating 010724
##we will check the overlapping with the CNS data and want to find something difference 
##here we used the database CNS to have a look 
ipt_dir <- 'opt3_check_H3K27me3_in_syntenic_CNS_122023///'
opt_dir <- 'opt3_check_H3K27me3_in_syntenic_CNS_results_122023//'

all_fl_list <- list.files(path= ipt_dir,pattern = 'txt')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_final_summary_file_add_CNS_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)
  
  target_species <- gsub('rice_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path)
  
  
  ##check the proportion under the H3K27me3
  ##how many H3K27me3 will be overlapped with the CNS
  ipt_CNS_dt <- ipt_dt[ipt_dt$OverlapCNS == 'yes',]
  length(unique(ipt_CNS_dt$rice_ACR))
  
  ipt_syn_dt <- ipt_dt[ipt_dt[[paste0(target_species,'_region')]] != 'none',]
  length(unique(ipt_syn_dt$rice_ACR))
  ipt_non_syn_dt <- ipt_dt[ipt_dt[[paste0(target_species,'_region')]] == 'none',]
  length(unique(ipt_non_syn_dt$rice_ACR))
  
  ipt_syn_CNS_dt <- ipt_syn_dt[ipt_syn_dt$OverlapCNS == 'yes',]
  length(unique(ipt_syn_CNS_dt$rice_ACR))
  
  ipt_non_syn_CNS_dt <- ipt_non_syn_dt[ipt_non_syn_dt$OverlapCNS == 'yes',]
  length(unique(ipt_non_syn_CNS_dt$rice_ACR))
  
  ipt_syn_CNS_H3K27me3_dt <- ipt_syn_CNS_dt[ipt_syn_CNS_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  length(unique(ipt_syn_CNS_H3K27me3_dt$rice_ACR))
  
  ipt_syn_CNS_noH3K27me3_dt <- ipt_syn_CNS_dt[ipt_syn_CNS_dt$riceOverH3K27Cate != 'OverlapH3K27me3',]
  length(unique(ipt_syn_CNS_noH3K27me3_dt$rice_ACR))
  
  ipt_non_syn_CNS_H3K27me3_dt <- ipt_non_syn_CNS_dt[ipt_non_syn_CNS_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  length(unique(ipt_non_syn_CNS_H3K27me3_dt$rice_ACR))
  
  ipt_non_syn_CNS_noH3K27me3_dt <- ipt_non_syn_CNS_dt[ipt_non_syn_CNS_dt$riceOverH3K27Cate != 'OverlapH3K27me3',]
  length(unique(ipt_non_syn_CNS_noH3K27me3_dt$rice_ACR))
  
  ##would it be the cns-acr preper
  ##20% of CNS-ACR overlapping H3K27me3 in both syn and non syn
  
  ##what percent of nonCNS-ACR overlapping H3K27me3
  ##in the syntenic region
  ##11% of nonCNS-ACR overlapping H3K27me3
  
  
  ipt_syn_noCNS_dt <- ipt_syn_dt[ipt_syn_dt$OverlapCNS == 'no',]
  length(unique(ipt_syn_noCNS_dt$rice_ACR))
  
  ipt_non_syn_noCNS_dt <- ipt_non_syn_dt[ipt_non_syn_dt$OverlapCNS == 'no',]
  length(unique(ipt_non_syn_noCNS_dt$rice_ACR))
  
  ipt_syn_noCNS_H3K27me3_dt <- ipt_syn_noCNS_dt[ipt_syn_noCNS_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  length(unique(ipt_syn_noCNS_H3K27me3_dt$rice_ACR))
  
  ipt_syn_noCNS_noH3K27me3_dt <- ipt_syn_noCNS_dt[ipt_syn_noCNS_dt$riceOverH3K27Cate != 'OverlapH3K27me3',]
  length(unique(ipt_syn_noCNS_noH3K27me3_dt$rice_ACR))
  
  ipt_non_syn_noCNS_H3K27me3_dt <- ipt_non_syn_noCNS_dt[ipt_non_syn_noCNS_dt$riceOverH3K27Cate == 'OverlapH3K27me3',]
  length(unique(ipt_non_syn_noCNS_H3K27me3_dt$rice_ACR))
  
  ipt_non_syn_noCNS_noH3K27me3_dt <- ipt_non_syn_noCNS_dt[ipt_non_syn_noCNS_dt$riceOverH3K27Cate != 'OverlapH3K27me3',]
  length(unique(ipt_non_syn_noCNS_noH3K27me3_dt$rice_ACR))
  
  ##within syntenic regions, we found H3K27me3 is more likely enriched in the CNS-ACR other than the nonCNS-ACRs
  syn_ACRcns_H3K27me3 <-  length(unique(ipt_syn_CNS_H3K27me3_dt$rice_ACR))
  syn_ACRcns_notH3K27me3 <- length(unique(ipt_syn_CNS_noH3K27me3_dt$rice_ACR))
  syn_ACRnoCNS_H3K27me3 <-  length(unique(ipt_syn_noCNS_H3K27me3_dt$rice_ACR))
  syn_ACRnoCNS_notH3K27me3 <-  length(unique(ipt_syn_noCNS_noH3K27me3_dt$rice_ACR))
  
  data <- matrix(c(syn_ACRcns_H3K27me3, syn_ACRnoCNS_H3K27me3, syn_ACRcns_notH3K27me3, syn_ACRnoCNS_notH3K27me3), nrow = 2)
  
  # Perform Fisher's exact test
  result <- fisher.test(data,alternative = 'greater')
  print(result)
  
  ACRcns_H3K27me3 <- syn_ACRcns_H3K27me3
  ACRcns <- length(unique(ipt_syn_CNS_dt$rice_ACR))
  prop_ACRcns_H3K27me3 <- ACRcns_H3K27me3/ACRcns
  
  ACRnoCNS_H3K27me3 <- syn_ACRnoCNS_H3K27me3
  ACRnoCNS <- length(unique(ipt_syn_noCNS_dt$rice_ACR))
  prop_ACRnoCNS_H3K27me3 <- ACRnoCNS_H3K27me3/ACRnoCNS
  
  ##The percentage of ACRs within H3K27me3 and not within H3K27me3  in genic, proximal, and distal manners for each species pair.
  ##The percentage for each class within H3K27me3 and H3K27me3 absent regions collectively sum to 100%.
  
  ##we will check the non-syntenic region
  nonsyn_ACRcns_H3K27me3 <-  length(unique(ipt_non_syn_CNS_H3K27me3_dt$rice_ACR))
  nonsyn_ACRcns <- length(unique(ipt_non_syn_CNS_dt$rice_ACR))
  prop_nonsyn_ACRcns_H3K27me3 <- nonsyn_ACRcns_H3K27me3/nonsyn_ACRcns
  
  nonsyn_ACRnoCNS_H3K27me3 <- length(unique(ipt_non_syn_noCNS_H3K27me3_dt$rice_ACR))
  nonsyn_ACRnoCNS <- length(unique(ipt_non_syn_noCNS_dt$rice_ACR))
  prop_nonsyn_ACRnoCNS_H3K27me3 <- nonsyn_ACRnoCNS_H3K27me3/nonsyn_ACRnoCNS
  
  
  nonsyn_ACRcns_H3K27me3 <-  length(unique(ipt_non_syn_CNS_H3K27me3_dt$rice_ACR))
  nonsyn_ACRcns_notH3K27me3 <- length(unique(ipt_non_syn_CNS_noH3K27me3_dt$rice_ACR))
  nonsyn_ACRnoCNS_H3K27me3 <-  length(unique(ipt_non_syn_noCNS_H3K27me3_dt$rice_ACR))
  nonsyn_ACRnoCNS_notH3K27me3 <-  length(unique(ipt_non_syn_noCNS_noH3K27me3_dt$rice_ACR))
  data <- matrix(c(nonsyn_ACRcns_H3K27me3, nonsyn_ACRnoCNS_H3K27me3, nonsyn_ACRcns_notH3K27me3, nonsyn_ACRnoCNS_notH3K27me3), nrow = 2)
  result <- fisher.test(data,alternative = 'greater')
  print(result)
  
  ##make a bar plot
  data_dt <- data.frame(
    spe <- c(target_species,target_species,target_species,target_species),
    synCate <- c('Syn','Syn','NonSyn','NonSyn'),
    cate <- c('ACRcns_H3K27me3','ACRnoCNS_H3K27me3','ACRcns_H3K27me3','ACRnoCNS_H3K27me3'),
    prop <- c(prop_ACRcns_H3K27me3,prop_ACRnoCNS_H3K27me3,prop_nonsyn_ACRcns_H3K27me3,prop_nonsyn_ACRnoCNS_H3K27me3)
    
  )
  colnames(data_dt) <- c('Species','synCate','Cate','Prop')
  
  return(data_dt)
  
})

combine_dt <- do.call(rbind,outs)

##make the plot
combine_dt$Species <- factor(combine_dt$Species, levels = c('maize','sorghum'))
combine_dt$synCate <- factor(combine_dt$synCate, levels = c('Syn','NonSyn'))
p <- ggplot(combine_dt, aes(fill=Cate, y=Prop, x=Species)) + 
  geom_bar(stat="identity", position=position_dodge())
#p <- p + geom_text(aes(label = Prop),position = position_stack(vjust = 0.5),vjust = 0.5,size=8)
p <- p+labs(x="\nCell type", y = "Number of feature / Total number of feature\n", cex.lab = 7) ## use \n to set the position
#p <- p+labs(x="\nCell type", y = "Number of feature\n", cex.lab = 7) ## use \n to set the position
p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             text = element_text(size = 30), ##change all the text size
             axis.text.y = element_text(size=40,colour = "black"),
             axis.text.x = element_text(size=40,colour = "black",angle = 45, hjust = 1),
             axis.title.x = element_text(size=40,colour = "black"),
             axis.title.y = element_text(size=40,colour = "black"),
             plot.margin = unit(c(1,1,1,1), "cm"))  ##change the margin more outer
p <- p + facet_wrap(~ synCate)

#p <- p + coord_cartesian(ylim = c(0, 100)) + ##threshold
#  scale_y_continuous(breaks=seq(0, 100, 20))

pdf(paste0(opt_dir,'/','opt_H3K27me3_overlap_cnsACR_nocnsACR_in_syntenic','.pdf'),width = 18,height = 12)
p
dev.off()






##########
##updating 122023
##we will check the overlapping with the rie atlas ACRs
ipt_dir <- 'opt3_check_H3K27me3_in_syntenic_capture_rice_atlas_122123////'
opt_dir <- 'opt3_check_H3K27me3_in_syntenic_capture_rice_atlas_results_122123///'

all_fl_list <- list.files(path= ipt_dir)

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_summary_H3K27me3_add_rice_atlas_acr.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_','',target_fl_path_nm)
  
  
  target_spe1 <- gsub('_to_.+','',target_fl_path_nm) ##maize
  target_spe2 <- gsub('.+_to_','',target_fl_path_nm) ##rice
  
  ipt_dt <- read.delim(target_fl_path)
  
  ##find the variable ACRs under the H3K27me3
  variable_ACR_broad_H3K27me3_dt <- ipt_dt[ipt_dt[[paste0(target_spe1,'_ACRID')]] != 'none' & 
                                            ipt_dt[[paste0(target_spe2,'_ACR')]] == 'none' & 
                                            ipt_dt[[paste0(target_spe1,'_InorNotH3K27m3')]] == 'BroadInflankH3K27me3peak',]
  variable_ACR_broad_H3K27me3_atlasACR_dt <- variable_ACR_broad_H3K27me3_dt[variable_ACR_broad_H3K27me3_dt$Atlas_ACR != 'none',] 
  
  
  ACR_num_variable_ACR_broad_H3K27me3_dt <- length(unique(variable_ACR_broad_H3K27me3_dt[[paste0(target_spe1,'_ACR')]]))
  
  ACR_num_variable_ACR_broad_H3K27me3_atlasACR_dt <- length(unique(variable_ACR_broad_H3K27me3_atlasACR_dt[[paste0(target_spe1,'_ACR')]]))
  
  ACR_num_variable_ACR_broad_H3K27me3_NOT_atlasACR_dt <- ACR_num_variable_ACR_broad_H3K27me3_dt - ACR_num_variable_ACR_broad_H3K27me3_atlasACR_dt
  
  
  variable_ACR_broad_H3K27me3_atlasACR_broadInroot_dt <- variable_ACR_broad_H3K27me3_atlasACR_dt[variable_ACR_broad_H3K27me3_atlasACR_dt$H3K27me3_organ_root_cate == 'BroadInflankH3K27me3peak',]
  length(variable_ACR_broad_H3K27me3_atlasACR_broadInroot_dt$maize_ACR)
  
  variable_ACR_broad_H3K27me3_atlasACR_broadInpanicle_dt <- variable_ACR_broad_H3K27me3_atlasACR_dt[variable_ACR_broad_H3K27me3_atlasACR_dt$H3K27me3_organ_panicle_cate == 'BroadInflankH3K27me3peak',]
  length(variable_ACR_broad_H3K27me3_atlasACR_broadInpanicle_dt$maize_ACR)
  
  variable_ACR_broad_H3K27me3_atlasACR_broadInpanicle_root_dt <- variable_ACR_broad_H3K27me3_atlasACR_dt[variable_ACR_broad_H3K27me3_atlasACR_dt$H3K27me3_organ_panicle_cate == 'BroadInflankH3K27me3peak' | 
                                                                                                           variable_ACR_broad_H3K27me3_atlasACR_dt$H3K27me3_organ_root_cate == 'BroadInflankH3K27me3peak',]
  length(variable_ACR_broad_H3K27me3_atlasACR_broadInpanicle_root_dt$maize_ACR)
  
  
  ##find the variable ACRs not overlap H3K27me3
  variable_ACR_broad_not_H3K27me3_dt <- ipt_dt[ipt_dt[[paste0(target_spe1,'_ACRID')]] != 'none' & 
                                       ipt_dt[[paste0(target_spe2,'_ACR')]] == 'none' & 
                                       ipt_dt[[paste0(target_spe1,'OverH3K27Cate')]] == 'none' &
                                         ipt_dt[[paste0(target_spe1,'_BroOrRes')]] == 'broad' ,]
  
  variable_ACR_broad_not_H3K27me3_atlasACR_dt <- variable_ACR_broad_not_H3K27me3_dt[variable_ACR_broad_not_H3K27me3_dt$Atlas_ACR != 'none',]
   
  ACR_num_variable_ACR_broad_not_H3K27me3_dt <- length(unique(variable_ACR_broad_not_H3K27me3_dt[[paste0(target_spe1,'_ACR')]]))
  
  ACR_num_variable_ACR_broad_not_H3K27me3_atlasACR_dt <- length(unique(variable_ACR_broad_not_H3K27me3_atlasACR_dt[[paste0(target_spe1,'_ACR')]]))

  ACR_num_variable_ACR_broad_not_H3K27me3_NOT_atlasACR_dt <- ACR_num_variable_ACR_broad_not_H3K27me3_dt - ACR_num_variable_ACR_broad_not_H3K27me3_atlasACR_dt
  
  
  variable_ACR_broad_not_H3K27me3_atlasACR_broadInroot_dt <- variable_ACR_broad_not_H3K27me3_atlasACR_dt[variable_ACR_broad_not_H3K27me3_atlasACR_dt$H3K27me3_organ_root_cate == 'BroadInflankH3K27me3peak',]
  length(unique(variable_ACR_broad_H3K27me3_atlasACR_broadInroot_dt$maize_ACR))
  
  variable_ACR_broad_not_H3K27me3_atlasACR_broadInpanicle_dt <- variable_ACR_broad_not_H3K27me3_atlasACR_dt[variable_ACR_broad_not_H3K27me3_atlasACR_dt$H3K27me3_organ_panicle_cate == 'BroadInflankH3K27me3peak',]
  length(unique(variable_ACR_broad_H3K27me3_atlasACR_broadInpanicle_dt$maize_ACR))
  
  variable_ACR_broad_not_H3K27me3_atlasACR_broadInpanicle_root_dt <- variable_ACR_broad_not_H3K27me3_atlasACR_dt[variable_ACR_broad_not_H3K27me3_atlasACR_dt$H3K27me3_organ_panicle_cate == 'BroadInflankH3K27me3peak' | 
                                                                                                                   variable_ACR_broad_not_H3K27me3_atlasACR_dt$H3K27me3_organ_root_cate == 'BroadInflankH3K27me3peak',]
  
  length(unique(variable_ACR_broad_not_H3K27me3_atlasACR_broadInpanicle_root_dt$maize_ACR))
  
  
  
  data <- matrix(c(ACR_num_variable_ACR_broad_H3K27me3_atlasACR_dt, ACR_num_variable_ACR_broad_H3K27me3_NOT_atlasACR_dt, ACR_num_variable_ACR_broad_not_H3K27me3_atlasACR_dt, ACR_num_variable_ACR_broad_not_H3K27me3_NOT_atlasACR_dt), nrow = 2)
  
  # Perform Fisher's exact test
  result <- fisher.test(data,alternative = 'greater')
  result
  
  
  data <- data.frame(
    Name = c("VACR_B_H3K27_atlas", "VACR_B_not_H3K27_atlas"),
    prop = c(ACR_num_variable_ACR_broad_H3K27me3_atlasACR_dt/ACR_num_variable_ACR_broad_H3K27me3_dt, 
             ACR_num_variable_ACR_broad_not_H3K27me3_atlasACR_dt/ACR_num_variable_ACR_broad_not_H3K27me3_dt),
    total = c(ACR_num_variable_ACR_broad_H3K27me3_dt,ACR_num_variable_ACR_broad_not_H3K27me3_dt),
    target = c(ACR_num_variable_ACR_broad_H3K27me3_atlasACR_dt,ACR_num_variable_ACR_broad_not_H3K27me3_atlasACR_dt),
    
    captureH3K27me3_root = c() 
    
    
  )
  
  data$spe <- target_spe1
  data$pval <- result$p.value
  
  return(data)
  
  
})

combine_dt <- do.call(rbind,outs)

p <- ggplot(data=combine_dt, aes(x=spe, y=prop, fill=Name)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_bar(stat="identity") +
  theme(
    plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
    axis.title = element_text(size =20, face="bold"),
    #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size = 15),
    strip.text = element_text(size=10),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  )  
#facet_wrap(~Rice_CellType,nrow = 1)

pdf(paste0(opt_dir,'/','opt_H3K27me3_and_noH3K27me3_overlap_rice_atlas_prop.pdf'),width = 8,height = 6 )
p
dev.off()





##########
##updating 020124
##we will check the exact p value of motif enrichment for the binomial test for different species

ipt_dir <- 'opt3_check_H3K27me3_motif_enrich_020124/////'
opt_dir <- 'opt3_check_H3K27me3_motif_enrich_results_020124////'

ipt_H3K27me3_dir <- 'opt3_check_H3K27me3_motif_enrich_020124/collect_all_H3K27me3_fl_dir_112023/'

all_fl_list <- list.files(path= ipt_dir, pattern = '.txt')

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('opt_binomial_test_plotR_fimopval0.0002_','',target_fl_path_nm)
  target_fl_path_nm <- gsub('.txt','',target_fl_path_nm)

  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  colnames(ipt_dt) <- c('motif','cate','prop')
  
  motif_nm_list <- unique(ipt_dt$motif)
    
  ipt_target_H3K27me3_dt <- read.delim(paste0(ipt_H3K27me3_dir,'/',target_fl_path_nm,'_H3K27me3_cate.txt'),header = F)
  
  head(ipt_target_H3K27me3_dt)
  
  ipt_target_H3K27me3_dt$ACRloc <- paste0(ipt_target_H3K27me3_dt$V1,'_',ipt_target_H3K27me3_dt$V2,'_',ipt_target_H3K27me3_dt$V3)
  ipt_target_H3K27me3_target_dt <- ipt_target_H3K27me3_dt[ipt_target_H3K27me3_dt$V10 == 'BroadInflankH3K27me3peak',]
  
  total_ACRnum <- length(unique(ipt_target_H3K27me3_target_dt$ACRloc))
  
  outs2 <- lapply(motif_nm_list, function(y){
    
    ipt_target_motif_dt <- ipt_dt[ipt_dt$motif == y,]
    
    control_ACR_mean_val <- mean(ipt_target_motif_dt[ipt_target_motif_dt$cate == 'control_ACR',]$prop)
    target_ACR_val <- ipt_target_motif_dt[ipt_target_motif_dt$cate == 'target_ACR',]$prop
    
    sucesses <- target_ACR_val*total_ACRnum
  
    result <- binom.test(as.integer(sucesses), total_ACRnum, p = control_ACR_mean_val)
    
    res <- c()
    res$motif <- y
    res$pval <- result$p.value
    res <- as.data.frame(res)
    
    return(res)
    
  })
  
    
  combine_dt <- do.call(rbind,outs2)
  combine_dt$spe <- target_fl_path_nm
  
  return(combine_dt)
})

combine_all_dt <- do.call(rbind,outs)

write.csv(combine_all_dt,paste0(opt_dir,'/opt_PRE_motif_enrichment_pval.csv'),quote = F)



##########
##updating 020824
##we will plot the len distribution and motif pie chart
ipt_dir <- c('add_check_blasted_region_len_020824/store_all_blasted_output_dir/')
opt_dir <- c('add_check_blasted_region_len_020824')

all_fl_list <- list.files(path= ipt_dir, pattern = '.txt')



##for the rice to maize
generate_pie <- function(df,cate) {
  #pie(df$number)
  #library(RColorBrewer)
  pct <- round(df$count/sum(df$count)*100)
  lbls <- paste(df[[cate]],pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(df$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  return(pie(df$count,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) )
}

library(RColorBrewer)
myPalette <- brewer.pal(10, "Set3") 

outs <- lapply(all_fl_list, function(x){
  
  target_fl_path <- paste0(ipt_dir,'/',x)
  
  target_fl_path_nm <- gsub('.+//','',target_fl_path)
  target_fl_path_nm <- gsub('_ACR_blasted_region_len_motif.txt','',target_fl_path_nm)
  target_fl_path_nm <- gsub('opt_final_','',target_fl_path_nm)
  
  ipt_dt <- read.delim(target_fl_path,header = F)
  head(ipt_dt)
  
  colnames(ipt_dt) <- c('ACR','rChr','rSt','rEd','cate','len')
  dim(ipt_dt)
  
  ipt_acr_blastlen_dt <- unique(ipt_dt[c('ACR','len','cate')])
  
  p <- ggplot(ipt_acr_blastlen_dt, aes(x = len)) +
    geom_histogram(color = "darkblue", fill = "lightblue", bins = 50)+
    coord_cartesian(xlim = c(0,500)) +
    scale_x_continuous(breaks=seq(0,500,10))+
    theme(text = element_text(size = 10)) + 
    theme(
      plot.title = element_text(face="bold.italic",size=15,hjust = 0.5),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=15,angle = 90, vjust = 0.5,hjust = 0.5,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=15,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=10),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )  
  
  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'_histgram.pdf'),width = 16, height = 6)
  print(p)
  dev.off()
  
  
  cate_dt <- as.data.frame(table(ipt_acr_blastlen_dt$cate))
  colnames(cate_dt) <- c('cate','count')

  pdf(paste0(opt_dir,'/opt_',target_fl_path_nm,'_motif_cate.pdf'),width = 6, height = 6)
  print(generate_pie(cate_dt,'cate'))
  dev.off()
  
  res <- c()
  res$spe_pair <- target_fl_path_nm
  res$ACRnum <- nrow(ipt_acr_blastlen_dt)
  
  res <- as.data.frame(res)
  return(res)
  

})

combine_dt <- do.call(rbind,outs)


























##########
##previous to check the evolutionary scores

ipt_broadTobroad_target_rice_dt$cate <- 'rice'
colnames(ipt_broadTobroad_target_rice_dt) <- c('ACR','dist','cate')

ipt_broadTobroad_target_maize_dt <- ipt_broadTobroad_dt[c('maize_ACR','maizeGeneDist')]
ipt_broadTobroad_target_maize_dt$cate <- 'maize'
colnames(ipt_broadTobroad_target_maize_dt) <- c('ACR','dist','cate')

ipt_broadTobroad_target_dt <- rbind(ipt_broadTobroad_target_rice_dt,ipt_broadTobroad_target_maize_dt)
ipt_broadTobroad_target_dt$dist <- as.numeric(ipt_broadTobroad_target_dt$dist)

##use the density plot
p<-ggplot(ipt_broadTobroad_target_dt, aes(x=dist, color=cate)) +
  geom_density(alpha=0.4)
#p <- p + coord_cartesian(xlim = c(0, 2))
#p <- p + scale_x_continuous(breaks=seq(0,2, 0.25))
p <- p + theme(
  plot.title = element_text(size=20,hjust = 0.5),
  axis.title = element_text(size =20),
  #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
  axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1),  ##change the text to italic
  axis.text.y = element_text(size=20,colour = "black"),
  axis.ticks = element_line(size = rel(2.5)),
  axis.ticks.length = unit(0.5, "cm"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  text = element_text(size = 15),
  strip.text = element_text(size=20),
  #legend.title = element_blank(),
  #legend.position = "none",
  plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
)
p
pdf('opt_CS_all_cates_density_addcate.pdf',width = 10,height = 6 )
p
dev.off()












