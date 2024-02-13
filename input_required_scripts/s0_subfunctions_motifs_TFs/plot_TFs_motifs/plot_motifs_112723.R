###################################################################################################
###################################################################################################
##                                                                                               ##
##                  functions for plotting accessibility of markers from cicero                  ##
##                                                                                               ##
###################################################################################################
###################################################################################################

##only plot the defined markers

# load libraries
library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
library(varistran)
library(edgeR)
library(parallel)
library(png)
library(phytools)


###################################################################################################
###################################################################################################
###################################################################################################

args <- commandArgs(trailingOnly=T)

all_gene_impute_acc_sparse_fl <- as.character(args[1])
meta_fl <- as.character(args[2])
markers <- as.character(args[3])


#GAobj_rds <- as.character(args[3]) ##GAobj.rds
#markers <- as.character(args[4])
output_dir <- as.character(args[4])

lim <- as.numeric(args[5])




plot.act.scores    <- function(meta_fl,output_dir,
                               act_fl=act_sparse_file, 
                               info_fl=NULL, 
                               top=NULL,
                               logT=F,
                               marker.dist=NULL,
                               outname="markerActivityScores.pdf", 
                               lim=0.95){
  
  
  
  if(file.exists(paste0(output_dir,'/temp_activity_motif_mtx.rds'))){
    acts <- readRDS(paste0(output_dir,'/temp_activity_motif_mtx.rds'))
  }else{
    ##transfer sparse to mtx
    activity <- read.table(act_fl,stringsAsFactors = T)
    acts <- sparseMatrix(i=as.numeric(activity$V1),
                         j=as.numeric(activity$V2),
                         x=as.numeric(activity$V3),
                         dimnames=list(levels(activity$V1), levels(activity$V2)))
    saveRDS(acts,paste0(output_dir,'/temp_activity_motif_mtx.rds'))
  }
  
  
 
  
  ##read the meta file
  df <- read.delim(meta_fl,row.names = 1)
  
  ##read the marker
  info <- read.delim(info_fl)

  
  ##corresponding the row and col of motifs, acc, and meta dt
  intersect_cells <- intersect(rownames(df),colnames(acts))
  df <- df[intersect_cells,]
  acts <- acts[,intersect_cells]
  
  ##change the row name of acts
  rownames(acts) <- gsub('_.+','',rownames(acts))
  

  # prep data
  df <- df[rownames(df) %in% colnames(acts),]
  acts <- acts[,which(rownames(df) %in% colnames(acts))]
  
  # reorder rows
  rownames(info) <- info$geneID
  info <- info[order(info$type),]
  info.genes <- rownames(info)
  act.genes <- rownames(acts)
  rd.cells <- rownames(df)

  
  # common genes
  common <- intersect(info.genes, act.genes)
  info <- info[which(rownames(info) %in% common),]
  info.ordered <- rownames(info)
  sub.scores <- acts[info.ordered,]
  gids <- info.ordered
  
  # setup plot size
  nrows <- ceiling(length(gids)/6)
  totals <- nrows*6
  ratio <- nrows/6
  
  # params
  png(file=paste0(output_dir,'/',outname), width=12, height=ratio*12, units="in", res=500, type="cairo")
  layout(matrix(c(1:totals), ncol=6, byrow=T))
  par(mar=c(2,2,1,1))
  
  # adjust cluster IDs
  message("begin plotting pre-defined markers...")
  for (i in 1:length(gids)){
    
    # copy meta data
    gene.index <- which(rownames(sub.scores) == gids[i])
    acv <- sub.scores[gene.index,]
    
    acv <- rescale(acv, c(-1, 1))
    
    # set up plot cols/sizes
    orderRow <- order(acv, decreasing=F)
    #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
    #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
    #cols <- inferno(100)
    #cols <- plasma(100)
    #cols <- colorRampPalette(rev(c(brewer.pal(11, "RdYlGn")[2:11])), bias=0.7)(100)

    #cols <- colorRampPalette(rev(c(brewer.pal(11, "RdBu")[1:11])), bias=1)(100)
    cols <- colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","gray94","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")),bias = 0.8)(100)

    #cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
    acv <- as.numeric(acv[orderRow])
    if(logT==T){
      acv <- log2(acv+1)
    }
    
    ##Note: important here, we will allow the df2 to have the same order as the acv that has acc value from small to large
    df2 <- df[orderRow,]
    ##change na to be -1 other than 0
    acv[is.na(acv)] <- -1
    acv[is.infinite(acv)] <- -1
    #upper.lim <- quantile(acv, lim)
    #acv[acv > upper.lim] <- upper.lim
    if(!is.null(marker.dist)){
      message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
      colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
    }else{
      min.acv <- min(acv) - (1e-6*min(acv))
      max.acv <- max(acv) + (1e-6*max(acv))
      message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
      if(min.acv == max.acv){
        next
      }
      
      ##Note: for each color it corresponds to one range of color: like color 1 (-1,-0.98]  and color 2 (-0.98,-0.96]
      ##Note: it allows the acv to be cut into different range, and allows each range to correspond one color
      ##Note: here we allow the acc value within different range 
      colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
    }
    colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
    colvec[is.na(colvec) & acv == -1] <- cols[1]
    ##change 0 to be -1
    #sizes <- rescale(acv, c(0.25, 0.3))
    
    # plot
    plot(df2$umap1, df2$umap2, col=colvec,
         #main=paste(info$name[i],info$type[i],sep="-"),
	       main=paste(info$common[i],info$type[i],sep="-"),
         #main=info$name[i],
         cex.main=0.6,
         xlab="", ylab="", bty="n",
         xaxt="n", yaxt="n", pch=16, cex=0.25)
    
    ##Note: we will generate a bar that shows exact 100 colors other than the range exactly corresponding to the real plotting.
    #add.color.bar(3, cols, title='',subtitle='', lims=NULL,prompt=F,x=min(df2$umap1+8),y=min(df2$umap2)+1,cex=0.5)
    
  }
  
  # turn device off
  dev.off()
  
}




plot.act.scores(meta_fl,output_dir,
                act_fl=all_gene_impute_acc_sparse_fl,
                info_fl=markers,
                logT=F,
                lim=lim,
                marker.dist=NULL,
                outname=paste0("combined.impute.known.Markers.motif.png"))

















