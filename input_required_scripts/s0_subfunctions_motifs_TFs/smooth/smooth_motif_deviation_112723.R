###################################################################################################
## project motif scores on umap embeddings
###################################################################################################
##updating 122122 we will allow the target_cluster not to be specified item
##add the target_cluster 
##add annotation
##change the LouvainClusters to LouvainClusters_afthm

# setup arguments
##provide the clusters file
##provide the motifs score
##we provide the same dimenstion file SVD.txt
##privide the prefix 
args <- commandArgs(T)
#if(length(args) != 4){
#	stop("Rscript project.motif.UMAP.R <reducedDims> <motif.scores> <pcs> <prefix>")
#}

# variables
input <- as.character(args[1])
motifs <- as.character(args[2])
pc.dat <- as.character(args[3])
prefix <- as.character(args[4])
target_cluster <- as.character(args[5])
output_dir <- as.character(args[6])

# load libs
library(RColorBrewer)
library(viridis)
library(scales)
library(pheatmap)
library(gplots)
library(RANN)
library(Matrix)

# functions
clusterMotifs <- function(x, prefix="all"){
    mat <- t(x)
    cor.mat <- as.matrix(cor(mat, method="spearman"))
    colnames(cor.mat) <- colnames(mat)
    rownames(cor.mat) <- colnames(mat)
    id <- data.frame(do.call(rbind, strsplit(colnames(mat), "_")))
    id$type <- substr(id$X2, 1, 3)
    rownames(id) <- rownames(cor.mat)
    id$X1 <- NULL
    id$X2 <- NULL
    
    # plot motif correlations
    cols <- magma(100)
    t.col <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(id$type)))
    pdf(paste0(prefix,".motif.correlations.pdf"), width=10, height=10)
    write.table(cor.mat, file=paste0(prefix,".motif.correlations.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    heatmap.2(cor.mat, col=cols, breaks=c(seq(-1,1,length.out=101)),  
              scale="none", trace="none",
              symm=T, useRaster=T,
              ColSideColors=t.col[factor(id$type)], 
              RowSideColors=t.col[factor(id$type)],
              distfun=function(x){as.dist(1-x)}, 
              hclustfun=function(x){hclust(x, method="ward.D2")},
              labRow=NA, labCol=NA)
    dev.off()
    
}
smoothByLouvain <- function(x, y, pc,target_cluster){

    # iterate over cluster
    #if (target_cluster == 'LouvainClusters_afthm'){
    #  clusts <- unique(y$LouvainClusters_afthm)
    #}
    #if (target_cluster == 'LouvainClusters'){
    #  clusts <- unique(y$LouvainClusters)
    #}
  
    clusts <- unique(y[[target_cluster]])
      
    motifs.ids <- rownames(x)
    out <- lapply(clusts, function(z){
        message(" - smoothing motifs for cluster ",z)
        
        #if (target_cluster == 'LouvainClusters_afthm'){
        #  lc <- subset(y, y$LouvainClusters_afthm==z)
        #}
        #if (target_cluster == 'LouvainClusters'){
        #  lc <- subset(y, y$LouvainClusters==z)
        #}
      
      
        lc <- subset(y, y[[target_cluster]]==z)
      
        #lc <- subset(y, y$LouvainClusters_afthm==z)
        sub.pc <- pc[rownames(lc),]
        sub.motif <- x[,rownames(lc)]
        sub.motif <- t(apply(sub.motif, 1, function(zz){
            if(sum(as.numeric(is.na(zz))) == length(zz)){
                zz <- rep(0, length(zz))
            }else{
                zz[is.na(zz)] <- mean(zz, na.rm=T)
            }
            return(zz)
        }))
        out.sub <- as.matrix(smoothData(sub.motif, k=10, step=3, npcs=ncol(sub.pc), rds=sub.pc, cl=z))
        out.sub <- out.sub[motifs.ids,]
        return(out.sub)
    })
    out <- do.call(cbind, out)
    out <- as.matrix(t(scale(t(out))))
    return(out)
}
smoothData <- function(x, k=15, step=3, npcs=15, cleanExp=F, df=NULL, rds=NULL, cl=NULL){
    
    # input
    data.use <- x
    
    # verbose
    if(!is.null(rds)){
        
        if(!is.null(df)){
            message("   * using UMAP manifold for smoothing ...")
            pcs <- df[,c("umap1","umap2")]
        }else{
            message("   * using prior PC space as manifold ...")
            pcs <- rds[colnames(x),c(1:npcs)]
        }
    }else{
        
        stop(" ! must pass PCs to smooth graph !")
    }
    
    # get KNN
    message("   * finding knn graph ...")
    knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
    j <- as.numeric(x = t(x = knn.graph))
    i <- ((1:length(x = j)) - 1) %/% k + 1
    edgeList = data.frame(i, j, 1)
    A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
    
    # Smooth graph
    message("   * smoothing graph ...")
    A = A + t(A)
    A = A / Matrix::rowSums(A)
    step.size = step
    if(step.size > 1){
        for(i in 1:step.size){
            message("     ~ step ",i)
            A = A %*% A
        }
    }
    
    # smooth data
    message("   * smoothing motif deviations ...")
    out <- t(A %*% t(data.use))
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    out <- as.matrix(out)
    its <- 0
    out <- t(apply(out, 1, function(z){
        its <<- its + 1
        num.na <- sum(as.numeric(is.na(z)))
        if(num.na == length(z)){
            message(" !! motif ",rownames(out)[its], " NA for Louvain cluster = ", cl)
            z <- rep(0, length(z))
        }else{
            z[is.na(z)] <- mean(z, na.rm=T)
        }
        return(z)}
    ))
    out[is.na(out)] <- mean(out, na.rm=T)
    
    # return sparse Matrix
    return(out)
}

# load data
message("loading data ...")
#a <- read.table(input)
a <- read.delim(input,row.names = 1)

b <- t(as.matrix(read.table(motifs)))
pcs <- read.table(pc.dat)
ids <- intersect(rownames(a), colnames(b))
ids <- intersect(ids, rownames(pcs))
a <- a[ids,]
b <- b[,ids]
pcs <- pcs[ids,]

# smooth scores
message("smoothing data ...")
b <- smoothByLouvain(b, a, pcs,target_cluster)
write.table(t(b), file=paste0(output_dir,'/',prefix,".smoothed_motifs.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# find # of correlated motifs
#clusterMotifs(b, prefix=prefix)

##we will close this function
# project motifs
#num_motifs <- nrow(b)
#message("collecting top ",num_motifs, " motifs ...")
#motif.scores <- apply(b, 1, function(x){var(x, na.rm=T)})
#motif.scores <- motif.scores[order(motif.scores, decreasing=T)]
#df <- as.data.frame(motif.scores)
#write.table(df, file=paste0(prefix,".variance.txt"), quote=F, row.names=T, col.names=T, sep="\t")
#topM <- head(motif.scores, n=num_motifs)
#motif.sub.scores <- b[rownames(b) %in% names(topM),]
#motif.sub.scores <- motif.sub.scores[,rownames(a)]

# plot
#message("plotting motif scores on umap embeddings ...")
#pdf(paste0(prefix,".motif.projections.UMAP.pdf"), width=15, height=40)
#nums <- ceiling(num_motifs/5)
#tnums <- nums*5
#layout(matrix(c(1:tnums), nrow=nums, ncol=5, byrow=T))
#par(mar=c(2,2,0.5,0.5))
#for(i in 1:num_motifs){
#    cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "goldenrod2", "firebrick3"), bias=1)(100)
#    scores <- as.numeric(motif.sub.scores[i,])
#    quants <- quantile(scores, c(0.01, 0.99))
#    scores[scores > quants[2]] <- quants[2]
#    scores[scores < quants[1]] <- quants[1]
#    min.q <- quants[1] - abs(quants[1]*0.1)
#    max.q <- quants[2] + abs(quants[2]*0.1)
#    col.vec <- cols[cut(scores, breaks=seq(from=min.q, to=max.q, length.out=101))]
#    plot(a$umap1, a$umap2, pch=16, col=col.vec, cex=0.3, 
#         xlab="", ylab="", xaxt="n", yaxt="n", main=rownames(motif.sub.scores)[i],
#         bty="none")
#}
#dev.off()

# plot by louvain cluster
#for (j in unique(a$LouvainClusters_afthm)){
#    sub.a <- subset(a, a$LouvainClusters_afthm==j)
#    motif.sub.Ascores <- motif.sub.scores[,rownames(sub.a)]
#    message("plotting motif scores on umap embeddings ...")
#    pdf(paste0(prefix,".Louvain",j,".motif.projections.UMAP.pdf"), width=15, height=40)
#    nums <- ceiling(num_motifs/5)
#    tnums <- nums*5
#    layout(matrix(c(1:tnums), nrow=nums, ncol=5, byrow=T))
#    par(mar=c(2,2,0.5,0.5))
#    for(i in 1:num_motifs){
#        cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "goldenrod2", "firebrick3"), bias=1)(100)
#        scores <- as.numeric(motif.sub.Ascores[i,])
#        quants <- quantile(scores, c(0.01, 0.99))
#        scores[scores > quants[2]] <- quants[2]
#        scores[scores < quants[1]] <- quants[1]
#        min.q <- quants[1] - abs(quants[1]*0.1)
#        max.q <- quants[2] + abs(quants[2]*0.1)
#        col.vec <- cols[cut(scores, breaks=seq(from=min.q, to=max.q, length.out=101))]
#        plot(sub.a$umapsub_1, sub.a$umapsub_2, pch=16, col=col.vec, cex=0.4, 
#             xlab="", ylab="", xaxt="n", yaxt="n", main=rownames(motif.sub.scores)[i],
#             bty="none")
#    }
#    dev.off()
#}
    









