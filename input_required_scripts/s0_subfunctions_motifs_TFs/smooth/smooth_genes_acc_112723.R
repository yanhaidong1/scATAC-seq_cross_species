###################################################################################################
## plot marker accessibility scores
###################################################################################################
##updating 062221 add the target cluster argument
##updating 062221 create a fold to store each type of markers in a seperate pdf that is suitable for the many markers within one cell type

# load arguments
args <- commandArgs(trailingOnly=T)
#if(length(args) != 5){stop("Rscript plot_marker_accessibility.R [meta] [gene_activity] [pcs.txt] [markers.bed] [threads]")}

#args    
meta <- as.character(args[1])
geneact <- as.character(args[2])
pcs <- as.character(args[3])
#mark <- as.character(args[4])
threads <- as.numeric(args[4])

##updating 062221
target_cluster <- as.character(args[5])
#plot_each_CT <- as.character(args[7]) ##use yes to initiate this argument
output_dir <- as.character(args[6])
function_script <- as.character(args[7])

##updating 071521
##if we update the markers, we need to re-run this script to obtain the new set of smooth marker file
#onlysmooth_marker <- as.character(args[10]) ##yes or no
#run_denovo <- as.character(args[11]) ##yes or no

# load functions
#source("functions.plot_marker_accessibility.R")
source(function_script)


# load data
if(file.exists(paste0(output_dir,"/GAobj.rds"))){
    dat <- readRDS(paste0(output_dir,"/GAobj.rds"))
}else{
    dat <- loadData(meta, pcs, geneact,target_cluster)
    ##updating 062221 save to the GAobj
    saveRDS(dat,paste0(output_dir,"/GAobj.rds"))
}
b.meta <- dat$b
activity.all <- dat$activity
h.pcs1 <- dat$h.pcs
#marker.info.dat <- dat$marker.info


run_smooth(b.meta, 
           activity.all,
           h.pcs1,
           output_dir,
           threads=threads,
           output="all")



# iterate over each major cluster
##rm smooth.markers=F
#out <- runMajorPriori(b.meta, activity.all, h.pcs1, marker.info.dat, plot_each_CT,marker_opt_dir,output_dir,onlysmooth_marker,threads=threads)
#saveRDS(out,paste0(output_dir,'/out.rds'))
##since the smooth function cannot handel big data we need to consider not to open this function
#if (run_denovo == 'yes'){
#  runMajorDeNovo(output_dir,out$b, out$activity, out$impute.activity, marker.info.dat,output="all")
#}


