##updating 080922 this script is to generate peak match using fimo results
##updating 080422 this script is to prepare the motif peak matrix
##it can be obtained from the matchMotifs

library(Matrix)


# arguments
args <- commandArgs(TRUE)

input_fimo_results_fl <- as.character(args[1])

motif_flt_cutoff <- as.numeric(args[2])

input.sp <- as.character(args[3])

organ_nm <- as.character(args[4])

output_dir <- as.character(args[5])

##check whether the motif matcheds file has been built
if (file.exists(paste0(output_dir,'/opt_all_motif_matches_mtx.rds'))){
  message (' - opt_all_motif_matches_mtx.rds is existed, and we do not need to rebuild this file')
}else{

  ##load files
  ACR_motif_dt <- read.delim(input_fimo_results_fl,header = F)
  
  ACR_motif_dt <- ACR_motif_dt[!grepl('scaffold',ACR_motif_dt$V1),]
  
  ACR_motif_dt$peaknm <- paste0('chr',ACR_motif_dt$V1,'_',ACR_motif_dt$V2,'_',ACR_motif_dt$V3)
  ACR_motif_dt$number <- 1
  
  ##filter motif
  if (organ_nm == 'rice'){
    ACR_motif_dt <- ACR_motif_dt[ACR_motif_dt$V11 < motif_flt_cutoff,]
  }else{
    ACR_motif_dt <- ACR_motif_dt[ACR_motif_dt$V10 < motif_flt_cutoff,]
  }
  
  ##aggregate the motif
  ACR_motif_match_dt <- as.data.frame(aggregate(number ~ peaknm + V9,data = ACR_motif_dt,sum))
  
  ##transfer to the matrix
  write.table(ACR_motif_match_dt,paste0(output_dir,'/temp_motif_dt.txt'),quote = F, sep = '\t')
  
  ACR_motif_match_dt <- read.table(paste0(output_dir,'/temp_motif_dt.txt'),stringsAsFactors = T)
  motif_mtx <- sparseMatrix(i=as.numeric(ACR_motif_match_dt$peaknm),
                              j=as.numeric(ACR_motif_match_dt$V9),
                              x=as.numeric(ACR_motif_match_dt$number),
                              dimnames=list(levels(ACR_motif_match_dt$peaknm), levels(ACR_motif_match_dt$V9)))
  mat <- as(motif_mtx, "dgCMatrix")
  
  ##change to the binary
  mat <- as.matrix((mat > 0) + 0)
  
  ##save to the output
  write.csv(mat,paste0(output_dir,'/opt_all_motif_matches_mtx.csv'),quote = F)
  saveRDS(mat,file=paste0(output_dir,'/opt_all_motif_matches_mtx.rds'))
}

##we also need to generate peak by cell mtx
# build counts matrix
message("Loading count matrix ...")
a <- read.table(input.sp,stringsAsFactors=T)
message("transfer to sparse Matrix")
a <- sparseMatrix(i=as.numeric(a$V1),
                  j=as.numeric(a$V2),
                  x=as.numeric(a$V3),
                  dimnames=list(levels(a$V1), levels(a$V2)))
#a <- a[as.character(peaks2$V4),]
message("save mtx")
saveRDS(a,paste0(output_dir,'/','opt_',organ_nm,'_peak_cell_mtx.rds'))
















