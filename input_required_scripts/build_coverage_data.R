##this script is to build a coverage dataframe that is used for the plotting in R
## calcTSSenrich.R


# load arguments
args <- commandArgs(T)


input_acr_mtx_fl <- as.character(args[1])
input_output_dir <- as.character(args[2])
prefix <- as.character(args[3])

acr_dt <- read.delim(input_acr_mtx_fl)
ave <- colMeans(acr_dt)
ave_dt <- as.data.frame(ave)
sd <- apply(acr_dt,2,sd)
n <- nrow(acr_dt)
error <- qnorm(0.975) * sd/sqrt(n)
ave_dt$error <- error 
write.csv(ave_dt,paste0(input_output_dir,'/opt_',prefix,'_ave_motif_cover_dt.csv'),quote=F)