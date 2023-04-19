# iterate arguments in the list and execute 2_4SGCCA_subsamples.R with the respective arg
setwd("/home/jarriaga/SGCCA_subsamples")
subtype = "TRU"
#pb <- txtProgressBar(min = start, max = end, style = 3)
start = 16
end = 30
for (arg in start:end) {
  #setTxtProgressBar(pb, arg) #update progress bar
  commands = paste(paste("Rscript 2_4SGCCA_subsamples_comm.R", subtype), arg)
  system(commands)
}
