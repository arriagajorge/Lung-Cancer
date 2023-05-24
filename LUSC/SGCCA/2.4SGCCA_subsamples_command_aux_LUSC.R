# iterate arguments in the list and execute 2_4SGCCA_subsamples.R with the respective arg
setwd("/home/mdiaz/workspace/LUSC")
subtype = "basal"
#pb <- txtProgressBar(min = start, max = end, style = 3)
start = 1
end = 100
for (arg in start:end) {
  print(arg)
  #setTxtProgressBar(pb, arg) #update progress bar
  commands = paste(paste("Rscript 2.4SGCCA_subsamples_comm_LUSC.R", subtype), arg)
  system(commands)
}