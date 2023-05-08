#setwd("/home/jarriaga/SGCCA_subsamples")
setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/")
subtypes = c("TRU", "prox.-prolif.", "prox.-inflam", "normal")

for (subtype in subtypes) {
  commands = paste("Rscript 2_3SGCCA_Final.R", subtype)
  system(commands)
}
