#setwd("/home/jarriaga/SGCCA_subsamples")
setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/")
subtypes = c("TRU", "prox.-prolif.", "prox.-inflam")

for (subtype in subtypes) {
  commands = paste("Rscript 1_6mfa.R", subtype)
  system(commands)
}