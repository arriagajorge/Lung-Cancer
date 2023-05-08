setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/")
subtypes = c("TRU", "prox.-prolif.", "prox.-inflam", "normal")

for (subtype in subtypes) {
  commands = paste("Rscript 2_5Joinsubsamples.R", subtype)
  system(commands)
}