setwd("/home/mdiaz/workspace/LUSC/selectedfeatures/")
#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
fun="GO:0035909"
subty="basal"


#get the components linked to the function
if(length(grep("GO",fun))>0){
  enrich=read_tsv("BP.enrichment")
}else{
  enrich=read_tsv("KEGG.enrichment")	
}
comp=enrich%>%filter(ID==fun&subtype==subty)%>%
  dplyr::select(component)%>%unlist

#get the features selected in those components
selected=read_tsv(paste(subty,"stable",sep='.'))
features=selected%>%filter(component==unlist(comp))%>%
  distinct(variable)%>%unlist
#print(paste("Function has",length(features),"features associated",sep=' '))

#get the data
data=data.table::fread(paste("/home/mdiaz/workspace/LUSC/",paste(subty,"eigenNormi",sep='.'), sep=""))
data=data[data$V1%in%features,]
write_tsv(data,paste(fun,subty,"mtrx",sep='.'),col_names=F)#needed for puma

write_tsv(data, "GO.tsv")
