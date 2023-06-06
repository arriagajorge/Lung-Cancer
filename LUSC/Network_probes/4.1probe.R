setwd("/home/mdiaz/workspace/LUSC/selectedfeatures/")
#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
###According to the subtype use 1 of the next samples
# fun="GO:0035909"
# subty="basal"
#fun="GO:1905205"
#subty="secretory"
# fun="GO:0002040"
# subty="classical"
# fun="GO:0051084"
# subty="primitive"

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

###CHOOSE ONE 
# write_tsv(data, "probe_basal.tsv")
# write_tsv(data, "probe_classical.tsv")
# write_tsv(data, "probe_secretory.tsv")
# write_tsv(data, "probe_primitive.tsv")
