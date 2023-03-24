#!/usr/bin/env Rscript
# for the moment not, we load data manually
# args=commandArgs(trailingOnly=TRUE)
# penalty_cpgs=as.numeric(args[1])
# penalty_transcripts=as.numeric(args[2])
# penalty_mir=as.numeric(args[3])

setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD")

# ' Read files in correct form
# '@param dataTable String name of the file. You add .csv, .tsv, etc.
# '@return data_table data_frame with the correct form
readMtrx <- function(dataTable){
  df <- data.table::fread(dataTable)
  return(as.matrix(df[,2:ncol(df)], rownames=df$V1))
}
penalty_cpgs <- readMtrx("CpGs.mtrx")
penalty_transcripts <- readMtrx("transcripts.mtrx")
penalty_mir <- readMtrx("miRNAs.mtrx")

head.matrix(penalty_cpgs)
head.matrix(penalty_transcripts)
head.matrix(penalty_mir)

subtype <- readMtrx("subtypeLUADord.tsv")
subtype_subtype <- as.factor(subtype[,"subtype"])
summary(subtype_subtype)
library(mixOmics)
#take model descriptors<-----------------recycled
describe=function(data,pc,pt,pm){
  #subsample observations
  #data=lapply(data,function(x) x[sample(1:n,subn),])
  resus=wrapper.sgcca(data,penalty=c(pc,pt,pm),scale=T,
                      scheme="centroid")
  #get results description
  description=as.data.frame(do.call(rbind,resus$AVE$AVE_X))
  description$nfeatures=sapply(resus$loadings,function(x) sum(x!=0))
  description$omic=rownames(description)
  description$penalty=resus$penalty
  colnames(description)[1]="AVE"
  return(description)
}

library(mixOmics)

######################DATA TO LIST OF MATRIXES PER MOLECULAR LEVEL
library(data.table)
library(parallel)

setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/multiomics/")
files=list.files()
# files=files[grep("eigenN",files)]
sizes = c(normal=5, 'prox.-inflam'=69 ,'prox.-prolif'=48,TRU=71)
sizes

cl <- parallel::makeCluster(4)
clusterExport(cl, c("describe","penalty_cpgs","penalty_transcripts",
                    "penalty_mir","wrapper.sgcca","files","sizes","fread"))

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
resus=do.call(rbind,parLapply(cl,1:10,function(x) {
  #1 more than the samples coz of the rownames 
  i=lapply(sizes,function(x) c(1,sample(2:x,5)))
  data=lapply(1:4,function(x) fread(files[x],select=i[[x]]))
  data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
  data=do.call(cbind,data)
  #separate omics
  data=apply(cbind(c(1,417032,427975),c(417031,427974,428180)),1,
             function(x) data[x[1]:x[2],])
  names(data)=c("CpGs","transcripts","miRNAs")
  describe(data,penalty_cpgs,penalty_transcripts,penalty_mir)}))
stopCluster(cl)
file=paste(penalty_cpgs,penalty_transcripts,penalty_mir,"tsv",sep='.')
write.table(resus,file,sep='\t',quote=F,row.names=F)