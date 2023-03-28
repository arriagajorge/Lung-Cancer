setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD")

penalty_cpgs <- 0.01
penalty_transcripts <- 0.01
penalty_mir <- 0.01

# ' Read files in correct form
# '@param dataTable String name of the file. You add .csv, .tsv, etc.
# '@return data_table data_frame with the correct form
readMtrx <- function(dataTable){
  df <- data.table::fread(dataTable)
  return(as.matrix(df[,2:ncol(df)], rownames=df$V1))
}

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

transposeData <- function(df){
  data = list()
  data[[1]] = t(df[[1]])
  data[[2]] = t(df[[2]])
  data[[3]] = t(df[[3]])
  names(data)=c("CpGs","transcripts","miRNAs")
  return(data)
}

library(data.table)
library(parallel)

setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/Normi/")
files=list.files()
files=files[grep("eigenN",files)]
sizes = c(normal=5, 'prox.-inflam'=69 ,'prox.-prolif'=48,TRU=71)
sizes

#funtion to transform data

# run
i=lapply(sizes,function(x) c(1,sample(2:x,4)))
data=lapply(1:4,function(x) fread(files[x],select=i[[x]]))
data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
data=do.call(cbind,data)

#separate omics
first_letter <- substr(rownames(data),1 ,1)
newData = list()
for (letter in unique(first_letter)) {
  newData[[letter]] <- data[first_letter == letter,] 
}
names(newData)=c("CpGs","transcripts","miRNAs")

# confirm separates well
# head(rownames(newData$CpGs)); tail(rownames(newData$CpGs))
# head(rownames(newData$transcripts)); tail(rownames(newData$transcripts))
# head(rownames(newData$miRNAs)); tail(rownames(newData$miRNAs))

# data=apply(cbind(c(1,417032,427975),c(417031,427974,428180)),1,
#            function(x) data[x[1]:x[2],])

describe(transposeData(newData),penalty_cpgs,penalty_transcripts,penalty_mir)
descr <- data.frame()
set.seed(11)
for (var in 1:10) {
  # run
  i=lapply(sizes,function(x) c(1,sample(2:x,4)))
  data=lapply(1:4,function(x) fread(files[x],select=i[[x]]))
  data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
  data=do.call(cbind,data)
  
  #separate omics
  first_letter <- substr(rownames(data),1 ,1)
  newData = list()
  for (letter in unique(first_letter)) {
    newData[[letter]] <- data[first_letter == letter,] 
  }
  names(newData)=c("CpGs","transcripts","miRNAs")
  
  descr <- rbind(descr, describe(transposeData(newData),
                                penalty_cpgs,penalty_transcripts,penalty_mir))
}

#lapply(1:10, function(x){})
file=paste(penalty_cpgs,penalty_transcripts,penalty_mir,"tsv",sep='.')
write.table(descr,file,sep='\t',quote=F,row.names=F)