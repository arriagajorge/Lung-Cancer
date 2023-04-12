#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
# subtype=args[1]
setwd("/run/media/jorge/Data/INMEGEN/Lung-Cancer/pruebas/SGCCA")

subtype = "TRU"
#ncomp=as.numeric(args[2])

library(igraph)
library(mixOmics)
library(data.table)
########################DATA
data=fread(paste("/run/media/jorge/Data/INMEGEN/Lung-Cancer/pruebas/SGCCA",
                 paste(subtype,"eigenNormi",sep='.'),sep='/'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
first_letter <- substr(rownames(data),1 ,1)
newData = list()
for (letter in unique(first_letter)) {
  newData[[letter]] <- data[first_letter == letter,] 
}
names(newData)=c("CpGs","transcripts","miRNAs")
data <- newData
penalty=c(CpGs=0.01,transcripts=0.01,miRNAs=0)#output of choose_penalty.R

transposeData <- function(df){
  data = list()
  data[[1]] = t(df[[1]])
  data[[2]] = t(df[[2]])
  data[[3]] = t(df[[3]])
  names(data)=c("CpGs","transcripts","miRNAs")
  return(data)
}
########################THE SGCCA
ncomp=nrow(data$miRNAs)-1#the last comp has all loadings>0
final=wrapper.sgcca(X=transposeData(data),penalty=penalty,scale=F,
                    scheme="centroid",ncomp=ncomp)#ncomp to explain 50% of transcripts matrix according to mfa.R
#get selected features
selected=lapply(final$loadings,function(y) 
  apply(y,2,function(x) x[x!=0]))
selected=as.data.frame(do.call(rbind,lapply(selected,function(y) 
  do.call(rbind,lapply(1:length(y),function(x) 
    cbind(names(y)[x],y[[x]],names(y[[x]])))))))
colnames(selected)=c("component","final","variable")

#####PLOT LOADINGS
library(ggplot2)
library(gridExtra)

selected$omic=substr(selected$variable,1,1)
selected$omic=gsub("E","transcripts",
                   gsub("h","miRNAs",gsub("c","CpGs",selected$omic)))
selected$final=as.numeric(as.character(selected$final))
png(paste(subtype,"loadings.png",sep='-'))
ggplot(selected,aes(x=omic,y=final))+
  geom_boxplot()+ylab("loading")+theme(text=element_text(size=18))
dev.off()

initial=wrapper.sgcca(X=data,penalty=rep(1,3),scale=F,
                      scheme="centroid",ncomp=ncomp)#ncomp to explain 50% of transcripts matrix according to mfa.R
rbind(rowSums(do.call(rbind,initial$AVE$AVE_X)),
      rowSums(do.call(rbind,final$AVE$AVE_X))) 
#          CpGs transcripts    miRNAs
#[1,] 0.7512414   0.5791811 0.5558149
#[2,] 0.6099801   0.5408933 0.5326386
selected=lapply(unique(selected$omic),function(x) 
  selected[selected$omic==x,])
temp=lapply(1:3,function(y) apply(selected[[y]],1,function(x) 
  initial$loadings[[y]][x[3],x[1]]))
selected=do.call(rbind,selected)
selected$initial=unlist(temp)
plots=lapply(unique(selected$omic),function(x)
  ggplot(selected[selected$omic==x,],
         aes(y=final,x=initial))+geom_point()+ggtitle(x)+
    theme(text=element_text(size=18)))
png(paste(subtype,"loadings_change.png",sep='-'))
grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

write.table(selected,paste(subtype,"selected",sep='.'),sep='\t',
            quote=F,row.names=F)

#net per enriched component
#source("function_networkAlt.R")
#g=network(final,comp=list(CpGs=7,transcripts=7,miRNAs=7),
#	blocks=1:3)$gR
#temp=as.data.frame(get.edgelist(g))
#temp$cor=E(g)$weight