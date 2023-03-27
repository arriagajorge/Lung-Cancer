setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/")
#!/usr/bin/env Rscript
library(data.table)
subtypeLUAD=read.table("subtypeLUAD.tsv",header=T,sep='\t')
# revert comments for unnormalized data
# expre=fread("RNAseqnormalized.tsv")
# miR=fread("miRNAseqNormi.tsv")
# methy=fread("methyM.tsv")
files=list.files("./multiomics",full.names=T)
#files=files[grep("eigenNormi",files)]
data=lapply(files,fread)
#expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
#miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
#methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)

data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(data)=gsub("./multiomics/","",files)
names(data)=gsub(".tsv","",names(data))
names(data)=gsub("RNAseqnormalized","transcripts",
                 gsub("miRNAseqNormi","miRNAs",gsub("methyM","CpGs",names(data))))
print(sapply(data,dim))
# CpGs miRNAseq transcripts
# [1,] 198709      206       10943
# [2,]    193      193         193
#print(sapply(data,function(x) head(rownames(x))))

#choose methy order
#subtype=subtype[order(match(subtype$samples,colnames(methy))),]
#expre=expre[,order(match(colnames(expre),subtype$samples))]
#miR=miR[,order(match(colnames(miR),subtype$samples))]
subtype=subtypeLUAD[order(match(subtypeLUAD$samples,colnames(data$CpGs))),]
data[2:3]=lapply(data[2:3],function(x)
  x[,order(match(colnames(x),subtype$samples))])
print(names(data))
#"CpGs"        "miRNAseq"    "transcripts"

subtypeLUAD$subtype <- as.factor(subtypeLUAD$subtype)
#data per subtype
concatenated=lapply(levels(subtypeLUAD$subtype),function(x) 
  list(CpGs=data$CpGs[,subtypeLUAD$subtype==x],
       transcripts=data$transcripts[,subtypeLUAD$subtype==x],
       miRNA=data$miRNAs[,subtypeLUAD$subtype==x]))

names(concatenated)=levels(subtypeLUAD$subtype)


##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
print(sapply(concatenated,dim))
# normal prox.-inflam prox.-prolif.    TRU
# [1,] 428180       428180        428180 428180
# [2,]      5           69            48     71
lapply(1:4,function(x) write.table(concatenated[[x]],
                                   paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))
##################################### SO FAAAAAAAAR
#to go back
#apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
#	function(x) data[x[1]:x[2],])

#next check PCs per subtype & data→→→→→→→→→→→→→→mfa.R
#next check inter-omic correlation→→→→→→→corrKernel.R<----optional
#drop near zero var features???????

###########################
#join omics
#!/usr/bin/env Rscript
#library(data.table)
files=list.files()
files=files[grep("mtrx",files)]
data=lapply(files,fread)
#data__ <- data[-1]
# data=lapply(data__,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
data=do.call(cbind,data)
data_ = data
data_=apply(cbind(c(1,417032,427975),c(417031,427974,428180)),1,
            function(x) data[x[1]:x[2],])
names(data_)=c("CpGs","transcripts","miRNAs")
write.table(data_[[1]],paste(names(data_)[1],"mtrx",sep='.'),sep='\t',quote=F)
################################ for 1.6
#run in terminal
# Rscript 1.6mfa.R normal
# Rscript 1.6mfa.R prox.-inflam
# Rscript 1.6mfa.R prox.-prolif.
# Rscript 1.6mfa.R TRU
