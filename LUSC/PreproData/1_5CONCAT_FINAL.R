setwd("/home/jvasquez/Documents/Lung-Cancer/LUSC/")
#!/usr/bin/env Rscript
library(data.table)
subtypeLUSC=read.table("subtypeLUSC.tsv",header=T,sep='\t')
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
subtype=subtypeLUSC[order(match(subtypeLUSC$samples,colnames(data$CpGs))),]
data[2:3]=lapply(data[2:3],function(x)
  x[,order(match(colnames(x),subtype$samples))])
print(names(data))
#"CpGs"        "miRNAseq"    "transcripts"

subtypeLUSC$subtype <- as.factor(subtypeLUSC$subtype)
#data per subtype
concatenated=lapply(levels(subtypeLUSC$subtype),function(x) 
  list(CpGs=data$CpGs[,subtypeLUSC$subtype==x],
       transcripts=data$transcripts[,subtypeLUSC$subtype==x],
       miRNA=data$miRNAs[,subtypeLUSC$subtype==x]))

names(concatenated)=levels(subtypeLUSC$subtype)


##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
print(sapply(concatenated,dim))
# normal prox.-inflam prox.-prolif.    TRU
# [1,] 428180       428180        428180 428180
# [2,]      5           69            48     71
lapply(1:5,function(x) write.table(concatenated[[x]],
                                   paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))

################################ for 1.6
#run in terminal
# Rscript 1.6mfa.R classical
# Rscript 1.6mfa.R secretory
# Rscript 1.6mfa.R basal
# Rscript 1.6mfa.R primitive
# Rscript 1.6mfa.R normal