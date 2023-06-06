#!/usr/bin/env Rscript
setwd("/home/mdiaz/workspace/LUSC/selectedfeatures")
########################PARAMETERS & PACKAGES
net="probe_primitive.filtered.alt"
library(igraph)
library(tidyverse)
library(biomaRt)
library(RCy3)

edges=read_tsv(net)
#known=read_tsv(gsub("tsv","mtrx.moti",net),col_names=F)

#pimp graph
g=graph.data.frame(edges[,1:2])
#node type maps if CpG, gene or miR
V(g)$type=substr(V(g)$name,1,1)
#node size maps degree
V(g)$size=degree(g)
#edge size maps corr
E(g)$width=abs(as.numeric(edges$MI))
#edge color maps if the interaction is known
#known=graph.data.frame(known[,1:2])2fdf38d7a222e4758b80fd147dd97fa4ac5cdeeb1a64a47d
#E(g)$color=edges$X3==1
#node color maps lfc
subtype="primitive"#unlist(strsplit(net,".",fixed=T))[2]
de=read_tsv("DE.genes.tsv")
dmir=read_tsv("DE.miR.tsv")
dm=read_tsv("DMcpgs-RUV.tsv")
da=list(cpgs=dm,genes=de,mir=dmir)
colnames(da$genes)[2]="id"
da2=lapply(da,function(x) 
  x[grep(subtype,x$contrast),]%>%filter(id%in%V(g)$name))
da3=data.frame(unique(mapply(rbind,do.call(rbind,lapply(da2[2:3],
                                                        function(x) x[,c("id","logFC","adj.P.Val")])),
                             da$cpgs[,c(2,7,4)])))#cpgs have different colnames
da4=da3[order(match(da3$id,V(g)$name)),]
V(g)$color=as.numeric(da4[1:7,]$logFC) ##Length of SECRETORY value must be 1 or 8, the number of target vertices, not 5 like in BASAL probe, WHY?
                                       ##Length of CLASSICAL value must be 1 or 5
                                       ##Length of BASAL value must be 1 or 5
                                       ##Length of PRIMITIVE value must be 1 or 7


#get readable names
methy=read_tsv("MapMethy.tsv")
methy=methy%>%filter(IlmnID%in%V(g)$name)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
                version=105)
myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
                              'entrezgene_id'),mart=mart)

#CpG names
temp=myannot%>%distinct(hgnc_symbol,entrezgene_id)%>%
  filter(!is.na(entrezgene_id))%>%
  merge(methy%>%dplyr::select(entrezgene_id,IlmnID))%>%
  distinct(hgnc_symbol,IlmnID)
nodes=data.frame("IlmnID"=V(g)$name)
nodes=merge(nodes,temp,by="IlmnID",all.x=T)
#gene names
temp=myannot%>%filter(ensembl_gene_id%in%V(g)$name)%>%
  dplyr::select(-entrezgene_id)
colnames(nodes)[1]="ensembl_gene_id"
nodes=merge(nodes,temp,by="ensembl_gene_id",all.x=T)
#merge hgnc_symbols for the 2 type of features
colnames(nodes)[1]="feature"
nodes=nodes%>%unite("name",c(hgnc_symbol.x,hgnc_symbol.y),na.rm=T)%>%
  #paste all the names linked to the same feature
  group_by(feature)%>%summarise(name=paste(name,collapse=','))
nodes$name[nodes$name==""]=nodes$feature[nodes$name==""]
nodes$name=gsub("hsa-mir","miR",nodes$name)
#order all
if(vcount(g)!=sum(V(g)$name%in%nodes$feature)){
  stop("names didn't match")
}
nodes=nodes[order(match(nodes$feature,V(g)$name)),]
V(g)$label=nodes$name

#start ctyoscape manually 1st, by system it fails randomly 
#system("Cytoscape",wait=F,ignore.stdout=T,ignore.stderr=T)
#plot the network
gc <- createNetworkFromIgraph(g, subtype)
setNodeLabelMapping(table.column = "label")
#layoutNetwork("hierarchical")
setNodeFontSizeDefault(18)
setEdgeLineWidthMapping(table.column="width",
                        widths=c(1,20))#pass widths or it'll make huge edges
setNodeShapeMapping(table.column="type",
                    table.column.values=list('c','E','h'),
                    shapes=list("ELLIPSE","RECTANGLE","HEXAGON"))#automatic mapping fails, so pass this 
#sets colors manually coz it fails with object 'res' not found
#setNodeColorMapping(table.column="color",
#	mapping.type="c",
#	table.column.values=list(min(V(g)$color),0,max(V(g)$color)),
#	colors=paletteColorBrewerRdBu(3),
#	network=getNetworkSuid())
setNodeSizeMapping(table.column="size",
                   mapping.type="continuous",
                   #pass values & sizes or it'll fail    
                   table.column.values=range(V(g)$size),
                   sizes=c(30,100))
#setEdgeColorMapping(table.column="color",
#    mapping.type="discrete",
#	table.column.values = c(TRUE,FALSE),
#	colors=c('#6600CC','#737373'))
#setEdgeTargetArrowShapeDefault('DELTA')
#setEdgeTargetArrowColorDefault('#737373')
clearEdgeBends()
saveSession(filename=paste(paste(unlist(strsplit(net,'.',fixed=T))[1:2],collapse='.'),"cys",sep='.'))
exportNetwork(gc, file = "probe_primitive.filtered.cys")