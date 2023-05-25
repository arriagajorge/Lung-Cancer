setwd("C:/Users/migue/workspace/LUSC/selectedfeatures")
library(tidyverse)
library(igraph)
library(RCy3)
library(ape)
library(phytools)

options("cyRestApi"="http://localhost:1234")
cytoscapePing()

bp=read_tsv("BP.enrichment")
k=read_tsv("KEGG.enrichment")
enriched=list(BP=bp,KEGG=k)

##############FIND "CROSSLINKED" FUNCTIONS
coenriched=do.call(rbind,enriched)%>%group_by(subtype)%>%
  group_map(~table(.x[,c("Description","component")]))
#matrix with the component intersections
intersection=lapply(coenriched,function(x) crossprod(t(x)))

# #matrix with the component unions
union3 <- function(vec) {
  n <- length(vec)
  sim_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      sim_matrix[i, j] <- sum(vec[c(i, j)] > 0)
    }
  }
  return(sim_matrix)
} 
un3 <- union3(coenriched[[3]])

# c2 <- coenriched[-3]
# union2 =lapply(c2,function(z) sapply(1:nrow(z),function(x)
#   sapply(1:nrow(z),function(y) sum(colSums(z[c(x,y),])>0))))

union = lapply(coenriched, function(z) {
  if (is.matrix(z) && nrow(z) >= 2) {
    sapply(1:nrow(z), function(x) {
      sapply(1:nrow(z), function(y) {
        tryCatch(
          {
            sum(colSums(z[c(x, y),]) > 0)
          },
          error = function(e) {
            0
          }
        )
      })
    })
  } else {
    matrix()
  }
})
union[[3]] <- un3

# ####Verify if Na, infinite or non numeric values exist in you matrix
# # lapply(coenriched, function(x) any(is.na(x)))
# # lapply(coenriched, function(x) any(is.infinite(x)))
# # lapply(coenriched, function(x) any(is.numeric(x)))
# 
# #####If exist, then use:
# which_inf <- which(sapply(coenriched, function(x) any(is.infinite(x))))
# names(coenriched)[which_inf]
# 
# which_na <- which(sapply(coenriched, function(x) any(is.na(x))))
# names(coenriched)[which_na]
# 
# which_numeric <- which(sapply(coenriched, function(x) any(is.numeric(x))))
# names(coenriched)[which_numeric]


#Jaccard index for the components
coenriched=lapply(1:4,function(x) intersection[[x]]/union[[x]])

coenriched[[3]][is.infinite(coenriched[[3]])] <- 1e+10

#don't forget 1-x to have identical sets together
trees=lapply(coenriched,function(x) hclust(as.dist(1-x)))

#get the groups that are enriched exactly in the same components
groups=lapply(trees,function(x) cutree(x,h=0))
groups=do.call(rbind,lapply(1:4,function(y) 
  data.frame(cbind("subtype"=unique(enriched$BP$subtype)[y],
                   "Description"=names(groups[[y]]),
                   "group"=groups[[y]]))))
write_tsv(groups,"Groups_per_component.tsv")

#copied from http://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Phylogenetic-trees.nb.html
names(trees) <- c("basal", "classical", "primitive", "secretory")
temp=ape::as.phylo(trees$secretory)
ig=ape::as.igraph.phylo(temp,FALSE)
ig=set_edge_attr(ig,"distance",value=temp$edge.length)
createNetworkFromIgraph(ig)
createColumnFilter('junctions', 'id', "^Node\\\\d+$", "REGEX")
junctions<-getSelectedNodes()
setNodeWidthBypass(junctions,1)
setNodeHeightBypass(junctions,1)
setNodeLabelBypass(junctions, "")

