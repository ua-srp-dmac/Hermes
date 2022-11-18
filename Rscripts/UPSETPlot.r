library(enrichplot)
library(cowplot)
library('org.Hs.eg.db')
library(clusterProfiler)
library(ggupset)
library(DOSE)
library(msigdbr)
library(ggplot2)
source("/var/www/Rfiles/Hermes/Rscripts/ChangeSVG.R")

args = commandArgs()

AresOut<- read.csv(file.path("/Hermes", args[7],"AresData.csv"))
AresOut <- AresOut[AresOut$pvalueadj<.05,]

ids <- mapIds(org.Hs.eg.db, as.character(AresOut$GeneName), 'ENTREZID', 'SYMBOL')
AresOut$ids <- ids
AresOut <- AresOut[!is.na(AresOut$ids),]

if( args[6] == "down"){

ids2 <- AresOut$ids[AresOut$Log2FoldChange <= as.numeric(args[4])]

} else if( args[6] == "up"){

ids2 <- AresOut$ids[AresOut$Log2FoldChange >= as.numeric(args[5])]

} else{

if(as.numeric(args[4]) == 0 && as.numeric(args[5])==0){

ids2 <- AresOut$ids

} else{

ids2 <- c(AresOut$ids[AresOut$Log2FoldChange <= as.numeric(args[4])],AresOut$ids[AresOut$Log2FoldChange >= as.numeric(args[5])])
}
}


if(args[8] == "DO"){

edo <- enrichDO(ids2, pvalueCutoff  = as.numeric(args[2]),  universe = AresOut$ids, minGSSize=5)

} else if(args[8] == "NCG"){

edo <- enrichNCG(ids2, pvalueCutoff  = as.numeric(args[2]),  universe = AresOut$ids, minGSSize=5)

} else if(args[8] == "DGN"){

edo <- enrichDGN(ids2, pvalueCutoff  = as.numeric(args[2]),  universe = AresOut$ids, minGSSize=5)

} else if(args[8] == "GO"){

edo <- enrichGO(ids2, OrgDb = org.Hs.eg.db, ont = args[9], pvalueCutoff  = as.numeric(args[2]),  universe = AresOut$ids, minGSSize=5)

} else if(args[8] == "Cell"){

cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))


edo <- enricher(ids2, TERM2GENE=cell_markers, minGSSize=5, pvalueCutoff  = as.numeric(args[2]))

} else if(args[8] == "Msi"){

m_t2g <- msigdbr(species = "Homo sapiens", category = args[10]) %>% 
  dplyr::select(gs_name, entrez_gene)


edo <- enricher(ids2, TERM2GENE=m_t2g, pvalueCutoff  = as.numeric(args[2]), minGSSize=5)

} else if(args[8] == "Wiki"){


edo <- enrichDO(ids2, pvalueCutoff  = as.numeric(args[2]),  universe = AresOut$ids)


} else if(args[8] == "Kegg"){

edo <- enrichKEGG(gene=ids2, organism = "hsa", pvalueCutoff  = as.numeric(args[2]),  universe = AresOut$ids, minGSSize=5)

}

if(!is.null(dim(edo)[1])){

if(dim(edo)[1] >0){

edo@result$Description <- gsub("_", " ", edo@result$Description)
edo@result$Description <- gsub("GO ", "", edo@result$Description)
edo@result$Description <- gsub("HALLMARK ", "", edo@result$Description)

temp <- edo@result$Description

for(i in 1:length(edo@result$Description)){
  
Lres <- length(str_split(temp[i]," ")[[1]])
if( Lres >= 3){
  edo@result$Description[i] <- paste0(str_split(temp[i]," ")[[1]][1],'..', str_split(temp[i]," ")[[1]][Lres],"_",i)
} else{ 
  edo@result$Description[i]  <-  temp[i]
  
}

}

temp <- cbind(temp,edo@result$Description)

p5 <- upsetplot(edo) + ggtitle("Upset Overlap Plot\n") + theme(plot.title = element_text(size=18,hjust = 0.5))
ggsave2(file=file.path("/Hermes",args[7], "plot5.svg"), plot=p5,width = 11,height = 8.5)

ChangeSVGUPSET(file.path("/Hermes", args[7] , "plot5.svg"),temp)

} else{

file.copy("/Hermes/error2.svg", file.path("/Hermes",args[7] , "plot5.svg"),overwrite = TRUE)

}

} else{

file.copy("/Hermes/error2.svg", file.path("/Hermes",args[7] , "plot5.svg"),overwrite = TRUE)

}
