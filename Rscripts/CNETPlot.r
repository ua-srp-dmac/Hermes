library(enrichplot)
library(cowplot)
library('org.Hs.eg.db')
library(clusterProfiler)
library(ggupset)
library(DOSE)
library(msigdbr)
library(ggplot2)

args = commandArgs()

AresOut<- read.csv(file.path("/Hermes", args[7],"AresData.csv"))
AresOut <- AresOut[AresOut$pvalueadj<.05,]
ids <- mapIds(org.Hs.eg.db, as.character(AresOut$GeneName), 'ENTREZID', 'SYMBOL')
AresOut$ids <- ids
AresOut <- AresOut[!is.na(AresOut$ids),]

if( args[6] == "down"){

ids2 <- AresOut$ids[AresOut$Log2FoldChange <= as.numeric(args[4])]
Fchange2 <- AresOut$Log2FoldChange[AresOut$Log2FoldChange <= as.numeric(args[4])]
names(Fchange2) <- ids2
Fchange2 <- sort(Fchange2, decreasing = TRUE)

} else if( args[6] == "up"){

ids2 <- AresOut$ids[AresOut$Log2FoldChange >= as.numeric(args[5])]
Fchange2 <- AresOut$Log2FoldChange[AresOut$Log2FoldChange >= as.numeric(args[5])]
names(Fchange2) <- ids2
Fchange2 <- sort(Fchange2, decreasing = TRUE)

} else{

if(as.numeric(args[4]) == 0 && as.numeric(args[5])==0){

ids2 <- AresOut$ids
Fchange2 <- AresOut$Log2FoldChange
names(Fchange2) <- ids2
Fchange2 <- sort(Fchange2, decreasing = TRUE)

} else{

ids2 <- c(AresOut$ids[AresOut$Log2FoldChange <= as.numeric(args[4])],AresOut$ids[AresOut$Log2FoldChange >= as.numeric(args[5])])
Fchange2 <- c(AresOut$Log2FoldChange[AresOut$Log2FoldChange <= as.numeric(args[4])], AresOut$Log2FoldChange[AresOut$Log2FoldChange >= as.numeric(args[5])])
names(Fchange2) <- ids2
Fchange2 <- sort(Fchange2, decreasing = TRUE)

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

if(dim(edo)[1] >1){

edoxB <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

p3A <- cnetplot(edoxB, foldChange=Fchange2, node_label= args[3]) + ggtitle("ORA Gene-Concept Network\n") + theme(plot.title = element_text(size=18,hjust = 0.5))

ggsave2(file=file.path("/Hermes",args[7] , "plot3.svg"), plot=p3A,width = 11,height = 8.5)

} else{

file.copy("/Hermes/error2.svg", file.path("/Hermes",args[7] , "plot3.svg"),overwrite = TRUE)

}


} else{

file.copy("/Hermes/error2.svg", file.path("/Hermes",args[7] , "plot3.svg"),overwrite = TRUE)

}
