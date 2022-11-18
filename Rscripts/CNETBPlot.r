library(enrichplot)
library(cowplot)
library('org.Hs.eg.db')
library(clusterProfiler)
library(ggplot2)

args = commandArgs()

Fchange2 <- readRDS(file = paste0("/Hermes/",args[4],"/Fchange2.rds"))
edo2 <- readRDS(file = paste0("/Hermes/",args[4],"/edo2.rds"))
edo2B <- edo2[edo2@result$p.adjust < as.numeric(args[2]), asis=T]


if(dim(edo2B)[1] >1){

edoxB <- setReadable(edo2B, 'org.Hs.eg.db', 'ENTREZID')


p3A <- cnetplot(edoxB, foldChange=Fchange2, node_label= args[3]) + ggtitle("GSEA Gene-Concept Network\n") + theme(plot.title = element_text(size=18,hjust = 0.5))
ggsave2(file=file.path("/Hermes",args[4], "plot3A.svg"), plot=p3A,width = 11,height = 8.5)

} else{

file.copy("/Hermes/error1.svg", file.path("/Hermes",args[4] , "plot3A.svg"),overwrite = TRUE)

}
