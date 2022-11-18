library(enrichplot)
library(cowplot)
library(clusterProfiler)
library(ggplot2)

args = commandArgs()

edo2 <- readRDS(file = paste0("/Hermes/",args[4],"/edo2.rds"))

if(!is.null(dim(edo2)[1])){

if(dim(edo2)[1] >1){


p7 <- eval(parse(text = paste0("gseaplot2(edo2, pvalue_table = TRUE, geneSetID = c(",args[2],"), subplots = c(1,2,3))")))+ ggtitle("GSEA Plot\n") + theme(plot.title = element_text(size=18,hjust = 0.5))

ggsave2(file=file.path("/Hermes",args[4], "plot7.svg"), plot=p7,width = 11,height = 8.5)


} else{

file.copy("/Hermes/error2.svg", file.path("/Hermes",args[4] , "plot1.svg"),overwrite = TRUE)

}



} else{

file.copy("/Hermes/error2.svg", file.path("/Hermes",args[4] , "plot1.svg"),overwrite = TRUE)

}

