library(enrichplot)
library(cowplot)
library(ggplot2)

source("/var/www/Rfiles/Hermes/Rscripts/ChangeSVG.R")

args = commandArgs()

edo2 <- readRDS(file = paste0("/Hermes/",args[3],"/edo2.rds"))
edo2B <- edo2[edo2@result$p.adjust < as.numeric(args[2]), asis=T]

if(dim(edo2B)[1] >1){

edo2B@result$Description <- gsub("_", " ", edo2B@result$Description)
edo2B@result$Description <- gsub("GO ", "", edo2B@result$Description)
edo2B@result$Description <- gsub("HALLMARK ", "", edo2B@result$Description)

temp <- edo2B@result$Description

for(i in 1:length(edo2B@result$Description)){
  
Lres <- length(str_split(temp[i]," ")[[1]])
if( Lres >= 3){
  edo2B@result$Description[i] <- paste0(str_split(temp[i]," ")[[1]][1],'..', str_split(temp[i]," ")[[1]][Lres],"_",i)
} else{ 
  edo2B@result$Description[i]  <-  temp[i]
  
}

}

temp <- cbind(temp,edo2B@result$Description)

p6 <- ridgeplot(edo2B) + ggtitle("Log2FoldChange Distribution\n") + labs(x = "Log2FoldChange") + theme(plot.title = element_text(size=18,hjust = 0.5)) 
ggsave2(file=file.path("/Hermes",args[3], "plot6.svg"), plot=p6,width = 11,height = 8.5)
ChangeSVGBasic(file.path("/Hermes",args[3], "plot6.svg"),temp)

} else{

file.copy("/Hermes/error1.svg", file.path("/Hermes",args[3] , "plot6.svg"),overwrite = TRUE)

}


