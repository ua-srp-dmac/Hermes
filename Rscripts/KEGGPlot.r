library('org.Hs.eg.db')
library(pathview)

args = commandArgs()

print(args)

Fchange2 <- readRDS(file = paste0("/Hermes/",args[3],"/Fchange2.rds"))


options(bitmapType='cairo')
pathview(gene.data  = Fchange2,
                     pathway.id = args[2],
		     kegg.dir=  file.path("/Hermes",args[3]) )
