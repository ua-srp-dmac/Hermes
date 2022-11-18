source("/var/www/Rfiles/Hermes/Rscripts/ChangeSVG.R")
library("RPostgres")

homedata <- "/var/www/data"
analysisdir =NULL
fpath = NULL
metafile = NULL
filename = NULL
setContentType("text/html")

cat('
<!DOCTYPE html>
<html>
<head>
<title>University of Arizona SRP-DMAC Olympus server</title>
<meta http-equiv="expires"; content="0"; charset=UTF-8">
<link rel="stylesheet" type="text/css" href="/AresC/semantic.min.css">

<style>
/* Chrome, Safari, Edge, Opera */
input::-webkit-outer-spin-button,
input::-webkit-inner-spin-button {
  -webkit-appearance: none;
  margin: 0;
}

/* Firefox */
input[type=number] {
  -moz-appearance: textfield;
}

#analysisoptiontext:hover{
font-size: 16px;
}

</style>
</head>


<body>

<div id="root">
<div class="App">
	<header class="App-header">
	<div class="ui inverted top fixed menu" style:"width:90%">
		<div class="ui container">
		<a class="header item" href="https://dmac.pharmacy.arizona.edu/Hermes/Home">
		<h1 >Hermes</h1>
		</a>
		<a class="item" href="https://dmac.pharmacy.arizona.edu">
		<h4 > Home</h4>
		</a>
		<a class="item" href="https://dmac.pharmacy.arizona.edu/Ares/redirect_uri?logout=https://dmac.pharmacy.arizona.edu">
		<h4>Logout</h4>
		</a>		
			<div class="ui float right dropdown link item">
   		 		<span class="text">Menu</span>
    					<div class="menu">    
      						<a class="item" href="https://dmac.pharmacy.arizona.edu"> Home </a>
      						<a class="item" href="https://dmac.pharmacy.arizona.edu/Hermes/Home"> Hermes </a>
						<a class="item" href="https://dmac.pharmacy.arizona.edu/Hermes/About/HermesAbout.html">About Hermes</a>    							
					</div>
			</div>
		</div>

	</div>
	</header>
<div class="ui container appBody" style = "padding-top: 140px;">


')

headconv <- read.csv("/var/www/Rfiles/headconversion.csv")

usremail <- SERVER$headers_in$OIDC_CLAIM_email

username <- headconv$Folder[headconv$Email==usremail]

Aperm <- headconv$Hermes[headconv$Email==usremail]

flaguser = 0

if(length(username) != 0){
if(Aperm == 1){
datadirectories <- list.dirs(path = homedata, full.names = FALSE, recursive = FALSE)

for(datadirs in datadirectories){

if(username == datadirs){

flaguser = 1
fpath = paste0("/var/www/data/", username, "/Hermes/Data")
}

}
}
}


filenames2 = list.dirs(path = "/Hermes", full.names = TRUE, recursive = FALSE)

if(length(filenames2) > 0){

for (i in filenames2){

if(difftime(Sys.time(), file.info(i)$ctime, units="hours") > 24){


unlink(i, recursive = TRUE)

}

}

}

flagmaintop = 0

if(flaguser == 1){

if(!is.null(GET[9])){

if(!is.na(names(GET[9]))){

if( names(GET[9])=="SubmitAnalysis"){

 flagmaintop = 1
}
}
}

if(flagmaintop == 0){

cat('	
	<div class="ui container" id="maintop">
		<p class="p-t-15" style = "font-size: large;">
		Select one file for analysis.</p>
		<div class="ui action input" > <form method="GET">
  		<input type="text" name="gid" id="filter" placeholder="Filter..." style = "height: 34px; border-radius: 5px;" maxlength="20" oninput ="filterchange()">
  		<button class="ui button">Filter</button>
		<button  name="reset" class="ui button" >Reset</button>
		</div>')

filenames = list.files(path = fpath)

if(!is.null(GET$gid)){
filenames = list.files(path = fpath ,pattern=GET$gid)
seh = GET$gid
}

if(!is.null(GET$reset)){
filenames = list.files(path = fpath)
GET$gid <- NULL
}
                                                   
cat('
	
			<div class="ui segment">	
			<div class="table-container" style="padding-top: 20px; height:350px;overflow-y: scroll;">
			<table class="ui very basic table">
				<thead class="">
				<tr class="">
				<th class="">
				Name</th>
				<th class="">
				Size</th>
				<th class="">
				Last Modified</th>
				<th class="" style="text-align: center;">
				Select
				</th>
				<th></th>
				<th></th>
				<th></th>
				<th></th>
				</tr>
				
				</thead>

				<tbody class="">')
if (length(list.files(path = fpath ,pattern=GET$gid)) != 0){				
	for(i in 1:length(filenames)){
	cat('
				<tr class="" >
				<td class="">')
	cat(filenames[i])
	cat('</td>
				<td class="">') 
	cat(paste0(round(file.size(file.path(fpath, filenames[i]))/1000))," KB")
	cat('</td>
				<td class="">') 
	cat(as.character(file.mtime(file.path(fpath, filenames[i]))))
	cat('</td>
	 <td class="" style="text-align: center;"> <input type="radio" name="frequency" checked="checked"') 
  cat(' value="')
  cat(filenames[i])
  cat('"></td>	
	<td></td>
	<td></td>
	<td></td>
	<td></td>			
				</tr>')
	}
}

				cat('</tbody>
				</table>
		</div>
</div>

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Select Database:</p>

<div class="ui fluid selection dropdown" >
 
<input type="hidden" name="analysisoptions" value="DO" id="analysisoptions" onchange="ShowOptions()">
<div class="default text" ></div>
	<div class="down menu">
		<div class="header">Enrichment</div>
		<div class="item" data-value="DO">Disease Ontology (DO) </div>
		<div class="item" data-value="NCG">Network of Cancer Gene (NCG) </div>
		<div class="item" data-value="DGN">Gene-disease Associations (DGN)</div>
		<div class="item" data-value="GO">Gene Ontology (GO)</div>
		<div class="item" data-value="Cell">Cell Marker</div>
		<div class="item" data-value="Msi">Molecular Signatures Database (MSigDb)</div>
		<div class="divider"></div>
		<div class="header">Pathways</div>
		<div class="item" data-value="Wiki">WikiPathways</div>
		<div class="item" data-value="Kegg">KEGG</div>

	</div>
</div>

</div>


<div class="ui container" style="padding-top: 20px; display:none;" id="GOdisplay">

<p class="p-t-15" style = "font-size: large;">
		Select Sub-Ontology:</p>


<div class="ui fluid selection dropdown" >
  <input type="hidden" value="ALL" name="GOont" id="GOont">
  <div class="default text"></div>
  <div class="down menu">
    <div class="item" data-value="ALL">All</div>
    <div class="item" data-value="BP">Biological Process</div>
    <div class="item" data-value="CC">Cellular Component</div>
    <div class="item" data-value="MF">Molecular Function</div>
  </div>
</div>
</div>

<div class="ui container" style="padding-top: 20px; display:none;" id="MSigdisplay">

<p class="p-t-15" style = "font-size: large;">
		Select Gene Set:</p>


<div class="ui fluid selection dropdown" >
  <input type="hidden" value="H" name="GOont" id="GOont">
  <div class="default text"></div>
  <div class="down menu">
    <div class="item" data-value="H">Hallmark Gene Sets</div>
    <div class="item" data-value="C1">Positional Gene Sets</div>
    <div class="item" data-value="C2">Curated Gene Sets</div>
    <div class="item" data-value="C3">Motif Gene Sets</div>
    <div class="item" data-value="C4">Computational Gene Sets</div>
    <div class="item" data-value="C5">GO Gene Sets</div>
    <div class="item" data-value="C6">Oncogenic Signatures</div>
    <div class="item" data-value="C7">Immunologic Signatures</div>

  </div>
</div>
</div>

<input type="hidden" value="both" name="FoldFilterValue" id="FoldFilterValue">

<div class="ui container" style="padding-top: 20px;">

<button type="button" class="ui button" id="downreg" onclick="foldchange(this.id)">
  Down-Regulated
</button>
<button type="button" class="ui button" id="both" onclick="foldchange(this.id)">
  Both
</button>
<button type="button" class="ui button" id="upreg" onclick="foldchange(this.id)">
  Up-regulated
</button>

</div>

<div class="ui container" style="padding-top: 20px; display: none;" id="foldtopT">
<div class="ui labeled input" >
<input type="number" placeholder="-1" style="text-align: center;" name="Fold1" id="Fold1T" value="-1">
  <div class="ui label" id="folddownT">
    &#8805; Log2FoldChange
  </div>
  <div class="ui label" id="foldupT">
     Log2FoldChange &#8805;
  </div>
  <div class="ui label" id="foldbothT">
    &#8805; Log2FoldChange &#8805;
  </div>

  <input type="number" placeholder="1" style="text-align: center;" name="Fold2" id="Fold2T" value="1">
</div>

</div>


		<div class="m-t-25" style = "padding-top: 20px;">
			<button class="ui black button" name=SubmitAnalysis onclick="maintop()">
			Submit Analysis</button>
		</div>
</form>

<br>
<br>

</div>

  <div id="Mainloader" class="ui disabled inverted dimmer" style ="height: 1050px;">
    <div class="ui text loader">Loading (this may take several minutes)</div>
  </div>

')
}


if(!is.null(GET[[2]]) & !is.null(GET[9])){

filename = GET[[2]]

if(file.exists(paste0("/var/www/data/", username, "/Hermes/Data/",filename)) & !is.na(names(GET[9]))){

if( names(GET[9])=="SubmitAnalysis")

{

if(GET[[7]] > 0|| GET[[8]] < 0){

cat('

<p class="p-t-15" style = "font-size: large; ">
		Please select Log2FoldChange cuttoff of less than or equal to 0 for down-regulated genes and greater than or equal to 0 for up-regulated genes.  </p>

')

} else {

usedb <- NULL
usedb$app <- "Hermes"
usedb$username <- SERVER$headers_in$OIDC_CLAIM_preferred_username

usedb <- as.data.frame(usedb)

con2 <- dbConnect(RPostgres::Postgres())

dbWriteTable(con2, "table1u", value = usedb, append = TRUE, row.names = FALSE)

dbDisconnect(con2)


metafile <- paste0("/var/www/data/", username, "/Hermes/Meta/",filename)
metadata <- read.csv(metafile)
metadata<-metadata[,2:3]

AresOut<- read.csv(file.path("/var/www/data/", username, "/Hermes/Data/",filename))

filename <- paste0(filename, round(rnorm(1,10000000, 1000000)),round(rnorm(1,10000000, 1000000)))

dir.create(file.path("/Hermes", filename))

write.csv(AresOut, paste0("/Hermes/", filename ,"/AresData.csv"))

AresOut <- AresOut[AresOut$pvalueadj<.05,]
ids <- mapIds(org.Hs.eg.db, as.character(AresOut$GeneName), 'ENTREZID', 'SYMBOL')
AresOut$ids <- ids
AresOut <- AresOut[!is.na(AresOut$ids),]

if(GET[[6]] == "both"){

if(GET[[7]] == 0 && GET[[8]]==0){

ids2 <- AresOut$ids
Fchange <- AresOut$Log2FoldChange
names(Fchange) <- ids2
Fchange <- sort(Fchange, decreasing = TRUE)

} else{

ids2 <- c(AresOut$ids[AresOut$Log2FoldChange <= as.numeric(GET[[7]])],AresOut$ids[AresOut$Log2FoldChange >= as.numeric(GET[[8]])])
Fchange <- c(AresOut$Log2FoldChange[AresOut$Log2FoldChange <= as.numeric(GET[[7]])], AresOut$Log2FoldChange[AresOut$Log2FoldChange >= as.numeric(GET[[8]])])
names(Fchange) <- ids2
Fchange <- sort(Fchange, decreasing = TRUE)
}

} else if(GET[[6]] == "downreg"){

ids2 <- AresOut$ids[AresOut$Log2FoldChange <= as.numeric(GET[[7]])]
Fchange <- AresOut$Log2FoldChange[AresOut$Log2FoldChange <= as.numeric(GET[[7]])]
names(Fchange) <- ids2
Fchange <- sort(Fchange, decreasing = TRUE)


} else{

ids2 <- AresOut$ids[AresOut$Log2FoldChange >= as.numeric(GET[[8]])]
Fchange <- AresOut$Log2FoldChange[AresOut$Log2FoldChange >= as.numeric(GET[[8]])]
names(Fchange) <- ids2
Fchange <- sort(Fchange, decreasing = TRUE)

}

Fchange2 <- Fchange


if(GET[[3]] == "DO"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Disease Ontology (DO) Analysis</p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

edo <- enrichDO(ids2, pvalueCutoff  = .05,  universe = AresOut$ids, minGSSize=5)
edo2 <- gseDO(Fchange2, pvalueCutoff = 1, minGSSize=5, seed = TRUE)
}

}
else if(GET[[3]] == "NCG"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Network of Cancer Gene (NCG) Analysis</p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

edo <- enrichNCG(ids2, pvalueCutoff  = .05,  universe = AresOut$ids, minGSSize=5)
edo2 <- gseNCG(Fchange2, pvalueCutoff = 1, minGSSize=5, seed = TRUE)
}

}

else if(GET[[3]] == "DGN"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Gene-disease Associations (DGN) Analysis</p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

edo <- enrichDGN(ids2, pvalueCutoff  = .05,  universe = AresOut$ids, minGSSize=5)
edo2 <- gseDGN(Fchange2, pvalueCutoff = 1, minGSSize=5, seed = TRUE)
}

}

else if(GET[[3]] == "GO"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Gene Ontology (GO) Analysis</p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

edo <- enrichGO(ids2, OrgDb = org.Hs.eg.db, ont = GET[[4]], pvalueCutoff  = .05,  universe = AresOut$ids, minGSSize=5)
edo2 <- gseGO(Fchange2, OrgDb = org.Hs.eg.db, ont = GET[[4]], pvalueCutoff = 1, minGSSize=5, seed = TRUE)
}

}

else if(GET[[3]] == "Cell"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Cell Marker Analysis</p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))

edo <- enricher(ids2, TERM2GENE=cell_markers, minGSSize=5, pvalueCutoff  = .05)
edo2 <- GSEA(Fchange2, TERM2GENE = cell_markers, minGSSize=5, pvalueCutoff = 1, seed = TRUE)
}
}

else if(GET[[3]] == "Msi"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Molecular Signatures Database (MSigDb) </p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

m_t2g <- msigdbr(species = "Homo sapiens", category = GET[[5]]) %>% 
  dplyr::select(gs_name, entrez_gene)

edo <- enricher(ids2, TERM2GENE=m_t2g, pvalueCutoff  = .05, minGSSize=5)
edo2 <- GSEA(Fchange2, TERM2GENE = m_t2g, pvalueCutoff = 1, minGSSize=5, seed = TRUE)
}

}

else if(GET[[3]] == "Wiki"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		WikiPathways Analysis Currently Unavaiable - Running Disease Ontology Instead</p>
</div>')

if(length(Fchange) == 0){

edo <- NA
edo2 <- NA

} else{

edo <- enrichDO(ids2, pvalueCutoff  = .05,  universe = AresOut$ids)
edo2 <- gseDO(Fchange2, pvalueCutoff = 1, seed = TRUE)
}

}

else if(GET[[3]] == "Kegg"){

cat('<div class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		KEGG Analysis</p>
</div>')

if(length(Fchange) == 0){
edo <- NA
edo2 <- NA

} else{

edo <- enrichKEGG(gene=ids2, organism = "hsa", pvalueCutoff  = .05,  universe = AresOut$ids, minGSSize=5)
edo2 <- gseKEGG(geneList = Fchange2, organism = "hsa", pvalueCutoff = 1, minGSSize=5, seed = TRUE)

if(dim(edo2)[1] > 0){

options(bitmapType='cairo')

pathview(gene.data  = Fchange2,
                     pathway.id = edo2@result$ID[1],
                     species    = "hsa",
		     kegg.dir=  file.path("/Hermes",filename) )

} 
}

}

saveRDS(Fchange2, file = paste0("/Hermes/",filename,"/Fchange2.rds"))
saveRDS(edo, file = paste0("/Hermes/",filename,"/edo.rds"))
saveRDS(edo2, file = paste0("/Hermes/",filename,"/edo2.rds"))

if(is.null(edo)){
edo <- NA
}

edo2F <-edo2

if(!is.na(edo2)){
edo2 <- edo2[edo2@result$p.adjust < .05, asis=T]
}


if(!is.na(edo2F)){

if(dim(edo2F)[1] > 0){

p7 <- gseaplot2(edo2F, geneSetID = 1,pvalue_table = TRUE,subplots = 1:3) + ggtitle("GSEA Plot\n") + theme(plot.title = element_text(size=18,hjust = 0.5))

ggsave(file=file.path("/Hermes",filename, "plot7.svg"), plot=p7,width = 11,height = 8.5)

} else{
file.copy("/Hermes/error3.svg", file.path("/Hermes",filename , "plot7.svg"),overwrite = TRUE)
}

} else{
file.copy("/Hermes/error3.svg", file.path("/Hermes",filename , "plot7.svg"),overwrite = TRUE)
}


if(!is.na(edo2)){

if(dim(edo2)[1] >1){

edo2@result$Description <- gsub("_", " ", edo2@result$Description)
edo2@result$Description <- gsub("GO ", "", edo2@result$Description)
edo2@result$Description <- gsub("HALLMARK ", "", edo2@result$Description)

temp <- edo2@result$Description

for(i in 1:length(edo2@result$Description)){
  
Lres <- length(str_split(temp[i]," ")[[1]])
if( Lres >= 3){
  edo2@result$Description[i] <- paste0(str_split(temp[i]," ")[[1]][1],'..', str_split(temp[i]," ")[[1]][Lres],"_",i)
} else{ 
  edo2@result$Description[i]  <-  temp[i]
  
}

}

temp <- cbind(temp,edo2@result$Description)

p6 <- ridgeplot(edo2) + ggtitle("Log2FoldChange Distribution\n") + labs(x = "Log2FoldChange") + theme(plot.title = element_text(size=18,hjust = 0.5)) 
ggsave(file=file.path("/Hermes",filename, "plot6.svg"), plot=p6,width = 11,height = 8.5)
ChangeSVGBasic(file.path("/Hermes",filename, "plot6.svg"),temp)

p2 <- dotplot(edo2, showCategory=20) + ggtitle("GSEA Dotplot\n") + theme(plot.title = element_text(size=18,hjust = 0.5)) 
ggsave(file=file.path("/Hermes",filename, "plot2.svg"), plot=p2,width = 11,height = 8.5)
ChangeSVGBasic(file.path("/Hermes",filename, "plot2.svg"),temp)

edoxB <- setReadable(edo2, 'org.Hs.eg.db', 'ENTREZID')
p3A <- cnetplot(edoxB, categorySize="pvalue", foldChange=Fchange2) + ggtitle("GSEA Gene-Concept Network\n") + theme(plot.title = element_text(size=18,hjust = 0.5))
ggsave(file=file.path("/Hermes",filename, "plot3A.svg"), plot=p3A,width = 11,height = 8.5)

} else{
file.copy("/Hermes/error3.svg", file.path("/Hermes",filename , "plot2.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot6.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot3A.svg"),overwrite = TRUE)
}

}else{
file.copy("/Hermes/error3.svg", file.path("/Hermes",filename , "plot2.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot6.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot3A.svg"),overwrite = TRUE)
}


if(!is.na(edo)){
if(dim(edo)[1] >0){

edo@result$Description <- gsub("_", " ", edo@result$Description)
edo@result$Description <- gsub("GO ", "", edo@result$Description)
edo@result$Description <- gsub("HALLMARK ", "", edo@result$Description)

tempB <- edo@result$Description

for(i in 1:length(edo@result$Description)){
  
Lres <- length(str_split(tempB[i]," ")[[1]])
if( Lres >= 3){
  edo@result$Description[i] <- paste0(str_split(tempB[i]," ")[[1]][1],'..', str_split(tempB[i]," ")[[1]][Lres],"_",i)
} else{ 
  edo@result$Description[i]  <-  tempB[i]
  
}
}

tempB <- cbind(tempB,edo@result$Description)

p1 <- dotplot(edo, showCategory=20) + ggtitle("ORA Dotplot\n") + theme(plot.title = element_text(size=18,hjust = 0.5)) 
ggsave(file=file.path("/Hermes",filename, "plot1.svg"), plot=p1,width = 11,height = 8.5)
ChangeSVGBasic(file.path("/Hermes",filename, "plot1.svg"),tempB)

p5 <- upsetplot(edo) + ggtitle("Upset Overlap Plot\n") + theme(plot.title = element_text(size=18,hjust = 0.5))
ggsave(file=file.path("/Hermes",filename, "plot5.svg"), plot=p5,width = 11,height = 8.5)
ChangeSVGUPSET(file.path("/Hermes",filename, "plot5.svg"),tempB)

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p3 <- cnetplot(edox, categorySize="pvalue", foldChange=Fchange) + ggtitle("ORA Gene-Concept Network\n") + theme(plot.title = element_text(size=18,hjust = 0.5))

ggsave(file=file.path("/Hermes",filename, "plot3.svg"), plot=p3,width = 11,height = 8.5)

p4 <- emapplot(edo, pie_scale=1.5) + ggtitle("Enrichment Map\n") + theme(plot.title = element_text(size=18,hjust = 0.5)) 
ggsave(file=file.path("/Hermes",filename, "plot4.svg"), plot=p4,width = 11,height = 8.5) 

} else{
file.copy("/Hermes/error3.svg", file.path("/Hermes",filename , "plot1.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot3.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot4.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot5.svg"),overwrite = TRUE)
}

} else{
file.copy("/Hermes/error3.svg", file.path("/Hermes",filename , "plot1.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot3.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot4.svg"),overwrite = TRUE)
file.copy("/Hermes/error3.svg", file.path("/Hermes", filename, "plot5.svg"),overwrite = TRUE)
}



cat('
<div id="context1" style = "padding-top: 20px;">
  <div class="ui secondary menu">
    <a class="item active" data-tab="firsttop">Gene Set Enrichment Analysis</a>
    <a class="item" data-tab="secondtop">Over Representation Enrichment Analysis</a>')

if(GET[[3]] == "Kegg"){
 cat('<a class="item" data-tab="thirdtop">Pathway Analysis</a>')
}


cat('
  </div>
<div class="ui tab segment active" data-tab="firsttop">
<div class="ui top attached tabular menu" >
  <a class="item active" data-tab="first">Files Selected</a>
  <a class="item" data-tab="second"><div class="ui raised segment"><img id = "gseaimagepre" src="') 
	cat(paste0(file.path("/HermesResults",filename, "plot7.svg")))
	cat('" alt="GSEA Plot" style="width:150px;height:150px;"> </img></div></a>
  <a class="item " data-tab="thirdA"><div class="ui raised segment"><img id = "cnetBimagepre" src="')
	cat(paste0(file.path("/HermesResults",filename, "plot3A.svg")))
	cat('" alt="Cnet Plot" style="width:150px;height:150px;"> </img></div></a>
  <a class="item " data-tab="sixth"><div class="ui raised segment"><img id = "ridgeimagepre" src="')
	cat(paste0(file.path("/HermesResults",filename, "plot6.svg")))
	cat('" alt="Ridge Plot" style="width:150px;height:150px;"> </img></div></a>
  <a class="item " data-tab="eighth"><div class="ui raised segment"><img id = "gseadotimagepre" src="')
	cat(paste0(file.path("/HermesResults",filename, "plot2.svg")))
	cat('" alt="Dot Plot GSEA" style="width:150px;height:150px;"></div></a>
</div>



<div class="ui bottom attached tab segment active" data-tab="first">
  
<div class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Selected File: ')
cat(filename)
cat('</p>
	
	<div class="ui segment">	
			<div class="table-container" style="padding-top: 20px; max-height:340px;overflow-y: scroll;">
			<table class="ui very basic table">
				<thead class="">
				<tr class="">
				<th style="text-align: center;">
				Group 1</th>
				<th style="text-align: center;">
				Group 2</th>
				</tr>	
				</thead>
				<tbody>
				')
	temp1 = NULL
	temp2 = NULL


	for(i in 1:length(metadata[,1])){

	if(metadata[i,1] == 'Treatment')
	{
	temp1=rbind(temp1,metadata[i,2])
	}

	if(metadata[i,1] == 'Control')
	{
	temp2=rbind(temp2,metadata[i,2])
	}

	
	}
	if(length(temp1) >= length(temp2)){
		for(i in 1:length(temp1)){

			cat('<tr>
			<td style="text-align: center;">')
			cat(temp1[i])
			cat('</td>
			<td style="text-align: center;">')
		if(!is.null(temp2[i]) &  !is.na(temp2[i])){
			cat(temp2[i])
			cat('</td>
			</tr>')	
			}
		else{
			cat('</td>
			</tr>')	
		}
		}

	}

	if(length(temp1) < length(temp2)){

		for(i in 1:length(temp2)){
			cat('<tr>
			<td style="text-align: center;">')

		if(!is.null(temp1[i])&  !is.na(temp1[i])){
			cat(temp1[i])
			cat('</td>')	
			}
		else{
			cat('</td>')	
		}

			cat('<td style="text-align: center;">')
			cat(temp2[i])
			cat('</td>
			</tr>')
		}

	}

		
	
			cat('	</tbody>
				</table>
		 	</div>
	</div>
	</div>

</div>

<div class="ui bottom attached tab segment" data-tab="second">


<div class="ui container" style = "padding-bottom: 20px;">

<div class="ui segment">
	<h4 class="ui header">Gene Set Enrichment Visualization</h4>
	<p class="p-t-15" style = "font-size: large;">
	Classic visualization of gene set enrichment analaysis (GSEA). Gene sets in the dropdown are ordered by adjusted p-value. Choose multiple gene sets to overlay them. 
	</p>
</div>
</div>


')

if(!is.na(edo2F)){

if(dim(edo2F)[1] > 0){

cat('


<div class="ui fluid multiple search selection dropdown" >

<input type="hidden" name="gseaset" id="gseaset">  
<div class="default text">Gene Sets</div>
<div class="menu" >
')

flag = 0

for(i in 1:length(edo2F@result$Description)){

if(edo2F@result$p.adjust[i] >= .05 && flag == 0){

cat('
<div class="divider"></div>
<div class="header" style="text-align: center;">Adjusted P-value Over 0.05</div>
')

flag = 1
}

cat('<div class="item" data-value="')
cat(i)
cat('">')
cat(edo2F@result$Description[i])
cat('</div>
')

}
cat('
</div>
</div>


		<div class="m-t-25" style = "padding-top: 20px;">
			<button class="ui black button" name=GSEAplot onclick="gsea()">
			Submit Analysis</button>
		</div>
<div class="ui divider"></div>

    <img  id = "gseaimage" style="height: 100%; width: 100%; object-fit: contain; padding-top: 20px;" src="')

cat(paste0(file.path("/HermesResults",filename, "plot7.svg")))
cat('">
<a href="')
cat(paste0(file.path("/HermesResults",filename, "plot7.svg")))
cat('" download="GSEAPlot"> 
<br>
<button class="ui button">Download Image </button>
</a>')
}else{

cat('

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits.</p>
</div>
')

}

cat('


  <div id="gsealoader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>

</div>')

cat('

<div class="ui bottom attached tab segment" data-tab="thirdA">

<div class="ui container" >

<div class="ui segment">
	<h4 class="ui header">Gene-Concept Network</h4>
	<p class="p-t-15" style = "font-size: large;">
	Network visualization of the links between genes and biological concepts.  
	</p>
</div>
</div>

')

if(!is.na(edo2)){

cat('

<div style = "padding-top: 20px;">

<div class="ui fluid selection dropdown"> 
	<input type="hidden" name="cnetBlabel" id="cnetBlabel">  
	<div class="default text">Label Options</div>
	<div class="menu" >
		<div class="item" data-value="all">All</div>
		<div class="item" data-value="category">Category</div>
		<div class="item" data-value="gene">Gene</div>
		<div class="item" data-value="none">None</div>
	</div>
</div>
</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="pcnetB" value="0.05" >
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnetB onclick="CnetB()">
	Submit Analysis</button>
</div>

 <div id="cnetBerror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Make sure to select all options and choose a p-value between 0 and 1.<p>
  </div> ')

if(dim(edo2)[1] >1){


cat('

<div id="warningCNETB" class="ui container" style="display: none;">
		<p class="p-t-15" style = "font-size: large;">
		Could not find enough enriched genes for visualization. You need at least two.</p>
</div>


<div class="ui divider"></div>


    <img id="cnetBimage" style="height: 100%; width: 100%; object-fit: contain" src="')

cat(paste0(file.path("/HermesResults",filename, "plot3A.svg")))
cat('">
<a id="downCNETB" href="')
cat(paste0(file.path("/HermesResults",filename, "plot3A.svg")))
cat('" download="CnetPlot"> 
<button class="ui button">Download Image </button>
</a>')

}else{

cat('

 <img  id = "cnetBimage" style="height: 100%; width: 100%; object-fit: contain; padding-top: 20px; display: none;" src="')

cat(paste0(file.path("/HermesResults",filename, "plot3A.svg")))
cat('">
<a id="downCNETB" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot3A.svg")))
cat('" download="CnetPlot"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningCNETB" class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Could not find enough enriched genes for visualization. You need at least two.</p>
</div>
')

}


}
else{

cat('

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large; padding-top: 20px;" >
		No genes found within log2foldchange limits.</p>
</div>
')

}


cat('

 <div id="cnetBloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>




</div>

<div class="ui bottom attached tab segment" data-tab="sixth">

<div class="ui container">

<div class="ui segment">
	<h4 class="ui header">Ridgeline Plot</h4>
	<p class="p-t-15" style = "font-size: large;">
	Visualization of fold change distriubtion for core enriched genes in the gene set enrichment analysis (GSEA). Useful for identifying up/down-regulated pathways.   
	</p>
</div>
</div> ')

if(!is.na(edo2)){

cat('
<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="pridge" value="0.05">
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnet onclick="ridge()">
	Submit Analysis</button>
</div>

 <div id="ridgeerror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Choose a p-value between 0 and 1.<p>
  </div>

')

if(dim(edo2)[1] >0){


cat('

<div id="warningRIDGE" class="ui container" style="display: none;">
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>


<div class="ui divider"></div>

<div id = "ridgeimage">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot6.svg")))
cat('"></object>
</div>
<a id="downRIDGE" href="')
cat(paste0(file.path("/HermesResults",filename, "plot6.svg")))
cat('" download="RidgePlot"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('

<div style= "display: none;" id = "ridgeimage">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot6.svg")))
cat('"></object>
</div>
<a style="display: none;" id="downRIDGE" href="')
cat(paste0(file.path("/HermesResults",filename, "plot6.svg")))
cat('" download="RidgePlot"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningRIDGE" class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('



<div class="ui container" >
		<p class="p-t-15" style = "font-size: large; padding-top: 20px;">
		No genes found within log2foldchange limits.</p>
</div>
')

}

cat('

 <div id="ridgeloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>


</div>


<div class="ui bottom attached tab segment" data-tab="eighth">

<div class="ui container">

<div class="ui segment">
	<h4 class="ui header">Gene Set Enrichment Analysis Dot Plot</h4>
	<p class="p-t-15" style = "font-size: large;">
	Dot plot visualization of gene set enrichment analysis.   
	</p>
</div>
</div>
')

if(!is.na(edo2)){

cat('

<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    Display Top:
  </div>
  <input type="number" placeholder="20" style="text-align: center;" id="topgseadot" value="20">
</div>

</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="pgseadot" value="0.05">
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnet onclick=" gseadot()">
	Submit Analysis</button>
</div>

 <div id="gseadoterror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Make sure Display Top is between 0 and 100. Make sure to choose a p-value between 0 and 1.<p>
  </div>


')


if(dim(edo2)[1] >0){

cat('

<div id="warningGSEADOT" class="ui container" style="padding-top: 20px;display: none;">
		<p class="p-t-15" style = "font-size: large; ">
		Could not find any enriched gene sets.</p>
</div>


<div class="ui divider"></div>

<div id = "gseadotimage">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot2.svg")))
cat('"></object>
</div>
<a id="downGSEADOT" href="')
cat(paste0(file.path("/HermesResults",filename, "plot2.svg")))
cat('" download="DotPlotGSEA"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('


<div id = "gseadotimage" style="display : none;">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot2.svg")))
cat('"></object>
</div>
<a style="display: none;" id="downGSEADOT" href="')
cat(paste0(file.path("/HermesResults",filename, "plot2.svg")))
cat('" download="DotPlotGSEA"> 
<button class="ui button">Download Image </button>
</a>


<div id="warningGSEADOT" class="ui container" style="padding-top: 20px;">
		<p class="p-t-15" style = "font-size: large; ">
		Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large; padding-top: 20px;">
		No genes found within log2foldchange limits.</p>
</div>
')

}

cat('

 <div id="gseadotloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>
</div>

</div>

<div class="ui tab segment" data-tab="secondtop">
    <div class="ui top attached tabular menu">
  <a class="item active" data-tab="third"><div class="ui raised segment"><img id = "cnetimagepre" src="') 
	cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
	cat('" alt="CNET Plot" style="width:150px;height:150px;"> </img></div></a>
  <a class="item " data-tab="fourth"><div class="ui raised segment"><img id = "emapimagepre" src="')
	cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
	cat('" alt="GSEA Plot" style="width:150px;height:150px;"> </img></div></a>
  <a class="item " data-tab="fifth"><div class="ui raised segment"><img id = "upsetimagepre" src="')
	cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
	cat('" alt="GSEA Plot" style="width:150px;height:150px;"> </img></div></a>
  <a class="item " data-tab="seventh"><div class="ui raised segment"><img id = "oradotimagepre" src="')
	cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
	cat('" alt="GSEA Plot" style="width:150px;height:150px;"> </img></div></a>   
 </div>')

cat('

<div class="ui bottom attached tab segment active" data-tab="third">

<div class="ui container" style = "padding-bottom: 20px;">

<div class="ui segment">
	<h4 class="ui header">Gene-Concept Network</h4>
	<p class="p-t-15" style = "font-size: large;">
	Network visualization of the links between genes and biological concepts.  
	</p>
</div>
</div>

')

cat('

<div class="ui container" style="padding-top: 20px;">

<button class="ui button" id="downreg" onclick="foldchange2(this.id)">
  Down-Regulated
</button>
<button class="ui button" id="both" onclick="foldchange2(this.id)">
  Both
</button>
<button class="ui button" id="upreg" onclick="foldchange2(this.id)">
  Up-regulated
</button>

</div>

<div class="ui container" style="padding-top: 20px; display: none;" id="foldtop">
<div class="ui labeled input" >
<input type="number" placeholder="-1" style="text-align: center;" name="Fold1" id="Fold1" value="-1">
  <div class="ui label" id="folddown">
    &#8805; Log2FoldChange
  </div>
  <div class="ui label" id="foldup">
     Log2FoldChange &#8805;
  </div>
  <div class="ui label" id="foldboth">
    &#8805; Log2FoldChange &#8805;
  </div>

  <input type="number" placeholder="1" style="text-align: center;" name="Fold2" id="Fold2" value="1">
</div>

</div>

<div style = "padding-top: 20px;">

<div class="ui fluid selection dropdown"> 
	<input type="hidden" name="cnetlabel" id="cnetlabel">  
	<div class="default text">Label Options</div>
	<div class="menu" >
		<div class="item" data-value="all">All</div>
		<div class="item" data-value="category">Category</div>
		<div class="item" data-value="gene">Gene</div>
		<div class="item" data-value="none">None</div>
	</div>
</div>
</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="pcnet" value="0.05">
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnet onclick="Cnet()">
	Submit Analysis</button>
</div>

 <div id="cneterror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Make sure to select all options and choose a p-value between 0 and 1.<p>
  </div>

<div class="ui divider"></div>

')

if(!is.na(edo)){

if(dim(edo)[1] >1){

cat('

<div id="warningCNET" class="ui container" style="display: none;">
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>

    <img id="cnetimage" style="height: 100%; width: 100%; object-fit: contain" src="')

cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
cat('">
<a id="downCNET" href="')
cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
cat('" download="CnetPlot"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('

    <img id="cnetimage" style="height: 100%; width: 100%; object-fit: contain; display: none;" src="')

cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
cat('">
<a id="downCNET" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
cat('" download="CnetPlot"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningCNET" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('

    <img id="cnetimage" style="height: 100%; width: 100%; object-fit: contain; display: none;" src="')

cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
cat('">
<a id="downCNET" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot3.svg")))
cat('" download="CnetPlot"> 
<button class="ui button">Download Image </button>
</a>


<div id="warningCNET" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits or could not find any enriched gene sets.</p>
</div>
')

}

cat('

 <div id="cnetloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>


</div>

<div class="ui bottom attached tab segment" data-tab="fourth">

<div class="ui container" style = "padding-bottom: 20px;">

<div class="ui segment">
	<h4 class="ui header">Enrichment Map</h4>
	<p class="p-t-15" style = "font-size: large;">
	Visualization that organizes enriched terms into a network with edges connecting overlapping gene sets. Mutually overlapping gene sets tend to cluster together, making it easy to identify functional modules.  
	</p>
</div>
</div>


')

cat('

<div class="ui container" style="padding-top: 20px;">

<button class="ui button" id="downreg" onclick="foldchange2(this.id)">
  Down-Regulated
</button>
<button class="ui button" id="both" onclick="foldchange2(this.id)">
  Both
</button>
<button class="ui button" id="upreg" onclick="foldchange2(this.id)">
  Up-regulated
</button>

</div>

<div class="ui container" style="padding-top: 20px; display: none;" id="foldtop2">
<div class="ui labeled input" >
<input type="number" placeholder="-1" style="text-align: center;" name="Fold1" id="Fold12" value="-1">
  <div class="ui label" id="folddown2">
    &#8805; Log2FoldChange
  </div>
  <div class="ui label" id="foldup2">
     Log2FoldChange &#8805;
  </div>
  <div class="ui label" id="foldboth2">
    &#8805; Log2FoldChange &#8805;
  </div>

  <input type="number" placeholder="1" style="text-align: center;" name="Fold2" id="Fold22" value="1">
</div>

</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="pemap" value="0.05">
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnet onclick="Emap()">
	Submit Analysis</button>
</div>

 <div id="emaperror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Choose a p-value between 0 and 1<p>
  </div>

<div class="ui divider"></div>

')

if(!is.na(edo)){

if(dim(edo)[1] >1){

cat('

<div id="warningEMAP" class="ui container" style="display: none;">
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits, or could not find any enriched gene sets.</p>
</div>

    <img id ="emapimage" style="height: 100%; width: 100%; object-fit: contain" src="')

cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
cat('">
<a id="downEMAP" href="')
cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
cat('" download="EMAPPlot"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('

  <img id ="emapimage" style="height: 100%; width: 100%; object-fit: contain; display: none;" src="')

cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
cat('">
<a id="downEMAP" style="display: none;"  href="')
cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
cat('" download="EMAPPlot"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningEMAP" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('

  <img id ="emapimage" style="height: 100%; width: 100%; object-fit: contain; display: none;" src="')

cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
cat('">
<a id="downEMAP" style="display: none;"  href="')
cat(paste0(file.path("/HermesResults",filename, "plot4.svg")))
cat('" download="EMAPPlot"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningEMAP" class="ui container">
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits, or could not find any enriched gene sets.</p>
</div>
')

}

cat('

 <div id="emaploader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>

</div>

<div class="ui bottom attached tab segment" data-tab="fifth">
<div class="ui container" style = "padding-bottom: 20px;">

<div class="ui segment">
	<h4 class="ui header">Gene-Concept Network</h4>
	<p class="p-t-15" style = "font-size: large;">
	 Visualization of complex associations between genes and gene sets, emphasizing the genes overallping among different gene sets.  
	</p>
</div>
</div>
')

cat('

<div class="ui container" style="padding-top: 20px;">

<button class="ui button" id="downreg" onclick="foldchange2(this.id)">
  Down-Regulated
</button>
<button class="ui button" id="both" onclick="foldchange2(this.id)">
  Both
</button>
<button class="ui button" id="upreg" onclick="foldchange2(this.id)">
  Up-regulated
</button>

</div>

<div class="ui container" style="padding-top: 20px; display: none;" id="foldtop3">
<div class="ui labeled input" >
<input type="number" placeholder="-1" style="text-align: center;" name="Fold13" id="Fold13" value="-1">
  <div class="ui label" id="folddown3">
    &#8805; Log2FoldChange
  </div>
  <div class="ui label" id="foldup3">
     Log2FoldChange &#8805;
  </div>
  <div class="ui label" id="foldboth3">
    &#8805; Log2FoldChange &#8805;
  </div>

  <input type="number" placeholder="1" style="text-align: center;" name="Fold23" id="Fold23" value="1">
</div>

</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="pupset" value="0.05">
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnet onclick="Upset()">
	Submit Analysis</button>
</div>

 <div id="upseterror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Choose a p-value between 0 and 1.<p>
  </div> 

<div class="ui divider"></div>

')

if(!is.na(edo)){
if(dim(edo)[1] >0){

cat('

<div id="warningUPSET" class="ui container"  style="display: none;">
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits, or could not find any enriched gene sets.</p>
</div>


<div id = "upsetimage">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
cat('"></object>
</div>
<a id="downUPSET" href="')
cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
cat('" download="Upset Plot"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('

<div id = "upsetimage" style="display: none;">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
cat('"></object>
</div>
<a id="downUPSET" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
cat('" download="Upset Plot"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningUPSET" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('

<div id = "upsetimage" style="display: none;">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
cat('"></object>
</div>
<a id="downUPSET" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot5.svg")))
cat('" download="Upset Plot"> 
<button class="ui button">Download Image </button>
</a>


<div id="warningUPSET" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits, or could not find any enriched gene sets.</p>
</div>
')

}


cat('

 <div id="upsetloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>

</div>

<div class="ui bottom attached tab segment" data-tab="seventh">

<div class="ui container" style = "padding-bottom: 20px;">

<div class="ui segment">
	<h4 class="ui header">Over Representation Analysis Dot Plot</h4>
	<p class="p-t-15" style = "font-size: large;">
	Dot plot visualization of over representation analysis.   
	</p>
</div>
</div>


')

cat('

<div class="ui container" style="padding-top: 20px;">

<button class="ui button" id="downreg" onclick="foldchange2(this.id)">
  Down-Regulated
</button>
<button class="ui button" id="both" onclick="foldchange2(this.id)">
  Both
</button>
<button class="ui button" id="upreg" onclick="foldchange2(this.id)">
  Up-regulated
</button>

</div>

<div class="ui container" style="padding-top: 20px; display: none;" id="foldtop4">
<div class="ui labeled input" >
<input type="number" placeholder="-1" style="text-align: center;" name="Fold14" id="Fold14" value="-1">
  <div class="ui label" id="folddown4">
    &#8805; Log2FoldChange
  </div>
  <div class="ui label" id="foldup4">
     Log2FoldChange &#8805;
  </div>
  <div class="ui label" id="foldboth4">
    &#8805; Log2FoldChange &#8805;
  </div>

  <input type="number" placeholder="1" style="text-align: center;" name="Fold2" id="Fold24" value="1">
</div>

</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    Display Top:
  </div>
  <input type="number" placeholder="20" style="text-align: center;" id="toporadot" value="20">
</div>

</div>


<div class="ui container" style="padding-top: 20px;">
<div class="ui labeled input">
  <div class="ui label">
    P-value Cutoff
  </div>
  <input type="number" placeholder="0.05" style="text-align: center;" id="poradot" value="0.05">
</div>

</div>

<div class="m-t-25" style = "padding-top: 20px;">
	<button class="ui black button" name=cnet onclick="Oradot()">
	Submit Analysis</button>
</div>

 <div id="oradoterror" style = "font-size: large; display:none; padding-top: 20px;">
    <p>Make sure Display Top is between 0 and 100. Make sure to choose a p-value between 0 and 1.<p>
  </div>

<div class="ui divider"></div>

')

if(!is.na(edo)){

if(dim(edo)[1] >0){

cat('

<div id="warningORADOT" class="ui container"  style="display: none;">
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits, or could not find any enriched gene sets.</p>
</div>


<div id = "oradotimage">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
cat('"></object>
</div>
<a id="downORADOT" href="')
cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
cat('" download="DotPlotORA"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('

<div id = "oradotimage" style="display: none;">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
cat('"></object>
</div>
<a id="downORADOT" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
cat('" download="DotPlotORA"> 
<button class="ui button">Download Image </button>
</a>

<div id="warningORADOT" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		 Could not find any enriched gene sets.</p>
</div>
')

}


}

else{

cat('

<div id = "oradotimage" style="display: none;">
    <object style="height: 100%; width: 100%; object-fit: contain" type="image/svg+xml" data="')
cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
cat('"></object>
</div>
<a id="downORADOT" style="display: none;" href="')
cat(paste0(file.path("/HermesResults",filename, "plot1.svg")))
cat('" download="DotPlotORA"> 
<button class="ui button">Download Image </button>
</a>


<div id="warningORADOT" class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits, or could not find any enriched gene sets.</p>
</div>
')

}



cat('
 <div id="oradotloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>

</div>
</div>
')

if(GET[[3]] == "Kegg"){

cat('
<div class="ui tab segment" data-tab="thirdtop">
<div class="ui top attached tabular menu" >
  <a class="item active" data-tab="nineth">Pathways</a>

</div>


<div class="ui bottom attached tab segment active" data-tab="nineth">

<div class="ui container" style = "padding-bottom: 20px;">

<div class="ui segment">
	<h4 class="ui header">Pathway</h4>
	<p class="p-t-15" style = "font-size: large;">
	Pathway visualization. Pathways in the dropdown are ordered by adjusted p-value.   
	</p>
</div>
</div>


')

if(!is.na(edo2F)){

if(dim(edo2F)[1] >0){

cat('
<div class="ui fluid search selection dropdown">

<input type="hidden" name="KEGGPath" id="KEGGPath">  
<div class="default text">Pathways</div>
<div class="menu" >
')

flag = 0

for(i in 1:length(edo2F@result$Description)){

if(edo2F@result$p.adjust[i] >= .05 && flag == 0){

cat('
<div class="divider"></div>
<div class="header" style="text-align: center;">Adjusted P-value Over 0.05</div>
')

flag = 1
}

cat('<div class="item" data-value="')
cat(edo2F@result$ID[i])
cat('">')
cat(edo2F@result$Description[i])
cat('</div>
')

}

cat('

</div>
</div>


		<div class="m-t-25" style = "padding-top: 20px;">
			<button class="ui black button" name=KEGGplot onclick="KEGG()">
			Submit Analysis</button>
		</div>
<div class="ui divider"></div>


    <img id="KEGGimg" style="height: 100%; width: 100%; object-fit: contain" src="')

cat(paste0(file.path("/HermesResults",filename, "newpath.png")))
cat('">
<a href="')
cat(paste0(file.path("/HermesResults",filename, "newpath.png")))
cat('" download="KeggPathway"> 
<button class="ui button">Download Image </button>
</a>')

} else{

cat('

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		Could not find any enriched pathways.</p>
</div>
')

}


}

else{

cat('

<div class="ui container" >
		<p class="p-t-15" style = "font-size: large;">
		No genes found within log2foldchange limits.</p>
</div>
')

}


}


cat('

  <div id="KEGGloader" class="ui disabled inverted dimmer">
    <div class="ui text loader">Loading</div>
  </div>

</div>
</div>
</div>
<br>
<br>
')

}

}

}

}
cat('

</div>

<script>
  function gsea() {
        $.ajax({
	  type: "POST",
          url: "/HermesResults/GSEAPlot.php",
          data: {"file" : document.getElementById("gseaset").value, "file2" : "NA", "file3" : "')
cat(filename)
cat('"},
beforeSend: function () {
    document.getElementById("gseaimage").style.visibility= "hidden";
    document.getElementById("gsealoader").className = "ui active inverted dimmer";
  },
  complete: function () {
   document.getElementById("gseaimage").src = document.getElementById("gseaimage").src +"?t=" + new Date().getTime();
setTimeout(function (){
     document.getElementById("gsealoader").className = "ui disabled inverted dimmer";
     document.getElementById("gseaimage").style.visibility= "visible";

     document.getElementById("gseaimagepre").src = document.getElementById("gseaimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },
});
    
  }
</script>

<script>
  function CnetB() {

if(document.getElementById("pcnetB").value > 1 || document.getElementById("pcnetB").value < 0 || document.getElementById("cnetBlabel").value == ""){

document.getElementById("cnetBerror").style.display = "block";

}

else{

document.getElementById("cnetBerror").style.display = "none";

        $.ajax({
	  type: "POST",
          url: "/HermesResults/CNETBPlot.php",
          data: {"file" : document.getElementById("pcnetB").value, "file2" : document.getElementById("cnetBlabel").value, "file3" : "')
cat(filename)
cat('"},

beforeSend: function () {
    document.getElementById("cnetBimage").style.visibility= "hidden";
    document.getElementById("cnetBloader").className = "ui active inverted dimmer";
  },
  complete: function () {

   document.getElementById("cnetBimage").src = document.getElementById("cnetBimage").src +"?t=" + new Date().getTime();

setTimeout(function (){

     document.getElementById("cnetBimage").style.display = "block";
     document.getElementById("downCNETB").style.display = "block";
     document.getElementById("warningCNETB").style.display = "none";

     document.getElementById("cnetBloader").className = "ui disabled inverted dimmer";
     document.getElementById("cnetBimage").style.visibility= "visible";

     document.getElementById("cnetBimagepre").src = document.getElementById("cnetBimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },

});

}  
  }
</script>


<script>
  function Cnet() {

if(document.getElementById("pcnet").value > 1 || document.getElementById("pcnet").value < 0 || document.getElementById("cnetlabel").value == ""){

document.getElementById("cneterror").style.display = "block";

}

else{
document.getElementById("cneterror").style.display = "none";

var regtype = "both";

if(document.getElementById("Fold1").style.visibility == "visible" && document.getElementById("Fold2").style.visibility == "hidden"){
regtype = "down";
}

else if(document.getElementById("Fold1").style.visibility == "hidden" && document.getElementById("Fold2").style.visibility == "visible"){
regtype = "up";
}

else{
regtype = "both";

}

        $.ajax({
	  type: "POST",
          url: "/HermesResults/CNETPlot.php",
          data: {"file" : document.getElementById("pcnet").value, "file2" : document.getElementById("cnetlabel").value, "file3" : document.getElementById("Fold1").value, "file4" : document.getElementById("Fold2").value, "file5" : regtype, "file6" : "')
cat(filename)
cat('", "file7" : "')

cat(GET[[3]])
cat('", "file8" : "')
cat(GET[[4]])
cat('", "file9" : "')
cat(GET[[5]])

cat('"},

beforeSend: function () {
    document.getElementById("cnetimage").style.visibility= "hidden";
    document.getElementById("cnetloader").className = "ui active inverted dimmer";
  },
  complete: function () {

   document.getElementById("cnetimage").src = document.getElementById("cnetimage").src +"?t=" + new Date().getTime();

setTimeout(function (){

     document.getElementById("cnetimage").style.display = "block";
     document.getElementById("downCNET").style.display = "block";
     document.getElementById("warningCNET").style.display = "none";


     document.getElementById("cnetloader").className = "ui disabled inverted dimmer";
     document.getElementById("cnetimage").style.visibility= "visible";

     document.getElementById("cnetimagepre").src = document.getElementById("cnetimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },

});
  }  
  }
</script>

<script>
  function Oradot() {

if(document.getElementById("poradot").value > 1 || document.getElementById("poradot").value < 0 || document.getElementById("toporadot").value < 0 || document.getElementById("toporadot").value > 100){

document.getElementById("oradoterror").style.display = "block";

}

else{

document.getElementById("oradoterror").style.display = "none";

var regtype = "both";

if(document.getElementById("Fold14").style.visibility == "visible" && document.getElementById("Fold24").style.visibility == "hidden"){
regtype = "down";
}

else if(document.getElementById("Fold14").style.visibility == "hidden" && document.getElementById("Fold24").style.visibility == "visible"){
regtype = "up";
}

else{
regtype = "both";

}

        $.ajax({
	  type: "POST",
          url: "/HermesResults/ORADOTPlot.php",
          data: {"file" : document.getElementById("poradot").value, "file2" : document.getElementById("toporadot").value, "file3" : document.getElementById("Fold14").value, "file4" : document.getElementById("Fold24").value, "file5" : regtype, "file6" : "')
cat(filename)
cat('", "file7" : "')

cat(GET[[3]])
cat('", "file8" : "')
cat(GET[[4]])
cat('", "file9" : "')
cat(GET[[5]])

cat('"},

beforeSend: function () {
    document.getElementById("oradotimage").style.visibility= "hidden";
    document.getElementById("oradotloader").className = "ui active inverted dimmer";
  },
  complete: function () {
document.getElementById("oradotimage").innerHTML += "";

setTimeout(function (){

     document.getElementById("oradotimage").style.display = "block";
     document.getElementById("downORADOT").style.display = "block";
     document.getElementById("warningORADOT").style.display = "none";

     document.getElementById("oradotloader").className = "ui disabled inverted dimmer";
     document.getElementById("oradotimage").style.visibility= "visible";

     document.getElementById("oradotimagepre").src = document.getElementById("oradotimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },

});
  }  
  }
</script>


<script>
  function Emap() {

if(document.getElementById("pemap").value > 1 || document.getElementById("pemap").value < 0){

document.getElementById("emaperror").style.display = "block";

}

else{

document.getElementById("emaperror").style.display = "none";


var regtype = "both";

if(document.getElementById("Fold12").style.visibility == "visible" && document.getElementById("Fold22").style.visibility == "hidden"){
regtype = "down";
}

else if(document.getElementById("Fold12").style.visibility == "hidden" && document.getElementById("Fold22").style.visibility == "visible"){
regtype = "up";
}

else{
regtype = "both";

}

        $.ajax({
	  type: "POST",
          url: "/HermesResults/EMAPPlot.php",
          data: {"file" : document.getElementById("pemap").value, "file2" : document.getElementById("cnetlabel").value, "file3" : document.getElementById("Fold12").value, "file4" : document.getElementById("Fold22").value, "file5" : regtype, "file6" : "')
cat(filename)
cat('", "file7" : "')

cat(GET[[3]])
cat('", "file8" : "')
cat(GET[[4]])
cat('", "file9" : "')
cat(GET[[5]])

cat('"},

beforeSend: function () {
    document.getElementById("emapimage").style.visibility= "hidden";
    document.getElementById("emaploader").className = "ui active inverted dimmer";
  },
  complete: function () {

   document.getElementById("emapimage").src = document.getElementById("emapimage").src +"?t=" + new Date().getTime();

setTimeout(function (){

     document.getElementById("emapimage").style.display = "block";
     document.getElementById("downEMAP").style.display = "block";
     document.getElementById("warningEMAP").style.display = "none";

     document.getElementById("emaploader").className = "ui disabled inverted dimmer";
     document.getElementById("emapimage").style.visibility= "visible";

     document.getElementById("emapimagepre").src = document.getElementById("emapimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },

});
}
  }
</script>

<script>
  function Upset() {

if(document.getElementById("pupset").value > 1 || document.getElementById("pupset").value < 0){

document.getElementById("upseterror").style.display = "block";

}

else{

document.getElementById("upseterror").style.display = "none";


var regtype = "both";

if(document.getElementById("Fold13").style.visibility == "visible" && document.getElementById("Fold23").style.visibility == "hidden"){
regtype = "down";
}

else if(document.getElementById("Fold13").style.visibility == "hidden" && document.getElementById("Fold23").style.visibility == "visible"){
regtype = "up";
}

else{
regtype = "both";

}

        $.ajax({
	  type: "POST",
          url: "/HermesResults/UPSETPlot.php",
          data: {"file" : document.getElementById("pupset").value, "file2" : "NA", "file3" : document.getElementById("Fold13").value, "file4" : document.getElementById("Fold23").value, "file5" : regtype, "file6" : "')
cat(filename)
cat('", "file7" : "')

cat(GET[[3]])
cat('", "file8" : "')
cat(GET[[4]])
cat('", "file9" : "')
cat(GET[[5]])

cat('"},

beforeSend: function () {
    document.getElementById("upsetimage").style.visibility= "hidden";
    document.getElementById("upsetloader").className = "ui active inverted dimmer";
  },
  complete: function () {
	document.getElementById("upsetimage").innerHTML += "";

setTimeout(function (){

     document.getElementById("upsetimage").style.display = "block";
     document.getElementById("downUPSET").style.display = "block";
     document.getElementById("warningUPSET").style.display = "none";


     document.getElementById("upsetloader").className = "ui disabled inverted dimmer";
     document.getElementById("upsetimage").style.visibility= "visible";
     document.getElementById("upsetimagepre").src = document.getElementById("upsetimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },

});
  }  
  }
</script>


<script>
  function ridge() {

if(document.getElementById("pridge").value > 1 || document.getElementById("pridge").value < 0){

document.getElementById("ridgeerror").style.display = "block";

}

else{

document.getElementById("ridgeerror").style.display = "none";


        $.ajax({
	  type: "POST",
          url: "/HermesResults/RIDGEPlot.php",
          data: {"file" : document.getElementById("pridge").value, "file3" : "')
cat(filename)
cat('"},

beforeSend: function () {
    document.getElementById("ridgeimage").style.visibility= "hidden";
    document.getElementById("ridgeloader").className = "ui active inverted dimmer";
  },
  complete: function () {

document.getElementById("ridgeimage").innerHTML += "";

setTimeout(function (){
     
     document.getElementById("ridgeimage").style.display = "block";
     document.getElementById("downRIDGE").style.display = "block";
     document.getElementById("warningRIDGE").style.display = "none";

     document.getElementById("ridgeloader").className = "ui disabled inverted dimmer";
     document.getElementById("ridgeimage").style.visibility= "visible";

     document.getElementById("ridgeimagepre").src = document.getElementById("ridgeimagepre").src +"?t=" + new Date().getTime();

}, 5000);
  },

});
  }  
  }
</script>

<script>
  function gseadot() {

if(document.getElementById("pgseadot").value > 1 || document.getElementById("pgseadot").value < 0 || document.getElementById("topgseadot").value < 0 || document.getElementById("topgseadot").value > 100){

document.getElementById("gseadoterror").style.display = "block";

}

else{

document.getElementById("gseadoterror").style.display = "none";

        $.ajax({
	  type: "POST",
          url: "/HermesResults/GSEADOTPlot.php",
          data: {"file" : document.getElementById("pgseadot").value, "file2" : document.getElementById("topgseadot").value, "file3" : "')
cat(filename)
cat('"},

beforeSend: function () {
    document.getElementById("gseadotimage").style.visibility= "hidden";
    document.getElementById("gseadotloader").className = "ui active inverted dimmer";
  },
  complete: function () {

document.getElementById("gseadotimage").innerHTML += "";

setTimeout(function (){

     document.getElementById("gseadotimage").style.display = "block";
     document.getElementById("downGSEADOT").style.display = "block";
     document.getElementById("warningGSEADOT").style.display = "none";


     document.getElementById("gseadotloader").className = "ui disabled inverted dimmer";
     document.getElementById("gseadotimage").style.visibility= "visible";

     document.getElementById("gseadotimagepre").src = document.getElementById("gseadotimagepre").src +"?t=" + new Date().getTime();
}, 5000);
  },

});
  }  
  }
</script>

<script>
  function KEGG() {

        $.ajax({
	  type: "POST",
          url: "/HermesResults/KEGGPlot.php",
          data: {"file" : document.getElementById("KEGGPath").value, "file3" : "')
cat(filename)
cat('"},

beforeSend: function () {
    document.getElementById("KEGGimg").style.visibility= "hidden";
    document.getElementById("KEGGloader").className = "ui active inverted dimmer";
  },
  complete: function () {

   document.getElementById("KEGGimg").src = document.getElementById("KEGGimg").src +"?t=" + new Date().getTime();

setTimeout(function (){
     document.getElementById("KEGGloader").className = "ui disabled inverted dimmer";
     document.getElementById("KEGGimg").style.visibility= "visible";

}, 5000);
  },

});
    
  }
</script>

<script>
  function maintop() {
document.getElementById("filter").value = document.getElementById("filter").value.replace(/\\\\/g, "");

document.getElementById("Mainloader").className = "ui active inverted dimmer";

  }
</script>

<script>
  function filterchange() 
{
document.getElementById("filter").value = document.getElementById("filter").value.replace(/\\\\/g, "");
}

</script>


<script>

function ShowOptions()
{
    if(document.getElementById("analysisoptions").value == "GO"){
		document.getElementById("GOdisplay").style.display = "block";
		document.getElementById("MSigdisplay").style.display = "none";
		}
     else if (document.getElementById("analysisoptions").value == "Msi"){
	document.getElementById("MSigdisplay").style.display = "block";
	document.getElementById("GOdisplay").style.display = "none";
       }
    else{
		document.getElementById("GOdisplay").style.display = "none";
		document.getElementById("MSigdisplay").style.display = "none";

		}

	
}

</script>

<script>

function foldchange(clickedid)
{
    if(clickedid == "downreg"){

		document.getElementById("FoldFilterValue").value = "downreg";

		document.getElementById("foldtopT").style.display = "block";
		document.getElementById("Fold1T").style.visibility = "visible";
		document.getElementById("Fold2T").style.visibility = "hidden";
		document.getElementById("folddownT").style.display = "block";
		document.getElementById("foldupT").style.display = "none";
		document.getElementById("foldbothT").style.display = "none";

		}
     else if (clickedid == "upreg"){

		document.getElementById("FoldFilterValue").value = "upreg";

		document.getElementById("foldtopT").style.display = "block";
		document.getElementById("Fold1T").style.visibility = "hidden";
		document.getElementById("Fold2T").style.visibility = "visible";
		document.getElementById("folddownT").style.display = "none";
		document.getElementById("foldupT").style.display = "block";
		document.getElementById("foldbothT").style.display = "none";


       }
    else if (clickedid == "both"){

		document.getElementById("FoldFilterValue").value = "both";

		document.getElementById("foldtopT").style.display = "block";
		document.getElementById("Fold1T").style.visibility = "visible";
		document.getElementById("Fold2T").style.visibility = "visible";
		document.getElementById("folddownT").style.display = "none";
		document.getElementById("foldupT").style.display = "none";
		document.getElementById("foldbothT").style.display = "block";
}
}

</script>

<script> 

function changepreanalysis()
{

if(document.getElementById("preanalysisoptions").style.display == "block")
{
document.getElementById("preanalysisoptions").style.display = "none";
}

else{
document.getElementById("preanalysisoptions").style.display = "block";
}

}

</script>

<script> 

function changepostanalysis()
{

if(document.getElementById("postanalysisoptions").style.display == "block")
{
document.getElementById("postanalysisoptions").style.display = "none";
}

else{
document.getElementById("postanalysisoptions").style.display = "block";
}

}

</script>

<script>
function foldchange2(clickedid)
{
    if(clickedid == "downreg"){
		document.getElementById("foldtop").style.display = "block";
		document.getElementById("Fold1").style.visibility = "visible";
		document.getElementById("Fold2").style.visibility = "hidden";
		document.getElementById("folddown").style.display = "block";
		document.getElementById("foldup").style.display = "none";
		document.getElementById("foldboth").style.display = "none";
		document.getElementById("foldtop2").style.display = "block";
		document.getElementById("Fold12").style.visibility = "visible";
		document.getElementById("Fold22").style.visibility = "hidden";
		document.getElementById("folddown2").style.display = "block";
		document.getElementById("foldup2").style.display = "none";
		document.getElementById("foldboth2").style.display = "none";
		document.getElementById("foldtop3").style.display = "block";
		document.getElementById("Fold13").style.visibility = "visible";
		document.getElementById("Fold23").style.visibility = "hidden";
		document.getElementById("folddown3").style.display = "block";
		document.getElementById("foldup3").style.display = "none";
		document.getElementById("foldboth3").style.display = "none";
		document.getElementById("foldtop4").style.display = "block";
		document.getElementById("Fold14").style.visibility = "visible";
		document.getElementById("Fold24").style.visibility = "hidden";
		document.getElementById("folddown4").style.display = "block";
		document.getElementById("foldup4").style.display = "none";
		document.getElementById("foldboth4").style.display = "none";
		}
     else if (clickedid == "upreg"){
		document.getElementById("foldtop").style.display = "block";
		document.getElementById("Fold1").style.visibility = "hidden";
		document.getElementById("Fold2").style.visibility = "visible";
		document.getElementById("folddown").style.display = "none";
		document.getElementById("foldup").style.display = "block";
		document.getElementById("foldboth").style.display = "none";
		document.getElementById("foldtop2").style.display = "block";
		document.getElementById("Fold12").style.visibility = "hidden";
		document.getElementById("Fold22").style.visibility = "visible";
		document.getElementById("folddown2").style.display = "none";
		document.getElementById("foldup2").style.display = "block";
		document.getElementById("foldboth2").style.display = "none";
		document.getElementById("foldtop3").style.display = "block";
		document.getElementById("Fold13").style.visibility = "hidden";
		document.getElementById("Fold23").style.visibility = "visible";
		document.getElementById("folddown3").style.display = "none";
		document.getElementById("foldup3").style.display = "block";
		document.getElementById("foldboth3").style.display = "none";
		document.getElementById("foldtop4").style.display = "block";
		document.getElementById("Fold14").style.visibility = "hidden";
		document.getElementById("Fold24").style.visibility = "visible";
		document.getElementById("folddown4").style.display = "none";
		document.getElementById("foldup4").style.display = "block";
		document.getElementById("foldboth4").style.display = "none";
       }
    else if (clickedid == "both"){
		document.getElementById("foldtop").style.display = "block";
		document.getElementById("Fold1").style.visibility = "visible";
		document.getElementById("Fold2").style.visibility = "visible";
		document.getElementById("folddown").style.display = "none";
		document.getElementById("foldup").style.display = "none";
		document.getElementById("foldboth").style.display = "block";
		document.getElementById("foldtop2").style.display = "block";
		document.getElementById("Fold12").style.visibility = "visible";
		document.getElementById("Fold22").style.visibility = "visible";
		document.getElementById("folddown2").style.display = "none";
		document.getElementById("foldup2").style.display = "none";
		document.getElementById("foldboth2").style.display = "block";
		document.getElementById("foldtop3").style.display = "block";
		document.getElementById("Fold13").style.visibility = "visible";
		document.getElementById("Fold23").style.visibility = "visible";
		document.getElementById("folddown3").style.display = "none";
		document.getElementById("foldup3").style.display = "none";
		document.getElementById("foldboth3").style.display = "block";
		document.getElementById("foldtop4").style.display = "block";
		document.getElementById("Fold14").style.visibility = "visible";
		document.getElementById("Fold24").style.visibility = "visible";
		document.getElementById("folddown4").style.display = "none";
		document.getElementById("foldup4").style.display = "none";
		document.getElementById("foldboth4").style.display = "block";
		}
	
}
</script>


<script type = "text/javascript" src="/AresC/jquery.min.js"></script>
<script type = "text/javascript" src="/AresC/semantic.min.js"></script>
<script type = "text/javascript"> 

$(".ui.dropdown")
  .dropdown()
;

</script>
<script type = "text/javascript">
$("#context1 .menu .item").tab({context: $("#context1")});
</script>
<script type = "text/javascript">
$(".tabular.menu .item").tab();
</script> ')
} else{

cat('

<p class="p-t-15" style = "font-size: large; ">
		There was a problem finding your account. Please contact us for assistance. </p>

')


}

cat('
</body>
</html>
')


