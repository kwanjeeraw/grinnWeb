#------Load required libraries
library(KEGGREST)
library(KEGGgraph)
library(RNeo4j)
library(UniProt.ws)
library(biomaRt)
library(XML)
library(RCurl)
library(jsonlite)

#------Connect to Neo4j DB
graph = startGraph("[URL]") #url of database; e.g. http://localhost:7474/db/data/
#clear(graph) #clear database

#------Upload global data
####add metabolite node
cat("... Loading compounds from KEGG\n")
CMP = keggList("compound") #load all kegg compounds
cat("... Returning ", length(CMP), " compounds from KEGG\n")
addConstraint(graph, "Metabolite", "mid")
for (i in 1:length(CMP)){
  comp = keggGet(names(CMP[i]))
  namelist = strsplit(comp[[1]]$NAME,split=";")
  if(i %% 500 == 0) cat("... Creating metabolite node: ",i,"\n")
  createNode(graph, "Metabolite", mid = comp[[1]]$ENTRY[[1]], 
            name = namelist[[1]], alias = if(length(namelist[-1]) == 0){c("")} else{unlist(namelist[-1])}, 
            xref = if(is.null(comp[[1]]$DBLINKS)){c("")} else{unlist(strsplit(comp[[1]]$DBLINKS,split="\t"))})
}

####add reaction node
cat("... Loading reactions from KEGG\n")
RX = keggList("reaction") #load all kegg reactions
cat("... Returning ", length(RX), " reactions from KEGG\n")
addConstraint(graph, "Reaction", "rid")
for (i in 1:length(RX)){
  reac = keggGet(names(RX[i]))
  if(i %% 100 == 0) cat("... Creating reaction node: ",i,"\n")
  createNode(graph, "Reaction", rid = reac[[1]]$ENTRY[[1]], 
            name = ifelse(is.null(reac[[1]]$NAME), "", reac[[1]]$NAME), equation = reac[[1]]$EQUATION, 
            definition = reac[[1]]$DEFINITION)
}

####add rpair node; output rpair-metabolite, reaction-rpair relations in csv files; relations can be duplicate!
cat("... Loading RPAIR from KEGG\n")
RP = keggList("rpair") #load all kegg reactions
cat("... Returning ", length(RP), " RPAIR from KEGG\n")
metRpair = data.frame()
rxRpair = data.frame()
addConstraint(graph, "Rpair", "rpid")

#function to generate reaction-rpair
genRxRpair <- function(x,y) {
	tmp=cbind(as.data.frame(strsplit(x,split=" ")),y)
	colnames(tmp) = c("c1","c2")
	tmp
}

for (i in 1:length(RP)){
  rpair = keggGet(names(RP[i]))
  createNode(graph, "Rpair", rpid = rpair[[1]]$ENTRY[[1]], name = rpair[[1]]$NAME, type = unlist(strsplit(rpair[[1]]$TYPE,split=" ")))
  metRpair = unique(rbind(metRpair,as.data.frame(cbind(rpair[[1]]$ENTRY[[1]],c(rpair[[1]]$ENTRY1$COMPOUND,rpair[[1]]$ENTRY2$COMPOUND))))) #genetate rpair-metabolite relations
  rxRpair = unique(rbind(rxRpair,do.call(rbind,lapply(rpair[[1]]$REACTION, FUN=genRxRpair, y=rpair[[1]]$ENTRY[[1]])))) #genetate reaction-rpair relations
  if(i %% 500 == 0) cat("... Creating and Merging RPAIR: ",i,"\n")
}
write.table(metRpair,file = "rpair_contain_met.csv", row.names=FALSE, sep=",", col.names=FALSE) #output relations
write.table(rxRpair,file = "rx_has_rpair.csv", row.names=FALSE, sep=",", col.names=FALSE) #output relations

####add regulator node; enhancer
addConstraint(graph, "Regulator", "rgid")
cat("... Loading enhancer from local file\n")
ENH = cbind(paste("ENH",seq(1,43011),sep=""),read.delim(file ="[PATH]", sep="\t", header=FALSE)) #path to local file; e.g. hg19_enhancers.bed; columns separated by \t; column#4 contains positions
cat("... Uploading", nrow(ENH), "enhancers\n")
for (i in 1:nrow(ENH)){
	if(i %% 500 == 0) cat("... Creating regulator node: ",i,"\n")
	createNode(graph, "Regulator", rgid = as.character(ENH[i,1]), 
            name = as.character(ENH[i,5]), type = "enhancer", 
            info = as.character(ENH[i,5]))
}
cat("... Done. \n")

#------Initialize global variables
addConstraint(graph, "Pathway", "pid")
addConstraint(graph, "Protein", "prid")
addConstraint(graph, "Gene", "gid")
uniprot = useMart("unimart") #use uniprot
uniprot = useMart("unimart", dataset="uniprot")
ensembl = useMart("ensembl") #use ensemble
genePrt = data.frame()
proteinRx = data.frame()
substRx = data.frame()
prodRx = data.frame()
log_con <- file("log.txt", open="a") #list of unmapped proteins

#------Set organism on ensembl
ensembl = useMart("ensembl", dataset="[ORGANISM]") #e.g. hsapiens_gene_ensembl; to show datasets listDatasets(ensembl)

#------Load model from SBML file
xSBML = xmlRoot(xmlTreeParse(file ="[PATH]")) #path to local file; e.g. models/hsa/BMID000000101155/BMID000000101155_url.xml
modelNode = xmlElementsByTagName(xSBML,"model")$model #model node

####add pathway node
cat("... Creating pathway node\n")
annoNode = xmlElementsByTagName(modelNode,"annotation")$annotation #annotation node
rdfDes = xmlToList(xmlElementsByTagName(xmlChildren(annoNode)[[1]],'Description')$Description) #list of rdf:Description
pid = as.list(xmlAttrs(modelNode))$id #get id
organism = lookupUniprotSpeciesFromTaxId(gsub("\\D","",rdfDes$occursIn[[1]])) #get organism
taxId(UniProt.ws) = gsub("\\D","",rdfDes$occursIn[[1]]) #set species
xref = gsub("http://identifiers.org/kegg.pathway/","",rdfDes$isDerivedFrom[[1]][1]) #get kegg ref
tryCatch(createNode(graph, "Pathway", pid=pid, name=as.list(xmlAttrs(modelNode))$name, organism=organism, xref=xref), finally = cat("\n")) #create pathway node

####add protein and gene node; output protein-reaction, gene-protein relations in csv files; relations can be duplicate!
lsSpecies = xmlElementsByTagName(modelNode,"listOfSpecies")$listOfSpecies #get list of species
df = as.data.frame(xmlSApply(lsSpecies, xmlAttrs)) #df of species
ind = which(apply(df, 2, function(x) any(grepl("SBO:0000252", x)))) #find index of enzymes with SBO term

##loop through list of enzyme species
cat("... Loading list of ", length(lsSpecies[ind])," enzyme species\n")
for (i in 1:length(lsSpecies[ind])){
#annSp = xmlElementsByTagName(lsSpecies[ind][[1]],"annotation")$annotation #annotation node; for testing
	annSp = xmlElementsByTagName(lsSpecies[ind][[i]],"annotation")$annotation #annotation node
	rdfSp = xmlToList(xmlElementsByTagName(xmlChildren(annSp)[[1]],'Description')$Description) #list of rdf:Description
	
	##if uniprot is not null
	if(length(grep("uniprot",rdfSp))){
		unpId = gsub("http://identifiers.org/uniprot/","",rdfSp[grep("uniprot",rdfSp)][[1]]) #get uniprot id
		rxId = gsub("http://identifiers.org/kegg.reaction/","",rdfSp[grep("kegg.reaction",rdfSp)][[1]]) #get reaction id
		
		##mapping uniprot id with biomart; organism specific
		attributes = c("accession","protein_name","ec_number","organism") #mapped columns
		unpBM = getBM(attributes=attributes, filters = "accession", values = unpId, mart = uniprot)
		
		columns = c("uniprot_swissprot_accession","entrezgene","external_gene_name","description","ensembl_gene_id", "ensembl_transcript_id","ensembl_peptide_id") #mapped columns
		ensBM = getBM(attributes=columns, filters = "uniprot_swissprot_accession", values = unpId, mart = ensembl)
		mappedId = merge(unpBM,ensBM,by.x = "accession", by.y ="uniprot_swissprot_accession")
		
		##loop through each uniprot id
		for (j in 1:length(unpId)){
		  prid = unpId[j]
		  unpInd = which(apply(mappedId, 1, function(x) any(grepl(prid, x)))) 
		  ##if uniprot is found in biomart
		  if(length(unpInd)){
			  pname = unique(mapply(mappedId[unpInd,],FUN=list)$protein_name)
			  ec = unique(mapply(mappedId[unpInd,],FUN=list)$ec_number)
			  pxref = unique(mapply(mappedId[unpInd,],FUN=list)$ensembl_peptide_id)
			  transcript = unique(mapply(mappedId[unpInd,],FUN=list)$ensembl_transcript_id)
			  organism = unique(mapply(mappedId[unpInd,],FUN=list)$organism)
			  
			  cat("... Creating and Merging protein and gene: ",prid,"\n")
			  tryCatch(createNode(graph, "Protein", prid=prid, name=pname, ec=ec, xref=pxref, transcript=transcript, organism=organism), error = function(e) e, finally = print("... Ignore error: protein already exists")) #create protein node
			  proteinRx = unique(rbind(proteinRx,cbind(rep(prid, times = length(rxId)), rxId))) #collect protein-reaction relations
			  
			  gid = paste0("entrez.",unique(mapply(mappedId[unpInd,],FUN=list)$entrezgene), sep="")
			  gname = unique(mapply(mappedId[unpInd,],FUN=list)$external_gene_name)
			  description = unique(mapply(mappedId[unpInd,],FUN=list)$description)
			  gxref = unique(mapply(mappedId[unpInd,],FUN=list)$ensembl_gene_id)
			  tryCatch(createNode(graph, "Gene", gid=gid, name=gname, description=description, xref=gxref, organism=organism), error = function(e) e, finally = print("... Ignore error: gene already exists")) #create gene node
			  genePrt = unique(rbind(genePrt,cbind(gid, prid))) #collect gene-protein relations
		  }
		  else{
				cat("ID: ",prid, file = log_con, sep="")
				print("Unable to upload: ")
				print(prid) #unmapped uniprot
				cat("\n", file = log_con)
		  }
		} #end for mapping
	} #end if grep uniprot

} #end for enzyme species

####update reaction type; output reaction-pathway, metabolite-reaction relations in csv files; relations can be duplicate!
lsReactions = xmlElementsByTagName(modelNode,"listOfReactions")$listOfReactions #get list of reactions
dfr = as.data.frame(xmlSApply(lsReactions, xmlAttrs)) #df of reactions

cat("... Writing out rx-ptw relation file\n")
write.table(cbind(apply(dfr["id",], 1,FUN=function(x) gsub("rn","", x)), rep(pid, times = length(lsReactions))),file = paste("rx_ispartof_ptw_",pid,".csv", sep=""), 
row.names=FALSE, append=TRUE, sep=",", col.names=FALSE) #output relations

##loop through list of reactions
cat("... Loading list of ", length(lsReactions)," reaction species\n")
for (i in 1:length(lsReactions)){
	if(i %% 10 == 0) cat("... Merging metabolite and reaction: ",i,"\n")
#rxSp = lsReactions[[1]] #for testing
	rxSp = lsReactions[[i]]
	rid = gsub("rn","",xmlAttrs(rxSp)["id"][[1]])
	
	rxNode = getUniqueNode(graph, "Reaction", rid = rid) #get reaction node
	updateProp(rxNode, reversible = xmlAttrs(rxSp)["reversible"][[1]]) #update reaction type
	
	lsReactants = xmlElementsByTagName(rxSp,"listOfReactants")$listOfReactants #get list of reactants
	attRt = xmlSApply(lsReactants, xmlAttrs)	
	dfrt = cbind(apply(as.data.frame(attRt["name",]), 1,FUN=function(x) gsub("[^C0-9]","", x)),rid)
	colnames(dfrt) = colnames(substRx)
	substRx = unique(rbind(substRx,dfrt)) #collection of substrates
	
	lsProducts = xmlElementsByTagName(rxSp,"listOfProducts")$listOfProducts #get list of products
	attPd = xmlSApply(lsProducts, xmlAttrs)
	dfpd = cbind(rid,apply(as.data.frame(attPd["name",]), 1,FUN=function(x) gsub("[^C0-9]","", x)))
	colnames(dfpd) = colnames(prodRx)
	prodRx = unique(rbind(prodRx,dfpd)) #collection of products
}

#------Relation output in csv files
write.table(proteinRx,file = "protein_catalyze_rx.csv", row.names=FALSE, sep=",", col.names=FALSE) #output relations
write.table(genePrt,file = "gene_encode_protein.csv", row.names=FALSE, sep=",", col.names=FALSE) #output relations
write.table(substRx,file = "met_substrateof_rx.csv", row.names=FALSE, sep=",", col.names=FALSE) #output relations
write.table(prodRx,file = "rx_produce_met.csv", row.names=FALSE, sep=",", col.names=FALSE) #output relations
close(log_con)

####add gene node; output regulator-gene relations in csv file; relations can be dubplicate!
#peakAnno = read.delim(file = "Documents/UCD/KEGG_DB/hg19.cage_peak_ann.txt")
  