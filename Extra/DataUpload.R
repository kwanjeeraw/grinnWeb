#------Load required libraries
library(KEGGREST)
library(KEGGgraph)
library(RNeo4j)
library(UniProt.ws)
library(biomaRt)
library(XML)
library(RCurl)
library(jsonlite)
library(stringr)
library (plyr)
library(RbioRXN)

##..Code to upload KEGG-based data
##..metabolite, reaction, pathway, gene = KEGG id
##..protein = uniprot id
#------Connect to Neo4j DB
graph = startGraph("[URL]") #url of database; e.g. http://localhost:7474/db/data/
#clear(graph) #clear database

#------Upload global data
####add metabolite node
####KEGG rest contains no description, no InChI
cat("... Loading compounds from KEGG\n")
CMP = keggList("compound") #load all kegg compounds
cat("... Returning ", length(CMP), " compounds from KEGG\n")
addConstraint(graph, "Metabolite", "GID")
for (i in 1:length(CMP)){
  comp = keggGet(names(CMP[i]))
  namelist = strsplit(comp[[1]]$NAME,split=";")
  keggid = paste0("KEGG:",comp[[1]]$ENTRY[[1]])
  if(i %% 500 == 0) cat("... Creating metabolite node: ",i,"\n")
  createNode(graph, "Metabolite", GID = comp[[1]]$ENTRY[[1]], 
            name = namelist[[1]], synonym = if(length(namelist[-1]) == 0){c("")} else{unlist(namelist[-1])}, 
            formula = ifelse(is.null(comp[[1]]$FORMULA[[1]]), "", comp[[1]]$FORMULA[[1]]), 
            molWeight = ifelse(is.null(comp[[1]]$MOL_WEIGHT[[1]]), 0, comp[[1]]$MOL_WEIGHT[[1]]), 
            xref = if(is.null(comp[[1]]$DBLINKS)){c("")} else{c(keggid,unlist(gsub(" ","",strsplit(comp[[1]]$DBLINKS,split="\t"))))})
}

####add kegg glycans to metabolite node
####KEGG rest contains no molWeight
GLY = keggList("glycan") #load all kegg glycans
cat("... Returning ", length(GLY), " glycans from KEGG\n")
for (i in 1:length(GLY)){
  comp = keggGet(names(GLY[i]))
  keggid = paste0("KEGG:",comp[[1]]$ENTRY[[1]])
  if(i %% 500 == 0) cat("... Creating metabolite node: ",i,"\n")
  createNode(graph, "Metabolite", GID = comp[[1]]$ENTRY[[1]], 
             name = ifelse(is.null(comp[[1]]$NAME[[1]]),"", comp[[1]]$NAME[[1]]), synonym = if(length(namelist[-1]) == 0){c("")} else{unlist(namelist[-1])}, 
             formula = ifelse(is.null(comp[[1]]$COMPOSITION[[1]]), "", comp[[1]]$COMPOSITION[[1]]), 
             molWeight = ifelse(is.null(comp[[1]]$MOL_WEIGHT[[1]]), 0, comp[[1]]$MOL_WEIGHT[[1]]), 
             xref = if(is.null(comp[[1]]$DBLINKS)){c("")} else{c(keggid,unlist(gsub(" ","",strsplit(comp[[1]]$DBLINKS,split="\t"))))})
}

####add InChI to metabolite node
chebi = get.ChEBI()
parsed_ChEBI = parse.ChEBI(chebi)
ind = which(parsed_ChEBI$KEGG == "")
kegg_ChEBI = parsed_ChEBI[-ind,]

for(i in 1:nrow(kegg_ChEBI)){
  querystring = paste0("MATCH (n { GID: '",kegg_ChEBI$KEGG[i],"' }) SET n.InChI = '",kegg_ChEBI$InChI[i], "' RETURN n")
  grinn::curlRequestCypher(querystring)
  print(cat(i,' loaded'))
}

####add reaction node
####KEGG rest contains no xref
####Reversible type gets updated from a model file
cat("... Loading reactions from KEGG\n")
RX = keggList("reaction") #load all kegg reactions
cat("... Returning ", length(RX), " reactions from KEGG\n")
addConstraint(graph, "Reaction", "GID")
for (i in 1:length(RX)){
  reac = keggGet(names(RX[i]))
  if(i %% 100 == 0) cat("... Creating reaction node: ",i,"\n")
  createNode(graph, "Reaction", GID = reac[[1]]$ENTRY[[1]], 
            name = ifelse(is.null(reac[[1]]$NAME), "", reac[[1]]$NAME),  
            description = reac[[1]]$DEFINITION, equation = reac[[1]]$EQUATION)
}

####add regulator node; enhancer
####future use
addConstraint(graph, "Regulator", "GID")
cat("... Loading enhancer from local file\n")
ENH = cbind(paste("ENH",seq(1,43011),sep=""),read.delim(file ="[PATH]", sep="\t", header=FALSE)) #path to local file; e.g. hg19_enhancers.bed; columns separated by \t; column#4 contains positions
cat("... Uploading", nrow(ENH), "enhancers\n")
for (i in 1:nrow(ENH)){
  if(i %% 500 == 0) cat("... Creating regulator node: ",i,"\n")
  createNode(graph, "Regulator", GID = as.character(ENH[i,1]), 
             name = as.character(ENH[i,5]), type = "enhancer", 
             info = as.character(ENH[i,5]))
}
cat("... Done. \n")

#------Initialize global variables
####GID = database id, used in every node type
addConstraint(graph, "Pathway", "GID")
addConstraint(graph, "Protein", "GID")
addConstraint(graph, "Gene", "GID")

uniprot = useMart("unimart") #use uniprot
uniprot = useMart("unimart", dataset="uniprot")
ensembl = useMart("ensembl") #use ensemble

ptwRx = data.frame()
genePrt = data.frame()
proteinRx = data.frame()
substRx = data.frame()
prodRx = data.frame()
log_con <- file("log.txt", open="a") #list of unmapped proteins

#------Set organism on ensembl
ensembl = useMart("ensembl", dataset="[ORGANISM]") #e.g. hsapiens_gene_ensembl; to show datasets listDatasets(ensembl)

#------Load model from SBML file
file = "[FILEPATH]"  #path to local file; e.g. models/hsa/BMID000000101155/BMID000000101155_url.xml
xSBML = xmlRoot(xmlTreeParse(file = file))
modelNode = xmlElementsByTagName(xSBML,"model")$model #model node

####add pathway node
####Biomodels contains no synonym, no description; assume manual curation
cat("... Creating pathway node\n")
annoNode = xmlElementsByTagName(modelNode,"annotation")$annotation #annotation node
rdfDes = xmlToList(xmlElementsByTagName(xmlChildren(annoNode)[[1]],'Description')$Description) #list of rdf:Description
ptwGID = gsub("meta_path_","",as.list(xmlAttrs(modelNode))$metaid) #use metaid|KEGG as GID

##organism is a global variable used for other node types
organism = lookupUniprotSpeciesFromTaxId(gsub("\\D","",rdfDes$occursIn[[1]])) #get organism

##!!!****slow, load first time running only****!!!##
taxId(UniProt.ws) = gsub("\\D","",rdfDes$occursIn[[1]]) #get species for UniProt.ws

xref = c(paste0("KEGG:",ptwGID),paste0("biomodels:",as.list(xmlAttrs(modelNode))$id)) #get xref
manCurate = TRUE
refId = keggGet(ptwGID)[[1]]$REFERENCE
describeIn = as.character(ldply(refId, data.frame)$REFERENCE)
if(length(describeIn) == 0){describeIn = {c("")}} #if there is no reference

#create pathway node
tryCatch(createNode(graph, "Pathway", GID=ptwGID, name=as.list(xmlAttrs(modelNode))$name,  
                    organism=organism, manCurate=manCurate, describeIn=describeIn, xref=xref), finally = cat("\n"))

####add protein and gene node; output protein-reaction, gene-protein relations in csv files; relations can be duplicate!, need to remove duplicate or set constrain 
lsSpecies = xmlElementsByTagName(modelNode,"listOfSpecies")$listOfSpecies #get list of species
df = as.data.frame(xmlSApply(lsSpecies, xmlAttrs)) #df of species
ind = which(apply(df, 2, function(x) any(grepl("SBO:0000252", x)))) #find index of enzymes with SBO term

##loop through list of enzyme species
cat("... Loading list of ", length(lsSpecies[ind])," enzyme species\n")

#start for enzyme species
for (i in 1:length(lsSpecies[ind])){
  cat("... Starting with enz ", i," \n")
#annSp = xmlElementsByTagName(lsSpecies[ind][[1]],"annotation")$annotation #annotation node; for testing
	annSp = xmlElementsByTagName(lsSpecies[ind][[i]],"annotation")$annotation #annotation node
	rdfSp = xmlToList(xmlElementsByTagName(xmlChildren(annSp)[[1]],'Description')$Description) #list of rdf:Description
	
	##if uniprot is not null
	if(length(grep("uniprot",rdfSp))){
		unpId = gsub("http://identifiers.org/uniprot/","",rdfSp[grep("uniprot",rdfSp)][[1]]) #get uniprot id
		rxId = gsub("http://identifiers.org/kegg.reaction/","",rdfSp[grep("kegg.reaction",rdfSp)][[1]]) #get reaction id
		
		##mapping uniprot id with biomart; biomart contains no synonym, no description, no km
		attributes = c("accession","protein_name","ec_number") #mapped columns
		unpBM = getBM(attributes=attributes, filters = "accession", values = unpId, mart = uniprot)
		
		##mapping uniprot id with biomart; organism specific; biomart contains no synonym
		columns = c("uniprot_swissprot_accession","entrezgene","external_gene_name","description","ensembl_gene_id", "ensembl_transcript_id","ensembl_peptide_id","start_position","end_position") #mapped columns
		ensBM = getBM(attributes=columns, filters = "uniprot_swissprot_accession", values = unpId, mart = ensembl)
    
		##mapping uniprot id with UniProt.ws; organism specific
		columns = c("UNIGENE","KEGG")
    uniWs <- select(UniProt.ws, unpId, columns, "UNIPROTKB")

		mappedId = merge(unpBM,ensBM,by.x = "accession", by.y ="uniprot_swissprot_accession")
		mappedId = merge(mappedId,uniWs,by.x = "accession", by.y ="UNIPROTKB")
    
		##loop through each uniprot id
		for (j in 1:length(unpId)){
		  prGID = unpId[j] #use Uniprot as GID
		  unpInd = which(apply(mappedId, 1, function(x) any(grepl(prGID, x))))
		  ##if uniprot is found in biomart
		  if(length(unpInd)){
			  pname = unique(mapply(mappedId[unpInd,],FUN=list)$protein_name)
			  ec = paste0("EC:",unique(mapply(mappedId[unpInd,],FUN=list)$ec_number))
			  peptide = paste0("ENSEMBL:",unique(mapply(mappedId[unpInd,],FUN=list)$ensembl_peptide_id))
        pxref = c(paste0("UniProt:",prGID),ec,peptide)
			  
			  cat("... Creating and Merging protein: ",prGID,"\n")
			  tryCatch(createNode(graph, "Protein", GID=prGID, name=pname, organism=organism, xref=pxref), error = function(e) e, finally = print("... Ignore error: protein already exists")) #create protein node
			  proteinRx = unique(rbind(proteinRx,cbind(rep(prGID, times = length(rxId)), rxId, 0))) #collect protein-reaction relations
			  
			  geGID = unique(mapply(mappedId[unpInd,],FUN=list)$KEGG)[1] #use kegg gene as GID; duplicate is removed; need to improve
			  gname = unique(mapply(mappedId[unpInd,],FUN=list)$external_gene_name)
			  description = unique(mapply(mappedId[unpInd,],FUN=list)$description)
        entrez = paste0("entrez:",unique(mapply(mappedId[unpInd,],FUN=list)$entrezgene))
			  transcript = paste0("ENSEMBL:",unique(mapply(mappedId[unpInd,],FUN=list)$ensembl_transcript_id))
			  unigene = paste0("UniGene:",unique(mapply(mappedId[unpInd,],FUN=list)$UNIGENE))
			  ensgene = paste0("ENSEMBL:",unique(mapply(mappedId[unpInd,],FUN=list)$ensembl_gene_id))
			  location = paste0(unique(mapply(mappedId[unpInd,],FUN=list)$start_position),"...",unique(mapply(mappedId[unpInd,],FUN=list)$end_position))
        gxref = c(paste0("KEGG:",geGID),ensgene,entrez,unigene,transcript)       
        
			  cat("... Creating and Merging gene: ",geGID,"\n")
			  tryCatch(createNode(graph, "Gene", GID=geGID, name=gname, description=description, organism=organism, location=location, xref=gxref), error = function(e) e, finally = print("... Ignore error: gene already exists")) #create gene node
			  genePrt = unique(rbind(genePrt,cbind(geGID, prGID))) #collect gene-protein relations
		  }
		  else{
				cat("ID: ",prGID, file = log_con, sep="")
				print("Unable to upload: ")
				print(prGID) #unmapped uniprot
				cat("\n", file = log_con)
		  }
		} #end for mapping
  } #end if grep uniprot

} #end for enzyme species

####update reaction type, reaction KEGG xref; output reaction-pathway, metabolite-reaction relations in csv files; relations can be duplicate!
####filter only main rpair
lsReactions = xmlElementsByTagName(modelNode,"listOfReactions")$listOfReactions #get list of reactions
dfr = as.data.frame(xmlSApply(lsReactions, xmlAttrs)) #df of reactions

##loop through list of reactions
cat("... Loading list of ", length(lsReactions)," reaction species\n")

#start reaction species
for (i in 1:length(lsReactions)){
  cat("... Starting with rx ", i," \n")
	if(i %% 10 == 0) cat("... Merging metabolite and reaction: ",i,"\n")
#rxSp = lsReactions[[1]] #for testing
	rxSp = lsReactions[[i]]
	rid = gsub("rn","",xmlAttrs(rxSp)["id"][[1]])
  reversible = xmlAttrs(rxSp)["reversible"][[1]]
	
	rxNode = getUniqueNode(graph, "Reaction", GID = rid) #get reaction node
	updateProp(rxNode, reversible=reversible, xref = c(paste0("KEGG:",rid))) #update reaction type
	
  ptwRx = unique(rbind(ptwRx,data.frame(ptwGID,rid))) #collection of pathway-reactions

  #keep main metabolites only; screen by rpair
  rxRpair = keggGet(rid)[[1]]$RPAIR
  mainInd = grep("main",rxRpair)
  findMain = function(x){
    unlist(strsplit(unlist(strsplit(x," "))[3],"_"))
  }
  mainList = unlist(lapply(rxRpair[mainInd], FUN=findMain))

	lsReactants = xmlElementsByTagName(rxSp,"listOfReactants")$listOfReactants #get list of reactants
	attRt = xmlSApply(lsReactants, xmlAttrs)	
  rtInd = which(gsub("cpd:","",attRt["name",]) %in% mainList)
  attRt = data.frame(attRt[,rtInd])
	dfrt = cbind(rid,apply(as.data.frame(attRt["name",]), 1,FUN=function(x) gsub("[^C0-9]","", x)),reversible)
#dfrt = cbind(rid,apply(as.data.frame(attRt["name",]), 1,FUN=function(x) gsub("gl:","", x)),reversible)
	colnames(dfrt) = colnames(substRx)
	substRx = unique(rbind(substRx,dfrt)) #collection of substrates
	
	lsProducts = xmlElementsByTagName(rxSp,"listOfProducts")$listOfProducts #get list of products
	attPd = xmlSApply(lsProducts, xmlAttrs)
  pdInd = which(gsub("cpd:","",attPd["name",]) %in% mainList)
  attPd = data.frame(attPd[,pdInd])
	dfpd = cbind(rid,apply(as.data.frame(attPd["name",]), 1,FUN=function(x) gsub("[^C0-9]","", x)),reversible)
#dfpd = cbind(rid,apply(as.data.frame(attPd["name",]), 1,FUN=function(x) gsub("gl:","", x)),reversible)
	colnames(dfpd) = colnames(prodRx)
	prodRx = unique(rbind(prodRx,dfpd)) #collection of products

} #end for reaction species

#------Relation output in csv files
write.table(ptwRx,file = "ptw_has_rx.csv", row.names=FALSE, sep=",", col.names=c("pathway","reaction")) #output relations
write.table(proteinRx,file = "protein_catalyze_rx.csv", row.names=FALSE, sep=",", col.names=c("protein","reaction","km")) #output relations
write.table(genePrt,file = "gene_encode_protein.csv", row.names=FALSE, sep=",", col.names=c("gene","protein")) #output relations
write.table(substRx,file = "rx_transform_met.csv", row.names=FALSE, sep=",", col.names=c("reaction","metabolite","reversible")) #output relations
write.table(prodRx,file = "rx_produce_met.csv", row.names=FALSE, sep=",", col.names=c("reaction","metabolite","reversible")) #output relations

close(log_con)

##...stat 01/30/2015...##
# > dim(proteinRx)
# [1] 3228    3
# > dim(genePrt)
# [1] 1397    2
# > dim(substRx)
# [1] 3499    3
# > dim(prodRx)
# [1] 3626    3
# > dim(ptwRx)
# [1] 2956    2

####add gene node; output regulator-gene relations in csv file; relations can be dubplicate!
#peakAnno = read.delim(file = "Documents/UCD/KEGG_DB/hg19.cage_peak_ann.txt")

####Don't need rpair
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

  