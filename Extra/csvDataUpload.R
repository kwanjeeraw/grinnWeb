library(RNeo4j)
library(stringr)

##..Code to upload data in csv file format
##..Data can be generated from somewhere
#------Connect to Neo4j DB
graph = startGraph("[URL]") #url of database; e.g. http://localhost:7474/db/data/

####Creating metabolite##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

#better way to deal with NA?
mydata$molWeight[is.na(mydata$molWeight)]<-0
mydata[is.na(mydata)]<-""

for (i in 1:nrow(mydata)){
	cat("... Creating metabolite node: ",i,"\n")
	GID = as.character(mydata$GID)[[i]]
	name = as.character(mydata$name)[[i]]

	synonym = strsplit(str_sub(as.character(mydata$synonym)[[i]],3,-3),split="\",\"")
	synonym  = if(length(synonym[[1]])==0){c("")} else{unlist(strsplit(synonym[[1]],split="\t"))}
	
	description = as.character(mydata$description)[[i]]
	formula = as.character(mydata$formula)[[i]]
	molWeight = mydata$molWeight[[i]]
	InChI = as.character(mydata$InChI)[[i]]
	
	xref = strsplit(str_sub(as.character(mydata$xref)[[i]],3,-3),split="\",\"")
	xref = if(length(xref[[1]])==0){c("")} else{unlist(strsplit(xref[[1]],split="\t"))}
	
	addConstraint(graph, "Metabolite", "GID")
	createNode(graph, "Metabolite", GID=GID, name=name, synonym=synonym, 
	           description=description, formula=formula, 
	           molWeight=molWeight, InChI=InChI, xref=xref)
}

####Creating reaction##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

mydata[is.na(mydata)]<-""

for (i in 1:nrow(mydata)){
	cat("... Creating reaction node: ",i,"\n")
	GID = as.character(mydata$GID)[[i]]
	name = as.character(mydata$name)[[i]]
	description = as.character(mydata$description)[[i]]
	equation = as.character(mydata$equation)[[i]]
	reversible = as.character(mydata$reversible)[[i]]

	xref = strsplit(str_sub(as.character(mydata$xref)[[i]],3,-3),split="\",\"")
	xref = if(length(xref[[1]])==0){c("")} else{unlist(strsplit(xref[[1]],split="\t"))}
  
	addConstraint(graph, "Reaction", "GID")
	createNode(graph, "Reaction", GID=GID, name=name, description=description, equation=equation, 
	           reversible=reversible, xref=xref)
}

####Creating protein##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

mydata[is.na(mydata)]<-""

for (i in 1:nrow(mydata)){
	cat("... Creating protein node: ",i,"\n")
	GID = as.character(mydata$GID)[[i]]
	name = as.character(mydata$name)[[i]]
	synonym = strsplit(str_sub(as.character(mydata$synonym)[[i]],3,-3),split="\",\"")
	synonym  = if(length(synonym[[1]])==0){c("")} else{unlist(strsplit(synonym[[1]],split="\t"))}
	
	description = as.character(mydata$description)[[i]]
	organism = as.character(mydata$organism)[[i]]
	
	xref = strsplit(str_sub(as.character(mydata$xref)[[i]],3,-3),split="\",\"")
	xref = if(length(xref[[1]])==0){c("")} else{unlist(strsplit(xref[[1]],split="\t"))}
	
	createNode(graph, "Protein", GID=GID, name=name, synonym=synonym, 
	           description=description, organism=organism, xref=xref)
}

####Creating gene##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

#!running!###
mydata = read.csv(file="gene.csv",header=T)
mydata[is.na(mydata)]<-""

for (i in 1:nrow(mydata)){
	cat("... Creating gene node: ",i,"\n")
	GID = as.character(mydata$GID)[[i]]
	name = as.character(mydata$name)[[i]]
	synonym = strsplit(str_sub(as.character(mydata$synonym)[[i]],3,-3),split="\",\"")
	synonym  = if(length(synonym[[1]])==0){c("")} else{unlist(strsplit(synonym[[1]],split="\t"))}
	
	description = as.character(mydata$description)[[i]]
	organism = as.character(mydata$organism)[[i]]
	location = as.character(mydata$location)[[i]]
  
	xref = strsplit(str_sub(as.character(mydata$xref)[[i]],3,-3),split="\",\"")
	xref = if(length(xref[[1]])==0){c("")} else{unlist(strsplit(xref[[1]],split="\t"))}
	
	createNode(graph, "Gene", GID=GID, name=name, synonym=synonym, 
	           description=description, organism=organism, location=location, xref=xref)
}

####Creating pathway##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

mydata[is.na(mydata)]<-""

for (i in 1:nrow(mydata)){
	cat("... Creating pathway node: ",i,"\n")
	GID = as.character(mydata$GID)[[i]]
	name = as.character(mydata$name)[[i]]
	synonym = strsplit(str_sub(as.character(mydata$synonym)[[i]],3,-3),split="\",\"")
	synonym  = if(length(synonym[[1]])==0){c("")} else{unlist(strsplit(synonym[[1]],split="\t"))}
	
	description = as.character(mydata$description)[[i]]
	organism = as.character(mydata$organism)[[i]]
	manCurate = as.character(mydata$manCurate)[[i]]
  
	describeIn = strsplit(str_sub(as.character(mydata$describeIn)[[i]],3,-3),split="\",\"")
	describeIn = if(length(describeIn[[1]])==0){c("")} else{unlist(strsplit(describeIn[[1]],split="\t"))}
	
	xref = strsplit(str_sub(as.character(mydata$xref)[[i]],3,-3),split="\",\"")
	xref = if(length(xref[[1]])==0){c("")} else{unlist(strsplit(xref[[1]],split="\t"))}
	
	createNode(graph, "Pathway", GID=GID, name=name, synonym=synonym, 
	           description=description, organism=organism, manCurate=manCurate, 
             describeIn=describeIn, xref=xref)
}

####Creating TRANSFORM##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

for (i in 1:nrow(mydata)){
  cat("... Creating TRANSFORM relation: ",i,"\n")
  met = as.character(mydata$metabolite)[[i]]
  rx = as.character(mydata$reaction)[[i]]
  to = getSingleNode(graph, "MATCH (x:Metabolite) WHERE x.GID = {gid} RETURN x", gid=met)
  from = getSingleNode(graph, "MATCH (x:Reaction) WHERE x.GID = {gid} RETURN x", gid=rx)
  reversible = as.character(mydata$reversible)[[i]]
  ##no need type
  ##type = strsplit(str_sub(as.character(mydata$type)[[i]],3,-3),split="\",\"")
  ##type = if(length(type[[1]])==0){c("")} else{unlist(strsplit(type[[1]],split="\t"))}
  
  createRel(from, "TRANSFORM", to, reversible=reversible)
}

####Creating PRODUCE##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

for (i in 1:nrow(mydata)){
  cat("... Creating PRODUCE relation: ",i,"\n")
  met = as.character(mydata$metabolite)[[i]]
  rx = as.character(mydata$reaction)[[i]]
  to = getSingleNode(graph, "MATCH (x:Metabolite) WHERE x.GID = {gid} RETURN x", gid=met)
  from = getSingleNode(graph, "MATCH (x:Reaction) WHERE x.GID = {gid} RETURN x", gid=rx)
  reversible = as.character(mydata$reversible)[[i]]
  ##no need type
  ##type = strsplit(str_sub(as.character(mydata$type)[[i]],3,-3),split="\",\"")
  ##type = if(length(type[[1]])==0){c("")} else{unlist(strsplit(type[[1]],split="\t"))}
  
  createRel(from, "PRODUCE", to, reversible=reversible)
}

####Creating HAS##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

for (i in 1:nrow(mydata)){
  cat("... Creating HAS relation: ",i,"\n")
  ptw = as.character(mydata$pathway)[[i]]
  rx = as.character(mydata$reaction)[[i]]
  to = getSingleNode(graph, "MATCH (x:Reaction) WHERE x.GID = {gid} RETURN x", gid=rx)
  from = getSingleNode(graph, "MATCH (x:Pathway) WHERE x.GID = {gid} RETURN x", gid=ptw)
  
  createRel(from, "HAS", to)
}

####Creating CATALYZE##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)
mydata$km[is.na(mydata$km)]<-0

for (i in 1:nrow(mydata)){
  cat("... Creating CATALYZE relation: ",i,"\n")
  pr = as.character(mydata$protein)[[i]]
  rx = as.character(mydata$reaction)[[i]]
  km = mydata$km[[i]]
  to = getSingleNode(graph, "MATCH (x:Reaction) WHERE x.GID = {gid} RETURN x", gid=rx)
  from = getSingleNode(graph, "MATCH (x:Protein) WHERE x.GID = {gid} RETURN x", gid=pr)

  createRel(from, "CATALYZE", to, km=km)
}

####Creating ENCODE##
file = "[FILEPATH]"
mydata = read.csv(file=file, header=T)

for (i in 1:nrow(mydata)){
  cat("... Creating ENCODE relation: ",i,"\n")
  pr = as.character(mydata$protein)[[i]]
  ge = as.character(mydata$gene)[[i]]
  to = getSingleNode(graph, "MATCH (x:Protein) WHERE x.GID = {gid} RETURN x", gid=pr)
  from = getSingleNode(graph, "MATCH (x:Gene) WHERE x.GID = {gid} RETURN x", gid=ge)
  
  createRel(from, "ENCODE", to)
}