#' \code{queryGrinn} queries through Neo4j REST API.
#'@import RCurl jsonlite igraph opencpu
#'@export
queryGrinn <- function(){}

#'@global list of cypher to query the path of different relationship types
relationList <- list(
  biochem = "MATCH (ptw:Pathway{organism:species})-[:HAS]->(rx:Reaction) WITH rx MATCH left<-[:TRANSFORM]-(rx)-[:PRODUCE]->right WHERE ANY(y IN keyword WHERE lower(y) = lower(left.GID)) AND ANY(y IN keyword WHERE lower(y) = lower(right.GID)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID",
  enzcatalyze = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Protein {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID",
  encgene = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Gene {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID",
  pathway = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Pathway {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID",
  pairwise = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:label {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, target.xref, source.GID, source.name, source.xref"
)

#'@export
#'@author K Wanich; kwanich@ucdavis.edu
#'@method perform HTTP request using RCurl
#'@param string of cypher
#'@return data frame of quried result
#'@note set url = db location
curlRequestCypher <- function(querystring){
  h = RCurl::basicTextGatherer()
  tryCatch({
    RCurl::curlPerform(url="localhost:7474/db/data/cypher",
                       postfields=paste('query',RCurl::curlEscape(querystring), sep='='),
                       writefunction = h$update,
                       verbose = FALSE
    ) 
    data <- jsonlite::fromJSON(h$value()) 
    result <- data$data #return as data.frame
    }, error = function(err) {
      print(err)
      result <- NULL #return NULL if not found
  }) # END tryCatch
}

#'@export
#'@author K Wanich; kwanich@ucdavis.edu
#'@method perform HTTP request using RCurl
#'@param string of URI
#'@return data frame of quried result
curlRequestUrlToDF <- function(url) {
  h = RCurl::basicTextGatherer()
  tryCatch({
    RCurl::curlPerform(url=url,
                       writefunction = h$update,
                       verbose = FALSE
    ) 
    result <- jsonlite::fromJSON(h$value()) #return as data.frame
    }, error = function(err) {
      print(err)
      result <- NULL #return NULL if not found
  }) # END tryCatch
}

#'@export
#'@author K Wanich; kwanich@ucdavis.edu
#'@method perform HTTP request using RCurl
#'@param string of URI
#'@return list of quried result
curlRequestUrlToList <- function(url) {
  h = RCurl::basicTextGatherer()
  tryCatch({
    RCurl::curlPerform(url=url,
                       writefunction = h$update,
                       verbose = FALSE
    ) 
    result <- jsonlite::fromJSON(h$value(), simplifyVector = FALSE) #return as list
    }, error = function(err) {
      print(err)
      result <- NULL #return NULL if not found
  }) # END tryCatch
}

#'@export
#'@author K Wanich; kwanich@ucdavis.edu
#'@method execute query for nodes
#'@description query for nodes from a cypher input,
#'quried results include node information and relationship information,
#'relationship information includes ID and name of adjacent nodes, name and additional attribute of relationships
#'@param string of cypher
#'@return data.frame of node information and repationship information
##need to improve query by cypher
queryNode <- function(querystring, querytype="node") {
  result <- curlRequestCypher(querystring)
  if(length(result) != 0){
    if(querytype =="node"){
      getNodeInfo(result) #call internal function
    }
  }
  else{#for cypher input
    result
  }
}

#'@method internal function to get complete node information
getNodeInfo <- function(data){
  for(i in 1:length(data)){
    outurl = data[[i]]$outgoing_relationships
    data[[i]] = RCurl::merge.list(data[[i]], lapply(getRelInfo(outurl,"out"),FUN=list))
    inurl = data[[i]]$incoming_relationships
    data[[i]] = RCurl::merge.list(data[[i]], lapply(getRelInfo(inurl,"in"),FUN=list))
  }
  result <- data #output node info + its relationships info
}

#'@method internal function to get all relationship info of a node
getRelInfo <- function(url,direction){
  result <- curlRequestUrlToList(url)
  if(direction=="out"){#outgoing rel
    rellist.names = c("end_name","end_GID","out_type","out_data") #set names
    rel <- sapply(rellist.names,function(x) NULL) #set default relationship value to null
    if(length(result)>0){
      for(i in 1:length(result)){
        dirurl = result[[i]]$end
        rel$end_name = c(rel$end_name, curlRequestUrlToList(dirurl)$data$name)
        rel$end_GID = c(rel$end_GID, curlRequestUrlToList(dirurl)$data$GID)
        rel$out_type = c(rel$out_type, result[[i]]$type)
        
        if(length(result[[i]]$data)>0){
          rel$out_data = c(rel$out_data, jsonlite::toJSON(result[[i]]$data))
        }
        else{
          rel$out_data = c(rel$out_data, jsonlite::toJSON(""))#no property default value
        }
      }
    }
  }#out
  else{#incoming rel
    rellist.names = c("start_name","start_GID","in_type","in_data") #set names
    rel <- sapply(rellist.names,function(x) NULL) #set default relationship value to null
    if(length(result)>0){
      for(i in 1:length(result)){
        dirurl = result[[i]]$start
        rel$start_name = c(rel$start_name, curlRequestUrlToList(dirurl)$data$name)
        rel$start_GID = c(rel$start_GID, curlRequestUrlToList(dirurl)$data$GID)
        rel$in_type = c(rel$in_type, result[[i]]$type)
  
        if(length(result[[i]]$data)>0){
          rel$in_data = c(rel$in_data, jsonlite::toJSON(result[[i]]$data))
        }
        else{
          rel$in_data = c(rel$in_data, jsonlite::toJSON(""))#no property default value
        }
      }
    }
  }#in
  rel #output relationships info
}

#'@export
#'@author K Wanich; kwanich@ucdavis.edu
#'@method generate a combined network
#'@description call internal function createBiochemNetwork, connectNodes to get different types of networks and merge them,
#'call internal function createCyNetwork to generate the combined network in json format recognized by cytoscapeJS,
#'select return type of result: json of nodes and edges or json of cytoscapreJS or all types
#'@param string of keywords, string of organism, string of return type
#'@example txtInput = "[\"C00152\", \"C01637\", \"C00135\", \"C01643\", \"C00079\", \"C03402\", \"C02988\", \"C03511\", \"C02163\", \"C02282\"]"
#'@example organism = "\"Homo sapiens\""
#'@example returnAs = "all" or "json" or "cytoscape", default is json of nodes and edges
#'@return list of nodes and edges in json or cytoscapeJS format
#'@note Grinn v1 provides 4 types of network: biochemical reaction, enzyme catalysis, encoding gene and metabolic pathway
integrateNetwork <- function(txtInput, organism, returnAs="json"){
  #generate 4 types of networks, add more types if db structure changed
  bcnw = createBiochemNetwork(txtInput,organism)
  enznw = connectNodes(txtInput,organism,"enzcatalyze")
  genenw = connectNodes(txtInput,organism,"encgene")
  ptwnw = connectNodes(txtInput,organism,"pathway")
   
  gu = igraph::graph.union(bcnw,enznw,genenw,ptwnw) #merge networks
  
  yield <- tryCatch({
    gudf = igraph::get.data.frame(gu,"both") #get data frame of generated network
    #format node name
    rownames(gudf$vertices) = NULL
    names(gudf) = c("nodes","edges")
    for(i in 1:nrow(gudf$nodes)){
      ind = which(!is.na(gudf$nodes[i,]))[1]
      gudf$nodes$sysname[i] = gudf$nodes[i,ind]
    }
    gudf$nodes = data.frame(id=gudf$nodes$name, name=gudf$nodes$sysname, href=paste0("node.html?GID=",gudf$nodes$name))
    
    #format edge attributes
    gudf$edges$reltype = paste(sep=",",gudf$edges$reltype_1,gudf$edges$reltype_2,
                               gudf$edges$reltype_3,gudf$edges$reltype_4) #more reltype if db structure changed
    gudf$edges = gudf$edges[,c("from","to","biochem","enzcatalyze","encgene","pathway","reltype")] #more intermediates if db structure changed
    colnames(gudf$edges)[1:2] = c("source","target")

    #network in cytoscapeJS format
    cynw = createCyNetwork(gudf$nodes, gudf$edges) 
    #network in json format
    pair = jsonlite::toJSON(gudf$edges)
    attb = jsonlite::toJSON(gudf$nodes)
    switch(returnAs,
           all = list(nodelist=attb, edgelist=pair, cynetwork=cynw),
           json = list(nodelist=attb, edgelist=pair),
           cytoscape = cynw,
           stop("data type not included"))
  }, error = function(err) {
    #on error return empty network
    list(nodes="", edges="")   
  }) # END tryCatch

}

#'@method internal function to create metabolite - RX - metabolite network
createBiochemNetwork <- function(txtInput, organism){
  #construct query string
  querystring = relationList$biochem
  querystring = gsub("keyword", txtInput, querystring)
  querystring = gsub("species", organism, querystring)

# querystring = "MATCH (ptw:Pathway{organism:\"Homo sapiens\"})-[:HAS]->(rx:Reaction) WITH rx MATCH left<-[:TRANSFORM]-(rx)-[:PRODUCE]->right WHERE ANY(y IN [\"C00024\",\"C00136\",\"C00010\",\"C05269\"] WHERE lower(y) = lower(left.GID)) AND ANY(y IN [\"C00024\",\"C00136\",\"C00010\",\"C05269\"] WHERE lower(y) = lower(right.GID)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID"
#print(querystring)
  data <- curlRequestCypher(querystring) #table of left.GID, left.name, right.GID, right.name, rx.GID, rx.name

  result <- tryCatch({   
      #create a temp graph from query result
      g <- igraph::graph.edgelist(matrix(data[,c(1,3)], ncol=2), directed = F)
      #set node name
      nodes = unique(rbind(data[,1:2], data[,3:4]))
      findNodeName = function(x){ which(nodes == x) } #internal function to find node name
      ind = lapply(igraph::V(g)$name,findNodeName)
      igraph::V(g)$sysname = nodes[unlist(ind),2]
      #set edge attribute
      igraph::E(g)$biochem = paste0(data[,5],'|',data[,6])
      #create network
      nw = igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list)
      if(igraph::ecount(nw)>0){# if graph is not empty
        igraph::E(nw)$biochem = lapply(igraph::E(nw)$biochem,unique) #remove duplicate
        combineAttb = function(x){ paste0(unlist(unlist(x)),collapse="||") } #combine reactions, format: GID|name||GID|name, each reaction seperated by ||
        igraph::E(nw)$biochem = lapply(igraph::E(nw)$biochem,combineAttb)
        cat("Returning... ",length(igraph::E(nw))," biochem relationships\n")
      }
      else{# if graph is empty
        igraph::E(nw)$biochem = ""
      }
      igraph::E(nw)$reltype = "biochem" #set relationship type as edge attribute
      result <- nw #return biochem network
    }, error = function(err) {
      #on error return empty network
      result <- igraph::graph.empty(n=0, directed=FALSE)
      igraph::V(result)$name = ""
      igraph::V(result)$sysname = ""
      result <- igraph::set.edge.attribute(result, name= "biochem", value= "")
      result <- igraph::set.edge.attribute(result, name= "reltype", value= "biochem")
    }) # END tryCatch
#..code to create Cy network..#
#   df = igraph::get.data.frame(nw,"both") #get data frame of generated network
#   #format name
#   rownames(df$vertices) = NULL
#   names(df) = c("nodes","edges")
#   names(df$nodes) = c("id","name")
#   colnames(df$edges)[1:2] = c("source","target")
#   network <- createCyNetwork(df$nodes, df$edges) 
}

#'@method internal function to create metabolite - PRT|GN|PTW - metabolite network
connectNodes <- function(txtInput, organism, reltype){
  #construct query string
  querystring = relationList[reltype]
  querystring = gsub("keyword", txtInput, querystring)
  querystring = gsub("species", organism, querystring)
  
# querystring = "UNWIND [\"met5\",\"met1\",\"met8\"] AS x WITH x MATCH (source:Protein{organism:\"Homo sapiens\"}), (target:Metabolite) WHERE lower(target.GID) = lower(x) 
# WITH target, source MATCH ptw = target<-[:TRANSFORM|PRODUCE]-()<-[:CATALYZE]-source 
#   RETURN target.GID, target.name, source.GID, source.name ORDER by source.GID"
#print(querystring)
  data <- curlRequestCypher(querystring)
  result <- tryCatch({ 
      if(typeof(data)=="list"){data = t(data.frame(data))} #check type of data, change to dataframe of characters if need
      data <- data[!duplicated(data), ] #table of target.GID, target.name, source.GID, source.name
      #create a temp graph from query result
      g <- igraph::graph.edgelist(matrix(data[,c(1,3)], ncol=2), directed = F)
      toNodes = unique(data[,1])

      #find if a node connects to a node from cocitation, connect nodes and create a network
      cocite = igraph::cocitation(g)[toNodes,toNodes] #get adj. matrix of nodes from cocitation
      nw <- igraph::graph.adjacency(cocite, mode = "undirected")
      nw <- igraph::simplify(nw, remove.multiple = T, remove.loops = T)
      el = igraph::get.edgelist(nw) #remove unconnected nodes
      nw <- igraph::graph.edgelist(el, directed = F)
      #set node name
      findNodeName = function(x){ which(data == x)[1] } #internal function to find node name
      ind = lapply(igraph::V(nw)$name,findNodeName)
      igraph::V(nw)$sysname = data[unlist(ind),2]
      #get intermediates
      gr = by(data,data[,3],data.frame)
      findRow = function(x){ which(nrow(x)>1) } #internal function to screen out path with no-intermediate node
      ind = lapply(gr, findRow)
      intmd = gr[names(unlist(ind))]
      findIntermediateNode = function(x){ #internal function to collect intermediates
        lsImd = ""
        for (i in 1:length(intmd)){
          if(all(x %in% intmd[i][[1]]$V1)){
            #format: GID|name||GID|name, each intermediate seperated by ||
            tmp = paste0(as.character(intmd[i][[1]]$V3[1]),'|',as.character(intmd[i][[1]]$V4[1]))
            lsImd = paste0(lsImd,tmp,sep="||")
    #lsImd = c(lsImd, paste0(as.character(intmd[i][[1]]$V3[1]),'|',as.character(intmd[i][[1]]$V4[1])))
          }  
        }
        lsImd = substr(lsImd, 1, nchar(lsImd)-2)
      }
    
      df = igraph::get.data.frame(nw) #get data frame of generated network
      nw = igraph::set.edge.attribute(nw, name= reltype, value= apply(df,1,findIntermediateNode)) #set intermediates as edge attribute, throw error if unfound
      igraph::E(nw)$reltype = reltype #set relationship type as edge attribute
      cat("Returning... ",length(igraph::E(nw))," ",reltype," relationships\n")
      result <- nw #return connected nodes
    }, error = function(err) {
    #on error return empty network
      result <- igraph::graph.empty(n=0, directed=FALSE)
      igraph::V(result)$name = ""
      igraph::V(result)$sysname = ""
      result <- igraph::set.edge.attribute(result, name= reltype, value= "")
      result <- igraph::set.edge.attribute(result, name= "reltype", value= reltype)
    }) # END tryCatch
#..code to create Cy network..#  
#   df = igraph::get.data.frame(nw,"both") #get data frame of generated network
#   #format name
#   rownames(df$vertices) = NULL
#   names(df) = c("nodes","edges")
#   names(df$nodes) = c("id","name")
#   colnames(df$edges)[1:2] = c("source","target")
#   network <- createCyNetwork(df$nodes, df$edges)
} 

#' @method internal function to format result network to cytoscapeJS input format
#' @references https://github.com/cytoscape/r-cytoscape.js/blob/master/cytoscapeJsSimpleNetwork.R
#' @description modify the R function from r-cytoscape.js (https://github.com/cytoscape/r-cytoscape.js/blob/master/cytoscapeJsSimpleNetwork.R), 
#' to return JSON output for generating cytoscape.js network 
createCyNetwork <- function(nodeData, edgeData, 
                                   nodeColor="#888888", nodeShape="ellipse") {  
  
  # There must be nodes and nodeData must have at least id and name columns
  if(nrow(nodeData) == 0 || !(all(c("id", "name") %in% names(nodeData)))) {
    return(list(nodes="", edges=""))
  }
  # There must be edges and edgeData must have at least source and target columns
  if(nrow(edgeData) == 0 || !(all(c("source", "target") %in% names(edgeData)))) {
    return(list(nodes="", edges=""))
  }

  # NODES
  ## Add color/shape columns if not present
  if(!("color" %in% colnames(nodeData))) {
    nodeData$color <- rep(nodeColor, nrow(nodeData))
  }
  if(!("shape" %in% colnames(nodeData))) {
    nodeData$shape <- rep(nodeShape, nrow(nodeData))
  }
  nodeEntries <- NULL
  for(i in 1:nrow(nodeData)) {   
    tmpEntries <- NULL
    for(col in colnames(nodeData)) {
      tmp2 <- paste0("\"", col, "\":\"", nodeData[i, col], "\"")
      tmpEntries <- c(tmpEntries, tmp2)
    }
    tmpEntries <- paste(tmpEntries, collapse=", ")
    tmp <- paste0("{ \"data\": { ", tmpEntries, "} }")  
    nodeEntries <- c(nodeEntries, tmp)
  }
  nodeEntries <- paste(nodeEntries, collapse=", ")
  
  # EDGES 
  edgeEntries <- NULL
  for(i in 1:nrow(edgeData)) {   
    tmpEntries <- NULL
    for(col in colnames(edgeData)) {
      tmp2 <- paste0("\"", col, "\":\"", edgeData[i, col], "\"")
      tmpEntries <- c(tmpEntries, tmp2)
    }
    tmpEntries <- paste(tmpEntries, collapse=", ")
    tmp <- paste0("{ \"data\": { ", tmpEntries, "} }")
    edgeEntries <- c(edgeEntries, tmp)
  }
  edgeEntries <- paste(edgeEntries, collapse=", ")
  network <- list(nodes=nodeEntries, edges=edgeEntries)

  #print(network)
  return(network)
}

#'@export
#'@author K Wanich; kwanich@ucdavis.edu
#'@method create a map of metabolites to the specified node type
#'@description connect metabolites to other node types,
#'call internal function createCyNetwork to generate the combined network in json format recognized by cytoscapeJS,
#'select return type of result: json of nodes and edges or json of cytoscapreJS or all types
#'@param string of keywords, string of organism, list of node type, string of return type
#'@example label = c("Protein","Gene","Pathway"), default is all 3 node types
#'@example returnAs = "all" or "json" or "cytoscape", default is json of nodes and edges
#'@return list of nodes and edges in json or cytoscapeJS format
#'@note Grinn v1 provides 3 node types for mapping: Protein, Gene, Pathway
pairToNode <- function(txtInput, organism, label=c("Protein","Gene","Pathway"), returnAs="json"){
  pair = data.frame() #list of mapped nodes
  attb = data.frame() #list of node attributes

  #for each label
  for(i in 1:length(label)){
    #construct query string
    querystring = relationList$pairwise
    querystring = gsub("keyword", txtInput, querystring)
    querystring = gsub("species", organism, querystring)
    querystring = gsub("label", label[i], querystring)
    data <- curlRequestCypher(querystring)
    
    if(length(data)>0){
      #collect list of mapped nodes and node attributes
      for(j in 1:length(data)){
        pair = rbind(pair,cbind(unlist(data[[j]][1]),unlist(data[[j]][4]),paste("PAIR_TO",label[i],sep = "_"))) #from-to-reltype
        attb = rbind(attb,cbind(unlist(data[[j]][1]),unlist(data[[j]][2]),(data[[j]][3]),"Metabolite")) #attribute from node
        attb = rbind(attb,cbind(unlist(data[[j]][4]),unlist(data[[j]][5]),(data[[j]][6]),label[i])) #attribute to node
      }
      cat("Returning... ",length(data)," pairs to ",label[i], "nodes\n")
    }
  
  }
  
  if(length(pair)>0){
    cat("**..TOTAL of ",nrow(unique(pair))," unique pairs found..**\n")
    colnames(pair) = c("source","target","reltype")
    colnames(attb) = c("GID","nodename","xref","nodetype")
    pair <- unique(pair)
    attb <- unique(attb)
    
    #create cytoscape network
    lsXref = unlist(lapply(attb[,3],FUN=function(x){paste(unlist(x),collapse = "||")}))    
    g <- igraph::graph.edgelist(as.matrix(pair[,-3]), directed = F)  
    g <- set.edge.attribute(g,name="reltype",value=as.character(pair[,3]))    
    g <- set.vertex.attribute(g, name="id", value=as.character(attb[,1]))
    g <- set.vertex.attribute(g, name="nodename", value=as.character(attb[,2]))
    g <- set.vertex.attribute(g, name="href", value=paste0("node.html?GID=",as.character(attb[,1])))
    g <- set.vertex.attribute(g, name="xref", value=lsXref)
    g <- set.vertex.attribute(g, name="nodetype", value=as.character(attb[,4]))
    gudf = igraph::get.data.frame(g,"both")
    rownames(gudf$vertices) = NULL
    names(gudf) = c("nodes","edges")
    colnames(gudf$edges)[1:2] = c("source","target")
    network = createCyNetwork(gudf$nodes, gudf$edges) #return integrated network
    
    #network in json format
    pair = jsonlite::toJSON(pair)
    attb = jsonlite::toJSON(attb)
  }
  else{# if no mapped node found
    print("Returning no data...")
    network = list(nodes="", edges="")
  }
  yield <- switch(returnAs,
             all = list(nodelist=attb, edgelist=pair, cynetwork=network),
             json = list(nodelist=attb, edgelist=pair),
             cytoscape = list(cynetwork=network),
             stop("data type not included"))
}

#'format txtInput from file
#'keywords seperated by new line
## txtInput = "[\"\",\"\"]"
formatTextInput <- function(){
  myinput = read.delim("tmp.txt", header=F, stringsAsFactors = F)
  txt = toString(unique(myinput))
  result = paste0("[",substr(txt,3,nchar(txt)-1),"]")
  result <- gsub("\n","",result)
}
#'format organism string
## organism = "\"\""
formatOrganismString <- function(org){
  result <- paste0("\"",org,"\"")
}
#write(unlist(intnw), "nw.txt", ncolumns = 1) #write network file
