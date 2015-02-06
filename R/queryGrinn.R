#' \code{queryGrinn} queries through Neo4j REST API.
#'@import RCurl jsonlite igraph opencpu
#'@export
queryGrinn <- function(){}

#'@global list of cypher to query the path of different relationship types
relationList <- list(
  biochem = "MATCH (ptw:Pathway{organism:species})-[:HAS]->(rx:Reaction) WITH rx MATCH left<-[:TRANSFORM]-(rx)-[:PRODUCE]->right WHERE ANY(y IN keyword WHERE lower(y) = lower(left.GID)) AND ANY(y IN keyword WHERE lower(y) = lower(right.GID)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID",
  enzcatalyze = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Protein {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID",
  encgene = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Gene {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID",
  biopathway = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:Pathway {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, source.GID, source.name ORDER by source.GID",
  pairwise = "UNWIND keyword AS x WITH x MATCH ptw = allShortestPaths((target:Metabolite)<-[*]-(source:label {organism:species})) WHERE lower(target.GID) = lower(x) RETURN DISTINCT target.GID, target.name, target.xref, source.GID, source.name, source.xref"
)

#'@author K Wanich; kwanich@ucdavis.edu
#'@method perform HTTP request using RCurl
#'@param string of cypher
#'@return data frame of given query
#'@note set url = db location
curlRequestCypher <- function(querystring){
  h = RCurl::basicTextGatherer()
  RCurl::curlPerform(url="localhost:7474/db/data/cypher",
                     postfields=paste('query',RCurl::curlEscape(querystring), sep='='),
                     writefunction = h$update,
                     verbose = FALSE
  ) 
  data <- jsonlite::fromJSON(h$value()) 
  result <- data$data #return as data.frame
}

#'@author K Wanich; kwanich@ucdavis.edu
#'@method perform HTTP request using RCurl
#'@param string of URI
#'@return list of given URI
curlRequestUrl <- function(url) {
  h = RCurl::basicTextGatherer()
  RCurl::curlPerform(url=url,
                     writefunction = h$update,
                     verbose = FALSE
  ) 
  result <- jsonlite::fromJSON(h$value()) #return as list
}

#'@author K Wanich; kwanich@ucdavis.edu
#'@method execute query for a node
#'@param string of cypher
#'@return list of node information including its repationships to other nodes
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
  result <- curlRequestUrl(url)
  if(direction=="out"){#outgoing rel
    rellist.names = c("end_name","end_GID","out_type","out_data") #set names
    rel <- sapply(rellist.names,function(x) NULL) #set default relationship value to null
    if(length(result)>0){
      for(i in 1:nrow(result)){
        dirurl = result$end[[i]]
        rel$end_name = c(rel$end_name, curlRequestUrl(dirurl)$data$name)
        rel$end_GID = c(rel$end_GID, curlRequestUrl(dirurl)$data$GID)
        rel$out_type = c(rel$out_type, result$type[[i]])
        
        if(ncol(result$data)>1){
          ind = which(lapply(result$data[i,], function(x) is.null(unlist(x))) == FALSE)
          rel$out_data = c(rel$out_data, jsonlite::toJSON(result$data[i,ind]))
        }
        else if(ncol(result$data)==0){
          rel$out_data = c(rel$out_data, jsonlite::toJSON(""))#no property default value
        }
        else{
          tmp = list(result$data[i,1])
          names(tmp) = names(result$data)
          rel$out_data = c(rel$out_data, jsonlite::toJSON(tmp))
        }
      }
    }
  }#out
  else{#incoming rel
    rellist.names = c("start_name","start_GID","in_type","in_data") #set names
    rel <- sapply(rellist.names,function(x) NULL) #set default relationship value to null
    if(length(result)>0){
      for(i in 1:nrow(result)){
        dirurl = result$start[[i]]
        rel$start_name = c(rel$start_name, curlRequestUrl(dirurl)$data$name)
        rel$start_GID = c(rel$start_GID, curlRequestUrl(dirurl)$data$GID)
        rel$in_type = c(rel$in_type, result$type[[i]])
  
        if(ncol(result$data)>1){
          ind = which(lapply(result$data[i,], function(x) is.null(unlist(x))) == FALSE)
          rel$in_data = c(rel$in_data, jsonlite::toJSON(result$data[i,ind]))
        }
        else if(ncol(result$data)==0){
          rel$in_data = c(rel$in_data, jsonlite::toJSON(""))#no property default value
        }
        else{
          tmp = list(result$data[i,1])
          names(tmp) = names(result$data)
          rel$in_data = c(rel$in_data, jsonlite::toJSON(tmp))
        }
      }
    }
  }#in
  rel #output relationships info
}

#'@author K Wanich; kwanich@ucdavis.edu
#'@method generate a combined network
#'@description call internal function createBiochemNetwork, connectNodes to get different types of networks and merge them
#'Grinn v1 provides 4 types of network: biochemical reaction, enzyme catalysis, encoding gene and metabolic pathway
#'call internal function createCyNetwork to generate the combined network in json format
#'@param string of keywords, string of organism
#'@example txtInput = "[\"C00152\", \"C01637\", \"C00135\", \"C01643\", \"C00079\", \"C03402\", \"C02988\", \"C03511\", \"C02163\", \"C02282\"]"
#'@example organism = "\"Homo sapiens\""
#'@return list of nodes and edges of the network, each in json format 
integrateNetwork <- function(txtInput, organism){
  #generate 4 types of networks, add more types if db structure changed
  bcnw = createBiochemNetwork(txtInput,organism)
  enznw = connectNodes(txtInput,organism,"enzcatalyze")
  genenw = connectNodes(txtInput,organism,"encgene")
  ptwnw = connectNodes(txtInput,organism,"biopathway")
  
  print("Merging network...")    
  gu = igraph::graph.union(bcnw,enznw,genenw,ptwnw) #merge networks
  
  network <- tryCatch({
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
    gudf$edges = gudf$edges[,c("from","to","biochem","enzcatalyze","encgene","biopathway","reltype")] #more intermediates if db structure changed
    colnames(gudf$edges)[1:2] = c("source","target")
    
    network <- createCyNetwork(gudf$nodes, gudf$edges) #return integrated network
  }, error = function(err) {
    #on error return empty network
    network <- list(nodes="", edges="")   
  }) # END tryCatch
}

#'@method internal function to create metabolite - RX - metabolite network
createBiochemNetwork <- function(txtInput, organism){
  #construct query string
  querystring = relationList$biochem
  querystring = gsub("keyword", txtInput, querystring)
  querystring = gsub("species", organism, querystring)

# querystring = "MATCH (ptw:Pathway{organism:\"Homo sapiens\"})-[:HAS]->(rx:Reaction) WITH rx MATCH left<-[:TRANSFORM]-(rx)-[:PRODUCE]->right WHERE ANY(y IN [\"C00024\",\"C00136\",\"C00010\",\"C05269\"] WHERE lower(y) = lower(left.GID)) AND ANY(y IN [\"C00024\",\"C00136\",\"C00010\",\"C05269\"] WHERE lower(y) = lower(right.GID)) RETURN left.GID, left.name, right.GID, right.name, rx.GID, rx.name ORDER BY left.GID"
  print("Querying bioChemNetwork...")
print(querystring)
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
      igraph::E(nw)$biochem = lapply(igraph::E(nw)$biochem,unique) #remove duplicate
      combineAttb = function(x){ paste0(unlist(unlist(x)),collapse="||") } #combine reactions, format: GID|name||GID|name, each reaction seperated by ||
      igraph::E(nw)$biochem = lapply(igraph::E(nw)$biochem,combineAttb) 
      igraph::E(nw)$reltype = "biochem" #set relationship type as edge attribute
print(length(igraph::E(nw)))
      result <- nw #return biochem network
    }, error = function(err) {
      #on error return empty network
      result <- igraph::graph.empty(n=0, directed=FALSE)
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
  print("Querying network...") 
print(querystring)
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
      nw = igraph::set.edge.attribute(nw, name= reltype, value= apply(df,1,findIntermediateNode)) #set intermediates as edge attribute
      igraph::E(nw)$reltype = reltype #set relationship type as edge attribute
print(length(igraph::E(nw)))
      result <- nw #return connected nodes
    }, error = function(err) {
    #on error return empty network
      result <- igraph::graph.empty(n=0, directed=FALSE) 
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
    return(NULL)
  }
  # There must be edges and edgeData must have at least source and target columns
  if(nrow(edgeData) == 0 || !(all(c("source", "target") %in% names(edgeData)))) {
    return(NULL)
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

#'@author K Wanich; kwanich@ucdavis.edu
#'@method create a map of metabolites to the specified node type
#'@param string of keywords, string of organism, string of node type
#'@return list of mapped nodes and node attributes in json format
#'@note Grinn v1 provides 3 node types for mapping: Protein, Gene, Pathway
mapToNode <- function(txtInput, organism, label){
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
        pair = rbind(pair,cbind(unlist(data[[j]][1]),unlist(data[[j]][4]))) #from-to
        attb = rbind(attb,cbind(unlist(data[[j]][1]),unlist(data[[j]][2]),(data[[j]][3]),"Metabolite")) #attribute from node
        attb = rbind(attb,cbind(unlist(data[[j]][4]),unlist(data[[j]][5]),(data[[j]][6]),label[i])) #attribute to node
      }
      cat("Returning... ",length(data)," pairs to ",label[i], "nodes\n")
    }
  
  }
  
  if(length(pair)>0){
    cat("\n**..TOTAL of ",nrow(unique(pair))," unique pairs found..**\n")
    colnames(pair) = c("from","to")
    colnames(attb) = c("GID","name","xref","type")
    pair <- jsonlite::toJSON(unique(pair))
    attb <- jsonlite::toJSON(unique(attb))
  }
  else{# if no mapped node found
    print("Returning no data...")
  }
  result <- list(pair,attb)
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
