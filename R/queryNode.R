#' \code{queryNode} query for nodes
#'@description query for nodes from a Cypher input. Quried results include node information and relationship information.
#'Relationship information includes ID and name of adjacent nodes, name and additional attributes of relationships
#'@usage queryNode(querystring)
#'@param querystring string of Cypher
#'@return list of data.frame containing node information and repationship information. Return empty list if found nothing.
#'@author Kwanjeera Wanichthanarak \email{kwanich@@ucdavis.edu}
#'@export
#'@seealso \code{\link{curlRequestCypher}}
#'@examples
#'# Query metabolites by name 
#'querystring = "UNWIND ['coa','co2','Malonyl-CoA'] AS x WITH x MATCH (n:Metabolite) WHERE lower(n.name) = lower(x) RETURN DISTINCT n" 
#'result = queryNode(querystring)

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