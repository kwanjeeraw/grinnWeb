% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/pairToNode.R
\name{pairToNode}
\alias{pairToNode}
\title{\code{pairToNode} generate an integrated network}
\usage{
pairToNode(txtInput, organism, label=c("Protein","Gene","Pathway"), returnAs="json", searchBy="grinn")
}
\arguments{
\item{txtInput}{string containing metabolite ids. Metabolite ids are from the specified database, see \code{searchBy}. Default is grinn id e.g. C00065.
Search using grinn ids, txtInput is in the following format: txtInput = "['id'1, 'id2']".
Search using other database ids, txtInput must be in the following format: txtInput = "['databaseName:id'1, 'databaseName:id2']".}

\item{organism}{string of species in the following format: organism = "'species'"}

\item{label}{character vector of node types to connect to metabolites, default is all three node types: c("Protein","Gene","Pathway")}

\item{returnAs}{string of output type. Specify the type of the returned network.
It can be "all", "json" or "cytoscape", but default is "json". "cytoscape" is the format used in Cytoscape.js}

\item{searchBy}{string of database name. Specify the database of which the ids used as the txtInput, see \code{txtInput}.
It can be "grinn" or "InChI" or "KEGG" or "PubChem" or "ChEBI". Default is "grinn".}
}
\value{
list of nodes and edges encapsulated in json format. Return empty list if found nothing.
}
\description{
from the list of metabolites, find and connect these metabolites to three type of nodes
including Protein, Gene and Pathway.
}
\examples{
# Query metabolites by PubChem ids and connect to Protein, Gene and Pathway
txtInput <- "['PubChem:3326','PubChem:3436','PubChem:7656']"
organism <- "'Homo sapiens'"
result <- pairToNode(txtInput, organism, returnAs="all", searchBy="PubChem")
# Query metabolites by grinn ids and and connect to Gene and Pathway
txtInput <- "['C00024','C00136','C05269']"
organism <- "'Homo sapiens'"
result <- pairToNode(txtInput, organism, label=c("Gene","Pathway"), returnAs="all", searchBy="grinn")
}
\author{
Kwanjeera W \email{kwanich@ucdavis.edu}
}
\references{
Modify the output generation for Cytoscape.js from \url{https://github.com/cytoscape/r-cytoscape.js/blob/master/cytoscapeJsSimpleNetwork.R}

Cytoscape.js, Network library for analysis and visualisation \url{http://js.cytoscape.org/}
}
\seealso{
\url{http://js.cytoscape.org/}
}

