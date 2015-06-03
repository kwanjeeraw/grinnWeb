grinnWeb
=========
An R web application for omics integration in the context of metabolic networks

Version: 1.0 (10 Mar 2015)

Description
=========
grinnWeb is a bioinformatics platform contains a graph database, and the R web application for metabolomics studies.
Based on metabolic networks from KEGG, grinnWeb can generate two types of networks: a network of connected metabolites and an integrated network of metabolites, proteins, genes and pathways, for network analyses and visualization.

Using online grinnWeb
=========
* For online usage under UC Davis network, click [here](http://grinn.genomecenter.ucdavis.edu/ocpu/user/kwanich/library/grinnWeb/www/)
 
Using local grinnWeb
=========
* <b>OR</b> for local usage, installation required [R](http://www.r-project.org/) and [OpenCPU](https://www.opencpu.org).
 1. Follow [the guidelines](https://www.opencpu.org/download.html) to install OpenCPU server locally.
 2. Require Neo4j 2.1.5 for the grinn internal database (local version), please send us an email for the grinn database files, currently available: Human database.

    - Download and then unzip [Neo4j server](http://neo4j.com/download/other-releases/)
    - Extract and move the grinn database files to the Neo4j server directory
    - Start the Neo4j server, for windows: Double-click on %NEO4J_HOME%\bin\Neo4j.bat, for linux: ./bin/neo4j start 
for more details see [here](http://neo4j.com/docs/stable/server-installation.html)
 3. Install grinn R package from github
 ```
    install.packages("devtools")
    library(devtools)
    install_github("kwanjeeraw/grinnWeb")
 ```
 
After installation, the following code is to run Grinn locally.
 ```
    library(grinnWeb)
    library(opencpu)
    opencpu$browse("library/grinnWeb/www")
 ```

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinnWeb/)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinnWeb/blob/master/LICENSE)
