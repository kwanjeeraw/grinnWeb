GrinnWeb
=========
An R web application to a graph database and R package for omics integration

Version: 1.0 (10 Mar 2015)

Description
=========
A bioinformatics platform contains a graph database, an R package and web application for metabolomics studies.
It can generate two types of networks: a network of connected metabolites and an integrated network of metabolites, proteins, genes and pathways, for network analyses and visualization.

Installation
=========
* <b>For online usage</b>, click [here](http://grinn.genomecenter.ucdavis.edu/ocpu/user/kwanich/library/grinn/www/)
* <b>For local usage</b>, install required [R](http://www.r-project.org/) and [OpenCPU](https://www.opencpu.org).
 1. Follow [the guidelines](https://www.opencpu.org/download.html) to install OpenCPU server locally.
 2. Install grinn R package from github
 ```
    install.packages("devtools")
    library(devtools)
    install_github("kwanjeeraw/grinnWeb")
 ```
 
Local running
=========
After installation, the following code is to run Grinn locally.
 ```
    library(grinnWeb)
    library(opencpu)
    opencpu$browse("library/grinnWeb/www")
 ```

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinn/)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinn/blob/master/LICENSE)
