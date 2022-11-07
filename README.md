# Introduction
ShinyCell is an extremely usefull R package generated by [Ouyang et al. ShinyCell: Simple and sharable visualisation of single-cell gene expression data. Bioinformatics, doi:10.1093/bioinformatics/btab209](http://dx.doi.org/10.1093/bioinformatics/btab209) that allows users to simply create interactive Shiny-based web applications to visualise single-cell data. Since ShinyCell designed for using singlecell dataset objects (such as h5ad, loom, Seurat,SingleCellExperiment) thus, in this repository we share required files to make it compatible with Debit-seq and spatial maps.

# How to run 
as [here](https://github.com/SGDDNB/ShinyCell#quick-start-guide) explains, to run the app locally, first make a folder with the favorite name and copy all of the shared files (that include .rds, .h5 and .txt files) plus two important .R files server.R and ui.R (that can be found in this repository) into the new created folder as shown below.

![alt text](https://github.com/atlasxomics/ShinyCell/blob/main/images/Screenshot%202022-11-07%20at%2010.28.29%20AM.png)

Then open either server.R or ui.R using RStudio and click on "Run App" in the top right corner. Now, you are all set to explore your data!

![alt text](https://github.com/atlasxomics/ShinyCell/blob/main/images/Screenshot%202022-11-07%20at%2010.28.53%20AM.png)

# Required Libraires
make sure your RStudio already installed the following librarires: shiny, shinyhelper, data.table, Matrix, DT, magrittr. Since all of them are CRAN packages you can install them by the following code:

install.packages(c('shiny', 'shinyhelper', 'data.table', 'Matrix', 'DT', 'magrittr'))

# Quick Start Guide
The shiny app contains seven tabs (highlighted in blue box), ...
