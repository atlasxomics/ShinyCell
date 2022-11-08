# Introduction
ShinyCell is an extremely usefull R package generated by [Ouyang et al. ShinyCell: Simple and sharable visualisation of single-cell gene expression data. Bioinformatics, doi:10.1093/bioinformatics/btab209](http://dx.doi.org/10.1093/bioinformatics/btab209) that allows users to simply create interactive Shiny-based web applications to visualise single-cell data. Since ShinyCell designed for using singlecell dataset objects (such as h5ad, loom, Seurat,SingleCellExperiment) thus, in this repository we share required files to make it compatible with Debit-seq and spatial maps.

# How to run 
as [here](https://github.com/SGDDNB/ShinyCell#quick-start-guide) explains, to run the app locally, first make a folder with the favorite name and copy all of the shared files (that include .rds, .h5 and .txt files) plus two important .R files server.R and ui.R (that can be found also [here](https://github.com/atlasxomics/ShinyCell/tree/main/R)) into the new created folder as shown below.

<p align="center">
<img src="https://github.com/atlasxomics/ShinyCell/blob/main/images/required_files.png" width="400" height="220">
  </p>
Next, you need to browse to the Rstudio working directory and change it to the new created folder as shown below in "R Sessions" secion.

<p align="center">
<img src="https://github.com/atlasxomics/ShinyCell/blob/main/images/rstudio_settings.png" width="280" height="260">
  </p>

Then open either server.R or ui.R using RStudio and click on "Run App" in the top right corner. Now, you are all set to explore your data!

<p align="center">
<img src="https://github.com/atlasxomics/ShinyCell/blob/main/images/how_to_run_app.png" width="900" height="245">
  </p>


# Required Libraires
make sure your RStudio already installed the following librarires: shiny, shinyhelper, data.table, Matrix, DT, magrittr, ggplot2, ggrepel, hdf5r, ggdendro, gridExtra, ggseqlogo and circlize. Since all of them are CRAN packages you can install them by the following code:

install.packages(c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "magrittr","ggplot2","ggrepel","hdf5r","ggdendro","gridExtra", "ggseqlogo", "circlize"))





# Quick Start Guide
Please first check [here](https://github.com/SGDDNB/ShinyCell#quick-start-guide) for a quick guide. The main differences between original and the version in this repository is the ability of visualization two different samples side by side.
The current version contains eleven tabs , and this might be developed to more tabs in future. 


<p align="center">
<img src="https://github.com/atlasxomics/ShinyCell/blob/main/images/app_feature.png" >
  </p>
