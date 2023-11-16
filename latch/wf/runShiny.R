library("ArchR")
library("GenomicRanges")
library('BSgenome')
library("dplyr")                                    # Load dplyr package
library("plyr")                                     # Load plyr package
library("readr")
library(qdap)
library('Seurat')
library(ShinyCell)
library('BSgenome.Mmusculus.UCSC.mm10')
library('BSgenome.Hsapiens.UCSC.hg38')
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library(data.table))


dir.create("/root/results", showWarnings = FALSE)
setwd("/root/results")

rawPath <- "/root/"
dataFiles <- dir(rawPath, "*.R$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

find_func <- function(tempdir,pattern){
    
  list.files(
  path = tempdir, # replace with the directory you want
  pattern = pattern, # has "test", followed by 0 or more characters,
                             # then ".csv", and then nothing else ($)
  full.names = TRUE # include the directory in the result
        , recursive = TRUE
)
    
}

args <- commandArgs(trailingOnly=TRUE)
tempdir <- args[1]

print(tempdir)
print(system(paste0("ls ", tempdir), intern = TRUE))
ArchRobj <- system(paste0("find ", tempdir, " -name '*_ArchRProject' -type d"), intern = TRUE)
print(ArchRobj)

proj3 <- loadArchRProject(path = ArchRobj, force = FALSE, showLogo = TRUE)


system(paste0("cp -r ", ArchRobj, " ."))

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tempdir <- args[1]

runs <- find_func(tempdir, "_SeuratObj*\\.rds$")

runs_spl <- strsplit(runs,"\\/|_")
runs_spl

n <- length(runs_spl[[1]]) - 3

fun1 <-  function(lst,nth){
    sapply(lst, `[`, nth)

}

names <- fun1(runs_spl,n+1)

names(runs) <- names

all <-  list()
for (i in seq_along(runs)){
        all[[i]] <- readRDS(runs[[i]])
}

for (i in seq_along(names)){
    
 all[[i]] <- RenameCells(all[[i]], new.names = paste0(unique(all[[i]]@meta.data$Sample)
                                               ,"#",colnames(all[[i]]),"-1"
                                                        ))
 
}

inputs1 <- all

# motif objects

runs <- find_func(tempdir, "_SeuratObjMotif*\\.rds$")

runs_spl <- strsplit(runs,"\\/|_")
runs_spl

n <- length(runs_spl[[1]]) - 3

fun1 <-  function(lst,nth){
  sapply(lst, `[`, nth)
  
}

names <- fun1(runs_spl,n+1)

names(runs) <- names

all2 <-  list()
for (i in seq_along(runs)){
  all2[[i]] <- readRDS(runs[[i]])
}


for (i in seq_along(names)){

  all2[[i]] <- RenameCells(all2[[i]], new.names = paste0(unique(all2[[i]]@meta.data$Sample)
                                                        ,"#",colnames(all2[[i]]),"-1"
  ))

}

inputs2 <- all2

inputs3 <- find_func(tempdir,"UMAPHarmony.csv")
req_genes1 <- find_func(tempdir,"req_genes1.csv")
req_genes2 <- find_func(tempdir, "req_genes2.csv")
req_genes3 <- find_func(tempdir,"req_genes3.csv")

req_motifs1 <- find_func(tempdir,"req_motifs1.csv")
req_motifs2 <- find_func(tempdir,"req_motifs2.csv")
req_motifs3 <- find_func(tempdir,"req_motifs3.csv")

project <- args[2]
groupBy <- args[3]



#===========================#===========================#===========================#=================

main_func <- function(seurat_lst){
find_samples_name <- function(seurat_lst){
    
    sapply(seq_along(seurat_lst), function(i) unique(seurat_lst[[i]]@meta.data$Sample))

    
}

samples <- find_samples_name(seurat_lst)
    
D00_fun <- function(seurat_lst){
   toRemove <- lapply(seurat_lst, function(x) {names(which(colSums(is.na(x@assays[[1]]@counts))>0))}) 
    mapply(function(x,y) x[,!colnames(x) %in% y],seurat_lst,toRemove)
}


D00 <- D00_fun(seurat_lst)
Spatial_D00_fun <- function(D00){
    
Spatial_D00 <- lapply(D00, function(x) as.data.frame(x@images[[1]]@coordinates[,c(5,4)]))
Spatial_D00 <- lapply(Spatial_D00, function(x) {colnames(x) <- paste0("Spatial_", 1:2)
                                               x
                                               })  
# Spatial_D00 <- lapply(Spatial_D00, setNames, nm = paste0("Spatial_", 1:2))

lapply(Spatial_D00, function(x) {x$Spatial_2 <- -(x$Spatial_2) 
                                x
                                })
#     return(Spatial_D00)
}
Spatial_D00 <- Spatial_D00_fun(D00)

                      
Spatial_D00_all_fun <- function(Spatial_D00){
    
    tmp <- lapply(seq_along(Spatial_D00), function(i) {bind_rows(Spatial_D00[-i])})
    
    tmp <- lapply(tmp, function(x) {x$Spatial_1<- 0
                             x
                             })
    tmp <- lapply(tmp, function(x) {x$Spatial_2<- 0
                             x
                             })
    
    tmp <- lapply(seq_along(Spatial_D00), function(i) {as.matrix(rbind(Spatial_D00[[i]],tmp[[i]]))
        
        
    })
    
}
           
Spatial_D00_all <- Spatial_D00_all_fun(Spatial_D00) 
                      
temp_fun <- function(D00){

temp <- lapply(D00,function(x) as.data.frame(x@assays[[1]]@counts))
temp <- lapply(temp, function(x)  {x$region <- rownames(x)
                                  x
                                  })  
         lapply(temp, function(x) {rownames(x) <- NULL
                                  x
                                  })      

    }
                
temp <- temp_fun(D00)                      
                      
                          
# merge seurat objects
combined_mat <-reduce(temp, full_join, by = "region");

rownames(combined_mat) <- combined_mat$region
combined_mat$region<- NULL
# removd extra cells
extra_cells <- setdiff(colnames(combined_mat),rownames(Spatial_D00_all[[1]]))
combined_mat <- combined_mat[,which(!colnames(combined_mat)%in%extra_cells)]
combined_mat <- as.matrix(combined_mat)
    
# clean columns of meta data per sample that attached sample's name before rbind
l <- D00
l <- lapply(l, function(x) { colnames(x@meta.data) <- gsub(paste0("_",Assays(x)),"",colnames(x@meta.data));x})
D00 <- l


# first get the list of meta data
list_of_metadata <- lapply(D00, function(x) x@meta.data)
# rbind meta data per samples
meta.data <- do.call("rbind", list_of_metadata)
write.csv(meta.data,'req_meta_data.csv', row.names = T)

    
combined <- CreateSeuratObject(
  counts = combined_mat,
  assay = "scATAC",
  meta.data = meta.data
)

combined@meta.data$Clusters <- factor(combined@meta.data$Clusters,levels = c(
paste0("C",seq_along(unique(combined@meta.data$Clusters)))))

Spatial_D00 <- list()
for (i in seq_along(samples)) {
    Spatial_D00[[i]] <- Spatial_D00_all[[i]][colnames(combined),]
    combined[[paste0(samples[i],"Spatial")]] <- CreateDimReducObject(embeddings = Spatial_D00[[i]]
                                              , key = paste0(samples[i],"Spatial_")
                                              , assay = DefaultAssay(combined))
    
}                           
                           
# we need to run Variable Features
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
                           
UMAPHarmony <- UMAPHarmony[match(colnames(combined),UMAPHarmony$X),]
rownames(UMAPHarmony) <- UMAPHarmony$X
UMAPHarmony$X <- NULL
UMAPHarmony <- as.matrix(UMAPHarmony)

combined[["UMAPHarmony"]] <- CreateDimReducObject(embeddings = UMAPHarmony
                                              , key = "UMAPHarmony_"
                                              , assay = DefaultAssay(combined))

return(combined)

}


#===========================#===========================#===========================#=================

spatial_in_tissue.obj <- inputs1
spatial_in_tissue_motif.obj <- inputs2

UMAPHarmony <- read.csv(inputs3)

combined <- main_func(spatial_in_tissue.obj)
combined_m <- main_func(spatial_in_tissue_motif.obj)



#===========================
scConf1 = createConfig(combined)
makeShinyFiles(combined, scConf1, gex.assay = "scATAC", gex.slot = "data",
               gene.mapping = TRUE, shiny.prefix = "sc1",
          #     shiny.dir = "../../../extra_analysis/ShinyApp/shinyAppMulti_Pieper_Lab/",
               default.gene1 = "Tiam1", default.gene2 = "Ccbe1",
          #       default.multigene = c("DCT","PDCD1","CD19","TIGIT","CD8A","SEMA4D","GZMA","CCL3","CD4","SDF4","B3GALT6"),
               default.dimred = c("UMAPHarmony_1", "UMAPHarmony_2"))




scConf2 = createConfig(combined_m)
makeShinyFiles(combined_m, scConf2, gex.assay = "scATAC", gex.slot = "counts",
              gene.mapping = TRUE, shiny.prefix = "sc2",
 # shiny.dir = "../../../extra_analysis/ShinyApp/shinyAppMulti_Pieper_Lab/",
              default.gene1 = "RFX3-1018", default.gene2 = "NEUROG2-1580",
  # default.multigene = c('WT1-266','SP6-275','SP5-279','SP4-180','TFAP2C-3','ZFX-158','SP2-232','SP1-267',),
              default.dimred = c("UMAPHarmony_1", "UMAPHarmony_2"))


citation = list(
  title   = paste0(project," Data Analysis")
                       )
makeShinyCodesMulti(
  shiny.title = paste0(project,"_Lab Data Analysis"), 
    shiny.footnotes = citation,
  shiny.prefix = c("sc1"
                , "sc2"
                  ),

  shiny.headers = c("Gene Accessibility"
                      , "Peak/Motifs"
                   ), 
  shiny.dir = "./shinyApp"
) 

#==============================================================

sc1def <- readRDS("/root/results/shinyApp/sc1def.rds")
sc2def <- readRDS("/root/results/shinyApp/sc2def.rds")


find_samples_name <- function(lst){

    sapply(seq_along(lst), function(i) unique(lst[[i]]@meta.data$Sample))
}   
samples <- find_samples_name(spatial_in_tissue.obj)


D00 <- list()
for (i in seq_along(samples)) {
D00[[i]] <- spatial_in_tissue.obj[[i]]
nal_cols <- which(colSums(is.na(D00[[i]]@assays[[1]]@counts))>0)
toRemove <- names(nal_cols)
D00[[i]] <- D00[[i]][,!colnames(D00[[i]]) %in% toRemove]
}
Spatial_D00 <- list()
for (i in seq_along(samples)) {
Spatial_D00[[i]] <- as.data.frame(D00[[i]]@images[[1]]@coordinates[,c(5,4)])
colnames(Spatial_D00[[i]]) <- paste0("Spatial_", 1:2)
Spatial_D00[[i]]$Spatial_2 <- -(Spatial_D00[[i]]$Spatial_2)
}
l <- Spatial_D00
xlim <- lapply(l, function(x) { xlim <- c(min(x[,1]),max(x[,1]));xlim})
ylim <- lapply(l, function(y) { ylim <- c(min(y[,2]),max(y[,2]));ylim})




sc1def$limits <- list()
for (i in seq_along(samples)) {
sc1def[['limits']][[paste0(samples[i],"Spatial1")]]<- c(
                                            min(xlim[[i]]),max(xlim[[i]])
                                           ,min(ylim[[i]]),max(ylim[[i]])
                                          )


}




sc2def$limits <- list()
for (i in seq_along(samples)) {
sc2def[['limits']][[paste0(samples[i],"Spatial1")]]<- c(
                                           min(xlim[[i]]),max(xlim[[i]])
                                          ,min(ylim[[i]]),max(ylim[[i]])
                                         )


}


#==================================



sc1def$meta1<- "Clusters"
sc1def$meta2<- "Condition"
sc1def$meta3<- "Sample"

sc2def$meta1<- "Clusters"
sc2def$meta2<- "Condition"
sc2def$meta3<- "Sample"


sc1def$Condition1 <- unique(combined@meta.data$Condition)[1]
sc1def$Condition2 <- unique(combined@meta.data$Condition)[2]

sc2def$Condition1 <- unique(combined_m@meta.data$Condition)[1]
sc2def$Condition2 <- unique(combined_m@meta.data$Condition)[2]



sc1def$Clusters <- read.csv(req_genes1)$x
sc1def$Condition <- read.csv(req_genes2)$x
sc1def$Sample <- read.csv(req_genes3)$x

sc2def$Clusters <- read.csv(req_motifs1)$x
sc2def$Condition <- read.csv(req_motifs2)$x
sc2def$Sample <- read.csv(req_motifs3)$x


sc1def$dimred[3] <- paste0(names(combined@reductions)[1],"1")
sc1def$dimred[4] <- paste0(names(combined@reductions)[1],"2")
sc1def$dimred[5] <- paste0(names(combined@reductions)[2],"1")
sc1def$dimred[6] <- paste0(names(combined@reductions)[2],"2")


sc2def$dimred[3] <- paste0(names(combined_m@reductions)[1],"1")
sc2def$dimred[4] <- paste0(names(combined_m@reductions)[1],"2")
sc2def$dimred[5] <- paste0(names(combined_m@reductions)[2],"1")
sc2def$dimred[6] <- paste0(names(combined_m@reductions)[2],"2")


saveRDS(sc1def,"/root/results/shinyApp/sc1def.rds")
saveRDS(sc2def,"/root/results/shinyApp/sc2def.rds")


sc1conf <- readRDS("/root/results/shinyApp/sc1conf.rds")
sc2conf <- readRDS("/root/results/shinyApp/sc2conf.rds")

fav <- which(sc1conf$ID == "Clusters" | sc1conf$ID == "Condition" | sc1conf$ID == "Sample")
rest <- which(!sc1conf$ID == "Clusters" & !sc1conf$ID == "Condition" & !sc1conf$ID == "Sample")

sc1conf <- sc1conf[c(fav,rest),]



fav <- which(sc2conf$ID == "Clusters" | sc2conf$ID == "Condition" | sc2conf$ID == "Sample")
rest <- which(!sc2conf$ID == "Clusters" & !sc2conf$ID == "Condition" & !sc2conf$ID == "Sample")

sc2conf <- sc2conf[c(fav,rest),]


saveRDS(sc1conf,"/root/results/shinyApp/sc1conf.rds")
saveRDS(sc2conf,"/root/results/shinyApp/sc2conf.rds")



# rawPath <- args[1]
rawPath <- "/root/results/"

dataFiles <- dir(rawPath, "*.rds$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

saveRDS(combined,"/root/results/shinyApp/combined.rds")
saveRDS(combined_m,"/root/results/shinyApp/combined_m.rds")


# 
# 
dataFiles <- dir(rawPath, "*.csv$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)


dataFiles <- dir(rawPath, "inpMarkers.txt", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

dataFiles <- dir(rawPath, "inpMarkers_motif.txt", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)


rawPath <- "/root/"
dataFiles <- dir(rawPath, "*.R$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

rawPath <- "/root/results/"
dataFiles <- dir(rawPath, "req_meta_data.csv", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)


system ("cp -r /root/results/*_ArchRProject /root/results/shinyApp")
system ("cp -r /root/www /root/results/shinyApp/www")

if (length(unique(proj3$Condition))<=1){
  
  file.remove("/root/results/shinyApp/ui.R")
  file.rename("/root/results/shinyApp/ui_v2.R","/root/results/shinyApp/ui.R")
} else {
  file.remove("/root/results/shinyApp/ui_v2.R")
  
}
 
 
system ("cp -r /root/results/shinyApp /root/shinyApp")

system ("rm -r /root/results")

setwd("/root")
system("ls -a |grep -v 'shinyApp' | xargs rm -r")

file.rename(list.files(pattern="shinyApp"), paste0(project,"_shinyApp"))


