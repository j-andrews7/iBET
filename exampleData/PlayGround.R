library(devtools)
devtools::load_all()
library(iBET)
library(DESeq2)
shh(library("scRNAseq"))
shh(library("scran"))
shh(library("scuttle"))
shh(library("scater")) 


#MBD3: 
res <- readRDS("/Volumes/JM/Development/iBET/exampleData/MBD3KO-res.rds")
dds <- readRDS("/Volumes/JM/Development/iBET/exampleData/MBD3KO-dds.rds")
res <- list("MBD3_res" = res)
dds <- list("MBD3_dds" = dds)


