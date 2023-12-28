
library(GenomicDataCommons)
library(magrittr)
library(TCGAutils)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

query <- GDCquery(
  project = "TCGA-PRAD", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Raw Simple Somatic Mutation",
  access = "controlled",
  workflow.type = "MuSE"
)
GDCdownload(query, method = "client", token.file = ".\\gdc-user-token.2023-08-21T15_57_01.259Z.txt")
