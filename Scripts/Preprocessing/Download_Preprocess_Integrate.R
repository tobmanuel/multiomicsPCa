

start_time <- Sys.time()

setwd("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\Raw_RNAseq") #location of the manifest downloaded from TCGA portal


# Change uuids to barcodes ---------------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   bioPkgs = c("GenomicDataCommons", "TCGAutils", "magrittr", "TCGAbiolinks", "DESeq2","tidyverse","airway")
# BiocManager::install(bioPkgs)

library(GenomicDataCommons)
library(magrittr)
library(TCGAutils)
manifest <- read.table("gdc_manifest.2023-08-08.txt", header = TRUE) #reads in the manifest
file_uuids <- manifest$id #selects the id column
head(file_uuids)

filebarcodes = UUIDtoBarcode(file_uuids, from_type = "file_id") #Converts the UUID to barcodes


## ----download Gene Expression Data------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)




query.exp.hg38 <- GDCquery(
  project = "TCGA-PRAD", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode =  filebarcodes$associated_entities.entity_submitter_id
)
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(
  query = query.exp.hg38,
  save = TRUE, 
  save.filename = "Genexp.rda"
)


PRADMatrix <- assay(expdat,"unstranded") #creates a matrix containing unstranded transcript counts

annot = data.frame(expdat@rowRanges@elementMetadata$gene_id,expdat@rowRanges@elementMetadata$gene_type)

annot = annot[which(annot$expdat.rowRanges.elementMetadata.gene_type == "protein_coding"),]

mRNAMatrix = PRADMatrix[annot$expdat.rowRanges.elementMetadata.gene_id,]



setwd("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data") #location of the phenotype data obtained at UCSC Xena portal
#Read phenotypic data
pheno <- read.csv("TCGA-PRAD.GDC_phenotype.tsv",sep = '\t') 

#create phenotype from chosen column
phenotype  <- data.frame(pheno [ ,c ("submitter_id.samples","primary_therapy_outcome_success")])

#Explore phenotype
table(phenotype$primary_therapy_outcome_success)

#Exclude samples with no Phenotype data
phenotype <- phenotype[-which(phenotype$primary_therapy_outcome_success == ""), ]



#Rename Samples (phenotype)
oldRowNames <- (phenotype$submitter_id.samples)
newNames <- substring((phenotype$submitter_id.samples),6,16)
phenotype$submitter_id.samples <- newNames


#Rename Samples (Gene expression)
colnames(mRNAMatrix )
oldfkpmColNames <- colnames(mRNAMatrix )
newColNames <- substring(colnames(mRNAMatrix ), 6, 16)
colnames(mRNAMatrix ) <- c(newColNames)


#find common samples
commonSamples <- intersect((colnames(mRNAMatrix )), phenotype$submitter_id.samples)
commonSamples


#create new gene expression data containing common samples
newPRAD  <- mRNAMatrix [,c(commonSamples)]

#create new phenotype data containing common samples
newPHENOTYPE <- phenotype[which(phenotype$submitter_id.samples %in% commonSamples), ]


#sort phenotype by gene expression data order
geneOrder=colnames(newPRAD)
geneOrder
sortedPhenotype <- newPHENOTYPE[match(geneOrder,newPHENOTYPE$submitter_id.samples), ]


all(colnames(newPRAD)==sortedPhenotype$submitter_id.samples)





# Differential Expression analysis ----------------------------------------

library(DESeq2)
library(tidyverse)
library(airway)

dds <- DESeqDataSetFromMatrix(countData = newPRAD,
                              colData = sortedPhenotype,
                              design = ~ primary_therapy_outcome_success)
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# set the factor level
table(newPHENOTYPE$primary_therapy_outcome_success)

dds$primary_therapy_outcome_success <- relevel(dds$primary_therapy_outcome_success, ref = "Complete Remission/Response")


#Run dds <- DESeq(dds)
dds <- DESeq(dds)
parRes<- results(dds, alpha = 0.01, contrast = c("primary_therapy_outcome_success", "Complete Remission/Response","Partial Remission/Response")) #checks the genes that differ between Complete and partial Remission

summary(parRes) #Visualize results and obtain the number of differentially expressed genes

proRes<- results(dds, alpha = 0.01, contrast = c("primary_therapy_outcome_success", "Complete Remission/Response","Progressive Disease")) #checks the genes that differ between Complete remission and progressive disease

summary(proRes) #Visualize results and obtain the number of differentially expressed genes

staRes<- results(dds, alpha = 0.01, contrast = c("primary_therapy_outcome_success", "Complete Remission/Response","Stable Disease")) #checks the genes that differ between Complete remission and stable disease

summary(staRes) #Visualize results and obtain the number of differentially expressed genes (DEG)


# MA plot for each constrast
plotMA(parRes) 
plotMA(proRes)
plotMA(staRes)

#Select DEGs
parRes= parRes[order(abs(parRes$log2FoldChange), decreasing = TRUE),] #sort results in respect to the absolute values of the log2FoldChange of each gene
parDEGs =rownames(parRes)[1:1434] #extract the DEGs


proRes= proRes[order(abs(proRes$log2FoldChange), decreasing = TRUE),] #sort results in respect to the absolute values of the log2FoldChange of each gene
proDEGs =rownames(proRes)[1:915] #extract the DEGs


staRes= staRes[order(abs(staRes$log2FoldChange), decreasing = TRUE),] #sort results in respect to the absolute values of the log2FoldChange of each gene
staDEGs =rownames(staRes)[1:850] #extract the DEGs



#Plot Venn diagram to reveal the distribution of differentially expressed genes
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")
library(ggvenn)

x <- list(partialResponse = parDEGs, progresiveDisease  = proDEGs,stableDisease=staDEGs) 

ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)


#select genes that are differentially expresssed exclusively between partial remission and other groups
nparDEGs = parDEGs [!(parDEGs %in% proDEGs)]
nparDEGs = nparDEGs [!(nparDEGs %in% staDEGs)]

#select genes that are differentially expresssed exclusively between progressive disease and other groups
nproDEGs = proDEGs [!(proDEGs %in% parDEGs)]
nproDEGs = nproDEGs [!(nproDEGs %in% staDEGs)]

#select genes that are differentially expresssed exclusively between stable and other groups
nstaDEGs = staDEGs [!(staDEGs %in% proDEGs)]
nstaDEGs = nstaDEGs [!(nstaDEGs %in% parDEGs)]

DEGs = union(union(nparDEGs,nproDEGs),nstaDEGs) #combine the exclusive DEGs


# # create new gene expression data containing differential expressed genes ----------------------
degMRNA = newPRAD[which(rownames(newPRAD) %in% DEGs),]

degMRNA = t(degMRNA)
degMRNA = data.frame(degMRNA)

colnames(degMRNA)<-paste(colnames(degMRNA),"mRNA",sep="_")


# 
# # Download miRNA data -----------------------------------------------------
setwd("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\Raw_miRNAseq")
manifest <- read.table("gdc_manifest.2023-08-08.txt", header = TRUE)
file_uuids <- manifest$id
head(file_uuids)

filebarcodes = UUIDtoBarcode(file_uuids, from_type = "file_id")

query.mirna <- GDCquery(
  project = "TCGA-PRAD",
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling",
  barcode = filebarcodes$associated_entities.entity_submitter_id,
  data.type = "miRNA Expression Quantification"
)
GDCdownload(query.mirna)
mirna <- GDCprepare(
  query = query.mirna,
  save = TRUE,
  save.filename = "mirna.rda"
)


miRNAMatrix <-mirna[ , grepl( "read_count" , colnames( mirna ) ) ]
rownames(miRNAMatrix)=mirna$miRNA_ID

#Rename Samples (miRNA)
colnames(miRNAMatrix)
oldfkpmColNames <- colnames(miRNAMatrix)
View(oldfkpmColNames)
newColNames <- substring(colnames(miRNAMatrix), 17, 27)
colnames(miRNAMatrix) <- c(newColNames)

#find common samples
commonSamples <- intersect((colnames(miRNAMatrix)), phenotype$submitter_id.samples)
commonSamples


#create new miRNA data containning common samples
newMIRNA  <- miRNAMatrix[,c(commonSamples)]




#create new phenotype data containning common samples
newPHENOTYPE <- phenotype[which(phenotype$submitter_id.samples %in% commonSamples), ]


#sort phenotype by gene expression data order
geneOrder=colnames(newMIRNA)
geneOrder
sortedPhenotype <- newPHENOTYPE[match(geneOrder,newPHENOTYPE$submitter_id.samples), ]
View(sortedPhenotype)

all(colnames(newMIRNA)==sortedPhenotype$submitter_id.samples)





# Differential Expression analysis ----------------------------------------

#complete Response

library(DESeq2)
library(tidyverse)
library(airway)

dds <- DESeqDataSetFromMatrix(countData = newMIRNA,
                              colData = sortedPhenotype,
                              design = ~ primary_therapy_outcome_success)
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# set the factor level
table(newPHENOTYPE$primary_therapy_outcome_success)

dds$primary_therapy_outcome_success <- relevel(dds$primary_therapy_outcome_success, ref = "Complete Remission/Response")


#Run dds <- DESeq(dds)
dds <- DESeq(dds)
parRes<- results(dds, alpha = 0.01, contrast = c("primary_therapy_outcome_success", "Complete Remission/Response","Partial Remission/Response"))

summary(parRes)

proRes<- results(dds, alpha = 0.01, contrast = c("primary_therapy_outcome_success", "Complete Remission/Response","Progressive Disease"))

summary(proRes)

staRes<- results(dds, alpha = 0.01, contrast = c("primary_therapy_outcome_success", "Complete Remission/Response","Stable Disease"))

summary(staRes)


# MA plot
plotMA(parRes)
plotMA(proRes)
plotMA(staRes)

#Select DEGs
parRes= parRes[order(abs(parRes$log2FoldChange), decreasing = TRUE),]
parDEGs =rownames(parRes)[1:19]


proRes= proRes[order(abs(proRes$log2FoldChange), decreasing = TRUE),]
proDEGs =rownames(proRes)[1:28]


staRes= staRes[order(abs(staRes$log2FoldChange), decreasing = TRUE),]
staDEGs =rownames(staRes)[1:26]



#Plot Venn diagram
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")
library(ggvenn)

x <- list(partialResponse = parDEGs, progresiveDisease  = proDEGs,stableDisease=staDEGs)

ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)

# 
nparDEGs = parDEGs [!(parDEGs %in% proDEGs)]
nparDEGs = nparDEGs [!(nparDEGs %in% staDEGs)]

nproDEGs = proDEGs [!(proDEGs %in% parDEGs)]
nproDEGs = nproDEGs [!(nproDEGs %in% staDEGs)]

nstaDEGs = staDEGs [!(staDEGs %in% proDEGs)]
nstaDEGs = nstaDEGs [!(nstaDEGs %in% parDEGs)]

DEGs = union(union(nparDEGs,nproDEGs),nstaDEGs)


# create new miRNA data containing differential expressed genes ----------------------
degMIRNA = newMIRNA[which(rownames(newMIRNA) %in% DEGs),]
degMIRNA = t(degMIRNA)
degMIRNA = data.frame(degMIRNA)

colnames(degMIRNA)<-paste(colnames(degMIRNA),"miRNA",sep="_")


# Download RPPA data -----------------------------------------------------
query.rppa <- GDCquery(
  project = "TCGA-PRAD", 
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)
GDCdownload(query.rppa) 
rppa <- GDCprepare(query.rppa)
rownames(rppa)=rppa$peptide_target
rppa = rppa[,6:357]



count_na_func <- function(x) sum(is.na(x)) 
rppa$count_na <- apply(rppa, 1, count_na_func)

rppa = rppa[(-which(rppa$count_na > 10)), ]

rppa = rppa[,1:352]
rppa [457,] = colSums(is.na(rppa)) 
rppa = rppa[1:456,c(rppa [457,] == 0)]

rppa = t(rppa)

rownames(rppa) <- substring(rownames(rppa), 6, 16)

rppa = data.frame(rppa)

colnames(rppa)<-paste(colnames(rppa),"rppa",sep="_")

# Download CNV data -------------------------------------------------------
setwd("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data")
cnv <- read.csv("TCGA-PRAD.gistic.tsv",sep = '\t')

# change CNV IDs
newColNews <- substring(names(cnv[,c(2:503)]), 6, 16)
oldfkpmColNames <- newColNews
gene <- "Gene.Symbol"
names(cnv[,c(2:503)]) <- newColNews
colnames(cnv) <- c(gene,newColNews)

names(cnv)<-gsub("\\.","-",names(cnv))

#CNV rows now represent samples
cnv <- t(cnv)
colnames(cnv) <- cnv[1,] #define column names
cnv<- cnv[-c(1),]
cnv = data.frame(cnv)

colnames(cnv)<-paste(colnames(cnv),"cnv",sep="_")

# Download Metylation Data --------------------------------------------------
query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "DNA Methylation",
  data.type = "Masked Intensities",
  platform = "Illumina Human Methylation 450"
)
GDCdownload(query, files.per.chunk=10)
betas <- GDCprepare(query)

meth = data.frame(betas@assays@data@listData[[1]])

meth$count_na <- apply(meth, 1, count_na_func)
table(meth$count_na)
meth = meth[(-which(meth$count_na > 0)), ]


meth[342904,] = colSums(is.na(meth)) 
meth = meth[1:342903,c(meth[342904,] == 0)]

names(meth)<-gsub("\\.","-",names(meth))

colnames(meth) <- substring(colnames(meth), 6, 16)

#find common samples
commonSamples <- intersect((colnames(meth)), phenotype$submitter_id.samples)
commonSamples


#create new miRNA data containning common samples
newMETH  <- meth[,c(commonSamples)]




#create new phenotype data containning common samples
newPHENOTYPE <- phenotype[which(phenotype$submitter_id.samples %in% commonSamples), ]


#sort phenotype by gene expression data order
geneOrder=colnames(newMETH)
geneOrder
sortedPhenotype <- newPHENOTYPE[match(geneOrder,newPHENOTYPE$submitter_id.samples), ]
sortedPhenotype$primary_therapy_outcome_success[which(sortedPhenotype$primary_therapy_outcome_success == "Complete Remission/Response")] = "CompleteRemission"
sortedPhenotype$primary_therapy_outcome_success[which(sortedPhenotype$primary_therapy_outcome_success == "Partial Remission/Response")] = "PartialRemission"
sortedPhenotype$primary_therapy_outcome_success[which(sortedPhenotype$primary_therapy_outcome_success == "Progressive Disease")] = "ProgressiveDisease"
sortedPhenotype$primary_therapy_outcome_success[which(sortedPhenotype$primary_therapy_outcome_success == "Stable Disease")] = "StableDisease"
View(sortedPhenotype)

all(colnames(newMETH)==sortedPhenotype$submitter_id.samples)



# Differential methylation ------------------------------------------------
#DE with limma https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
library(edgeR)


d0 <- DGEList(newMETH)
d0 <- calcNormFactors(d0,method="none")

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

group <- as.factor(sortedPhenotype$primary_therapy_outcome_success)
group

mm <- model.matrix(~0 + group)

y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

par_contr <- makeContrasts(groupCompleteRemission - groupPartialRemission, levels = colnames(coef(fit)))
par_contr

tmp <- contrasts.fit(fit, par_contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)

length(which(top.table$adj.P.Val < 0.01))

parDMGs = rownames(top.table)[1:length(which(top.table$adj.P.Val < 0.01))]


pro_contr <- makeContrasts(groupCompleteRemission - groupProgressiveDisease, levels = colnames(coef(fit)))
pro_contr

tmp <- contrasts.fit(fit, pro_contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)

length(which(top.table$adj.P.Val < 0.01))

proDMGs = rownames(top.table)[1:length(which(top.table$adj.P.Val < 0.01))]

sta_contr <- makeContrasts(groupCompleteRemission - groupStableDisease, levels = colnames(coef(fit)))
sta_contr

tmp <- contrasts.fit(fit, sta_contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)

length(which(top.table$adj.P.Val < 0.01))

staDMGs = rownames(top.table)[1:length(which(top.table$adj.P.Val < 0.01))]


#Plot Venn diagram
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")
library(ggvenn)

x <- list(partialResponse = parDMGs, progresiveDisease  = proDMGs,stableDisease=staDMGs)

ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)

# 
nparDMGs = parDMGs [!(parDMGs %in% proDMGs)]
nparDMGs = nparDMGs [!(nparDMGs %in% staDMGs)]

nproDMGs = proDMGs [!(proDMGs %in% parDMGs)]
nproDMGs = nproDMGs [!(nproDMGs %in% staDMGs)]

nstaDMGs = staDMGs [!(staDMGs %in% proDMGs)]
nstaDMGs = nstaDMGs [!(nstaDMGs %in% parDMGs)]

DMGs = union(union(nparDMGs,nproDMGs),nstaDMGs)


# create new methylation data containing differential expressed genes ----------------------
degMETH = newMETH[which(rownames(newMETH) %in% DMGs),]
degMETH = t(degMETH)
degMETH = data.frame(degMETH)

colnames(degMETH)<-paste(colnames(degMETH),"meth",sep="_")



# Import SNV --------------------------------------------------------------
setwd("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\SNV")

library("BEDMatrix")

m = BEDMatrix("snv_nomiss.bed") #imports BED file containing snv data

SnpMartrix = data.frame(m[1:nrow(m), 1:ncol(m)]) #coverts to dataframe
names(SnpMartrix)<-gsub("\\.","_",names(SnpMartrix)) #adjust naming convention

#remove normal tissue sample data from the dataframe (retaining only tumor sample data)
remove = as.vector(seq(1, 1006, by= 2)) 
SnpMartrix = SnpMartrix[-c(remove),]

fileUUID <- read.table("subset_vcfsaa", quote="\"", comment.char="") #read in sample order used during merging
prad.snv.manifest <- read.delim("gdc_manifest.2023-08-24.txt") #read in snv manifest
fileUUID$V2 = prad.snv.manifest$id[match(fileUUID$V1,prad.snv.manifest$filename)] #create database for matching IDs 


fileBarcode= UUIDtoBarcode(fileUUID$V2, from_type = "file_id") #convert to UUID

index = substr(fileBarcode$associated_entities.entity_submitter_id,14,14)== '0' #create index of only tumor samples

fileBarcode = fileBarcode[index,] #use index to subset fileBarcodes

SnpMartrix = data.frame(t(SnpMartrix))

colnames(SnpMartrix) = fileBarcode$associated_entities.entity_submitter_id #assign Barcodes to individuals

colnames(SnpMartrix) <- substring(colnames(SnpMartrix), 6, 16)
names(SnpMartrix)<-gsub("\\.","-",names(SnpMartrix))

snv = data.frame(t(SnpMartrix))
colnames(snv)<-paste(colnames(snv),"snv",sep="_")
row.names(snv)
# Collect Common samples --------------------------------------------------
commonSamples <- intersect(intersect(intersect(intersect(intersect((rownames(degMRNA)), (rownames(degMIRNA))),(rownames(cnv))),(rownames(rppa))),(rownames(degMETH))),(rownames(snv)))
commonSamples

setwd("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data")


# Create New Datasets from Common Samples ---------------------------------

#Gene Expression

degMRNA <- degMRNA[c(commonSamples),]
View (degMRNA)


#miRNA
degMIRNA <- degMIRNA[c(commonSamples),]
View (degMIRNA)

#CNV
cnv <- cnv[c(commonSamples),]
View (cnv)


#RPPA
rppa <- rppa[c(commonSamples),]
View (rppa)

#RPPA
degMETH <- degMETH[c(commonSamples),]
View (degMETH)

snv <- snv[c(commonSamples),]
View (snv)


# Homogenize rows and columns ----------------------------------------------


#Make Order the same
degMIRNA <- degMIRNA[match(rownames(degMRNA),rownames(degMIRNA)), ]
cnv <- cnv[match(rownames(degMRNA),rownames(cnv)), ]
rppa <-rppa[match(rownames(degMRNA),rownames(rppa)),]
degMETH <- degMETH[match(rownames(degMRNA),rownames(degMETH)),]
snv <- snv[match(rownames(degMRNA),rownames(snv)),]


# verification ------------------------------------------------------------

all(rownames(degMRNA)==rownames(degMIRNA))
all(rownames(degMRNA)==rownames(cnv))
all(rownames(degMRNA)==rownames(rppa))
all(rownames(degMRNA)==rownames(degMETH))
all(rownames(degMRNA)==rownames(snv))

# Merge files -------------------------------------------------------------

cnv_snv = data.frame(cbind(cnv, snv))

cnv_meth = data.frame(cbind(cnv, degMETH))

cnv_mRNA = data.frame(cbind(cnv, degMRNA))

cnv_miRNA = data.frame(cbind(cnv, degMIRNA))

cnv_rppa = data.frame(cbind(cnv, rppa))

snv_meth = data.frame(cbind(snv, degMETH))

snv_mRNA = data.frame(cbind(snv, degMRNA))

gc() #reduces memory usage

snv_miRNA = data.frame(cbind(snv, degMIRNA))

snv_rppa = data.frame(cbind(snv, rppa))

meth_mRNA = data.frame(cbind(degMETH, degMRNA))

meth_miRNA = data.frame(cbind(degMETH, degMIRNA))

meth_rppa = data.frame(cbind(degMETH, rppa))

mRNA_miRNA = data.frame(cbind(degMRNA, degMIRNA))

mRNA_rppa = data.frame(cbind(degMRNA, rppa))

miRNA_rppa = data.frame(cbind(degMIRNA, rppa))

cnv_snv_meth = data.frame(cbind(cnv,snv, degMETH))

cnv_snv_mRNA = data.frame(cbind(cnv,snv, degMRNA))

cnv_snv_miRNA = data.frame(cbind(cnv,snv, degMIRNA))

cnv_snv_rppa = data.frame(cbind(cnv,snv, rppa))

gc()

cnv_meth_mRNA = data.frame(cbind(cnv,degMETH, degMRNA))

cnv_meth_miRNA = data.frame(cbind(cnv,degMETH, degMIRNA))

cnv_meth_rppa = data.frame(cbind(cnv,degMETH, rppa))

cnv_mRNA_miRNA = data.frame(cbind(cnv,degMRNA, degMIRNA))

cnv_mRNA_rppa = data.frame(cbind(cnv,degMRNA, rppa))

cnv_miRNA_rppa = data.frame(cbind(cnv,degMIRNA, rppa))

snv_meth_mRNA = data.frame(cbind(snv,degMETH, degMRNA))

snv_meth_miRNA = data.frame(cbind(snv,degMETH, degMIRNA))

snv_meth_rppa = data.frame(cbind(snv,degMETH, rppa))

snv_mRNA_miRNA = data.frame(cbind(snv,degMRNA, degMIRNA))

snv_mRNA_rppa = data.frame(cbind(snv,degMRNA, rppa))

snv_miRNA_rppa = data.frame(cbind(snv,degMIRNA, rppa))

meth_mRNA_miRNA = data.frame(cbind(degMETH,degMRNA, degMIRNA))

meth_mRNA_rppa = data.frame(cbind(degMETH,degMRNA, rppa))

meth_miRNA_rppa = data.frame(cbind(degMETH,degMIRNA, rppa))

mRNA_miRNA_rppa = data.frame(cbind(degMRNA,degMIRNA, rppa))

cnv_snv_meth_mRNA = data.frame(cbind(cnv,snv,degMETH, degMRNA))

cnv_snv_meth_miRNA = data.frame(cbind(cnv,snv,degMETH, degMIRNA))

cnv_snv_meth_rppa = data.frame(cbind(cnv,snv,degMETH, rppa))

cnv_snv_mRNA_miRNA = data.frame(cbind(cnv,snv,degMRNA, degMIRNA))

cnv_snv_mRNA_rppa = data.frame(cbind(cnv,snv,degMRNA, rppa))

cnv_snv_miRNA_rppa = data.frame(cbind(cnv,snv,degMIRNA, rppa))

cnv_meth_mRNA_miRNA = data.frame(cbind(cnv,degMETH,degMRNA, degMIRNA))

cnv_meth_mRNA_rppa = data.frame(cbind(cnv,degMETH,degMRNA, rppa))

cnv_meth_miRNA_rppa = data.frame(cbind(cnv,degMETH,degMIRNA, rppa))

cnv_mRNA_miRNA_rppa = data.frame(cbind(cnv,degMRNA,degMIRNA, rppa))

snv_meth_mRNA_miRNA = data.frame(cbind(snv,degMETH,degMRNA, degMIRNA))

snv_meth_mRNA_rppa = data.frame(cbind(snv,degMETH,degMRNA, rppa))

snv_meth_miRNA_rppa = data.frame(cbind(snv,degMETH,degMIRNA, rppa))

snv_mRNA_miRNA_rppa = data.frame(cbind(snv,degMRNA,degMIRNA, rppa))

meth_mRNA_miRNA_rppa = data.frame(cbind(degMETH,degMRNA,degMIRNA, rppa))

cnv_snv_meth_mRNA_miRNA = data.frame(cbind(cnv,snv,degMETH,degMRNA, degMIRNA))

cnv_snv_meth_mRNA_rppa = data.frame(cbind(cnv,snv,degMETH,degMRNA, rppa))

cnv_snv_meth_miRNA_rppa = data.frame(cbind(cnv,snv,degMETH,degMIRNA, rppa))

cnv_snv_mRNA_miRNA_rppa = data.frame(cbind(cnv,snv,degMRNA,degMIRNA, rppa))

cnv_meth_mRNA_miRNA_rppa = data.frame(cbind(cnv,degMETH,degMRNA,degMIRNA, rppa))

snv_meth_mRNA_miRNA_rppa = data.frame(cbind(snv,degMETH,degMRNA,degMIRNA, rppa))

cnv_snv_meth_mRNA_miRNA_rppa = data.frame(cbind(cnv,snv,degMETH,degMRNA,degMIRNA, rppa))

# check proportions -------------------------------------------------------

beforCORR = c()
omics = list(cnv,snv, degMETH,degMRNA, degMIRNA, rppa ,cnv_snv, cnv_meth,
             cnv_mRNA, cnv_miRNA,cnv_rppa, snv_meth, snv_mRNA,snv_miRNA, snv_rppa,
             meth_mRNA,meth_miRNA, meth_rppa, mRNA_miRNA, mRNA_rppa,
             miRNA_rppa,cnv_snv_meth, cnv_snv_mRNA,cnv_snv_miRNA,cnv_snv_rppa,
             cnv_meth_mRNA, cnv_meth_miRNA, cnv_meth_rppa, cnv_mRNA_miRNA,
             cnv_mRNA_rppa,cnv_miRNA_rppa, snv_meth_mRNA, snv_meth_miRNA,
             snv_meth_rppa, snv_mRNA_miRNA , snv_mRNA_rppa ,snv_miRNA_rppa,
             meth_mRNA_miRNA, meth_mRNA_rppa, meth_miRNA_rppa, mRNA_miRNA_rppa,
             cnv_snv_meth_mRNA, cnv_snv_meth_miRNA, cnv_snv_meth_rppa, 
             cnv_snv_mRNA_miRNA, cnv_snv_mRNA_rppa, cnv_snv_miRNA_rppa,
             cnv_meth_mRNA_miRNA, cnv_meth_mRNA_rppa, cnv_meth_miRNA_rppa,
             cnv_mRNA_miRNA_rppa, snv_meth_mRNA_miRNA, snv_meth_mRNA_rppa,
             snv_meth_miRNA_rppa, snv_mRNA_miRNA_rppa, meth_mRNA_miRNA_rppa, 
             cnv_snv_meth_mRNA_miRNA, cnv_snv_meth_mRNA_rppa, cnv_snv_meth_miRNA_rppa,
             cnv_snv_mRNA_miRNA_rppa, cnv_meth_mRNA_miRNA_rppa, snv_meth_mRNA_miRNA_rppa,
             cnv_snv_meth_mRNA_miRNA_rppa)
length(omics)

omics_names = c("cnv","snv", "meth","mRNA", "miRNA", "rppa" ,"cnv_snv", "cnv_meth",
             "cnv_mRNA", "cnv_miRNA","cnv_rppa", "snv_meth", "snv_mRNA","snv_miRNA", "snv_rppa",
             "meth_mRNA","meth_miRNA", "meth_rppa", "mRNA_miRNA", "mRNA_rppa",
             "miRNA_rppa","cnv_snv_meth", "cnv_snv_mRNA","cnv_snv_miRNA","cnv_snv_rppa",
             "cnv_meth_mRNA", "cnv_meth_miRNA", "cnv_meth_rppa", "cnv_mRNA_miRNA",
             "cnv_mRNA_rppa","cnv_miRNA_rppa", "snv_meth_mRNA", "snv_meth_miRNA",
             "snv_meth_rppa", "snv_mRNA_miRNA" , "snv_mRNA_rppa" ,"snv_miRNA_rppa",
             "meth_mRNA_miRNA", "meth_mRNA_rppa", "meth_miRNA_rppa", "mRNA_miRNA_rppa",
             "cnv_snv_meth_mRNA", "cnv_snv_meth_miRNA", "cnv_snv_meth_rppa", 
             "cnv_snv_mRNA_miRNA", "cnv_snv_mRNA_rppa", "cnv_snv_miRNA_rppa",
             "cnv_meth_mRNA_miRNA", "cnv_meth_mRNA_rppa", "cnv_meth_miRNA_rppa",
             "cnv_mRNA_miRNA_rppa", "snv_meth_mRNA_miRNA", "snv_meth_mRNA_rppa",
             "snv_meth_miRNA_rppa", "snv_mRNA_miRNA_rppa", "meth_mRNA_miRNA_rppa", 
             "cnv_snv_meth_mRNA_miRNA", "cnv_snv_meth_mRNA_rppa", "cnv_snv_meth_miRNA_rppa",
             "cnv_snv_mRNA_miRNA_rppa", "cnv_meth_mRNA_miRNA_rppa", "snv_meth_mRNA_miRNA_rppa",
             "cnv_snv_meth_mRNA_miRNA_rppa")

for (omic in omics) {
  
  nCNV = length(omic[grepl("_cnv",colnames(omic))])
  nSNV = length(omic[grepl("_snv",colnames(omic))])
  nMETH = length(omic[grepl("_meth",colnames(omic))])
  nMRNA = length(omic[grepl("_mRNA",colnames(omic))])
  nMIRNA = length(omic[grepl("_miRNA",colnames(omic))])
  nRPPA = length(omic[grepl("_rppa",colnames(omic))])
  beforCORR = rbind(beforCORR,c(nCNV,nSNV,nMETH,nMRNA,nMIRNA,nRPPA))
}

row.names(beforCORR) = omics_names

write.csv(beforCORR,"BeforeCorrelation.csv", row.names = TRUE)


# Remove Highly correlated features -------------------------------------------------------------
rm_corr <- function(i){
  y = i[,arrayInd(which(apply(i, 2, var, na.rm=TRUE) != 0),dim(i)[2])] #removes constant columns
  
  tmp <- cor(as.data.frame(sapply(y, as.numeric))) #calculates correlation
  tmp[!lower.tri(tmp)] <- 0
  
  
  z <- 
    y[, !apply(tmp, 2, function(x) any(abs(x) > 0.7, na.rm = TRUE))] #removes highly correlated
  
  return(z)
}


cnv = rm_corr(cnv)
snv = rm_corr(snv)

gc()

degMETH = rm_corr(degMETH)
degMRNA = rm_corr(degMRNA)
degMIRNA = rm_corr(degMIRNA)
rppa = rm_corr(rppa)

gc()

cnv_snv = rm_corr(cnv_snv)
cnv_meth = rm_corr(cnv_meth)
cnv_mRNA = rm_corr(cnv_mRNA)
cnv_miRNA = rm_corr(cnv_miRNA)
cnv_rppa = rm_corr(cnv_rppa)

gc()

snv_meth = rm_corr(snv_meth)
snv_mRNA = rm_corr(snv_mRNA)

gc()

snv_miRNA = rm_corr(snv_miRNA)
snv_rppa = rm_corr(snv_rppa)

gc()

meth_mRNA = rm_corr(meth_mRNA)
meth_miRNA = rm_corr(meth_miRNA)
meth_rppa = rm_corr(meth_rppa)
mRNA_miRNA = rm_corr(mRNA_miRNA)
mRNA_rppa = rm_corr(mRNA_rppa)
miRNA_rppa = rm_corr(miRNA_rppa)
gc()


cnv_snv_meth = rm_corr(cnv_snv_meth)
cnv_snv_mRNA = rm_corr(cnv_snv_mRNA)
gc()


cnv_snv_miRNA = rm_corr(cnv_snv_miRNA)
cnv_snv_rppa = rm_corr(cnv_snv_rppa)

gc()

cnv_meth_mRNA = rm_corr(cnv_meth_mRNA)
cnv_meth_miRNA = rm_corr(cnv_meth_miRNA)
cnv_meth_rppa = rm_corr(cnv_meth_rppa)
cnv_mRNA_miRNA = rm_corr(cnv_mRNA_miRNA)
cnv_mRNA_rppa = rm_corr(cnv_mRNA_rppa)
cnv_miRNA_rppa = rm_corr(cnv_miRNA_rppa)
gc()


snv_meth_mRNA = rm_corr(snv_meth_mRNA)
snv_meth_miRNA = rm_corr(snv_meth_miRNA)
gc()


snv_meth_rppa = rm_corr(snv_meth_rppa)
snv_mRNA_miRNA = rm_corr(snv_mRNA_miRNA)
snv_mRNA_rppa = rm_corr(snv_mRNA_rppa)


gc()


snv_miRNA_rppa = rm_corr(snv_miRNA_rppa)
meth_mRNA_miRNA = rm_corr(meth_mRNA_miRNA)
gc()


meth_mRNA_rppa = rm_corr(meth_mRNA_rppa)
meth_miRNA_rppa = rm_corr(meth_miRNA_rppa)
mRNA_miRNA_rppa = rm_corr(mRNA_miRNA_rppa)
gc()


cnv_snv_meth_mRNA = rm_corr(cnv_snv_meth_mRNA)
cnv_snv_meth_miRNA = rm_corr(cnv_snv_meth_miRNA)

gc()

cnv_snv_meth_rppa = rm_corr(cnv_snv_meth_rppa)
cnv_snv_mRNA_miRNA = rm_corr(cnv_snv_mRNA_miRNA)
cnv_snv_mRNA_rppa = rm_corr(cnv_snv_mRNA_rppa)
gc()


cnv_snv_miRNA_rppa = rm_corr(cnv_snv_miRNA_rppa)
cnv_meth_mRNA_miRNA = rm_corr(cnv_meth_mRNA_miRNA)
cnv_meth_mRNA_rppa = rm_corr(cnv_meth_mRNA_rppa)

gc()
save.image("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\TCGA-PRAD-workspace3.RData") #save workspace as a fail-safe in case of system failure

cnv_meth_miRNA_rppa = rm_corr(cnv_meth_miRNA_rppa)
cnv_mRNA_miRNA_rppa = rm_corr(cnv_mRNA_miRNA_rppa)
snv_meth_mRNA_miRNA = rm_corr(snv_meth_mRNA_miRNA)
gc()


save.image("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\TCGA-PRAD-workspace4.RData")

snv_meth_mRNA_rppa = rm_corr(snv_meth_mRNA_rppa)
snv_meth_miRNA_rppa = rm_corr(snv_meth_miRNA_rppa)
gc()

save.image("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\TCGA-PRAD-workspace5.RData")
snv_mRNA_miRNA_rppa = rm_corr(snv_mRNA_miRNA_rppa)
meth_mRNA_miRNA_rppa = rm_corr(meth_mRNA_miRNA_rppa)
cnv_snv_meth_mRNA_miRNA = rm_corr(cnv_snv_meth_mRNA_miRNA)

gc()
save.image("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\TCGA-PRAD-workspace6.RData")




cnv_snv_meth_mRNA_rppa = rm_corr(cnv_snv_meth_mRNA_rppa)
gc()


cnv_snv_meth_miRNA_rppa = rm_corr(cnv_snv_meth_miRNA_rppa)
gc()


cnv_snv_mRNA_miRNA_rppa = rm_corr(cnv_snv_mRNA_miRNA_rppa)
gc()
save.image("C:\\Users\\FEDGEN-LAB1-SYS7\\Desktop\\multiomics\\Data\\TCGA-PRAD-workspace7.RData")

cnv_meth_mRNA_miRNA_rppa = rm_corr(cnv_meth_mRNA_miRNA_rppa)
gc()


snv_meth_mRNA_miRNA_rppa = rm_corr(snv_meth_mRNA_miRNA_rppa)
gc()


cnv_snv_meth_mRNA_miRNA_rppa = rm_corr(cnv_snv_meth_mRNA_miRNA_rppa)
gc()

# check proportions -------------------------------------------------------

afterCORR = c()
omics = list(cnv,snv, degMETH,degMRNA, degMIRNA, rppa ,cnv_snv, cnv_meth,
             cnv_mRNA, cnv_miRNA,cnv_rppa, snv_meth, snv_mRNA,snv_miRNA, snv_rppa,
             meth_mRNA,meth_miRNA, meth_rppa, mRNA_miRNA, mRNA_rppa,
             miRNA_rppa,cnv_snv_meth, cnv_snv_mRNA,cnv_snv_miRNA,cnv_snv_rppa,
             cnv_meth_mRNA, cnv_meth_miRNA, cnv_meth_rppa, cnv_mRNA_miRNA,
             cnv_mRNA_rppa,cnv_miRNA_rppa, snv_meth_mRNA, snv_meth_miRNA,
             snv_meth_rppa, snv_mRNA_miRNA , snv_mRNA_rppa ,snv_miRNA_rppa,
             meth_mRNA_miRNA, meth_mRNA_rppa, meth_miRNA_rppa, mRNA_miRNA_rppa,
             cnv_snv_meth_mRNA, cnv_snv_meth_miRNA, cnv_snv_meth_rppa, 
             cnv_snv_mRNA_miRNA, cnv_snv_mRNA_rppa, cnv_snv_miRNA_rppa,
             cnv_meth_mRNA_miRNA, cnv_meth_mRNA_rppa, cnv_meth_miRNA_rppa,
             cnv_mRNA_miRNA_rppa, snv_meth_mRNA_miRNA, snv_meth_mRNA_rppa,
             snv_meth_miRNA_rppa, snv_mRNA_miRNA_rppa, meth_mRNA_miRNA_rppa, 
             cnv_snv_meth_mRNA_miRNA, cnv_snv_meth_mRNA_rppa, cnv_snv_meth_miRNA_rppa,
             cnv_snv_mRNA_miRNA_rppa, cnv_meth_mRNA_miRNA_rppa, snv_meth_mRNA_miRNA_rppa,
             cnv_snv_meth_mRNA_miRNA_rppa)
length(omics)

for (omic in omics) {
  nCNV = length(omic[grepl("_cnv",colnames(omic))])
  nSNV = length(omic[grepl("_snv",colnames(omic))])
  nMETH = length(omic[grepl("_meth",colnames(omic))])
  nMRNA = length(omic[grepl("_mRNA",colnames(omic))])
  nMIRNA = length(omic[grepl("_miRNA",colnames(omic))])
  nRPPA = length(omic[grepl("_rppa",colnames(omic))])
  afterCORR = rbind(afterCORR,c(nCNV,nSNV,nMETH,nMRNA,nMIRNA,nRPPA))
}

row.names(afterCORR) = omics_names

write.csv(afterCORR,"AfterCorrelation.csv", row.names = TRUE)




# Add Phenotypes ----------------------------------------------------------
commonsamplePhenotypes = phenotype[which(phenotype$submitter_id.samples %in% commonSamples), ]
commonsamplePhenotypes <-commonsamplePhenotypes[match(rownames(degMRNA),commonsamplePhenotypes$submitter_id.samples), ]
commonsamplePhenotypes$primary_therapy_outcome_success
table(commonsamplePhenotypes$primary_therapy_outcome_success)


cnv$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv), commonsamplePhenotypes$submitter_id.samples)]
snv$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv), commonsamplePhenotypes$submitter_id.samples)]
degMETH$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(degMETH), commonsamplePhenotypes$submitter_id.samples)]
degMRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(degMRNA), commonsamplePhenotypes$submitter_id.samples)]
degMIRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(degMIRNA), commonsamplePhenotypes$submitter_id.samples)]
rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth), commonsamplePhenotypes$submitter_id.samples)]
cnv_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_mRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_meth$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth), commonsamplePhenotypes$submitter_id.samples)]
snv_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_mRNA), commonsamplePhenotypes$submitter_id.samples)]
snv_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_miRNA), commonsamplePhenotypes$submitter_id.samples)]
snv_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_rppa), commonsamplePhenotypes$submitter_id.samples)]
meth_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_mRNA), commonsamplePhenotypes$submitter_id.samples)]
meth_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_miRNA), commonsamplePhenotypes$submitter_id.samples)]
meth_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_rppa), commonsamplePhenotypes$submitter_id.samples)]
mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_mRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_mRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_mRNA), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_miRNA), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
snv_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
meth_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
meth_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
meth_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_mRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_mRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_rppa), commonsamplePhenotypes$submitter_id.samples)] 
cnv_snv_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
meth_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(meth_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_mRNA_miRNA$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_mRNA_miRNA), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_mRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_mRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_meth_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_meth_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
snv_meth_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(snv_meth_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]
cnv_snv_meth_mRNA_miRNA_rppa$Outcome = commonsamplePhenotypes$primary_therapy_outcome_success[match(rownames(cnv_snv_meth_mRNA_miRNA_rppa), commonsamplePhenotypes$submitter_id.samples)]

setwd("C:/Users/FEDGEN-LAB1-SYS7/Desktop/multiomics/Data")

# Save Processed Files ----------------------------------------------------
file_names<-paste(omics_names,"csv",sep=".") 
write.csv(cnv,paste(".\\processed\\" , file_names[1]), row.names = TRUE)
write.csv(snv,paste(".\\processed\\" , file_names[2]), row.names = TRUE)
write.csv(degMETH,paste(".\\processed\\" , file_names[3]), row.names = TRUE)
write.csv(degMRNA,paste(".\\processed\\" , file_names[4]), row.names = TRUE)
write.csv(degMIRNA,paste(".\\processed\\" , file_names[5]), row.names = TRUE)
write.csv(rppa,paste(".\\processed\\" , file_names[6]), row.names = TRUE)
write.csv(cnv_snv,paste(".\\processed\\" , file_names[7]), row.names = TRUE)
write.csv(cnv_meth,paste(".\\processed\\" , file_names[8]), row.names = TRUE)
write.csv(cnv_mRNA,paste(".\\processed\\" , file_names[9]), row.names = TRUE)
write.csv(cnv_miRNA,paste(".\\processed\\" , file_names[10]), row.names = TRUE)
write.csv(cnv_rppa,paste(".\\processed\\" , file_names[11]), row.names = TRUE)
write.csv(snv_meth,paste(".\\processed\\" , file_names[12]), row.names = TRUE)
write.csv(snv_mRNA,paste(".\\processed\\" , file_names[13]), row.names = TRUE)
write.csv(snv_miRNA,paste(".\\processed\\" , file_names[14]), row.names = TRUE)
write.csv(snv_rppa,paste(".\\processed\\" , file_names[15]), row.names = TRUE)
write.csv(meth_mRNA,paste(".\\processed\\" , file_names[16]), row.names = TRUE)
write.csv(meth_miRNA,paste(".\\processed\\" , file_names[17]), row.names = TRUE)
write.csv(meth_rppa,paste(".\\processed\\" , file_names[17]), row.names = TRUE)
write.csv(mRNA_miRNA,paste(".\\processed\\" , file_names[19]), row.names = TRUE)
write.csv(mRNA_rppa,paste(".\\processed\\" , file_names[20]), row.names = TRUE)
write.csv(miRNA_rppa,paste(".\\processed\\" , file_names[21]), row.names = TRUE)
write.csv(cnv_snv_meth,paste(".\\processed\\" , file_names[22]), row.names = TRUE)
write.csv(cnv_snv_mRNA,paste(".\\processed\\" , file_names[23]), row.names = TRUE)
write.csv(cnv_snv_miRNA,paste(".\\processed\\" , file_names[24]), row.names = TRUE)
write.csv(cnv_snv_rppa,paste(".\\processed\\" , file_names[25]), row.names = TRUE)
write.csv(cnv_meth_mRNA,paste(".\\processed\\" , file_names[26]), row.names = TRUE)
write.csv(cnv_meth_miRNA,paste(".\\processed\\" , file_names[27]), row.names = TRUE)
write.csv(cnv_meth_rppa,paste(".\\processed\\" , file_names[28]), row.names = TRUE)
write.csv(cnv_mRNA_miRNA,paste(".\\processed\\" , file_names[29]), row.names = TRUE)
write.csv(cnv_mRNA_rppa,paste(".\\processed\\" , file_names[30]), row.names = TRUE)
write.csv(cnv_miRNA_rppa,paste(".\\processed\\" , file_names[31]), row.names = TRUE)
write.csv(snv_meth_mRNA,paste(".\\processed\\" , file_names[32]), row.names = TRUE)
write.csv(snv_meth_miRNA,paste(".\\processed\\" , file_names[33]), row.names = TRUE)
write.csv(snv_meth_rppa,paste(".\\processed\\" , file_names[34]), row.names = TRUE)
write.csv(snv_mRNA_miRNA,paste(".\\processed\\" , file_names[35]), row.names = TRUE)
write.csv(snv_mRNA_rppa,paste(".\\processed\\" , file_names[36]), row.names = TRUE)
write.csv(snv_miRNA_rppa,paste(".\\processed\\" , file_names[37]), row.names = TRUE)
write.csv(meth_mRNA_miRNA,paste(".\\processed\\" , file_names[38]), row.names = TRUE)
write.csv(meth_mRNA_rppa,paste(".\\processed\\" , file_names[39]), row.names = TRUE)
write.csv(meth_miRNA_rppa,paste(".\\processed\\" , file_names[40]), row.names = TRUE)
write.csv(mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[41]), row.names = TRUE)
write.csv(cnv_snv_meth_mRNA,paste(".\\processed\\" , file_names[42]), row.names = TRUE)
write.csv(cnv_snv_meth_miRNA,paste(".\\processed\\" , file_names[43]), row.names = TRUE)
write.csv(cnv_snv_meth_rppa,paste(".\\processed\\" , file_names[44]), row.names = TRUE)
write.csv(cnv_snv_mRNA_miRNA,paste(".\\processed\\" , file_names[45]), row.names = TRUE)
write.csv(cnv_snv_mRNA_rppa,paste(".\\processed\\" , file_names[46]), row.names = TRUE)
write.csv(cnv_snv_miRNA_rppa,paste(".\\processed\\" , file_names[47]), row.names = TRUE)
write.csv(cnv_meth_mRNA_miRNA,paste(".\\processed\\" , file_names[48]), row.names = TRUE)
write.csv(cnv_meth_mRNA_rppa,paste(".\\processed\\" , file_names[49]), row.names = TRUE)
write.csv(cnv_meth_miRNA_rppa,paste(".\\processed\\" , file_names[50]), row.names = TRUE)
write.csv(cnv_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[51]), row.names = TRUE)
write.csv(snv_meth_mRNA_miRNA,paste(".\\processed\\" , file_names[52]), row.names = TRUE)
write.csv(snv_meth_mRNA_rppa,paste(".\\processed\\" , file_names[53]), row.names = TRUE)
write.csv(snv_meth_miRNA_rppa,paste(".\\processed\\" , file_names[54]), row.names = TRUE)
write.csv(snv_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[55]), row.names = TRUE)
write.csv(meth_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[56]), row.names = TRUE)
write.csv(cnv_snv_meth_mRNA_miRNA,paste(".\\processed\\" , file_names[57]), row.names = TRUE)
write.csv(cnv_snv_meth_mRNA_rppa,paste(".\\processed\\" , file_names[58]), row.names = TRUE)
write.csv(cnv_snv_meth_miRNA_rppa,paste(".\\processed\\" , file_names[59]), row.names = TRUE)
write.csv(cnv_snv_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[60]), row.names = TRUE)
write.csv(cnv_meth_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[61]), row.names = TRUE)
write.csv(snv_meth_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[62]), row.names = TRUE)
write.csv(cnv_snv_meth_mRNA_miRNA_rppa,paste(".\\processed\\" , file_names[63]), row.names = TRUE)

end_time <- Sys.time()

# Print the execution times
cat("Parallel execution time:",
    end_time - start_time, "\n")

