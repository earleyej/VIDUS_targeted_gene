# Eric J. Earley, PhD
# RTI International
# April 2024

# This code will run differential transcript usage analysis using the R package 'satuRn' (https://github.com/statOmics/satuRn)
# vignette: https://statomics.github.io/satuRn/articles/Vignette.html


# This code runs in 3 chunks
# 1. fitDTU - fit a quasi-binomial generalized linear model for transcript usage profiles
# 2. testDTU - tests for differential transcript usage between groups
# 3. plotDTU - plotting


# I made a docker image that pre-loads all of these packages: 
# docker pull rtibiocloud/saturn:v1.4.0_2e333c2
#	docker run -v /shared/eearley/vidus:/vidus/ -it rtibiocloud/saturn:v1.4.0_2e333c2

suppressMessages(library(satuRn))
suppressMessages(library(AnnotationHub))
suppressMessages(library(ensembldb))
suppressMessages(library(edgeR))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(ggplot2))
suppressMessages(library(DEXSeq))
suppressMessages(library(stageR))



#### Load Data ####
# This object contains the filtered.vidus.tx.dds data
load("/vidus/VIDUS_HIV_DTE_deseq2_2023_11_28.RData")
pheno = read.table("/vidus/master.pheno.drugs.5cellTypes.txt", sep="\t", header=T, stringsAsFactors=F)


# get count data
tx = assay(filtered.vidus.tx.dds)

# this has 61k rows
# it has already been filtered for low-count transcripts
# see code 03a_DTE_prepare_data.R



# load transcript - gene table that I made in 03a_DTE_prepare_data.R
txInfo = readRDS("/vidus/gencode_v28_tx_id_name_map.rds")
# remove version number from transcript ID
txInfo = data.frame(txInfo)
txInfo$TRANSCRIPTID = gsub("\\..*$","",txInfo$TRANSCRIPTID)
rownames(tx) = gsub("\\..*$","",rownames(tx))


# Remove genes that only have 1 transcript
txInfo <- txInfo[txInfo$TRANSCRIPTID %in% rownames(tx), ]
dim(txInfo) #61,523
txInfo <- subset(txInfo, 
                 duplicated(GENENAME) | duplicated(GENENAME, fromLast = TRUE))
dim(txInfo) #53,406
tx <- tx[which(
  rownames(tx) %in% txInfo$TRANSCRIPTID), ]
dim(tx) #53,320

# make the txInfo and tx colinear
txInfo <- txInfo[match(rownames(tx), txInfo$TRANSCRIPTID), ]

# Subset the pheno file 
# convert age to 5 year bins
pheno$age.bin = cut(pheno$ageatint,breaks=5)
# Age, sex, RNA quality (RIN), 5 deconvolution-derived cell class proportions, and the top 5 Principal Components
pheno.analysis = pheno[,c("iid","hiv","female","age.bin","RNA_Quality_Score","cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8",paste0("PC",c(1:5)))]

# create design matrix from the pheno file
# All three main functions of satuRn require a SummarizedExperiment object as an input class
# also rename the columns in txInfo
colnames(txInfo) = c("isoform_id","gene_id")
sumExp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = tx),
    colData = pheno.analysis,
    rowData = txInfo
)
metadata(sumExp)$formula <- ~ 0 + hiv + female + as.factor(colData(sumExp)$age.bin) + RNA_Quality_Score + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8 + PC1 + PC2 + PC3 + PC4 + PC5




# Fit quasi-binomial generalized linear model
# ~5 min.
system.time({
sumExp <- satuRn::fitDTU(
    object = sumExp,
    formula = ~ 0 + hiv + female + as.factor(colData(sumExp)$age.bin) + RNA_Quality_Score + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8 + PC1 + PC2 + PC3 + PC4 + PC5,
    parallel = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    verbose = TRUE
)
})



# set up the contrast matrix
design <- model.matrix(~ 0 + hiv + female + as.factor(colData(sumExp)$age.bin) + RNA_Quality_Score + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8 + PC1 + PC2 + PC3 + PC4 + PC5, 
                       data=pheno.analysis) 



