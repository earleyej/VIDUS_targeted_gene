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
library(BiocParallel)
options(MulticoreParam=MulticoreParam(workers=6))


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

# restrict to just the targeted genes for the HIV paper
pull <- c("CD44",
          "CX3CR1",
          "EPSTI1",
          "HCP5",
          "IFI44L",
          "IFIT3",
          "IFITM1",
          "MX1",
          "NLRC5",
          "PARP9",
          "PLSCR1",
          "RASSF3",
          "RIN2",
          "RSAD2",
          "TAP1",
          "TAP2",
          "TNF",
          "TNIP3")
txInfo = txInfo[txInfo$GENENAME %in% pull,]
dim(tx) #60K
tx = tx[rownames(tx) %in% txInfo$TRANSCRIPTID,]
dim(tx) #119

# Remove genes that only have 1 transcript
txInfo <- txInfo[txInfo$TRANSCRIPTID %in% rownames(tx), ]
dim(txInfo) #119
txInfo <- subset(txInfo, 
                 duplicated(GENENAME) | duplicated(GENENAME, fromLast = TRUE))
dim(txInfo) #116
tx <- tx[which(
  rownames(tx) %in% txInfo$TRANSCRIPTID), ]
dim(tx) #116

# make the txInfo and tx colinear
txInfo <- txInfo[match(rownames(tx), txInfo$TRANSCRIPTID), ]



# Subset the pheno file columns
# convert age to 5 year bins
pheno$age.bin = cut(pheno$ageatint,breaks=5)
# Age, sex, RNA quality (RIN), 5 deconvolution-derived cell class proportions, and the top 5 Principal Components
pheno.analysis = pheno[,c("iid","hiv","female","age.bin","RNA_Quality_Score","cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8",paste0("PC",c(1:5)))]
#### TESTING ###
pheno.analysis$hiv = ifelse(pheno.analysis$hiv == 0,"N","Y")
#### END ####


# create design matrix from the pheno file
# All three main functions of satuRn require a SummarizedExperiment object as an input class
# also rename the columns in txInfo
colnames(txInfo) = c("isoform_id","gene_id")
sumExp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = tx),
    colData = pheno.analysis,
    rowData = txInfo
)
#metadata(sumExp)$formula <- ~ 0 + hiv + female + as.factor(colData(sumExp)$age.bin) + RNA_Quality_Score + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8 + PC1 + PC2 + PC3 + PC4 + PC5
metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$hiv)


# Fit quasi-binomial generalized linear model
# ~5 min. for the full transcriptome
system.time({
sumExp <- satuRn::fitDTU(
    object = sumExp,   
    parallel = FALSE,
    formula = ~ 0 + hiv,
    BPPARAM = BiocParallel::bpparam(),
    verbose = TRUE
)
})
#formula = ~ 0 + as.factor(hiv) + female + as.factor(colData(sumExp)$age.bin) + RNA_Quality_Score + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8 + PC1 + PC2 + PC3 + PC4 + PC5,


# set up the contrast matrix
# using limma
group <- as.factor(pheno.analysis$hiv)
design <- model.matrix(~ 0 + group) # construct design matrix
#colnames(design) <- c(0,1)
L <- limma::makeContrasts(
    Contrast1 = "groupN - groupY",
    levels = design
)

#manually
#L = data.frame("Contrast1" = c(-1,1))
#rownames(L) = c(0,1)


# Perform the contrast
sumExp <- satuRn::testDTU(
    object = sumExp,
    contrasts = L,
    diagplot1 = TRUE,
    diagplot2 = TRUE,
    sort = FALSE
)


# save spot
save.image("/vidus/DTU_saturn_2024_05_02.RData")



# Results
head(rowData(sumExp)[["fitDTUResult_Contrast1"]])

#












## Testing below using smaller subsets of data
# 1,000 rows = 37 seconds
# 5,000 rows = 165 seconds
# 20,000 rows = 642 seconds (10 min)
tx.sub = tx[c(1:20000),]
txInfo.sub = txInfo[c(1:20000),]
# remove lonely transcripts
txInfo.sub <- subset(txInfo.sub, 
                 duplicated(gene_id) | duplicated(gene_id, fromLast = TRUE))
tx.sub <- tx[which(
  rownames(tx.sub) %in% txInfo.sub$isoform_id), ]

# create object
sumExp.sub <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = tx.sub),
    colData = pheno.analysis,
    rowData = txInfo.sub
)
metadata(sumExp.sub)$formula <- ~ 0 + hiv 
# fit simple model
sumExp.sub <- satuRn::fitDTU(
    object = sumExp.sub,
    formula = ~ 0 + hiv,
    parallel = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    verbose = TRUE
)
# perform contrast
L = contr.treatment(pheno.analysis$hiv)
system.time({
sumExp.sub <- satuRn::testDTU(
    object = sumExp.sub,
    contrasts = L,
    diagplot1 = TRUE,
    diagplot2 = TRUE,
    sort = FALSE
)
})
