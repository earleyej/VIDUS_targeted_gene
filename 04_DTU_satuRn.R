# Eric J. Earley, PhD
# RTI International
# April 2024

# This code will run differential transcript usage analysis using the R package 'satuRn' (https://github.com/statOmics/satuRn)
# vignette: https://statomics.github.io/satuRn/articles/Vignette.html


# This code runs in 3 chunks
# 1. fitDTU - fit a quasi-binomial generalized linear model for transcript usage profiles
# 2. testDTU - tests for differential transcript usage between groups
# 3. plotDTU - plotting


#### Environment Setup ####
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("satuRn")


library(satuRn)
library(AnnotationHub)
library(ensembldb)
library(edgeR)
library(SummarizedExperiment)
library(ggplot2)
library(DEXSeq)
library(stageR)



#### Load Data ####


