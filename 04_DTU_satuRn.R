# Eric J. Earley, PhD
# RTI International
# April 2024

# This code will run differential transcript usage analysis using the R package 'satuRn' (https://github.com/statOmics/satuRn)
# vignette: https://statomics.github.io/satuRn/articles/Vignette.html


# This code runs in 3 chunks
# 1. fitDTU - fit a quasi-binomial generalized linear model for transcript usage profiles
# 2. testDTU - tests for differential transcript usage between groups
# 3. plotDTU - plotting


# I made a docker image: 
# 	docker pull rtibiocloud/saturn:v1.4.0_2e333c2
#	docker run -v /shared/eearley/vidus:/vidus/ -it rtibiocloud/saturn:v1.4.0_2e333c2


library(satuRn)
library(AnnotationHub)
library(ensembldb)
library(edgeR)
library(SummarizedExperiment)
library(ggplot2)
library(DEXSeq)
library(stageR)



#### Load Data ####
# This object contains the filtered.vidus.tx.dds data
load("/vidus/VIDUS_HIV_DTE_deseq2_2023_11_28.RData")

# get count data
tx = assay(filtered.vidus.tx.dds)

# this has 61k rows




