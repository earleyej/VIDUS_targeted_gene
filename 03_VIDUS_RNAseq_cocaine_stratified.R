# prepare libraries and load
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

#docker run -v /rti-01/eearley/vidus:/vidus/ -i -t rtibiocloud/deseq2:1.22.2
setwd("vidus")


load("./VIDUS_HIV_DGE_deseq2_dds_2023_11_13.RData")




######################################
#### Differential Gene Expression ####
######################################

# DGE model - BASE + absolute.cell.proportions
common.vars <- c('female', 'ageatint',
    "RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
    "cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")
contrast.var <- "hiv"
vidus.full.formula <- paste0("~", paste0(c(contrast.var, common.vars), collapse=" + "))
design(filtered.vidus.gene.dds) <- as.formula(vidus.full.formula)
design(filtered.vidus.gene.dds)

# perform the regression
vidus.fit <- DESeq(filtered.vidus.gene.dds, test = "Wald",
        fitType = "parametric", sfType = "ratio", betaPrior = F,
        parallel = F)

# Access results
resultsNames(vidus.fit)
hiv.results <- DESeq2::results(vidus.fit, name = "hiv_1_vs_0", alpha = 0.05, cooksCutoff = Inf)



##########################
#### apeGLM shrinkage ####
##########################
#apply apeGLM shrinkage to fold changes
# this takes a long time, consider parallelizing
hiv.shrunk.results <- hiv.results
hiv.shrunk.results <- lfcShrink(vidus.fit, res = hiv.shrunk.results, 
    coef = "hiv_1_vs_0", type = "apeglm", parallel = F)



#### save final output ####
save("hiv.shrunk.results",file="./hiv.shrunk.dge.cocaineONLY.results.rda")


