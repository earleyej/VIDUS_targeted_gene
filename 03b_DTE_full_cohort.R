# prepare libraries and load
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(4))

parallel=T


#docker run -v /shared/eearley/vidus:/vidus/ -i -t rtibiocloud/deseq2:1.22.2


# INPUTS
#filtered.vidus.gene.dds
#



setwd("vidus")

print("loading RData file")
load("./VIDUS_HIV_DTE_deseq2_2023_11_28.RData")
print("finished loading RData file")



######################################
#### Differential Gene Expression ####
######################################


common.vars <- c('female', 'age.bin',
    "RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
    "cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")
contrast.var <- "hiv"
vidus.full.formula <- paste0("~", paste0(c(contrast.var, common.vars), collapse=" + "))


print("Model fitting")
design(filtered.vidus.tx.dds) <- as.formula(vidus.full.formula)
design(filtered.vidus.tx.dds)
# perform the regression
vidus.fit <- DESeq(filtered.vidus.tx.dds, test = "Wald",
        fitType = "parametric", sfType = "ratio", betaPrior = F,
        parallel = parallel)
resultsNames(vidus.fit)
hiv.results <- DESeq2::results(vidus.fit, name = "hiv_1_vs_0", alpha = 0.05, cooksCutoff = Inf)
save.image(file="./model.fit.dte.2023_11_28.RData")

##########################
#### apeGLM shrinkage ####
##########################
#apply apeGLM shrinkage to fold changes
# this takes a long time, consider parallelizing
print("apeGLM shrinkage")
hiv.shrunk.results <- lfcShrink(vidus.fit, res = hiv.results, 
    coef = "hiv_1_vs_0", type = "apeglm", parallel = parallel)


#### save final output ####
save("hiv.shrunk.results",file="./hiv.shrunk.dte.results.2023_11_28.rda")
