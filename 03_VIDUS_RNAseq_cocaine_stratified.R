# prepare libraries and load
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(4))

parallel=T


#docker run -v /rti-01/eearley/vidus:/vidus/ -i -t rtibiocloud/deseq2:1.22.2
setwd("vidus")

print("loading RData file")
load("./VIDUS_HIV_DGE_deseq2_cocaineStratified_2023_11_14.RData")
print("finished loading RData file")



######################################
#### Differential Gene Expression ####
######################################


common.vars <- c('female', 'ageatint',
    "RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
    "cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")
contrast.var <- "hiv"
vidus.full.formula <- paste0("~", paste0(c(contrast.var, common.vars), collapse=" + "))


# Cocaine only
print("Model fitting cocaine subset")
design(filtered.vidus.coc.gene.dds) <- as.formula(vidus.full.formula)
design(filtered.vidus.coc.gene.dds)
# perform the regression
vidus.fit.coc <- DESeq(filtered.vidus.coc.gene.dds, test = "Wald",
        fitType = "parametric", sfType = "ratio", betaPrior = F,
        parallel = parallel)
resultsNames(vidus.fit.coc)
hiv.results.coc <- DESeq2::results(vidus.fit.coc, name = "hiv_1_vs_0", alpha = 0.05, cooksCutoff = Inf)

# Non-cocaine only
print("Model fitting non-cocaine subset")
design(filtered.vidus.noncoc.gene.dds) <- as.formula(vidus.full.formula)
design(filtered.vidus.noncoc.gene.dds)
# perform the regression
vidus.fit.noncoc <- DESeq(filtered.vidus.noncoc.gene.dds, test = "Wald",
        fitType = "parametric", sfType = "ratio", betaPrior = F,
        parallel = parallel)
resultsNames(vidus.fit.noncoc)
hiv.results.noncoc <- DESeq2::results(vidus.fit.noncoc, name = "hiv_1_vs_0", alpha = 0.05, cooksCutoff = Inf)

##########################
#### apeGLM shrinkage ####
##########################
#apply apeGLM shrinkage to fold changes
# this takes a long time, consider parallelizing
print("apeGLM shrinkage for cocaine subset")
hiv.shrunk.results <- lfcShrink(vidus.fit.coc, res = hiv.results.coc, 
    coef = "hiv_1_vs_0", type = "apeglm", parallel = parallel)

print("apeGLM shrinkage for non-cocaine subset")
hiv.shrunk.results <- lfcShrink(vidus.fit.noncoc, res = hiv.results.noncoc, 
    coef = "hiv_1_vs_0", type = "apeglm", parallel = parallel)

#### save final output ####
save("hiv.shrunk.results",file="./hiv.shrunk.dge.cocaineONLY.results.rda")


