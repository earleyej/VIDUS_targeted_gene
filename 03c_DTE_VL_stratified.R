# prepare libraries and load
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(4))

parallel=T


#docker run -v /shared/eearley/vidus:/vidus/ -i -t rtibiocloud/deseq2:1.22.2
setwd("vidus")

print("loading RData file")
load("./VIDUS_HIV_DTE_vl_deseq2_2023_11_28.RData")
load("./VIDUS_HIV_DTE_novl_deseq2_2023_11_28.RData")

print("finished loading RData file")

######################################
#### Differential Gene Expression ####
######################################
common.vars <- c('female', 'age.bin',
    "RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
    "cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")
contrast.var <- "hiv"
vidus.full.formula <- paste0("~", paste0(c(contrast.var, common.vars), collapse=" + "))



# VL only
print("Model fitting")
design(filtered.vidus.vl.tx.dds) <- as.formula(vidus.full.formula)
design(filtered.vidus.vl.tx.dds)
vidus.fit.vl <- DESeq(filtered.vidus.vl.tx.dds, test = "Wald",
         fitType = "parametric", sfType = "ratio", betaPrior = F,
         parallel = parallel)
resultsNames(vidus.fit.vl)
save(vidus.fit.vl, file="./model.fit.dte.VL.2023_11_28.RData")


# no VL only
print("Model fitting")
design(filtered.vidus.novl.tx.dds) <- as.formula(vidus.full.formula)
vidus.fit.novl <- DESeq(filtered.vidus.novl.tx.dds, test = "Wald",
         fitType = "parametric", sfType = "ratio", betaPrior = F,
         parallel = parallel)
resultsNames(vidus.fit.novl)
save(vidus.fit.novl, file="./model.fit.dte.noVL.RData")





##########################
#### apeGLM shrinkage ####
##########################

# VL
#apply apeGLM shrinkage to fold changes
# this takes a long time, consider parallelizing
print("apeGLM shrinkage")
hiv.results.vl <- DESeq2::results(vidus.fit.vl, 
                               name = "hiv_1_vs_0", 
                               alpha = 0.05, 
                               cooksCutoff = Inf)
hiv.shrunk.results.vl <- lfcShrink(vidus.fit.vl, 
                                res = hiv.results.vl, 
                                coef = "hiv_1_vs_0", 
                                type = "apeglm", 
                                parallel = parallel)
# save final output
save("hiv.shrunk.results.vl",file="./hiv.shrunk.dte.results.vl.2023_11_28.rda")




# no VL
print("apeGLM shrinkage")
hiv.results.novl <- DESeq2::results(vidus.fit.novl, 
                               name = "hiv_1_vs_0", 
                               alpha = 0.05, 
                               cooksCutoff = Inf)
hiv.shrunk.results.novl <- lfcShrink(vidus.fit.novl, 
                                res = hiv.results.novl, 
                                coef = "hiv_1_vs_0", 
                                type = "apeglm", 
                                parallel = parallel)
# save final output
save("hiv.shrunk.results.novl",file="./hiv.shrunk.dte.results.novl.2023_11_28.rda")
