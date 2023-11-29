# prepare libraries and load
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(2))

parallel=T


#docker run -v /shared/eearley/vidus:/vidus/ -i -t rtibiocloud/deseq2:1.22.2
setwd("vidus")

print("loading RData file")
load("./VIDUS_HIV_DTE_coc_deseq2_2023_11_28.RData")
load("./VIDUS_HIV_DTE_noncoc_deseq2_2023_11_28.RData")
print("finished loading RData file")



######################################
#### Differential Gene Expression ####
######################################


######################################
#### Differential Gene Expression ####
######################################
common.vars <- c('female', 'age.bin',
    "RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
    "cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")
contrast.var <- "hiv"
vidus.full.formula <- paste0("~", paste0(c(contrast.var, common.vars), collapse=" + "))



# cocaine only
print("Model fitting")
design(filtered.vidus.coc.tx.dds) <- as.formula(vidus.full.formula)
design(filtered.vidus.coc.tx.dds)
vidus.fit.coc <- DESeq(filtered.vidus.coc.tx.dds, test = "Wald",
         fitType = "parametric", sfType = "ratio", betaPrior = F,
         parallel = parallel)
resultsNames(vidus.fit.coc)
save(vidus.fit.coc, file="./model.fit.dte.coc.2023_11_29.RData")


# no cocaine only
print("Model fitting")
design(filtered.vidus.nococ.tx.dds) <- as.formula(vidus.full.formula)
vidus.fit.nococ <- DESeq(filtered.vidus.nococ.tx.dds, test = "Wald",
         fitType = "parametric", sfType = "ratio", betaPrior = F,
         parallel = parallel)
resultsNames(vidus.fit.nococ)
save(vidus.fit.novl, file="./model.fit.dte.nococ.2023_11_29.RData")





##########################
#### apeGLM shrinkage ####
##########################

# cocaine only
#apply apeGLM shrinkage to fold changes
print("apeGLM shrinkage")
hiv.results.coc <- DESeq2::results(vidus.fit.coc, 
                               name = "hiv_1_vs_0", 
                               alpha = 0.05, 
                               cooksCutoff = Inf)
hiv.shrunk.results.coc <- lfcShrink(vidus.fit.coc, 
                                res = hiv.results.coc, 
                                coef = "hiv_1_vs_0", 
                                type = "apeglm", 
                                parallel = parallel)
# save final output
save("hiv.shrunk.results.coc",file="./hiv.shrunk.dte.results.coc.2023_11_29.rda")




# no VL
print("apeGLM shrinkage")
hiv.results.nococ <- DESeq2::results(vidus.fit.nococ, 
                               name = "hiv_1_vs_0", 
                               alpha = 0.05, 
                               cooksCutoff = Inf)
hiv.shrunk.results.nococ <- lfcShrink(vidus.fit.nococ, 
                                res = hiv.results.nococ, 
                                coef = "hiv_1_vs_0", 
                                type = "apeglm", 
                                parallel = parallel)
# save final output
save("hiv.shrunk.results.novl",file="./hiv.shrunk.dte.results.nococ.2023_11_28.rda")
