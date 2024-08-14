library(DESeq2)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(stringr)
setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/Eric Johnson/VIDUS/candidate genes study/Results/")
#load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_input.rda")
load("R_data_files/VIDUS_HIV_DGE_deseq2_2023_12_13.RData")


# load phenotype
pheno = read.table("../Data/master.pheno.drugs.5cellTypes.txt", sep="\t", header=T, stringsAsFactors = F)
pheno$age.bin = cut(pheno$ageatint,breaks=5)

# load genes and ENSG IDs
gtf <- read.table(gzfile("../Data/gencode.v28.annotation.gtf.gz"), sep="\t", header=F, stringsAsFactors = F, quote="\"")
gtf$gene_name = str_match(gtf$V9, "gene_name (.*?);")[,2]
gtf$transcript = str_match(gtf$V9, "transcript_id (.*?);")[,2] # first row for a gene does not have a transcriptID. Think of the 1st row like a gene-specific annotation and all subsequent rows as transcript
gtf$transcript_type = str_match(gtf$V9, "transcript_type (.*?);")[,2]
gtf$ensg = str_match(gtf$V9, "gene_id (.*?);")[,2]
pull <- c("CD44","CX3CR1","EPSTI1","HCP5","IFI44L","IFIT3","IFITM1","MX1","NLRC5","PARP9","PLSCR1","RASSF3","RIN2","RSAD2","TAP1","TAP2","TNF","TNIP3")
vidus.gtf <- gtf[gtf$gene_name %in% pull,c("gene_name","ensg","transcript","transcript_type")]
vidus.gtf <- vidus.gtf[!duplicated(vidus.gtf),]
unique(vidus.gtf$gene_name)





#################################
#### PLOTTING w/ correlation ####
#################################


### LOOP HERE
p <- list()
corr.res = NULL
mod.res = NULL
for (i in 1:length(pull)) {
  g = unique(vidus.gtf[vidus.gtf$gene_name == pull[i],"gene_name"])
  print(g)
  ensg = unique(vidus.gtf[vidus.gtf$gene_name == pull[i],"ensg"])
  x<-plotCounts(filtered.vidus.gene.dds, 
                gene=ensg, 
                intgroup="hiv",
                returnData=T)
  rownames(x) = gsub("93-0*","",rownames(x))
  x.vl = x[rownames(x) %in% pheno[pheno$viral_suppressed == 0,"gwas_code"],]
  x.vl$gwas_code = rownames(x.vl)
  x.vl=merge(x.vl, 
             pheno[,c("gwas_code","proximal_log10_vl",'female', 'age.bin',"RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',"cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")],
             by="gwas_code")
  
  corr.res = rbind(corr.res,
                   data.frame("gene" = g,
                              "corr" = cor(x.vl$count,x.vl$proximal_log10_vl))
  )
  
  # model viral load
  #mod = lm(proximal_log10_vl ~ count + female + age.bin + RNA_Quality_Score + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8, data=x.vl)
  mod = lm(proximal_log10_vl ~ count, data=x.vl)
  print(summary(mod))
  mod.res = rbind(mod.res,
                  data.frame("gene" = g,
                             "beta" = summary(mod)$coefficient[2,1],
                             "P" = summary(mod)$coefficient[2,4]))
  
  
  # save predictions of the model in the new data frame 
  # together with variable you want to plot against
  predicted_df <- data.frame(pred = predict(mod, x.vl), 
                             count=x.vl$count)
  
  
  p[[g]] = ggplot(x.vl, aes(x=proximal_log10_vl, y=count )) +
    geom_point() +
    ggtitle(g) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_smooth(method = "lm", se = T) +
    geom_text(x=5, 
             y=0.8*max(x.vl$count), 
             label= paste0("P = ", round(mod.res$P[length(mod.res$P)], digits=2) ))
    #geom_line(color='red',data = predicted_df, aes(x=pred, y=count))
    # geom_text(x=5, 
    #           y=0.8*max(x.vl$count), 
    #           label= paste0("cor = ",format(round(cor(x.vl$count,x.vl$proximal_log10_vl),2),nsmall=2)))
  # stat_cor(aes(label=after_stat(rr.label)), 
  #          label.x=5, 
  #          label.y=0.8*max(x.vl$count)) +
  # 
  
}


corr.res
mod.res
n <- length(p)
nCol <- floor(sqrt(n))
pdf(file="vidus_counts_vs_VL_N33_2023_11_14.pdf",
    height=10,width=13)
do.call("grid.arrange", c(p, ncol=nCol))
dev.off()









# ######################
# #### Linear model ####
# ######################
# 
# i=1
# out = NULL
# 
# # UNIVARIATE
# for (i in 1:length(pull)) {
#   g = unique(vidus.gtf[vidus.gtf$gene_name == pull[i],"gene_name"])
#   print(g)
#   ensg = unique(vidus.gtf[vidus.gtf$gene_name == pull[i],"ensg"])
#   x<-plotCounts(filtered.vidus.gene.dds, 
#                 gene=ensg, 
#                 intgroup="hiv",
#                 returnData=T)
#   rownames(x) = gsub("93-0*","",rownames(x))
#   x$gwas_code = rownames(x)
#   x=merge(x, 
#           pheno[,c("gwas_code","proximal_log10_vl",'female', 'age.bin',"RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',"cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")],
#           by="gwas_code")
#   
#   mod.res = summary(lm(proximal_log10_vl ~ count, data=x))
#   res = mod.res$coefficients[rownames(mod.res$coefficients) == "count",]
#   out = rbind(out,
#               t(data.frame(c(g,res))))
# }
# out.univ = out
# 
# 
# # MULTIVARIATE
# out = NULL
# for (i in 1:length(pull)) {
#   g = unique(vidus.gtf[vidus.gtf$gene_name == pull[i],"gene_name"])
#   print(g)
#   ensg = unique(vidus.gtf[vidus.gtf$gene_name == pull[i],"ensg"])
#   x<-plotCounts(filtered.vidus.gene.dds, 
#                 gene=ensg, 
#                 intgroup="hiv",
#                 returnData=T)
#   rownames(x) = gsub("93-0*","",rownames(x))
#   x$gwas_code = rownames(x)
#   x=merge(x, 
#           pheno[,c("gwas_code","proximal_log10_vl",'female', 'age.bin',"RNA_Quality_Score", 'PC1', 'PC2', 'PC3', 'PC4', 'PC5',"cd4Tcells","Bcells","granulocytes","Monocytes","T.cells.CD8")],
#           by="gwas_code")
#   mod = lm(proximal_log10_vl ~ count + female + age.bin + RNA_Quality_Score + PC1 + PC2 + PC3 + PC4 + PC5 + cd4Tcells + Bcells + granulocytes + Monocytes + T.cells.CD8, data=x)
#   mod.res = summary(mod)
#   res = mod.res$coefficients[rownames(mod.res$coefficients) == "count",]
#   out = rbind(out,
#               t(data.frame(c(g,res))))
#   
#   # plot and add the regression line
#   mod.res = rbind(mod.res,
#                   data.frame("gene" = g,
#                              "beta" = summary(mod)$coefficient[2],
#                              "P" = summary(mod)$coefficient[8]))
#   
# }
# out.multiv = out
# out.multiv








