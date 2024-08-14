# this code will check for co-expression of downstream genes to the 18 candidates


library(ggplot2)
library(DESeq2)
library(gridExtra)
library(cowplot)
library(stringr)


setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/Eric Johnson/VIDUS/candidate genes study/Results/")

# load phenotype
pheno = read.table("../Data/master.pheno.drugs.5cellTypes.txt", sep="\t", header=T, stringsAsFactors = F)

# load genes and ENSG IDs
gtf <- read.table(gzfile("../Data/gencode.v28.annotation.gtf.gz"), sep="\t", header=F, stringsAsFactors = F, quote="\"")
gtf$gene_name = str_match(gtf$V9, "gene_name (.*?);")[,2]
gtf$transcript = str_match(gtf$V9, "transcript_id (.*?);")[,2] # first row for a gene does not have a transcriptID. Think of the 1st row like a gene-specific annotation and all subsequent rows as transcript
gtf$transcript_type = str_match(gtf$V9, "transcript_type (.*?);")[,2]
gtf$ensg = str_match(gtf$V9, "gene_id (.*?);")[,2]

# limiting to candidate genes
pull <- c("CD44","CX3CR1","EPSTI1","HCP5","IFI44L","IFIT3","IFITM1","MX1","NLRC5","PARP9","PLSCR1","RASSF3","RIN2","RSAD2","TAP1","TAP2","TNF","TNIP3")
vidus.gtf <- gtf[gtf$gene_name %in% pull,c("gene_name","ensg","transcript","transcript_type")]
vidus.gtf <- vidus.gtf[!duplicated(vidus.gtf),]
unique(vidus.gtf$gene_name)




# gene-level results
#load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_modelFit.rda") # contains 'vidus.fit'
load("R_data_files/hiv.shrunk.dge.results.2023_11_27.rda")



#### MHC Class 1 genes ####
mhc.table = read.table("../pathway_analysis/MHC_class1_genes.txt", sep="\t", header=T, stringsAsFactors = F)
mhc.class1 = mhc.table[mhc.table$Class == "class I",]
# remove pseudogenes
mhc.class1 = mhc.class1[!grepl("pseudogene",mhc.class1$Details),]
mhc.gtf = gtf[gtf$gene_name %in% mhc.class1$Gene, c("gene_name","ensg","transcript","transcript_type")]

get.results<-function(df=hiv.results,
                      prefix = "",
                      gtf = vidus.gft) {
  res<-as.data.frame(df[rownames(df) %in% unique(gtf$ensg),])
  res$ensg <- rownames(res)
  res <- merge(res,
               gtf[,c("ensg","gene_name")],by="ensg")
  row.names(res)<-NULL
  res <- res[!duplicated(res),]
  res<-res[,c("gene_name","ensg","log2FoldChange","pvalue")]
  #res$fdr <- p.adjust(res$p,method="fdr")
  res$bon <- p.adjust(res$p,method="bonferroni")
  colnames(res)[c(3,4,5)]<-c(paste0(prefix,"_LFC"),
                               paste0(prefix,"_p"),
                               paste0(prefix,"_bon"))
  colnames(res)[1]<-"gene"
  return(res)
}

res.base.cells <- get.results(prefix="base-cells",
                              df=hiv.shrunk.results,
                              gtf=mhc.gtf)
#add baseMean - this will be the same for all experiments
hiv.shrunk.results$ensg <- rownames(hiv.shrunk.results)
res.base.cells <- merge(res.base.cells, 
                        as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% mhc.gtf$ensg,c("ensg","baseMean")]), by="ensg")
mhc.results = res.base.cells[order(res.base.cells$`base-cells_bon`),]
head(mhc.results)
write.table(mhc.results, file="vidus_hiv_dge_coexpression_mhc.txt", sep="\t", col.names=T, row.names=F, quote=F)


#### NF-KB targets - TNF activation ####
nfkb = read.table("../Data/HALLMARK_TNFA_SIGNALING_VIA_NFKB_gene_list.txt", sep="\t", header=F, stringsAsFactors = F)
nfkb.gtf = gtf[gtf$gene_name %in% nfkb$V1, c("gene_name","ensg","transcript","transcript_type")]
dim(nfkb.gtf)
nfkb.gtf = nfkb.gtf[!duplicated(nfkb.gtf),]
dim(nfkb.gtf)

nfkb.res = get.results(prefix="base-cells",
                       df=hiv.shrunk.results,
                       gtf=nfkb.gtf)
nfkb.res = merge(nfkb.res,
                 as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% nfkb.gtf$ensg, c("ensg","baseMean")]), 
                 by="ensg")
dim(nfkb)
nfkb.res = nfkb.res[order(nfkb.res$`base-cells_bon`),]
write.table(nfkb.res, file="vidus_hiv_dge_coexpression_nfkb_targets_2023_12_08.txt",
            sep="\t", col.names=T, row.names=F, quote=F)
#



#### Interferon alpha activated genes ####
ifna = read.table("../Data/HALLMARK_INTERFERON_ALPHA_RESPONSE_gene_list.txt")
ifna.gtf = gtf[gtf$gene_name %in% ifna$V1, c("gene_name","ensg","transcript","transcript_type")]
dim(ifna.gtf)
ifna.gtf = ifna.gtf[!duplicated(ifna.gtf),]
dim(ifna.gtf)

ifna.res = get.results(prefix="base-cells",
                       df=hiv.shrunk.results,
                       gtf=ifna.gtf)
# remove the candidate genes from this list. Many of the candidates are activated by IFN-alpha signaling
dim(ifna.res)
candidates <- c("CD44","CX3CR1","EPSTI1","HCP5","IFI44L","IFIT3","IFITM1","MX1","NLRC5","PARP9","PLSCR1","RASSF3","RIN2","RSAD2","TAP1","TAP2","TNF","TNIP3")
ifna.res = ifna.res[!ifna.res$gene %in% candidates,]
dim(ifna.res)

hiv.shrunk.results$ensg <- rownames(hiv.shrunk.results)
ifna.res = merge(ifna.res,
                 as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% ifna.gtf$ensg, c("ensg","baseMean")]), 
                 by="ensg")
ifna.res = ifna.res[order(ifna.res$`base-cells_bon`),]
# reorder columns
ifna.res = ifna.res[,c("ensg","gene","baseMean","base-cells_LFC","base-cells_p","base-cells_bon")]
write.table(ifna.res, file="vidus_hiv_dge_coexpression_ifna_targets_2023_12_12.txt",
            sep="\t", col.names=T, row.names=F, quote=F)






#### NF-KB pathway ####
#### NOTE: This is UPSTREAM of NLRC5 and maybe the other candidates, as well. 
# https://www.gsea-msigdb.org/gsea/msigdb/cards/BIOCARTA_NFKB_PATHWAY
# from BIOCARTA_NFKB_PATHWAY from c2 mSigDB
# nfkb = read.table("../Data/nfkb_pathway_c2_mSigDB.txt",sep="\t",header=T, stringsAsFactors = F)
# nfkb.gtf = gtf[gtf$gene_name %in% nfkb$Gene, c("gene_name","ensg","transcript","transcript_type")]
# dim(nfkb.gtf)
# nfkb.res = get.results(prefix="base-cells",
#                        df=hiv.shrunk.results,
#                        gtf=nfkb.gtf)
# nfkb.res = merge(nfkb.res,
#                  as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% nfkb.gtf$ensg, c("ensg","baseMean")]), by="ensg")
# dim(nfkb)
# nfkb.res = nfkb.res[order(nfkb.res$`base-cells_bon`),]
# 
# write.table(nfkb.res, file="vidus_hiv_dge_coexpression_nfkb.txt", sep="\t", col.names=T, row.names=F, quote=F)
