#### Environment Setup ####

setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/Eric Johnson/VIDUS/candidate genes study/Results/")
library(ggplot2)
library(DESeq2)
library(gridExtra)
library(cowplot)
library(grid)
library(stringr)






#### Create gene/transcript annotation table ####
# gtf <- read.table(gzfile("C:/Users/eearley/Documents/Eric Johnson/VIDUS/candidate genes study/Data/gencode.v28.annotation.gtf.gz"),sep=";",fill=T,stringsAsFactors = F)
# head(gtf)
# gtf$gene_name <- gsub("gene_name ","",gtf$V4)
# gtf$gene_name <- gsub(" ","",gtf$gene_name)
# gtf$transcript <- gsub("transcript_id ","",gtf$V2)
# gtf$transcript <- gsub(" ","",gtf$transcript)
# gtf$transcript_type <- gsub("transcript_type ","",gtf$V5)
# gtf$transcript_type <- gsub(" ","",gtf$transcript_type)
# gtf$ensg <- gsub("^.+ ","",gtf$V1)
gtf <- read.table(gzfile("../Data/gencode.v28.annotation.gtf.gz"), sep="\t", header=F, stringsAsFactors = F, quote="\"")
head(gtf)
dim(gtf)
gtf$gene_name = str_match(gtf$V9, "gene_name (.*?);")[,2]
gtf$transcript = str_match(gtf$V9, "transcript_id (.*?);")[,2] # first row for a gene does not have a transcriptID. Think of the 1st row like a gene-specific annotation and all subsequent rows as transcript
gtf$transcript_type = str_match(gtf$V9, "transcript_type (.*?);")[,2]
gtf$ensg = str_match(gtf$V9, "gene_id (.*?);")[,2]



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


vidus.gtf <- gtf[gtf$gene_name %in% pull,c("gene_name","ensg","transcript","transcript_type")]
vidus.gtf <- vidus.gtf[!duplicated(vidus.gtf),]
unique(vidus.gtf$gene_name)


#### Fuctions ####

violin.plot <- function(gene=gene,
                        ens=ens,
                        dds=vidus.fit) { #filtered.vidus.gene.dds) {
  #custom margins (top,right,bottom,left)
  custom.margins<-unit(c(0.2,0.5,-1.2,0),"cm") 
  
  x<-plotCounts(dds, gene=ens, intgroup="hiv",returnData=T)
  x$count<-x$count/1000
  p1 <- ggplot(x, aes(x=hiv, y=count, fill=hiv)) + 
    geom_violin() + 
    #geom_jitter(shape=16, position=position_jitter(0.2))
    #stat_summary(fun=median, geom="point") +
    geom_boxplot(width=0.05) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.margin=custom.margins,
          axis.text = element_text(size=12),
          axis.text.x = element_text(size=16,face="bold"),
          axis.text.y = element_text(size=14),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size=16,face="bold"),
          legend.title = element_text(size=16,face="bold",hjust=0),
          legend.position = "none",
          axis.ticks.x = element_blank()) +
    scale_x_discrete(labels=c("","")) +
    scale_y_continuous(limits=c(0,(1.2*max(x$count)))) +
    ggtitle(gene) +
    ylab("") +
    xlab("") +
    scale_fill_discrete(name = "HIV Status", labels = c("-", "+"))  
  
  p1
  
  return(p1)
  
}

violin.plot.wVL.woVL <- function(gene=gene,
                        ens=ens,
                        dds=vidus.fit) { 
  #custom margins (top,right,bottom,left)
  custom.margins<-unit(c(0.2,0.7,-0.5,0),"cm") 
  
  #extract vl counts
  x<-plotCounts(dds, gene=ens, intgroup="hiv",returnData=T)
  x$count<-x$count/1000
  x <- cbind(x,samp.vl[,c("class")])
  colnames(x)[3]<-"class"
  
  p1 <- ggplot(x, aes(x=class, y=count, fill=class)) + 
    geom_violin() + 
    #geom_jitter(shape=16, position=position_jitter(0.2))
    #stat_summary(fun=median, geom="point") +
    geom_boxplot(width=0.05) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.margin=custom.margins,
          axis.text = element_text(size=12),
          axis.text.x = element_text(size=16,face="bold"),
          axis.text.y = element_text(size=14),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size=16,face="bold"),
          legend.title = element_text(size=16,face="bold",hjust=0),
          axis.ticks.x = element_blank()) +
    #legend.position = "none") +
    scale_x_discrete(labels=c("","","")) +
    #scale_y_continuous(labels=)
    ggtitle(gene) +
    ylab("") +
    xlab("") +
    scale_fill_discrete(name = "HIV Status", labels = c("Control", "No VL","Yes VL"))  
  
  p1
  
  return(p1)
  
}


#### GENE LEVEL ####
#control vs. HIV+ (N=588) all samples - 
# Results from model:
# SNP ~ HIV + sex + age + RIN + PC1-5 + 5-cell-type-prop.
load("R_data_files/hiv.shrunk.dge.results.2023_11_27.rda") # contains 'hiv.shrunk.results'

#load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_input.rda")
#load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_modelFit.rda") # contains 'vidus.fit'
load("R_data_files/VIDUS_HIV_DGE_deseq2_2023_12_13.RData") #filtered.vidus.gene.dds

p <- list()

# GENES - base count plots 
x <- vidus.gtf[,c("gene_name","ensg")]
gene.list <- x[!duplicated(x),]
for (i in 1:nrow(gene.list)) {
  g<-gene.list$gene[i]
  print(g)
  
  p[[g]]<-violin.plot(dds=filtered.vidus.gene.dds,#filtered.vidus.gene.dds,
                      gene=gene.list$gene[i],
                      ens=gene.list$ensg[i])
  legend<-get_legend(p[[g]])
  hiv.shrunk.dge.5cells.results.rda<-p[[g]] + theme(legend.position = "none")
}

pdf(file="../manuscript/figures/fig1.violin_plots_DGE.pdf",width=10,height=8)
grid.arrange(p$CD44,p$CX3CR1,p$EPSTI1,
             p$HCP5,p$IFI44L,
             p$IFIT3,p$IFITM1,p$MX1,
             p$NLRC5,p$PARP9,p$PLSCR1,
             p$RASSF3,p$RIN2,p$RSAD2,
             p$TAP1,p$TAP2,p$TNF,p$TNIP3,
             legend,
             ncol=3,
             left=textGrob("Normalized Counts (x1,000 reads)",
                           gp=gpar(fontsize=16),
                           rot=90))
#bottom=textGrob("HIV status",gp=gpar(fontsize=16)),
dev.off()






#### TRANSCRIPT LEVEL ####

gene.list<-read.table("vidus_hiv_dte_5_cell_types.txt",sep="\t",header=T,stringsAsFactors = F) # aka Supplemental table 3 - DTE results (N=588)


#restrict to bonferonni < 0.05
gene.list<-gene.list[gene.list$Bonferroni < 0.05,]


#5cell-types
load("R_data_files/hiv.shrunk.dte.results.2023_11_28.rda")
load("R_data_files/VIDUS_HIV_DTE_deseq2_2023_11_28.RData") #filtered.vidus.tx.dds

#find ensenbl ID 
#for (e in gene.list$transcript) {
#  x<-rownames(hiv.shrunk.results)[grep(e,rownames(hiv.shrunk.results))]
#  gene.list[gene.list$transcript == e,"transcript"]<-x
#}

# keep only protein coding
dim(gene.list)
protein_coding = vidus.gtf[vidus.gtf$transcript_type == "protein_coding",]
gene.list = gene.list[gene.list$transcript %in% protein_coding$transcript,]
dim(gene.list)




violin.plot.dte <- function(gene=gene,
                            ens=ens,
                            dds=filtered.vidus.tx.dds) {
  #custom margins (top,right,bottom,left)
  custom.margins<-unit(c(0.5,0.5,0,0.5),"cm")
  
  
  
  if (length(ens) > 1) {
    x<-plotCounts(dds, gene=ens[1], intgroup="hiv",returnData=T)
    x$transcript <-ens[1]
    rownames(x)<-NULL
    for (e in ens[2:length(ens)]) {
      #for (e in ens[2:3]) {
      tmp<-plotCounts(dds, gene=e, intgroup="hiv",returnData=T)
      tmp$transcript<-e
      rownames(tmp)<-NULL
      x<-rbind(x,tmp)
    }
  } else {
    x<-plotCounts(dds, gene=ens[1], intgroup="hiv",returnData=T)
    x$transcript <-ens[1]
    rownames(x)<-NULL
  }
  
  x$count<-x$count/1000
  
  # strip away the version number from the transcript ID
  x$transcript = gsub("\\..*$","",x$transcript)
  
  
  p1 <- ggplot(x, aes(x=transcript, y=count, fill=hiv)) +
    geom_violin() +
    #geom_jitter(shape=16, position=position_jitter(0.2)) +
    #stat_summary(fun=median, geom="point") +
    #geom_boxplot(outlier.size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size=22),
          plot.margin=custom.margins,
          axis.text = element_text(size=18),
          #axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1),
          axis.text.x = element_text(size=12, vjust = 1),
          axis.text.y = element_text(size=15),
          #axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size=16,face="bold"),
          legend.title = element_text(size=16,face="bold",hjust=0)) +
    scale_fill_discrete(name = "HIV Status", labels = c("-", "+"))  +
    #scale_x_discrete(labels=c(rep("",length(unique(x$transcript))))) +
    ggtitle(gene) +
    ylab("") +
    xlab("")
  
  p1
  return(p1)
  
}

p <- list()

#unique(full.gene.list$Name)

for (i in 1:length(unique(gene.list$gene_name))) {
  g<-unique(gene.list$gene_name)[i]
  print(g)
  #gene.list <- full.gene.list[full.gene.list$Name == g,]
  
  p[[g]]<-violin.plot.dte(gene = g,
                          dds = filtered.vidus.tx.dds,
                          ens=c(unique(gene.list[gene.list$gene_name == g, "transcript"])))
  legend<-get_legend(p[[g]])
  p[[g]]<-p[[g]] + theme(legend.position = "none")
}


pdf(file="../manuscript/figures/fig2.violinplots_DTE.pdf",width=15,height=10)
grid.arrange(p$EPSTI1,
             p$IFIT3,
             p$MX1,
             p$PARP9,
             p$PLSCR1,
             p$RIN2,
             ncol=2,
             bottom=textGrob("Transcript ID",gp=gpar(fontsize=24,face="bold")),
             left=textGrob("Normalized Counts (x1,000 reads)",
                           gp=gpar(fontsize=24,face="bold"),
                           rot=90))
dev.off()



# make table of transcript results
# combined HIV (N=588)
load("R_data_files/hiv.shrunk.dte.5cells.results.rda")
res.dte.cell<-as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% vidus.gtf$transcript,])
res.dte.cell$transcript<-rownames(res.dte.cell)
row.names(res.dte.cell)<-NULL
res.dte.cell <- merge(res.dte.cell,vidus.gtf,by="transcript")
res.dte.cell<-res.dte.cell[,c("gene_name","transcript","transcript_type","baseMean","log2FoldChange","pvalue")]
res.dte.cell$bonferonni<- res.dte.cell$p * dim(res.dte.cell)[1]
res.dte.cell$bonferonni<-ifelse(res.dte.cell$bonferonni > 1,1,res.dte.cell$bonferonni)
res.dte.cell = res.dte.cell[order(res.dte.cell$gene_name, res.dte.cell$transcript),]

# no VL (N=555)
load("R_data_files/vidus_hiv_acquisition_dte_deseq2_5cell-typeAbs_woVL_shrunk.rda") #
tmp = as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% vidus.gtf$transcript,])
tmp$transcript<-rownames(tmp)
row.names(tmp)<-NULL
tmp <- merge(tmp,vidus.gtf,by="transcript")
tmp<-tmp[,c("gene_name","transcript","transcript_type","baseMean","log2FoldChange","pvalue")]
tmp$bonferonni<- tmp$p * dim(tmp)[1]
tmp$bonferonni<-ifelse(tmp$bonferonni > 1,1,tmp$bonferonni)
tmp = tmp[order(tmp$gene_name, tmp$transcript),]
tmp = tmp[,c("transcript","baseMean","log2FoldChange","pvalue","bonferonni")]
colnames(tmp)[c(2:5)] = c("baseMean_noVL","L2FC_noVL","p_noVL","bon_noVL")
res.dte.cell = merge(res.dte.cell, tmp, by="transcript")



# only VL (N=394)
load("R_data_files/vidus_hiv_acquisition_dte_deseq2_5cell-typeAbs_wVL_shrunk.rda") #
tmp = as.data.frame(hiv.shrunk.results[rownames(hiv.shrunk.results) %in% vidus.gtf$transcript,])
tmp$transcript<-rownames(tmp)
row.names(tmp)<-NULL
tmp <- merge(tmp,vidus.gtf,by="transcript")
tmp<-tmp[,c("gene_name","transcript","transcript_type","baseMean","log2FoldChange","pvalue")]
tmp$bonferonni<- tmp$p * dim(tmp)[1]
tmp$bonferonni<-ifelse(tmp$bonferonni > 1,1,tmp$bonferonni)
tmp = tmp[order(tmp$gene_name, tmp$transcript),]
tmp = tmp[,c("transcript","baseMean","log2FoldChange","pvalue","bonferonni")]
colnames(tmp)[c(2:5)] = c("baseMean_yesVL","L2FC_yesVL","p_yesVL","bon_yesVL")
res.dte.cell = merge(res.dte.cell, tmp, by="transcript")


res.dte.cell = res.dte.cell[order(res.dte.cell$gene_name, res.dte.cell$transcript),]

write.table(res.dte.cell, file="vidus_dte_results_2023_09_26.txt", sep="\t", col.names=T, row.names=F, quote=F)

#
























#### GENE LEVEL - control vs. HIV+ w/ and w/o VL ####
#load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_wVL_shrunk.rda")
#vidus.fit.vl<-vidus.fit
#load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_woVL_shrunk.rda")
#vidus.fit.wovl<-vidus.fit
load("R_data_files/vidus_hiv_acquisition_dge_deseq2_5cell-typeAbs_modelFit.rda")

#classify by control, noVL, yesVL
tmp <- read.table("C:/Users/eearley/Documents/Eric Johnson/VIDUS/candidate genes study/Data/master.pheno.drugs.5cellTypes.txt",sep="\t",header=T,stringsAsFactors = F)
samp.vl <- tmp[,c("iid","hiv","viral_suppressed")]
samp.vl$class <- ifelse(samp.vl$hiv == 0,"control",
                        ifelse(samp.vl$hiv == 1 & samp.vl$viral_suppressed == 1,"noVL","yesVL"))

#load in gene list
x <- vidus.gtf[,c("gene_name","ensg")]
gene.list <- x[!duplicated(x),]

for (i in 1:nrow(gene.list)) {
  g<-gene.list$gene[i]
  print(g)
  
  p[[g]]<-violin.plot.wVL.woVL(dds=vidus.fit,#filtered.vidus.gene.dds,
                      gene=gene.list$gene[i],
                      ens=gene.list$ensg[i])
  legend<-get_legend(p[[g]])
  p[[g]]<-p[[g]] + theme(legend.position = "none")
}

# export
pdf(file="../manuscript/fig1.violin_plots_DGE_VLclasses.pdf",width=20,height=16)
grid.arrange(p$CX3CR1,p$EPSTI1,p$IFI44L,
             p$IFIT3,p$IFITM1,p$MX1,
             p$NLRC5,p$PARP9,p$PLSCR1,
             p$RASSF3,p$RIN2,p$RSAD2,
             p$TAP1,p$TAP2,p$TNF,
             legend,
             ncol=3,
             bottom=textGrob("HIV status",gp=gpar(fontsize=16)),
             left=textGrob("Normalized Counts (x1,000 reads)",
                           gp=gpar(fontsize=16),
                           rot=90))
dev.off()












































#### GENES - residual plots ####
# see this thread: https://support.bioconductor.org/p/60567/#92847
# "it's important to note that "residuals" make less sense in a GLM context than in an LM context. 
#For example, the likelihood ratio test does not involve the sum of squared residuals 
#but instead the log likelihood of the data given the coefficient estimates."


# so what should I be plotting instead of residuals?
# I want to visualize the LFC based on the adjustment of the model coefficients

# try residuals anyway.
# key gene is CX3CR1 which shows up-reg with base counts
# but LFC is negative
fitted.common.scale = assays(vidus.fit)[["mu"]]/normalizationFactors(vidus.fit)
residuals <- counts(filtered.vidus.gene.dds,normalized=TRUE) - fitted.common.scale

custom.margins<-unit(c(0.2,0.7,-0.5,0),"cm") 

#x<-plotCounts(dds, gene=ens, intgroup="hiv",returnData=T)
#x$count<-x$count/1000
p1 <- ggplot(residuals, aes(x=hiv, y=count, fill=hiv)) + 
  geom_violin() + 
  #geom_jitter(shape=16, position=position_jitter(0.2))
  #stat_summary(fun=median, geom="point") +
  geom_boxplot(width=0.05) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.margin=custom.margins,
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=16,face="bold"),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=16,face="bold"),
        legend.title = element_text(size=16,face="bold",hjust=0)) +
  #legend.position = "none") +
  scale_x_discrete(labels=c("","")) +
  #scale_y_continuous(labels=)
  ggtitle(gene) +
  ylab("") +
  xlab("") +
  scale_fill_discrete(name = "HIV Status", labels = c("-", "+"))  





