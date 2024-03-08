library(IsoformSwitchAnalyzeR)
library("edgeR",quiet = T,warn.conflicts = F)
library("RColorBrewer",quiet = T,warn.conflicts = F)
library("genefilter",quiet = T,warn.conflicts = F)
library("ggrepel",quiet = T,warn.conflicts = F)
library("ggfortify",quiet = T,warn.conflicts = F)
library("cluster",quiet = T,warn.conflicts = F)
library("factoextra",quiet = T,warn.conflicts = F)
library("ggplot2",quiet = T,warn.conflicts = F)
library("M3C",quiet = T,warn.conflicts = F)
library("AnnotationDbi",quiet = T,warn.conflicts = F)
#library("xlsx",quiet = T,warn.conflicts = F)
library("clusterProfiler",quiet = T,warn.conflicts = F)
library("png",quiet = T,warn.conflicts = F)
library("curl",quiet = T,warn.conflicts = F)
library("corrplot",quiet = T,warn.conflicts = F)
library("ggpubr",quiet = T,warn.conflicts = F)
library("ggpmisc",quiet = T,warn.conflicts = F)
#library("CDSeq",quiet = T,warn.conflicts = F)
library("plotly",quiet = T,warn.conflicts = F)
library("dplyr",quiet = T,warn.conflicts = F)
library("NormqPCR",quiet = T,warn.conflicts = F)

######## Import the cuantifications done by salmon

salmonquant_MG= importIsoformExpression(sampleVector=c( "/Midguts/Midgut_Infected_7d_M2.fastq.dsrc.fastq//quant.sf",
                                                        "/Midguts/Midgut_Non-Infected_7d_M2.fastq.dsrc.fastq/quant.sf",
                                                        "/Midguts/Midgut_Infected_7d_M3.fastq.dsrc.fastq/quant.sf",
                                                        "/Midguts/Midgut_Non-Infected_7d_M3.fastq.dsrc.fastq/quant.sf"))

salmonquant_SG= importIsoformExpression(sampleVector=c( "/Salivary_glands/Salivary_glands_Infected_14d_M2.fastq.dsrc.fastq/quant.sf",
                                                        "/Salivary_glands/Salivary_glands_Non-Infected_14d_M2.fastq.dsrc.fastq/quant.sf",
                                                        "/Salivary_glands/Salivary_glands_Infected_14d_M3.fastq.dsrc.fastq/quant.sf",
                                                        "/Salivary_glands/Salivary_glands_Non-Infected_14d_M1.fastq.dsrc.fastq/quant.sf"))

################## DUI

### Create the matrix design

Design_MG <- data.frame(sampleID = colnames(salmonquant_MG$abundance)[-1],
                        condition = c("Infected","Non-Infected","Infected","Non-Infected")
)
## try to correct the batch effect
Design_MG$batch <- factor(c('a','a','b','b'))

Design_SG <-data.frame(sampleID = colnames(salmonquant_SG$abundance)[-1],
                                    condition = c("Infected","Non-Infected","Infected","Non-Infected")
)
## try to correct the batch effect
Design_SG$batch <- factor(c('a','a','b','c')) # same as add infection in the desing matrix


### import the desing matrix, the annotation and the counts
#I ignore the message that returns when I execute this option.
dataList_MG <- importRdata(
  isoformCountMatrix   = salmonquant_MG$counts,
  isoformRepExpression = salmonquant_MG$abundance,
  designMatrix         = Design_MG,
  isoformExonAnnoation = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="Infected",condition_2="Non-Infected"),
  showProgress = FALSE)
dataList_MG <- preFilter(dataList_MG)
dataList_SG <- importRdata(
  isoformCountMatrix   = salmonquant_SG$counts,
  isoformRepExpression = salmonquant_SG$abundance,
  designMatrix         = Design_SG,
  isoformExonAnnoation = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame( condition_1 ="Non-Infected",condition_2="Infected"),
  showProgress = FALSE)
dataList_SG <- preFilter(dataList_SG)



# PART 1
## this part includes gene and low isoforms expression filtering, statistical analysis to identify isoform switches and annotate these switches with ORF and write the nucleotide and amino acids secuences in fasta files
analycedSwitchList_MG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_MG,
  pathToGTF = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  pathToOutput = '/ctrol_vs_inf/',
  outputSequences      = F, # change to TRUE when analyzing your own data
  prepareForWebServers = F , # change to TRUE if you will use webservers for external sequence analysis
)

#                  Comparison nrIsoforms nrSwitches nrGenes
#1 Infected vs Non_Infected          6          4      4

analycedSwitchList_SG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_SG,
  pathToGTF = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  pathToOutput = '/test_DUI/',
  outputSequences      = F, # change to TRUE when analyzing your own data
  prepareForWebServers = F  # change to TRUE if you will use webservers for external sequence analysis
)
#                  Comparison nrIsoforms nrSwitches nrGenes
#1 Infected vs Non_Infected         22         13      13


extractSwitchSummary( analycedSwitchList_MG )
extractSwitchSummary( analycedSwitchList_SG )


MG<-analycedSwitchList_MG$isoformSwitchAnalysis

pval_mg<-MG[which(MG$padj < 0.05),]
dif_pval_mg<-pval_mg[which(abs(pval_mg$dIF) > 0.1),] #7

write.table(dif_pval_mg, file = "/ctrol_vs_inf/Isoforms_DUI_ctrol_vs_inf_MG_padj_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


SG<-analycedSwitchList_SG$isoformSwitchAnalysis

pval_sg<-SG[which(SG$padj < 0.05),]
dif_pval_sg<-pval_sg[which(abs(pval_sg$dIF) > 0.1),] #22
write.table(dif_pval_sg, file = "/ctrol_vs_inf/Isoforms_DUI_ctrol_vs_inf_SG_padj_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

#############################################################################################################################################3
#############################################################################################################################################

########################    DESEQ2    #########################################################################
####################################### SG
counts_SG<-salmonquant_SG$counts
colnames(counts_SG)<- c("ID", "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M1")

#### I will filter for genes with more than one isoform.
counts_SG$gene<-gsub("-R.","",counts_SG$ID)
a<-gsub("AGAP[0-9]{6}-RA","FALSE",counts_SG$ID) ### make -RA isoforms FALSE
pos<-which(counts_SG$ID == a) # we take out RB's positions and subtract 1 from his position to get his RA as well.
pos_ra<- pos-1

pos_def<-c(pos,pos_ra)
pos_def<-sort(pos_def)# I order from lowest to highest
pos_def<-unique(pos_def)# I remove duplicates
counts_SG_def<-counts_SG[pos_def,]#3160
## these are genes with more than one isoform.
##############################


#### From the information found on this page, filter by a CPM value that equals 10 counts http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
#### so let's calculate the cpm of our samples.

reads_i2<-26351420
cpm_i2 <- ((counts_SG_def$Infected_M2*10^6)/reads_i2)

reads_ni2<-31528786
cpm_ni2 <- ((counts_SG_def$Non_Infected_M2*10^6)/reads_ni2)

reads_i3<-27146350
cpm_i3 <- ((counts_SG_def$Infected_M3*10^6)/reads_i3)

reads_ni3<-35573386
cpm_ni3 <- ((counts_SG_def$Non_Infected_M1*10^6)/reads_ni3)

counts_SG_def <- as.data.frame(cbind(counts_SG_def,cpm_i2,cpm_ni2,cpm_i3,cpm_ni3))


#### I now calculate for each sample how many cpm equals 10 counts.
# infected m2 10*10^6/reads_i2 = 0.379
# no infected m2 10*10^6/reads_ni2 = 0.317
# infected m3 10*10^6/reads_i3 = 0.368
# no infected m1 / m3 10*10^6/reads_ni3 = 0.281
#### being very restrictive any sample can`t have less tan 10 counts

iso_filter<-filter(counts_SG_def, cpm_i2>=0.379  & cpm_ni2>=0.317  & cpm_i3>=0.368 & cpm_ni3>0.281109) #1679 isoformas quedan

write.table(iso_filter, file = "/ctrol_vs_inf/All_isoforms_filtered_by_cpm_more_than_one_transcript_ctrl_inf_SG.txt",quote = FALSE, sep="\t", row.names = FALSE, col.names = T)

library(readr)
library(DESeq2)
library(data.table)
library(dplyr)

Counts_I2 <- data.frame(iso_filter$ID,iso_filter$Infected_M2)
Counts_NI2 <- data.frame(iso_filter$ID,iso_filter$Non_Infected_M2)
Counts_I3<-data.frame(iso_filter$ID,iso_filter$Infected_M3)
Counts_NI3 <- data.frame(iso_filter$ID,iso_filter$Non_Infected_M1)

identical(Counts_I2$ID,Counts_NI3$ID) # T
identical(Counts_I2$ID,Counts_NI2$ID) # T
identical(Counts_I2$ID,Counts_I3$ID) # T

All_counts <- cbind(Counts_I2, Counts_NI2$iso_filter.Non_Infected_M2, Counts_I3$iso_filter.Infected_M3, Counts_NI3$iso_filter.Non_Infected_M1)
colnames(All_counts)=c("ID", "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")
All_counts <- unique(All_counts)
head(All_counts)

Counts_final <- apply(All_counts[,-1],2,as.integer)
rownames(Counts_final) <- All_counts$ID
rownames(Counts_final)
head(Counts_final)

# load the conditions matrix
sampleFiles <- c( "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")
sampleCondition <- c("Infection","Control","Infection","Control")
sampleTable <- data.frame(sampleName = sampleFiles,
                          condition = sampleCondition)
rownames(sampleTable) <- sampleTable$sampleName
sampleTable
sampleTable$condition <- factor(sampleTable$condition,levels=c("Control","Infection"))

#Check that the iso_names of the columns match
all(rownames(sampleTable) %in% colnames(Counts_final)) #T
all(rownames(sampleTable) == colnames(Counts_final))#T
#si el orden no coincide
#all_counts <- all_counts[, rownames(sampleTable)]
#all(rownames(sampleTable) == colnames(all_counts))

## execute DESeq2
dataset2<- DESeqDataSetFromMatrix(countData = Counts_final,
                                  colData = sampleTable,
                                  design = ~ condition)
dataset2


######## Pre-filtering ################# I only keep readings that add up to at least 1 in total.

keep <- rowSums(counts(dataset2)) >= 0
dataset3 <- dataset2[keep,]
######################### Normalization

dataset3 <- estimateSizeFactors(dataset3)
data_counts <-counts(dataset3,normalized=TRUE)
dim(data_counts) #1678
dim(dataset3) #1678

# PCAs:
vsd<- vst(dataset3, blind=FALSE)
# vsd is now the normalized log2-transformed data
colData(dataset3)
pdf(file = "/ctrol_vs_inf/PCA_ctrol_vs_inf_SG_R.pdf")
plotPCA(vsd, intgroup=c("condition"))
dev.off()
# Correlation between the counts:
cor(assay(vsd))
a <- pheatmap::pheatmap(cor(assay(vsd)))
pdf(file = "/ctrol_vs_inf/Heatmap_ctrol_vs_inf_SG_R.pdf")
a
dev.off()


        #################################### QC
        counts_SG<-counts(dataset2)

        counts_SG_log<-edgeR::cpm((counts_SG), log=TRUE)
        colnames(counts_SG_log)<-c( "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")

        counts_SG_no_log<-edgeR::cpm((counts_SG), log=F)
        colnames(counts_SG_no_log)<-c( "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")


        ### 5.3. MDS_norm == MDS
        z <-  limma::plotMDS(counts_SG_no_log, col=c("black","red","black","red"), gene.selection = "pairwise", plot=F)
        edge <- sd(z$x)
        plotMDS(counts_SG_no_log,  col=c("black","red","black","red"), gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
        title(main="MDS-PCoA Sample Names Norm")


        ### 5.4. MDS_log_norm == MDS_Log
        par(mfrow=c(1,1))
        z <- limma::plotMDS(counts_SG_log, gene.selection = "pairwise")
        z
        edge <- sd(z$x)
        #cat(edge)
        pdf(file = "/ctrol_vs_inf/PCOA_ctrol_vs_inf_SG_R.pdf")
        plotMDS(counts_SG_log, gene.selection = "pairwise",col=c("black","red","black","red"), xlim=c(min(z$x)-edge,max(z$x) + edge))
        title(main="MDS-PCoA log2 Sample Names Norm")
        dev.off()


        ### 6. PCA by type
        data_pca <- as.matrix(counts_SG)
        data_pca <- as.data.frame(t(data_pca))
        data_pca.PC = prcomp(data_pca)
        data_pca$Type <- c("Infected", "Non_Infected", "Infected", "Non_Infected")
        data_pca$Name <- colnames(counts_SG)
        pdf(file = "/ctrol_vs_inf/PCA_by_type_ctrol_vs_inf_SG_R.pdf")
        plot(autoplot(data_pca.PC,  label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
        dev.off()

        ### 8.3. Dendogram cluster raw norm no log
        par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
        pr.hc.c <- hclust(na.omit(dist(t(edgeR::cpm(counts_SG,log=F)),method = "euclidean")))
        pdf(file = "/ctrol_vs_inf/Dendogram_ctrol_vs_inf_SG_R.pdf")
        plot(pr.hc.c, xlab="Sample Distance ",main=paste("Hierarchical Clustering of normalized counts ", sep=""), labels=colnames(counts_SG), cex=0.5)
        dev.off()

        ### 8.3. Dendogram cluster raw norm log
        par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
        pr.hc.c <- hclust(na.omit(dist(t(edgeR::cpm(counts_SG,log=T)),method = "euclidean")))
        pdf(file = "/ctrol_vs_inf/Dendogram_log_ctrol_vs_inf_SG_R.pdf")
        plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts " , sep=""), labels=colnames(counts_SG), cex=0.5)
        dev.off()



############## Differential expression analysis
dataset3 <- DESeq(dataset3)
res_iso <- results(dataset3, alpha= 0.05)
resultsNames(dataset3)
res_iso <- results(dataset3, name="condition_Infection_vs_Control",alpha=0.05)
res_iso
summary(res_iso)

sum(res_iso$pvalue < 0.05,na.rm = T) # 65
sum(res_iso$pvalue < 0.01,na.rm = T) # 20
sum(res_iso$padj < 0.05,na.rm = T) # 2
sum(res_iso$padj < 0.01,na.rm = T) # 0

####### I have to get the ones with a lFC greater than 1 and -1 (a change of at least 10%).
res_iso1<-as.data.frame(res_iso)

library(data.table)
res_iso2<-res_iso1[!(res_iso1$log2FoldChange %between% c(-1,1)),]
sum(res_iso2$pvalue < 0.05,na.rm = T) # 55
sum(res_iso2$pvalue < 0.01,na.rm = T) # 20
sum(res_iso2$padj < 0.05,na.rm = T) # 2
sum(res_iso2$padj < 0.01,na.rm = T) # 0

isoforms<-as.data.frame(res_iso2)
isoforms_DEMG_SG <- isoforms[isoforms$padj <= 0.05,c(1,2,3,4,5,6)]
dim(isoforms_DEMG_SG)#2
isoforms_DEMG_SG

write.table(isoforms_DEMG_SG, file = "/ctrol_vs_inf/Isoforms_DEMG_ctrol_vs_inf_SG_padj_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


####################################################################################################################################################
#####################################3 MG


counts_MG<-salmonquant_MG$counts
colnames(counts_MG)<- c("ID", "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")

#### I will filter for genes with more than one isoform.
counts_MG$gene<-gsub("-R.","",counts_MG$ID)
a<-gsub("AGAP[0-9]{6}-RA","FALSE",counts_MG$ID) ### make -RA isoforms FALSE
pos<-which(counts_MG$ID == a) # we take out RB's positions and subtract 1 from his position to get his RA as well.
pos_ra<- pos-1

pos_def<-c(pos,pos_ra)
pos_def<-sort(pos_def)# I order from lowest to highest
pos_def<-unique(pos_def)# I remove duplicates
counts_MG_def<-counts_MG[pos_def,]#3160
## these are genes with more than one isoform.
##############################


#### From the information found on this page, filter by a CPM value that equals 10 counts http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
#### so let's calculate the cpm of our samples.


reads_i2<-69607948
cpm_i2 <- ((counts_MG_def$Infected_M2*10^6)/reads_i2)

reads_ni2<-28196140
cpm_ni2 <- ((counts_MG_def$Non_Infected_M2*10^6)/reads_ni2)

reads_i3<-81735510
cpm_i3 <- ((counts_MG_def$Infected_M3*10^6)/reads_i3)

reads_ni3<-35501460
cpm_ni3 <- ((counts_MG_def$Non_Infected_M3*10^6)/reads_ni3)

counts_MG_def <- as.data.frame(cbind(counts_MG_def,cpm_i2,cpm_ni2,cpm_i3,cpm_ni3))


#### I now calculate for each sample how many cpm equals 10 counts.
# infected m2   10*10^6/reads_i2 = 0.143
# no infected m2 10*10^6/reads_ni2 = 0.354
# infected m3 10*10^6/reads_i3 = 0.122
# no infected m1 / m3 10*10^6/reads_ni3 = 0.281

iso_filter<-filter(counts_MG_def, cpm_i2>0.143  & cpm_ni2>0.354  & cpm_i3>0.122 & cpm_ni3>0.281 ) #1394 isoformas quedan

write.table(iso_filter, file = "/ctrol_vs_inf/All_isoforms_filtered_by_cpm_more_than_one_transcript_ctrl_inf_MG.txt",quote = FALSE, sep="\t", row.names = FALSE, col.names = T)

#### isoforms that for one condition have at least 10 counts and the other condition have both replicates 0 or more than 0 and less than 10 counts
# iso_filter<-filter(counts_MG_def, cpm_i2>0.143  & cpm_ni2>0.354  & cpm_i3>0.122 & cpm_ni3>0.281)
# iso_filter<-filter(iso_filter, !(cpm_i2 > 0 & cpm_i3== 0 & cpm_i2 < 10 | cpm_i2 == 0 & cpm_i3 > 0 & cpm_i3 < 10 ))
# iso_filter<-filter(iso_filter, !(cpm_ni2 > 0 & cpm_ni3== 0 & cpm_ni2 < 10 | cpm_ni2 == 0 & cpm_ni3 > 0 & cpm_ni3 < 10 ))

library(readr)
library(DESeq2)
library(data.table)
library(dplyr)


Counts_I2 <- data.frame(iso_filter$ID,iso_filter$Infected_M2)
Counts_NI2 <- data.frame(iso_filter$ID,iso_filter$Non_Infected_M2)
Counts_I3<-data.frame(iso_filter$ID,iso_filter$Infected_M3)
Counts_NI3 <- data.frame(iso_filter$ID,iso_filter$Non_Infected_M3)

identical(Counts_I2$ID,Counts_NI3$ID) # T
identical(Counts_I2$ID,Counts_NI2$ID) # T
identical(Counts_I2$ID,Counts_I3$ID) # T

All_counts <- cbind(Counts_I2, Counts_NI2$iso_filter.Non_Infected_M2, Counts_I3$iso_filter.Infected_M3, Counts_NI3$iso_filter.Non_Infected_M3)
colnames(All_counts)=c("ID", "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")
All_counts <- unique(All_counts)
head(All_counts)

Counts_final <- apply(All_counts[,-1],2,as.integer)
rownames(Counts_final) <- All_counts$ID
rownames(Counts_final)
head(Counts_final)

# load the conditions matrix
sampleFiles <- c( "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")
sampleCondition <- c("Infection","Control","Infection","Control")
sampleTable <- data.frame(sampleName = sampleFiles,
                          condition = sampleCondition)
rownames(sampleTable) <- sampleTable$sampleName
sampleTable
sampleTable$condition <- factor(sampleTable$condition,levels=c("Control","Infection"))

#Check that the iso_names of the columns match
all(rownames(sampleTable) %in% colnames(Counts_final)) #T
all(rownames(sampleTable) == colnames(Counts_final))#T
#si el orden no coincide
#all_counts <- all_counts[, rownames(sampleTable)]
#all(rownames(sampleTable) == colnames(all_counts))

### execute DESeq2
dataset2<- DESeqDataSetFromMatrix(countData = Counts_final,
                                  colData = sampleTable,
                                  design = ~ condition)
dataset2


######## Pre-filtering ################# I only keep readings that add up to at least 1 in total.

keep <- rowSums(counts(dataset2)) >= 0
dataset3 <- dataset2[keep,]
######################### Normalization

dataset3 <- estimateSizeFactors(dataset3)
data_counts <-counts(dataset3,normalized=TRUE)
dim(data_counts) #1394
dim(dataset3) #1394

# PCAs:
vsd<- vst(dataset3, blind=FALSE)
colData(dataset3)
pdf(file = "/ctrol_vs_inf/PCA_ctrol_vs_inf_MG_R.pdf")
plotPCA(vsd, intgroup=c("condition"))
dev.off()
# Correlation between the counts:
cor(assay(vsd))
a <- pheatmap::pheatmap(cor(assay(vsd)))
pdf(file = "/ctrol_vs_inf/Heatmap_ctrol_vs_inf_MG_R.pdf")
a
dev.off()

            #################################### QC
            counts_MG<-counts(dataset2)

            counts_MG_log<-edgeR::cpm((counts_MG), log=TRUE)
            colnames(counts_MG_log)<-c( "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")

            counts_MG_no_log<-edgeR::cpm((counts_MG), log=F)
            colnames(counts_MG_no_log)<-c( "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M3")



            ### 5.3. MDS_norm == MDS
            z <-  limma::plotMDS(counts_MG_no_log, col=c("black","red","black","red"), gene.selection = "pairwise", plot=F)
            edge <- sd(z$x)
            plotMDS(counts_MG_no_log,  col=c("black","red","black","red"), gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
            title(main="MDS-PCoA Sample Names Norm")


            ### 5.4. MDS_log_norm == MDS_Log

            z <- limma::plotMDS(counts_MG_log, gene.selection = "pairwise")
            z
            edge <- sd(z$x)
            #cat(edge)
            pdf(file = "/ctrol_vs_inf/PCOA_ctrol_vs_inf_MG_R.pdf")
            plotMDS(counts_MG_log, gene.selection = "pairwise",col=c("black","red","black","red"), xlim=c(min(z$x)-edge,max(z$x) + edge))
            title(main="MDS-PCoA log2 Sample Names Norm")
            dev.off()


            ### 6. PCA by type
            data_pca <- as.matrix(counts_MG)
            data_pca <- as.data.frame(t(data_pca))
            data_pca.PC = prcomp(data_pca)
            data_pca$Type <- c("Infected", "Non_Infected", "Infected", "Non_Infected")
            data_pca$Name <- colnames(counts_MG)
            pdf(file = "/ctrol_vs_inf/PCA_by_type_ctrol_vs_inf_MG_R.pdf")
            plot(autoplot(data_pca.PC,  label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
            dev.off()


            ### 8.3. Dendogram cluster raw norm no log
            par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
            pr.hc.c <- hclust(na.omit(dist(t(edgeR::cpm(counts_MG,log=F)),method = "euclidean")))
            pdf(file = "/ctrol_vs_inf/Dendogram_ctrol_vs_inf_MG_R.pdf")
            plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of normalized counts ", sep=""), labels=colnames(counts_MG), cex=0.5)
            dev.off()


            ### 8.3. Dendogram cluster raw norm log
            par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
            pr.hc.c <- hclust(na.omit(dist(t(edgeR::cpm(counts_MG,log=T)),method = "euclidean")))
            pdf(file = "/ctrol_vs_inf/Dendogram_log_ctrol_vs_inf_MG_R.pdf")
            plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts " , sep=""), labels=colnames(counts_MG), cex=0.5)
            dev.off()





############## Differential expresison analysis
dataset3 <- DESeq(dataset3)
res_iso <- results(dataset3, alpha= 0.05)
resultsNames(dataset3)
res_iso <- results(dataset3, name="condition_Infection_vs_Control",alpha=0.05)
res_iso
summary(res_iso)
sum(res_iso$pvalue < 0.05,na.rm = T) # 52
sum(res_iso$pvalue < 0.01,na.rm = T) # 16
sum(res_iso$padj < 0.05,na.rm = T) # 4
sum(res_iso$padj < 0.01,na.rm = T) # 4

####### I have to get the ones with a lFC greater than 1 and -1 (a change of at least 10%).
res_iso1<-as.data.frame(res_iso)

library(data.table)
res_iso2<-res_iso1[!(res_iso1$log2FoldChange %between% c(-1,1)),]
sum(res_iso2$pvalue < 0.05,na.rm = T) # 19
sum(res_iso2$pvalue < 0.01,na.rm = T) # 13
sum(res_iso2$padj < 0.05,na.rm = T) # 4
sum(res_iso2$padj < 0.01,na.rm = T) # 4

isoforms<-as.data.frame(res_iso2)
isoforms_DEMG_MG <- isoforms[isoforms$padj <= 0.05,c(1,2,3,4,5,6)]
dim(isoforms_DEMG_MG)#4
isoforms_DEMG_MG

write.table(isoforms_DEMG_MG, file = "/ctrol_vs_inf/Isoforms_DEMG_ctrol_vs_inf_MG_padj_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


######### comparison of our isoforms with the results of (https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-12-216)
## SG
dui<- dif_pval_sg
dui$gene<-gsub("-R.*","",dui$isoform_id)

deg<- isoforms_DEMG_SG
deg$isoform_ID<-rownames(deg)
deg$gene<-gsub("-R.*","",deg$isoform_ID)
malar<- read.delim("/Malar2013/Genes_paper_malar_2013.txt", header=F)


malar[which(malar$V2 %in% dui$gene),] #0
malar[which(malar$V2 %in% deg$gene),] #0


### MG
dui<- dif_pval_mg
dui$gene<-gsub("-R.*","",dui$isoform_id)

deg<- isoforms_DEMG_MG
deg$isoform_ID<-rownames(deg)
deg$gene<-gsub("-R.*","",deg$isoform_ID)
malar<- read.delim("/Malar2013/Genes_paper_malar_2013.txt", header=F)


malar[which(malar$V2 %in% dui$gene),] #1
malar[which(malar$V2 %in% deg$gene),] #0
