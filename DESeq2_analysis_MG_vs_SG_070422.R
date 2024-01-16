library(readr)
library(DESeq2)
library(Rsubread)
library(data.table)
library(dplyr)
#### DEseq2 ####

iso_filter<- read.delim(" /DESeq2/All_isoforms_filtered_by_cpm_more_than_one_transcript.txt", header=T)

########

Counts_MG1 <- data.frame(iso_filter$ID,iso_filter$I1d7)
Counts_MG3 <- data.frame(iso_filter$ID,iso_filter$I3d7)
Counts_SG1<-data.frame(iso_filter$ID,iso_filter$I1d14)
Counts_SG3 <- data.frame(iso_filter$ID,iso_filter$I3d14)

identical(Counts_MG1$ID,Counts_MG3$ID) # T
identical(Counts_MG1$ID,Counts_SG1$ID) # T
identical(Counts_MG1$ID,Counts_SG3$ID) # T

All_counts <- cbind(Counts_MG1, Counts_MG3$iso_filter.I3d7, Counts_SG1$iso_filter.I1d14, Counts_SG3$iso_filter.I3d14)
colnames(All_counts)=c("ID", "I1d7", "I3d7", "I1d14", "I3d14")
All_counts <- unique(All_counts)
head(All_counts)

Counts_final <- apply(All_counts[,-1],2,as.integer)
rownames(Counts_final) <- All_counts$ID
rownames(Counts_final)
head(Counts_final)

# load the condition matrix
sampleFiles <- c( "I1d7", "I3d7", "I1d14", "I3d14")
sampleCondition <- c("Midgut","Midgut","Salivary_Gland","Salivary_Gland")
sampleTable <- data.frame(sampleName = sampleFiles,
                          condition = sampleCondition)
rownames(sampleTable) <- sampleTable$sampleName
sampleTable
sampleTable$condition <- factor(sampleTable$condition,levels=c("Midgut","Salivary_Gland"))

#Check that the iso_names of the columns match

all(rownames(sampleTable) %in% colnames(Counts_final)) #T
all(rownames(sampleTable) == colnames(Counts_final))#T
# if the order does not match
#all_counts <- all_counts[, rownames(sampleTable)]
#all(rownames(sampleTable) == colnames(all_counts))

#  before executing DESeq2 I have to check that the count variables are integer and the condition factor variables are integer.
str(Counts_final)
str(sampleTable)
sampleTable$sampleName <- as.factor(sampleTable$sampleName)
class(sampleTable$sampleName)
#  Other options include
#class(counts$infec1)
#class(counts[,c("infec1")])

# load the package and run the programme
############ DESeq2 ###################

dataset2<- DESeqDataSetFromMatrix(countData = Counts_final,
                                  colData = sampleTable,
                                  design = ~ condition)
dataset2


######## Pre-filtering ################# I only keep readings that add up to at least 1 in total.

keep <- rowSums(counts(dataset2)) >= 0
dataset3 <- dataset2[keep,]

######################### Normalization ######################

dataset3 <- estimateSizeFactors(dataset3)
data_counts <-counts(dataset3,normalized=TRUE)
dim(data_counts) #1561
dim(dataset3) #1561
write.table(data_counts," /DESeq2/DESEQ2_norm_counts_MG_vs_SG.txt", quote = F,row.names = T, col.names = T,sep = "\t")

############## Differential expression analysis ############
dataset3 <- DESeq(dataset3)
res_iso <- results(dataset3, alpha= 0.05)
resultsNames(dataset3)
res_iso <- results(dataset3, name="condition_Salivary_Gland_vs_Midgut",alpha=0.05)
res_iso
summary(res_iso)
sum(res_iso$pvalue < 0.05,na.rm = T) # 3993 a 503
sum(res_iso$pvalue < 0.01,na.rm = T) # 2959 a 378
sum(res_iso$padj < 0.05,na.rm = T) # 3176 a 397
sum(res_iso$padj < 0.01,na.rm = T) # 2418 a 298

####### I have to get the ones with a lFC greater than 1 and -1 (a change of at least 10%).

res_iso1<-as.data.frame(res_iso)

res_iso2<-res_iso1[!(res_iso1$log2FoldChange %between% c(-1,1)),]
sum(res_iso2$pvalue < 0.05,na.rm = T) # 463
sum(res_iso2$pvalue < 0.01,na.rm = T) # 374
sum(res_iso2$padj < 0.05,na.rm = T) # 392
sum(res_iso2$padj < 0.01,na.rm = T) # 298

isoforms<-as.data.frame(res_iso2)
isoforms_DEI <- isoforms[isoforms$pvalue <= 0.05,c(1,2,3,4,5,6)]
dim(isoforms_DEI)#463
isoforms_DEI
write.table(rownames(isoforms_DEI), file = " /DESeq2/Isoforms_DEI_MG_vs_SG_0.05_id.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

sg<- isoforms_DEI[isoforms_DEI$log2FoldChange < 0,]
mg<- isoforms_DEI[isoforms_DEI$log2FoldChange > 0,]
write.table(rownames(sg), file = " /DESeq2/Isoforms_DEI_MG_vs_SG_0.05_up_SG.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table(rownames(mg), file = " /DESeq2/Isoforms_DEI_MG_vs_SG_0.05_up_MG.txt",quote = FALSE, sep="\t", row.names = F, col.names = F)
