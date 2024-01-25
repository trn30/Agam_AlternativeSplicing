library(dplyr)
#load("/isoform_2/Counts.RData")

####### I will load the counts of all the isoforms, and make a representation of how the counts are distributed. This way I get the limit over which I have to cut so that I don't have so many 0s.
######## Midguts vs. Salivary glands

#### RNA ####

counts_p<- read.delim(" /Isoform_DEI/Isoforms_DEI_switch_DESeq2/Counts_SG_switchisoform.txt", header=T)
counts_p<- counts_p[,c(6:10,12:15)]
colnames(counts_p)<- c("ID", "Infected_M2", "Non_Infected_M2", "Infected_M3", "Non_Infected_M1","Infected_M2_lenght", "Non_Infected_M2_lenght", "Infected_M3_lenght", "Non_Infected_M1_lenght")

####  I will filter for genes with more than one isoform.
counts_p$gene<-gsub("-R.","",counts_p$ID)
a<-gsub("AGAP[0-9]{6}-RA","FALSE",counts_p$ID) ### make -RA isoforms FALSE
pos<-which(counts_p$ID == a) #  we take out RB's positions and subtract 1 from his position to get his RA as well.
pos_ra<- pos-1

pos_def<-c(pos,pos_ra)
pos_def<-sort(pos_def)# I order from lowest to highest
pos_def<-unique(pos_def)# I remove duplicates
counts_p_def<-counts_p[pos_def,]#3160
## these are the genes with more than one isoform
##############################


#### From the information found on this page, we will filter by a CPM value that equals 10 counts ( http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf )
#### Therefore, we will calculate the cpm of our samples.

reads_i2<-26351420
cpm_i2 <- ((counts_p_def$Infected_M2*10^6)/reads_i2)

reads_ni2<-31528786
cpm_ni2 <- ((counts_p_def$Non_Infected_M2*10^6)/reads_ni2)

reads_i3<-27146350
cpm_i3 <- ((counts_p_def$Infected_M3*10^6)/reads_i3)

reads_ni3<-35573386
cpm_ni3 <- ((counts_p_def$Non_Infected_M1*10^6)/reads_ni3)

counts_p_def <- as.data.frame(cbind(counts_p_def,cpm_i2,cpm_ni2,cpm_i3,cpm_ni3))


#### Now I calculate for each sample how many cpm equals 10 counts.
# infected m2 10*10^6/reads_i2 = 0.379
# no infected m2 10*10^6/reads_ni2 = 0.317
# infected m3 10*10^6/reads_i3 = 0.368
# no infected m1 / m3 10*10^6/reads_ni3 = 0.281

c_clean<-filter(counts_p_def, cpm_i2>0.379  & cpm_ni2>0.317  & cpm_i3>0.368 & cpm_ni3>0.281 ) #1679 isoformas quedan


write.table(c_clean, file = " /DESeq2/All_isoforms_filtered_by_cpm_more_than_one_transcript_ctrl_inf_SG.txt",quote = FALSE, sep="\t", row.names = FALSE, col.names = T)

########################################################################################################
iso_filter<- read.delim(" /DESeq2/All_isoforms_filtered_by_cpm_more_than_one_transcript_ctrl_inf_SG.txt", header=T)


#### DESeq2 Analysis ####

library(readr)
library(DESeq2)
library(Rsubread)
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

# I check that the iso_names of the columns are the same
all(rownames(sampleTable) %in% colnames(Counts_final)) #T
all(rownames(sampleTable) == colnames(Counts_final))#T
# if the order does not match
#all_counts <- all_counts[, rownames(sampleTable)]
#all(rownames(sampleTable) == colnames(all_counts))


# Before executing DESeq2 I have to check that the count variables are integer and the condition factor variables are integer.
# str(Counts_final)
# str(sampleTable)
# sampleTable$sampleName <- as.factor(sampleTable$sampleName)
# class(sampleTable$sampleName)
# Other options include
# class(counts$infec1)
# class(counts[,c("infec1")])


# load the package and run the programme

dataset2<- DESeqDataSetFromMatrix(countData = Counts_final,
                                  colData = sampleTable,
                                  design = ~ condition)
dataset2


######## Pre-filtering ################# I will only keep readings that add up to at least 1 in total

keep <- rowSums(counts(dataset2)) >= 0
dataset3 <- dataset2[keep,]
######################### Normalization

dataset3 <- estimateSizeFactors(dataset3)
data_counts <-counts(dataset3,normalized=TRUE)
dim(data_counts) #1679
dim(dataset3) #1679
write.table(data_counts," /DESeq2/DESEQ2_norm_counts_SG_Ctrl_vs_Infec.txt", quote = F,row.names = T, col.names = T,sep = "\t")

############## Differential expression analysis
dataset3 <- DESeq(dataset3)
res_iso <- results(dataset3, alpha= 0.05)
resultsNames(dataset3)
res_iso <- results(dataset3, name="condition_Infection_vs_Control",alpha=0.05)
res_iso
summary(res_iso)
sum(res_iso$pvalue < 0.05,na.rm = T) # 66
sum(res_iso$pvalue < 0.01,na.rm = T) # 20
sum(res_iso$padj < 0.05,na.rm = T) # 2
sum(res_iso$padj < 0.01,na.rm = T) # 0

####### I have to get the ones with a lFC greater than 1 and -1 (a change of at least 10%).
res_iso1<-as.data.frame(res_iso)

library(data.table)
res_iso2<-res_iso1[!(res_iso1$log2FoldChange %between% c(-1,1)),]
sum(res_iso2$pvalue < 0.05,na.rm = T) # 56
sum(res_iso2$pvalue < 0.01,na.rm = T) # 20
sum(res_iso2$padj < 0.05,na.rm = T) # 2
sum(res_iso2$padj < 0.01,na.rm = T) # 0

isoforms<-as.data.frame(res_iso2)
isoforms_DEI <- isoforms[isoforms$pvalue <= 0.05,c(1,2,3,4,5,6)]
dim(isoforms_DEI)#56
isoforms_DEI
write.table((isoforms_DEI), file = " /DESeq2/Isoforms_DEI_Ctrol_vs_Infec_SG_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)



########################################################################################################


#### Isoform usage analysis ####
#
load("/isoform_2/Counts.RData")

### Create the matrix design

Design_SG <- data.frame(sampleID = colnames(salmonquant_SG$abundance)[-1],
                        condition = c("Infected","Non-Infected","Infected","Non-Infected")
)


### import the desing matrix, the annotation and the counts
# I ignore the table that returns when executing this option.

dataList_SG2 <- importRdata(
  isoformCountMatrix   = salmonquant_SG$counts,
  isoformRepExpression = salmonquant_SG$abundance,
  designMatrix         = Design_SG,
  isoformExonAnnoation = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="Infected",condition_2="Non-Infected"),
  showProgress = FALSE)


# PART 1
## this part includes gene and low isoforms expression filtering, statistical analysis to identify isoform switches and annotate these switches with ORF and write the nucleotide and amino acids secuences in fasta files
analycedSwitchList_SG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_SG2,
  genomeObject = Agam,
  pathToOutput = ' /Iso_usage/',
  outputSequences      = TRUE, # change to TRUE when analyzing your own data
  prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)
#                  Comparison nrIsoforms nrSwitches nrGenes
#1 Infected vs Non_Infected          22         13     13

extractSwitchSummary( analycedSwitchList_SG )
DUI_sg<-analycedSwitchList_SG$isoformSwitchAnalysis ### 35 isoformas son DUI correspondientes a 13 genes

write.table((DUI_sg), file = " /Iso_usage/Isoforms_DUI_Ctrol_vs_Infec_SG_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
