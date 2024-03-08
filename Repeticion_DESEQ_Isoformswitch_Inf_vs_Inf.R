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

salmonquant_MGSG= importIsoformExpression(sampleVector=c( "/Paired/EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf",
                                                        "/Paired/EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf",
                                                        "/Paired/EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf",
                                                        "/Paired/EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf"))


################## DUI

### Create the matrix design

Design_MGSG <- data.frame(sampleID = colnames(salmonquant_MGSG$abundance)[-1],
                        condition = c("Midguts","Salivary_glands","Midguts","Salivary_glands"),
                        Infection= c("I1","I1","I3","I3")
)


### import the desing matrix, the annotation and the counts
#ignoro el mesaque que me devuelve al ejecutar esta opcion
dataList_MGSG <- importRdata(
  isoformCountMatrix   = salmonquant_MGSG$counts,
  isoformRepExpression = salmonquant_MGSG$abundance,
  designMatrix         = Design_MGSG,
  isoformExonAnnoation = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="Midguts",condition_2="Salivary_glands"),
  showProgress = FALSE)
dataList_MGSG <- preFilter(dataList_MGSG) 

# PART 1
## this part includes gene and low isoforms expression filtering, statistical analysis to identify isoform switches and annotate these switches with ORF and write the nucleotide and amino acids secuences in fasta files
analycedSwitchList_MGSG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_MGSG,
  pathToGTF = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  pathToOutput = '/test_DUI/',
  outputSequences      = F, # change to TRUE when analyzing your own data
  prepareForWebServers = F  # change to TRUE if you will use webservers for external sequence analysis
  )
#                  Comparison nrIsoforms nrSwitches nrGenes
# Midguts vs Salivary_glands        247        175     145


extractSwitchSummary(analycedSwitchList_MGSG )

MGSG<-analycedSwitchList_MGSG$isoformSwitchAnalysis

pval_MGSG<-MGSG[which(MGSG$padj < 0.05),]
dif_pval_MGSG<-pval_MGSG[which(abs(pval_MGSG$dIF) > 0.1),] ##247

write.table(dif_pval_MGSG, file = "/ctrol_vs_inf/Isoforms_DUI_inf_Mg_vs_inf_SG_padj_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


mg1<- dif_pval_MGSG[dif_pval_MGSG$dIF < 0,]#123
sg1<- dif_pval_MGSG[dif_pval_MGSG$dIF > 0,]#124
write.table(mg1, file = "/ctrol_vs_inf/Isoforms_DUI_inf_MG_vs_inf_SG_padj_0.05_up_MG.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
write.table(sg1, file = "/ctrol_vs_inf/Isoforms_DUI_inf_MG_vs_inf_SG_padj_0.05_up_SG.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

#############################################################################################################################################3
#############################################################################################################################################

########################    DESEQ2    #########################################################################
####################################### MGSG
counts_MGSG<-salmonquant_MGSG$counts
colnames(counts_MGSG)<- c("ID","I1d7", "I1d14", "I3d7", "I3d14")

#### I will filter for genes with more than one isoform.
counts_MGSG$gene<-gsub("-R.","",counts_MGSG$ID)
a<-gsub("AGAP[0-9]{6}-RA","FALSE",counts_MGSG$ID) ### hago que las isofromas -RA pongan FALSE
pos<-which(counts_MGSG$ID == a) # sacamos las posiciones de RB y restamos 1 a su posicion para tener su RA tambien
pos_ra<- pos-1

pos_def<-c(pos,pos_ra)
pos_def<-sort(pos_def)# ordeno de menor a mayor
pos_def<-unique(pos_def)# elimino los duplicados
counts_MGSG_def<-counts_MGSG[pos_def,]#3160
## these are genes with more than one isoform.


##############################


#### From the information found on this page, we will filter by a CPM value that equals 10 counts http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
#### so let's calculate the cpm of our samples.

reads_mg1<-28594227
cpm_mg1 <- ((counts_MGSG_def$I1d7*10^6)/reads_mg1)

reads_sg1<-22100970
cpm_sg1 <- ((counts_MGSG_def$I1d14*10^6)/reads_sg1)

reads_mg3<-18984545
cpm_mg3 <- ((counts_MGSG_def$I3d7*10^6)/reads_mg3)

reads_sg3<-45308711
cpm_sg3 <- ((counts_MGSG_def$I3d14*10^6)/reads_sg3)

counts_MGSG_def <- as.data.frame(cbind(counts_MGSG_def,cpm_mg1,cpm_sg1,cpm_mg3,cpm_sg3))

#### I now calculate for each sample how many cpm equals 10 counts.
# I1D7   10*10^6/reads_mg1 = 0.35
# I1d14  10*10^6/reads_sg1 = 0.45
# I3d7   10*10^6/reads_mg3 = 0.52
# I3d14  10*10^6/reads_sg3 = 0.22
iso_filter<-filter(counts_MGSG_def, cpm_mg1>0.35  & cpm_sg1>0.45  & cpm_mg3>0.52  & cpm_sg3>0.22 ) #1561 isoformas quedan
write.table(iso_filter, file = "/ctrol_vs_inf/All_isoforms_filtered_by_cpm_more_than_one_transcript_inf_MG_vs_inf_SG.txt",quote = FALSE, sep="\t", row.names = FALSE, col.names = T)

library(readr)
library(DESeq2)
library(data.table)
library(dplyr)


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
#si el orden no coincide
#all_counts <- all_counts[, rownames(sampleTable)]
#all(rownames(sampleTable) == colnames(all_counts))

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
#write.table(data_counts,"/DESeq2/DESEQ2_norm_counts_MG_vs_SG.txt", quote = F,row.names = T, col.names = T,sep = "\t")

# PCAs:
vsd<- vst(dataset3, blind=FALSE)
# vsd is now the normalized log2-transformed data
colData(dataset3)
pdf(file = "/ctrol_vs_inf/PCA_Inf_MG_vs_inf_SG_R.pdf")
plotPCA(vsd, intgroup=c("condition"))
dev.off()

# Correlation between the counts:
cor(assay(vsd))
a <- pheatmap::pheatmap(cor(assay(vsd)))
pdf(file = "/ctrol_vs_inf/Heatmap_Inf_MG_vs_inf_SG_R.pdf")
a
dev.off()

#################################### QC
####################################
            counts_MGSG<-counts(dataset2)

            counts_MGSG_log<-edgeR::cpm((counts_MGSG), log=TRUE)
            colnames(counts_MGSG_log)<-c( "I1d7", "I3d7", "I1d14", "I3d14")

            counts_MGSG_no_log<-edgeR::cpm((counts_MGSG), log=F)
            colnames(counts_MGSG_no_log)<-c( "I1d7", "I3d7", "I1d14", "I3d14")


            ### 5.3. MDS_norm == MDS
            z <-  limma::plotMDS(counts_MGSG_no_log, col=c("black","red","black","red"), gene.selection = "pairwise", plot=F)
            edge <- sd(z$x)
            plotMDS(counts_MGSG_no_log,  col=c("black","red","black","red"), gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
            title(main="MDS-PCoA Sample Names Norm")


            ### 5.4. MDS_log_norm == MDS_Log
            par(mfrow=c(1,1))
            z <- limma::plotMDS(counts_MGSG_log, gene.selection = "pairwise")
            z
            edge <- sd(z$x)
            #cat(edge)
            pdf(file = "/ctrol_vs_inf/PCOA_Inf_MG_vs_inf_SG_R.pdf")
            plotMDS(counts_MGSG_log, gene.selection = "pairwise",col=c("black","red","black","red"), xlim=c(min(z$x)-edge,max(z$x) + edge))
            title(main="MDS-PCoA log2 Sample Names Norm")
            dev.off()


            ### 6. PCA by type
            data_pca <- as.matrix(counts_MGSG)
            data_pca <- as.data.frame(t(data_pca))
            data_pca.PC = prcomp(data_pca)
            data_pca$Type <- c("Midgut", "Midgut", "Salivary_glands", "Salivary_glands")
            data_pca$Name <- colnames(counts_MGSG)
            pdf(file = "/ctrol_vs_inf/PCA_by_type_inf_MG_vs_inf_SG_R.pdf")

            plot(autoplot(data_pca.PC,  label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
            dev.off()

            ### 8.3. Dendogram cluster raw norm no log
            par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
            pr.hc.c <- hclust(na.omit(dist(t(edgeR::cpm(counts_MGSG,log=F)),method = "euclidean")))
            pdf(file = "/ctrol_vs_inf/Dendogram_inf_MG_vs_inf_SG_R.pdf")
            plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of normalized counts  ", sep=""), labels=colnames(counts_MGSG), cex=0.5)
            dev.off()
            ### 8.3. Dendogram cluster raw norm log
            par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
            pr.hc.c <- hclust(na.omit(dist(t(edgeR::cpm(counts_MGSG,log=T)),method = "euclidean")))
            pdf(file = "/ctrol_vs_inf/Dendogram_log_inf_MG_vs_inf_SG_R.pdf")
            plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts " , sep=""), labels=colnames(counts_MGSG), cex=0.5)
            dev.off()


############## Differential expression analysis ############
dataset3 <- DESeq(dataset3)
res_iso <- results(dataset3, alpha= 0.05)
resultsNames(dataset3)
res_iso <- results(dataset3, name="condition_Salivary_Gland_vs_Midgut",alpha=0.05)
res_iso
summary(res_iso)
sum(res_iso$pvalue < 0.05,na.rm = T) #  503
sum(res_iso$pvalue < 0.01,na.rm = T) #   377
sum(res_iso$padj < 0.05,na.rm = T) #   397
sum(res_iso$padj < 0.01,na.rm = T) #   298

####### I have to get the ones with a lFC greater than 1 and -1 (a change of at least 10%).
res_iso1<-as.data.frame(res_iso)

res_iso2<-res_iso1[!(res_iso1$log2FoldChange %between% c(-1,1)),]
sum(res_iso2$pvalue < 0.05,na.rm = T) # 463
sum(res_iso2$pvalue < 0.01,na.rm = T) # 373
sum(res_iso2$padj < 0.05,na.rm = T) # 392
sum(res_iso2$padj < 0.01,na.rm = T) # 298

isoforms<-as.data.frame(res_iso2)
isoforms_DEI <- isoforms[isoforms$padj <= 0.05,c(1,2,3,4,5,6)]
dim(isoforms_DEI)#399
isoforms_DEI<- na.omit(isoforms_DEI) # 392
write.table(isoforms_DEI, file = "/ctrol_vs_inf/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

mg<- isoforms_DEI[isoforms_DEI$log2FoldChange < 0,]#182
sg<- isoforms_DEI[isoforms_DEI$log2FoldChange > 0,]#210
write.table(mg, file = "/ctrol_vs_inf/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05_up_MG.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
write.table(sg, file = "/ctrol_vs_inf/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05_up_SG.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

######### comparison of our isoforms with the results of (https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-12-216)

dui<- dif_pval_MGSG
dui$gene<-gsub("-R.*","",dui$isoform_id)

deg<- isoforms_DEI
deg$isoform_ID<-rownames(deg)
deg$gene<-gsub("-R.*","",deg$isoform_ID)
malar<- read.delim("/Malar2013/Genes_paper_malar_2013.txt", header=F)


malar[which(malar$V2 %in% dui$gene),] #4 3 unique
malar[which(malar$V2 %in% deg$gene),] #9 8 unique
