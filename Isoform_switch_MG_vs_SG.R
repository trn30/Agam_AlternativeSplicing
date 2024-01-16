library(IsoformSwitchAnalyzeR)
library(BSgenome.Agam.VectorBase.P412)
library(BSgenome.Agam.VectorBase.P410)
load(" /RNA-seq/RData/isoform_switch.RData")

######## Import the counts done by salmon

salmonquant_single= importIsoformExpression(sampleVector=c(" /RNA-seq/out_salmon/single/EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_single_mapped_gomis/quant.sf",
                                                           " /RNA-seq/out_salmon/single/EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_single_mapped_gomis/quant.sf",
                                                           " /RNA-seq/out_salmon/single/EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_single_mapped_gomis/quant.sf",
                                                           " /RNA-seq/out_salmon/single/EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_single_mapped_gomis/quant.sf"))

salmonquant_paired= importIsoformExpression(sampleVector=c( " /RNA-seq/out_salmon/paired/EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_paired_mapped_gomis/quant.sf",
                                                            " /RNA-seq/out_salmon/paired/EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_paired_mapped_gomis/quant.sf",
                                                            " /RNA-seq/out_salmon/paired/EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_paired_mapped_gomis/quant.sf",
                                                            " /RNA-seq/out_salmon/paired/EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA_paired_mapped_gomis/quant.sf"))

head (salmonquant$TPM, 2)
head(salmonquant$counts, 21)
head(salmonquant$abundance)

write.table(salmonquant_single$counts, file = " /tappAS/Expresion_Matrix_SG_vs_MG_DESeq2_2_RAW.txt", col.names = T, row.names = T, quote = F, sep = "\t")


            ### Create the matrix design

Design <- data.frame(sampleID = colnames(salmonquant$abundance)[-1],
                     condition = c("MG","SG","MG","SG"),
                     Infection= c("I1","I1","I3","I3")
)
Design <- as.factor(Design)


### import the desing matrix, the annotation and the counts
dataList_single <- importRdata(
  isoformCountMatrix   = salmonquant_single$counts,
  isoformRepExpression = salmonquant_single$abundance,
  designMatrix         = Design,
  isoformExonAnnoation = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="MG",condition_2="SG"),
  showProgress = FALSE)
dataList_single
#I ignore the message I get back when executing this option.

dataList_paired <- importRdata(
  isoformCountMatrix   = salmonquant_paired$counts,
  isoformRepExpression = salmonquant_paired$abundance,
  designMatrix         = Design,
  isoformExonAnnoation = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="MG",condition_2="SG"),
  showProgress = FALSE)
dataList_paired


# PART 1
## this part includes gene and low isoforms expression filtering, statistical analysis to identify isoform switches and annotate these switches with ORF and write the nucleotide and amino acids secuences in fasta files
analycedSwitchList_single <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_single,
  alpha = 0.05,
  dIFcutoff = 0.1,
  switchTestMethod='DEXSeq',
  genomeObject = Agam,
  pathToOutput = ' /Iso_usage/single_end',
  outputSequences      = TRUE, # change to TRUE whan analyzing your own data
  prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)

analycedSwitchList_paired <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_paired,
  dIFcutoff = 0.1,
  switchTestMethod='DEXSeq',
  genomeObject = Agam,
  pathToOutput = ' /Iso_usage/paired_end',
  outputSequences      = TRUE, # change to TRUE whan analyzing your own data
  prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)
# The number of isoform switches found were:
# Comparison nrIsoforms nrSwitches nrGenes
# 1   MG vs SG        221        160     133 -> SINGLE
# 1   MG vs SG        222        169     138 -> PAIRED
extractSwitchSummary( analycedSwitchList_single )

r<-analycedSwitchList_single$isoformFeatures

#### PART 2
# This part consists of importing and incorporating the results of all external sequence analysis,
# analyse alternative splicing, predict the functional consequences of isoform switches, and map out
# (i) individual genes with isoform switches, and
# ii) genome-wide patterns in the consequences of isoform switching.

analycedSwitchList2_single <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = analycedSwitchList_single,
  alpha                     = 0.05,
  dIFcutoff                 = 0.1,   # Cutoff for defining switch size - set high for short runtime in example data
  n                         = NA,
  removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
  pathToCPATresultFile      = (" /RNA-seq/External_analysis/single_end/CPAT_result.txt"),
  pathToPFAMresultFile      = (" /RNA-seq/External_analysis/single_end/Pfam_result.txt"),
  pathToIUPred2AresultFile  = (" /RNA-seq/External_analysis/single_end/IUPred2_result"),
  pathToSignalPresultFile   = (" /RNA-seq/External_analysis/single_end/SignalP_result.txt"),
  codingCutoff              = 0.39, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4131051/
  outputPlots               = FALSE
)

analycedSwitchList2_paired <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = analycedSwitchList_paired,
  dIFcutoff                 = 0.05,   # Cutoff for defining switch size - set high for short runtime in example data
  n                         = NA,
  removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
  pathToCPATresultFile      = (" /RNA-seq/External_analysis/paired_end/CPAT_result.txt"),
  pathToPFAMresultFile      = (" /RNA-seq/External_analysis/paired_end/Pfam_result.txt"),
  pathToIUPred2AresultFile  = (" /RNA-seq/External_analysis/paired_end/IUPred2_result"),
  pathToSignalPresultFile   = (" /RNA-seq/External_analysis/paired_end/SignalP_result.txt"),
  codingCutoff              = 0.39, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4131051/
  outputPlots               = FALSE
)
# Comparison nrIsoforms nrSwitches nrGenes
# 1   MG vs SG        104         79      64 -> SINGLE
# 1   MG vs SG        102         76      60-> PAIRED


# PAIRED
# gene_ref    gene_id gene_name condition_1 condition_2
# 1 geneComp_00004577 AGAP005563        NA          MG          SG
# 2 geneComp_00008185 AGAP010147        NA          MG          SG
# 3 geneComp_00006515 AGAP008036        NA          MG          SG
# gene_switch_q_value switchConsequencesGene Rank
# 1        2.682164e-12                   TRUE    1
# 2        8.415466e-11                   TRUE    2
# 3        3.800292e-10                   TRUE    3


#VOLCANO PLOT PAIRED
volcano<- ggplot(data=analycedSwitchList2$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()
pdf(file =' /RNA-seq/out_isoform/paired_end/volcano_plot.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()


#########################################################################################################################
### achieve usage ratio

rat<-analycedSwitchList2$isoformSwitchAnalysis
rat$ratio<-rat$dIF*100
dei<- read.delim(" /DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T)
dei$isoform_id<-rownames(dei)
mix<-merge(rat,dei,by="isoform_id")#72

#VOLCANO %usage vs DEI (232 son DIU y DEI )
volcano<- ggplot(data=mix, aes(x=log2FoldChange, y=abs(ratio))) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  labs(x=' LogFC( DEI )', y='% usage change') +
  theme_bw()
x11()
volcano

pdf(file =' /Iso_usage/volcano_plot_usage_DEI.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()



##### I add a usage of 0 to those dei which do not have

posicion1<- which(dei$isoform_id %in% rat$isoform_id)
posicion2<- which(rat$isoform_id %in% dei$isoform_id)
dei$usage<-0
dei[posicion1,8]<-rat[posicion2,11] #### ahora tengo las isoformas que tienen usage con el usage correspondiente

mix2<- dei


volcano<- ggplot(data=mix2, aes(x=log2FoldChange, y=abs(usage))) +
  geom_point(aes( color=abs(usage) > 0  ), # default cutoff
             size=1) +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('Isoform with Switch', values = c('black','red')) +
  labs(x=' LogFC( DEI )', y='% usage change')  +
  theme_bw()

x11()
volcano

pdf(file =' /Iso_usage/volcano_plot_usage_DEI_all.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()

write.table(dei$isoform_id, file = " /Iso_usage/Isoforms_DEI_MG_vs_SG_0.05.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

### I get the gene symbols

immunoDB<- read.delim(" /immunodb_genes_agam.txt", header=T)
immuno_id<-immunoDB$UNIQID
immunoDB$gene_id<-immunoDB$UNIQID
#### VEO QUE dei SON IMMUNES
dei$gene_id<-gsub("-R.","",dei$isoform_id)
imu<-merge(immunoDB,dei,by="gene_id")
### represento el volcano con genes immunes

imu$delabel <- NA
imu$delabel[abs(imu$usage) > 10] <- imu$NAME[abs(imu$usage) > 10]

volcano<- ggplot(data=imu, aes(x=log2FoldChange, y=abs(usage),label=delabel)) +
  geom_point(aes( color=abs(usage) > 0  ), # default cutoff
             size=1
  ) +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('Isoform with Switch', values = c('black','red')) +
  labs(x=' LogFC( DEI )', y='% usage change') +
  geom_text(nudge_x = -0.55,nudge_y = 0.5,size=4) +
  theme_bw()

x11()
volcano

pdf(file =' /Iso_usage/volcano_plot_usage_DEI_immuno_genes.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()
