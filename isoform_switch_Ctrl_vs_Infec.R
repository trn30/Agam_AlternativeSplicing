library(IsoformSwitchAnalyzeR)
library(BSgenome.Agam.VectorBase.P412)
######## Import the counts done by salmon

salmonquant_MG= importIsoformExpression(sampleVector=c( " /RNA-seq/out_salmon/Midguts/Midgut_Infected_7d_M2/quant.sf",
                                                        " /RNA-seq/out_salmon/Midguts/Midgut_Non-Infected_7d_M2/quant.sf",
                                                        " /RNA-seq/out_salmon/Midguts/Midgut_Infected_7d_M3/quant.sf",
                                                        " /RNA-seq/out_salmon/Midguts/Midgut_Non-Infected_7d_M3/quant.sf"))

salmonquant_SG= importIsoformExpression(sampleVector=c( " /RNA-seq/out_salmon/Salivary_glands/Salivary_glands_Infected_14d_M2/quant.sf",
                                                        " /RNA-seq/out_salmon/Salivary_glands/Salivary_glands_Non-Infected_14d_M2/quant.sf",
                                                        " /RNA-seq/out_salmon/Salivary_glands/Salivary_glands_Infected_14d_M3/quant.sf",
                                                        " /RNA-seq/out_salmon/Salivary_glands/Salivary_glands_Non-Infected_14d_M1/quant.sf"))



### Create the matrix design

Design_MG <- data.frame(sampleID = colnames(salmonquant_MG$abundance)[-1],
                        condition = c("Infected","Non-Infected","Infected","Non-Infected")
)
Design_SG <-Design_MG <- data.frame(sampleID = colnames(salmonquant_SG$abundance)[-1],
                                    condition = c("Infected","Non-Infected","Infected","Non-Infected")
)


### import the desing matrix, the annotation and the counts
# I ignore the table that returns when executing this option.
dataList_MG2 <- importRdata(
  isoformCountMatrix   = salmonquant_MG$counts,
  isoformRepExpression = salmonquant_MG$abundance,
  designMatrix         = Design_MG,
  isoformExonAnnoation = "/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf.gz",
  isoformNtFasta       = "/genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="Infected",condition_2="Non-Infected"),
  showProgress = FALSE)

dataList_SG <- importRdata(
  isoformCountMatrix   = salmonquant_SG$counts,
  isoformRepExpression = salmonquant_SG$abundance,
  designMatrix         = Design_SG,
  isoformExonAnnoation = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf.gz",
  isoformNtFasta       = " /genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame( condition_1 ="Non_Infected",condition_2="Infected"),
  showProgress = FALSE)

load(" /RData/analisis_2.RData")

save.image(" /genomic_data_AgamP4/analisis_2_parte_2.RData")

# PART 1
## this part includes gene and low isoforms expression filtering, statistical analysis to identify isoform switches and annotate these switches with ORF and write the nucleotide and amino acids secuences in fasta files
analycedSwitchList_MG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_MG,
  dIFcutoff = 0.1,
  switchTestMethod='DEXSeq',
  genomeObject = Agam,
  pathToOutput = ' /out_isoform/Midguts',
  outputSequences      = TRUE, # change to TRUE when analyzing your own data
  prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)
#                  Comparison nrIsoforms nrSwitches nrGenes
#1 Infected vs Non_Infected          6          4      4




analycedSwitchList_SG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_SG,
  dIFcutoff = 0.1,
  switchTestMethod='DEXSeq',
  genomeObject = Agam,
  pathToOutput = ' /out_isoform/Salivary_glands',
  outputSequences      = TRUE, # change to TRUE when analyzing your own data
  prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)
#                  Comparison nrIsoforms nrSwitches nrGenes
#1 Infected vs Non_Infected         22         13      13


extractSwitchSummary( analycedSwitchList_MG )
extractSwitchSummary( analycedSwitchList_SG )

# The nucleotide and amino acid sequences of these isoforms have been submitted to the directory provided.
# These sequences allow external analysis of protein domes (Pfam), coding potential (CPAT/CPC2) or signal peptides (SignalIP).
# See ?analyzeCPAT, ?analyzeCPC2, ?analyzePFAM or ?analyzeSignalIP (under details) for suggested ways to run these three tools.
## netsurf2 analysis must be done with files of 100 sequences maximum and 100000 amino acids.
# To split the file use the pyfasta program with this command " pyfasta split -n 5 isoformSwitchAnalyzeR_isoform_AA_subset_1_of_1.fasta"
# https://www.biostars.org/p/13270/

#### PART 2
# This part consists of importing and incorporating the results of all external sequence analysis,
# analyse alternative splicing, predict the functional consequences of isoform switches, and map out
# (i) individual genes with isoform switches, and
# ii) genome-wide patterns in the consequences of isoform switching.

analycedSwitchList2_MG <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = analycedSwitchList_MG,
  dIFcutoff                 = 0.1,   # Cutoff for defining switch size - set high for short runtime in example data
  n                         = NA,
  removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
  pathToCPATresultFile      = (" /External_analysis/Midguts/CPAT_MG.txt"),
  pathToPFAMresultFile      = (" /External_analysis/Midguts/Pfam_MG"),
  pathToIUPred2AresultFile  = (" /External_analysis/Midguts/MG.result"),
  pathToSignalPresultFile   = NULL,
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'IDR_identified',
    'IDR_type'),
  codingCutoff              = 0.39, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4131051/
  outputPlots               = FALSE
)


analycedSwitchList2_SG <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = analycedSwitchList_SG,
  dIFcutoff                 = 0.1,   # Cutoff for defining switch size - set high for short runtime in example data
  n                         = NA,
  removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
  pathToCPATresultFile      = (" /External_analysis/Salivary_glands/CPAT_SG.txt"),
  pathToPFAMresultFile      = (" /External_analysis/Salivary_glands/Pfam_SG"),
  pathToIUPred2AresultFile  = (" /External_analysis/Salivary_glands/SG.result"),
  pathToSignalPresultFile   = (" /External_analysis/Salivary_glands/SignalP_SG.txt"),
  codingCutoff              = 0.39, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4131051/
  outputPlots               = FALSE
)


### Volcano plot
volcano_MG<- ggplot(data=analycedSwitchList_MG$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
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

volcano_SG<- ggplot(data=analycedSwitchList_MG$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
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
