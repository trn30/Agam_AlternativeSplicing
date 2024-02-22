#### volcano plot representation of DEMG vs. % change of usage
library(IsoformSwitchAnalyzeR)


salmonquant_MGSG= importIsoformExpression(sampleVector=c( "/Paired/EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf",
                                                          "/Paired/EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf",
                                                          "/Paired/EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf",
                                                          "/Paired/EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq/quant.sf"))


Design_MGSG <- data.frame(sampleID = colnames(salmonquant_MGSG$abundance)[-1],
                          condition = c("Midguts","Salivary_glands","Midguts","Salivary_glands"),
                          Infection= c("I1","I1","I3","I3")
)



### import the desing matrix, the annotation and the counts
#I don't know the message I get back when I execute this option.
dataList_MGSG <- importRdata(
  isoformCountMatrix   = salmonquant_MGSG$counts,
  isoformRepExpression = salmonquant_MGSG$abundance,
  designMatrix         = Design_MGSG,
  isoformExonAnnoation = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  isoformNtFasta       = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa",
  comparisonsToMake = data.frame(condition_1="Midguts",condition_2="Salivary_glands"),
  showProgress = FALSE)

# PART 1
## this part includes gene and low isoforms expression filtering, statistical analysis to identify isoform switches and annotate these switches with ORF and write the nucleotide and amino acids secuences in fasta files
analycedSwitchList_MGSG <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = dataList_MGSG,
  #dIFcutoff = 0.1,
  #switchTestMethod='DEXSeq',
  #genomeObject = Agam,
  pathToGTF = "/Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gtf",
  pathToOutput = '/test_DUI/',
  outputSequences      = F, # change to TRUE when analyzing your own data
  prepareForWebServers = F  # change to TRUE if you will use webservers for external sequence analysis
)
#                  Comparison nrIsoforms nrSwitches nrGenes
# Midguts vs Salivary_glands        247        175     145


MGSG<-analycedSwitchList_MGSG$isoformSwitchAnalysis

pval_MGSG<-MGSG[which(MGSG$padj < 0.05),]
dif_pval_MGSG<-pval_MGSG[which(abs(pval_MGSG$dIF) > 0.1),] ##247

rat<-dif_pval_MGSG
rat$ratio<-rat$dIF*100

### I extract the DEMG genes

demg<- read.delim("/DEMG/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05.txt", header=T)
demg$id<-rownames(demg)
g<-demg$id
demg$id<-gsub("-R.*","",demg$id)
g<-gsub("-R.*","",g)

#now I relate to each DUI with the FoldChange value of the DEMG isoform.

rat_g<-rat
rat_g$id<-gsub("-R.*","",rat$isoform_id)
i<-rat_g$id
i<- unique(i)

mer<-merge(rat_g,demg, by="id")

#### I have to add to the gene dataset those genes that are DUI

demg_f<- demg
demg_f$DUI<-g %in% i
demg_f$usage<-ifelse(g %in% i, rat_g$ratio, 0)
sum(
  demg_f$DUI == TRUE
)

# VOLCANO %usage vs demg
volcano<- ggplot(data=demg_f, aes(x=log2FoldChange, y=abs(usage))) +
  geom_point(aes( color=abs(usage) > 0  ),   size=1 )+
  geom_vline(xintercept = c(-1, 1), linetype='dashed') +
  scale_color_manual('DEMG ', values = c('black','red')) + scale_fill_manual('demg', values = c('green','blue')) +
  labs(x=' LogFC( DEMG )', y='% usage change')  +
  theme_bw()
volcano

pdf(file ='/volcano_plot_usage_DEMG.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()

### immuno genes
# merge of diffbind peaked isoforms with immunodb immune genes
immunoDB<- read.delim("/Alternative_splicing/immunodb_genes_agam.txt", header=T)
immuno_id<-immunoDB$UNIQID
demg_f$isoform_id<-rownames(demg_f)

immunoDB$id<-immunoDB$UNIQID

pp<-merge(immunoDB,demg_f,by="id")

### add the symbol
sim <- read.delim("/Iso_usage/Isoforms_vectorbase_information.txt")
sim$isoform_id<- sim$source_id
mix21<- merge(sim, pp, by="isoform_id")
mix2$Gene.Name.or.Symbol[mix2$Gene.Name.or.Symbol == "N/A"]<-NA
mix2$gene_id<-mix2$Gene.ID
mix2$delabel <- NA
mix2$delabel[abs(mix2$usage) > 0] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 0]
mix2$delabel[abs(mix2$usage) > 50] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 50]
mix2$delabel[abs(mix2$usage) > 50] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 50]
mix2$delabel[abs(mix2$usage) > 50] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 50]


### volcano with isoforms that are demg + dui and immunes
volcano<- ggplot(data=mix2, aes(x=log2FoldChange, y=abs(usage),label=delabel)) +
  geom_point(aes( color=abs(usage) > 0  ), # default cutoff
             size=1
  ) +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('DEMG + DUI', values = c('black','red')) +
  labs(x=' LogFC(DEMG)', y='% usage change') +
  geom_label_repel() +
  theme_bw()

volcano

### I do the same but for all the genes not only the immune ones.

mix2<- merge(sim, demg_f, by="isoform_id")
mix2$Gene.Name.or.Symbol[mix2$Gene.Name.or.Symbol == "N/A"]<-NA
mix2$gene_id<-mix2$Gene.ID
mix2$delabel <- NA
mix2$delabel[abs(mix2$usage) > 0] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 0]
mix2$delabel2[mix2$Gene.ID %in% mix21$Gene.ID] <- mix2$Gene.Name.or.Symbol[mix2$Gene.ID %in% mix21$Gene.ID]


volcano<- ggplot(data=mix2, aes(x=log2FoldChange, y=abs(usage),label=delabel2)) +
  geom_point(aes( color=abs(usage) > 0  ), # default cutoff
             size=1
  ) +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('DEMG + DUI', values = c('black','red')) +
  labs(x=' LogFC( DEMG )', y='% usage change') +
  geom_label_repel() +
  theme_bw()

volcano

### I will add to the volcano that only the immunes genes have symbolo

########

mix2<- merge(sim, demg_f, by="isoform_id")
mix2$Gene.Name.or.Symbol[mix2$Gene.Name.or.Symbol == "N/A"]<-NA
mix2$gene_id<-mix2$Gene.ID
mix2$delabel <- NA
mix2$delabel[abs(mix2$usage) > 0] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 0]
mix2$delabel2[mix2$Gene.ID %in% mix21$Gene.ID] <- mix2$Gene.Name.or.Symbol[mix2$Gene.ID %in% mix21$Gene.ID]


volcano<- ggplot(data=mix2, aes(x=log2FoldChange, y=abs(usage),label=delabel2)) +
  geom_point(aes( color=abs(usage) > 0  ), # default cutoff
             size=1
  ) +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('DEMG + DUI ', values = c('black','red')) +
  labs(x=' LogFC( DEMG )', y='% usage change') +
  geom_label_repel() +
  theme_bw()

volcano


pdf(file ='/volcano_plot_usage_DEMG_130224.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()
