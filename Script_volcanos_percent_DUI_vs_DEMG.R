####  volcano plot representation of SDR vd. % change of usage
library(IsoformSwitchAnalyzeR)
library(BSgenome.Agam.VectorBase.P412)

counts_p<- read.delim(" /Isoform_DEI/Isoforms_DEI_switch_DESeq2/Counts_MG_vs_SG_switchisoform.txt", header=T)


load(" /Iso_usage/env_DESeq2_SG_vs_MG_comparison.R.RData")
load(" /RNA-seq/RData/isoform_switch.RData")

rat<-analycedSwitchList2_paired$isoformSwitchAnalysis
rat$ratio<-rat$dIF*100

### I extract the DEG genes

deg<- read.delim(" /DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T)
deg$id<-rownames(deg)
g<-deg$id
deg$id<-gsub("-R.*","",deg$id)
g<-gsub("-R.*","",g)

# now I relate to each isoform the FoldChange value of the gene

rat_g<-rat
rat_g$id<-gsub("-R.*","",rat$isoform_id)
i<-rat_g$id
i<- unique(i)

mer<-merge(rat_g,deg, by="id")

#### I have to add to the gene dataset those genes that are DUI

deg_f<- deg
deg_f$DUI<-g %in% i
deg_f$usage<-ifelse(g %in% i, rat_g$ratio, 0)


# VOLCANO %usage vs SDR (114 are DUI and DEI )
volcano<- ggplot(data=deg_f, aes(x=log2FoldChange, y=abs(usage))) +
  geom_point(aes( color=abs(usage) > 0  ),   size=1 )+
  geom_vline(xintercept = c(-1, 1), linetype='dashed') +
  scale_color_manual('DEMG ', values = c('black','red')) + scale_fill_manual('DEG', values = c('green','blue')) +
  labs(x=' LogFC( DEMG )', y='% usage change')  +
  theme_bw()
x11()
volcano

pdf(file =' /Iso_usage/volcano_plot_usage_DEMG.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()

### immuno genes
# merge of diffbind peaks isoforms with immunodb immune genes
immunoDB<- read.delim(" /immunodb_genes_agam.txt", header=T)
immuno_id<-immunoDB$UNIQID
deg_f$isoform_id<-rownames(deg_f)

immunoDB$id<-immunoDB$UNIQID

pp<-merge(immunoDB,deg_f,by="id")

### add symbol
sim <- read.delim(" /Iso_usage/Isoforms_vectorbase_information.txt")
sim$isoform_id<- sim$source_id
mix21<- merge(sim, pp, by="isoform_id")
mix2$Gene.Name.or.Symbol[mix2$Gene.Name.or.Symbol == "N/A"]<-NA
mix2$gene_id<-mix2$Gene.ID
mix2$delabel <- NA
mix2$delabel[abs(mix2$usage) > 0] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 0]
mix2$delabel[abs(mix2$usage) > 50] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 50]
mix2$delabel[abs(mix2$usage) > 50] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 50]
mix2$delabel[abs(mix2$usage) > 50] <- mix2$Gene.Name.or.Symbol[abs(mix2$usage) > 50]


### volcano with isoforms that are deg + dui and immunes
volcano<- ggplot(data=mix2, aes(x=log2FoldChange, y=abs(usage),label=delabel)) +
  geom_point(aes( color=abs(usage) > 0  ), # default cutoff
             size=1
  ) +
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('Isoform with Switch', values = c('black','red')) +
  labs(x=' LogFC(DEMG)', y='% usage change') +
  geom_text() +
  theme_bw()
x11()
volcano


###  I do the same but for the genes not just the images.

mix2<- merge(sim, deg_f, by="isoform_id")
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
  scale_color_manual('Isoform with Switch', values = c('black','red')) +
  labs(x=' LogFC( DEMG )', y='% usage change') +
  geom_text_repel() +
  theme_bw()
x11()
volcano



###  I will add to the volcano that only the immunes genes have symbolo

########

mix2<- merge(sim, deg_f, by="isoform_id")
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
  scale_color_manual('Isoform with Switch', values = c('black','red')) +
  labs(x=' LogFC( DEMG )', y='% usage change') +
  geom_text_repel() +
  theme_bw()
x11()
volcano


pdf(file =' /volcano_plot_usage_DEMG_081122.pdf', onefile = FALSE, height=6, width = 9)
volcano
dev.off()
