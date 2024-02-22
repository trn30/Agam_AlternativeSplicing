#### check for overlapping zones
# Linux terminal
# bedtools intersect -a /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -b /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -wa -wb -f 0.25 > /overlapping/overlapping_genes_with_genes.bed
#
# bedtools intersect -a /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -wa -wb -f 0.25 > /overlapping/overlapping_isoforms_body.bed
# bedtools intersect -a /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -wa -wb -f 0.25 > /overlapping/overlapping_isoform_promoters.bed

### overlap of promoters with isoform promoters
deg<- read.delim("/DEMG/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05.txt", header=T)
deg$ID<-rownames(deg)
dui<- read.delim("/DUI/Isoforms_DUI_inf_Mg_vs_inf_SG_padj_0.05.txt", header=T)

prom_over<- read.delim("/overlapping/overlapping_isoform_promoters.bed", header=F)
prom_over$V4<-gsub("-R*.","",prom_over$V4)
prom_over$V10<-gsub("-R*.","",prom_over$V10)
prom_over2<-ifelse(prom_over$V4 != prom_over$V10,"yes","")
sum(prom_over2=="yes")#4824
prom_over3<-prom_over[prom_over$V4 != prom_over$V10,]
prom_over4<- unique(prom_over3$V4)#3754

deg$gene<-gsub("-R*.","",deg$ID)
deg_prom_over<-deg[which(!(deg$ID %in% prom_over3$V4)),]#392

dui$gene<-gsub("-R*.","",dui$isoform_id)
dui_prom_over<-dui[which(!(dui$gene%in% prom_over3$V4)),]#200

### overlap of genes with genes
gene_over<- read.delim("/overlapping/overlapping_genes_with_genes.bed", header=F) #

##### I create two vectors with the two columns chained together
v4<-c(gene_over$V4,gene_over$V10)
v10<-c(gene_over$V10,gene_over$V4)
### I create a new double dataframe and change the columns to the previously created vectors.
gene_over_<-rbind(gene_over,gene_over)
gene_over_$V4<-v4
gene_over_$V10<-v10

gene_over2<-ifelse(gene_over_$V4 != gene_over_$V10,"yes","")
sum(gene_over2=="yes")#1752
gene_over3<-gene_over_[gene_over_$V4 != gene_over_$V10,]
gene_over4<- unique(gene_over3$V4)#1262


deg$gene<-gsub("-R*.","",deg$ID)
deg_gene_over<-deg[which(!(deg$gene  %in% gene_over3$V4)),]#344 of392
write.table((deg_gene_over$ID), file = "/DEMG_without_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
deg_gene_over2<-deg[which((deg$gene  %in% gene_over3$V4)),]#48 of392
write.table((deg_gene_over2$ID), file = "/DEMG_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
t<-unique(deg_gene_over2$gene)#38 genes overlapping

dui_gene_over<-dui[which(!(dui$gene %in% gene_over3$V4)),]#200 of 247
write.table((dui_gene_over$isoform_id), file = "/DUI_without_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
dui_gene_over2<-dui[which((dui$gene %in% gene_over3$V4)),]#47 of 227
write.table((dui_gene_over2$isoform_id), file = "/DUI_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
t<-unique(dui_gene_over2$gene)#31 genes overlapping

###############################################################################################################################################
