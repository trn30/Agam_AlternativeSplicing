#### check for overlapping zones
bedtools intersect -a  genomic_data_AgamP4/Genes_AgamP4_release_54.bed -b  genomic_data_AgamP4/Genes_AgamP4_release_54.bed -wa -wb -f 0.25 >  /overlapping/overlapping_genes_with_genes.bed

bedtools intersect -a  genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b  genomic_data_AgamP4/Genes_AgamP4_release_54.bed -wa -wb -f 0.25 >  /overlapping/overlapping_isoforms_body.bed
bedtools intersect -a  genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b  genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -wa -wb -f 0.25 >  /overlapping/overlapping_isoform_promoters.bed

### overlap of promoters with isoform promoters
deg<- read.delim(" /DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T) #4933 DEG
deg$ID<-rownames(deg)
dui<- read.delim(" /Iso_usage/Isoforms_DUI_MG_vs_SG_0.05 (1).txt", header=T) #211 dui

prom_over<- read.delim(" /overlapping/overlapping_isoform_promoters.bed", header=F) #
prom_over$V4<-gsub("-R*.","",prom_over$V4)
prom_over$V10<-gsub("-R*.","",prom_over$V10)
prom_over2<-ifelse(prom_over$V4 != prom_over$V10,"yes","")
sum(prom_over2=="yes")#4824
prom_over3<-prom_over[prom_over$V4 != prom_over$V10,]
prom_over4<- unique(prom_over3$V4)#3754

deg$gene<-gsub("-R*.","",deg$ID)
deg_prom_over<-deg[which(!(deg$ID %in% prom_over3$V4)),]#368


dui_prom_over<-dui[which(!(dui$gene%in% prom_over3$V4)),]#155

######## overlap of promoters with genes
body_over<- read.delim(" /overlapping/overlapping_isoforms_body.bed", header=F) #
body_over$V4<-gsub("-R*.","",body_over$V4)
body_over$V10<-gsub("-R*.","",body_over$V10)
##### creo dos vectoces con las dos columnas encadenadas
v4<-c(body_over$V4,body_over$V10)
v10<-c(body_over$V10,body_over$V4)
### creo un nuevo dataframe doble y cambio las columnas por los vectores anteriormente creadas
body_over_<-rbind(body_over,body_over)
body_over_$V4<-v4
body_over_$V10<-v10


body_over2<-ifelse(body_over_$V4 != body_over_$V10,"yes","")
sum(body_over2=="yes")#9834
body_over3<-body_over_[body_over_$V4 != body_over_$V10,]
body_over4<- unique(body_over3$V4)#5626

deg$gene<-gsub("-R*.","",deg$ID)
deg_body_over<-deg[which(!(deg$gene  %in% body_over3$V4)),]#159 de 456 por lo que quedan 304
write.table((deg_body_over$ID), file = " /overlapping/DEMG_without_overlapping_genes_with_promoters_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
deg_body_over2<-deg[which((deg$gene  %in% body_over3$V4)),]#159 de 456 por lo que quedan 304
write.table((deg_body_over2$ID), file = " /overlapping/DEMG_overlapping_genes_with_promoters_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
t<-unique(deg_body_over2$gene)#131 genes overlapping

dui_body_over<-dui[which(!(dui$gene %in% body_over3$V4)),]#80 de 211 por lo que quedan 131
write.table((dui_body_over$ID), file = " /overlapping/DUI_without_overlapping_genes_with_promoters_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
dui_body_over2<-dui[which((dui$gene %in% body_over3$V4)),]#80 de 211 por lo que quedan 131
write.table((dui_body_over2$ID), file = " /overlapping/DUI_overlapping_genes_with_promoters_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
t<-unique(dui_body_over2$gene)#48 genes overlapping

### gene overlap with genes
gene_over<- read.delim(" /overlapping/overlapping_genes_with_genes.bed", header=F) #

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
deg_gene_over<-deg[which(!(deg$gene  %in% gene_over3$V4)),]#56 de 456 por lo que quedan 407
write.table((deg_gene_over$ID), file = " /overlapping/DEMG_without_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
deg_gene_over2<-deg[which((deg$gene  %in% gene_over3$V4)),]#56 de 456 por lo que quedan 407
write.table((deg_gene_over2$ID), file = " /overlapping/DEMG_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
t<-unique(deg_gene_over2$gene)#44 genes overlapping

dui_gene_over<-dui[which(!(dui$gene %in% gene_over3$V4)),]#31 de 211 por lo que quedan 180
write.table((dui_gene_over$ID), file = " /overlapping/DUI_without_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
dui_gene_over2<-dui[which((dui$gene %in% gene_over3$V4)),]#31 de 211 por lo que quedan 180
write.table((dui_gene_over2$ID), file = " /overlapping/DUI_overlapping_genes_ID_definitive.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
t<-unique(dui_gene_over2$gene)#19 genes overlapping

###############################################################################################################################################
