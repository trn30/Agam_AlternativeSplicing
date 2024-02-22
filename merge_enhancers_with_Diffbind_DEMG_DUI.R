####################### I will try to see how many enhacers match diffbind peaks and DEMG and DUI isoforms.

### add to the DEMG and DUI the positions of the whole genomic region
deg<- read.delim("/DEMG/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05.txt", header=T)
deg$ID<-rownames(deg)
dui<- read.delim("/DUI/Isoforms_DUI_inf_Mg_vs_inf_SG_padj_0.05.txt", header=T) #211 dui
colnames(dui)[3]<-"ID"

region_pro<- read.delim("/Genomes/genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed", header=F)
colnames(region_pro)<-c("Chr","Start","End","ID","Strand","Score")
region_iso<- read.delim("/Genomes/genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed", header=F)
colnames(region_iso)<-c("Chr","Start_iso","End_iso","ID","Strand","Score")

deg_all<-merge(deg,region_pro,by="ID")
deg_all<-merge(deg_all,region_iso,by="ID")
deg_all<-deg_all[,c(8,9,15,1,16,17,3)]
deg_all$gene<-gsub("-R.*","",deg_all$ID)
unique(deg_all$gene) #302 genes

dui_all<-merge(dui,region_pro,by="ID")
dui_all<-merge(dui_all,region_iso,by="ID")
dui_all<-dui_all[,c(11,12,18,1,19:20)]
dui_all$gene<-gsub("-R.*","",dui_all$ID)
unique(dui_all$gene) #354 genes

write.table((deg_all), file = "/Enhancers/DEMG_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_all), file = "/Enhancers/DUI_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### A.GAMBIAE ENHANCERS OBTAINED BY PREVIOUS AUTHORES

enh_a<- read.delim("_def/enhancers_JL/S11_novel_regulatory_regions_a_JL", header=T)
enh_a<-enh_a[,c(1:4,7,9:12)]
enh_a$THS.Peak.ID<-ifelse(enh_a$THS.Peak.ID == "",NA,enh_a$THS.Peak.ID)
enh_a<-na.omit(enh_a)#1883 omito las que no tiene pico asociado
library(tidyr)
enh_a<-separate_rows(enh_a,Target.by.others, sep = ",")
enh_a$Target.coincidence <-enh_a$Annotated.Gene.ID == enh_a$Target.by.others
enh_a$Chromosome<-paste("AgamP4",enh_a$Chromosome,sep = "_")
unique(enh_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)

# now I look at how many coincide in each column of possibilities and take to the gene.coincidence column the ones in the deg group.

deg_a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% deg_all$gene),]
a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% deg_all$gene),7]
a1<-dplyr::pull(a1, Annotated.Gene.ID)
deg_a1$Gene.coincience<-a1

deg_a2<-enh_a[which(enh_a$Target.by.others %in% deg_all$gene),]
a2<-enh_a[which(enh_a$Target.by.others %in% deg_all$gene),8]
a2<-dplyr::pull(a2, Target.by.others)
deg_a2$Gene.coincience<-a2


deg_a_1<-rbind(deg_a1,deg_a2)
deg_a<-deg_a_1[,c(1:5,10)]
deg_a<-deg_a[!duplicated(deg_a),]
unique(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#88 enhancers present as defined by previous authors
unique(deg_a$Gene.coincience)# corresponding to 44 genes

dui_a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% dui_all$gene),]
a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% dui_all$gene),7]
a1<-dplyr::pull(a1, Annotated.Gene.ID)
dui_a1$Gene.coincience<-a1


dui_a2<-enh_a[which(enh_a$Target.by.others %in% dui_all$gene),]
a2<-enh_a[which(enh_a$Target.by.others %in% dui_all$gene),8]
a2<-dplyr::pull(a2, Target.by.others)
dui_a2$Gene.coincience<-a2


dui_a_1<-rbind(dui_a1,dui_a2) #240 enhancers
dui_a<-dui_a_1[,c(1:5,10)]
dui_a<-dui_a[!duplicated(dui_a),]
unique(dui_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#33 enhancers present as defined by previous authors
unique(dui_a$Gene.coincience)# corresponding to 17 genes

#### ENHANCERS A.GAMBIAE HOMOLOGOUS TO THOSE OF DROSOPHILA

enh_b<- read.delim("_def/enhancers_JL/S11_novel_regulatory_regions_b_JL", header=T)
enh_b<-enh_b[,c(1:4,7,9:12)]
enh_b$THS.Peak.ID<-ifelse(enh_b$THS.Peak.ID == "",NA,enh_b$THS.Peak.ID)
enh_b<-na.omit(enh_b)#383 omito las que no tiene pico asociado
library(tidyr)
enh_b<-separate_rows(enh_b,Target.by.others..An..gambiae.ortholog., sep = ",")
enh_b$Target.coincidence <-enh_b$Annotated.Gene.ID == enh_b$Target.by.others..An..gambiae.ortholog.
enh_b$Chromosome<-paste("AgamP4",enh_b$Chromosome,sep = "_")
unique(enh_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)

deg_b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% deg_all$gene),]
b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% deg_all$gene),7]
b1<-dplyr::pull(b1, Annotated.Gene.ID)
deg_b1$Gene.coincience<-b1

deg_b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% deg_all$gene),]
b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% deg_all$gene),8]
b2<-dplyr::pull(b2, Target.by.others..An..gambiae.ortholog. )
deg_b2$Gene.coincience<-b2

deg_b_1<-rbind(deg_b1,deg_b2)
deg_b<-deg_b_1[,c(1:5,10)]
deg_b<-deg_b[!duplicated(deg_b),]
unique(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#21enhancers present defined by homology with drosophila
unique(deg_b$Gene.coincience)#corresponding to 22 genes


dui_b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% dui_all$gene),]
b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% dui_all$gene),7]
b1<-dplyr::pull(b1, Annotated.Gene.ID)
dui_b1$Gene.coincience<-b1

dui_b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% dui_all$gene),]
b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% dui_all$gene),8]
b2<-dplyr::pull(b2, Target.by.others..An..gambiae.ortholog.)
dui_b2$Gene.coincience<-b2

dui_b_1<-rbind(dui_b1,dui_b2)
dui_b<-dui_b_1[,c(1:5,10)]
dui_b<-dui_b[!duplicated(dui_b),]
unique(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#13 enhancers present defined by homology with drosophila
unique(dui_b$Gene.coincience)#corresponding to 12 genes


### now with bedtools I will see how many diffbind peaks match with DEMG and DUI enhancers.
### linux terminal

 # bedtools intersect -a /Enhancers/merge_DEMG_enhancers_a.bed -b /ATAC/DiffBind/Diffbind_peaks_JL.bed -wa -wb -f 0.51 > /Enhancers/merge_DEMG_enhancers_a_diffbind_peaks.bed
 # bedtools intersect -a /Enhancers/merge_DUI_enhancers_a.bed -b /ATAC/DiffBind/Diffbind_peaks_JL.bed -wa -wb -f 0.51 > /Enhancers/merge_DUI_enhancers_a_diffbind_peaks.bed
 #
 # bedtools intersect -a /Enhancers/merge_DEMG_enhancers_b.bed -b /ATAC/DiffBind/Diffbind_peaks_JL.bed -wa -wb -f 0.51 > /Enhancers/merge_DEMG_enhancers_b_diffbind_peaks.bed
 # bedtools intersect -a /Enhancers/merge_DUI_enhancers_b.bed -b /ATAC/DiffBind/Diffbind_peaks_JL.bed -wa -wb -f 0.51 > /Enhancers/merge_DUI_enhancers_b_diffbind_peaks.bed



mer_demg_b<- read.delim("/Enhancers/merge_DEMG_enhancers_b_diffbind_peaks.bed", header=F)
unique(mer_demg_b$V4)#1 ENHANCERS
unique(mer_demg_b$V6)#1 GENES

mer_demg_a<- read.delim("/Enhancers/merge_DEMG_enhancers_a_diffbind_peaks.bed", header=F)
p<- unique(mer_demg_a$V4)#21 ENHANCERS
unique(mer_demg_a$V6)#15 GENES

mer_dui_b<- read.delim("/Enhancers/merge_DUI_enhancers_b_diffbind_peaks.bed", header=F)
unique(mer_dui_b$V4)#0 ENHANCERS
unique(mer_dui_b$V6)#0GENES

mer_dui_a<- read.delim("/Enhancers/merge_DUI_enhancers_a_diffbind_peaks.bed", header=F)
pz<-unique(mer_dui_a$V4)#12 ENHANCERS
unique(mer_dui_a$V6)#7 GENES
