

deg<- read.delim(" /DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T) #4933 DEG
deg$ID<-rownames(deg)
dui<- read.delim(" /Iso_usage/Isoforms_DUI_MG_vs_SG_0.05 (1).txt", header=T) #211 dui

#### add to the DEMG and DUI the positions of the whole genomic region
region_pro<- read.delim(" genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed", header=F)
colnames(region_pro)<-c("Chr","Start","End","ID","Strand","Score")
region_iso<- read.delim(" genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed", header=F)
colnames(region_iso)<-c("Chr","Start_iso","End_iso","ID","Strand","Score")

deg_all<-merge(deg,region_pro,by="ID")
deg_all<-merge(deg_all,region_iso,by="ID")
deg_all<-deg_all[,c(8,9,15,1,16,17,3)]
deg_all$gene<-gsub("-R.*","",deg_all$ID)
unique(deg_all$gene) #354 genes

dui_all<-merge(dui,region_pro,by="ID")
dui_all<-merge(dui_all,region_iso,by="ID")
dui_all<-dui_all[,c(24,25,31,1,32,33,6)]
dui_all$gene<-gsub("-R.*","",dui_all$ID)
unique(dui_all$gene) #354 genes

write.table((deg_all), file = " /enhancers_JL/DEMG_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_all), file = " /enhancers_JL/DUI_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)


#### A.GAMBIAE ENHANCERS OBTAINED BY PREVIOUS AUTHORS

enh_a<- read.delim(" /enhancers_JL/S11_novel_regulatory_regions_a_JL", header=T)
enh_a<-enh_a[,c(1:4,7,9:12)]
enh_a$THS.Peak.ID<-ifelse(enh_a$THS.Peak.ID == "",NA,enh_a$THS.Peak.ID)
enh_a<-na.omit(enh_a)#1883 I omit those with no associated spout
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
unique(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#100 enhancers presentes definidos por previos autores
unique(deg_a$Gene.coincience)#correspondientes a 48 genes

write.table((deg_a), file = " /enhancers_JL/merge_DEMG_enhancers_a.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((deg_a[,6]), file = " /enhancers_JL/merge_DEMG_enhancers_a_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

dui_a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% dui_all$gene),]
a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% dui_all$gene),7]
a1<-dplyr::pull(a1, Annotated.Gene.ID)
dui_a1$Gene.coincience<-a1


dui_a2<-enh_a[which(enh_a$Target.by.others %in% dui_all$gene),]
a2<-enh_a[which(enh_a$Target.by.others %in% dui_all$gene),8]
a2<-dplyr::pull(a2, Target.by.others)
dui_a2$Gene.coincience<-a2


dui_a_1<-rbind(dui_a1,dui_a2) #233 enhancers
dui_a<-dui_a_1[,c(1:5,10)]
dui_a<-dui_a[!duplicated(dui_a),]
unique(dui_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#31 enhancers presentes definidos por previos autores
unique(dui_a$Gene.coincience)#correspondientes a 15 genes

write.table((dui_a), file = " /enhancers_JL/merge_DUI_enhancers_a.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_a[,6]), file = " /enhancers_JL/merge_DUI_enhancers_a_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### ENHANCERS A.GAMBIAE HOMOLOGUES OF THOSE OF DROSOPHILA

enh_b<- read.delim(" /enhancers_JL/S11_novel_regulatory_regions_b_JL", header=T)
enh_b<-enh_b[,c(1:4,7,9:12)]
enh_b$THS.Peak.ID<-ifelse(enh_b$THS.Peak.ID == "",NA,enh_b$THS.Peak.ID)
enh_b<-na.omit(enh_b)#383 I omit those with no associated spout
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
unique(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#25 enhancers presentes definidos por homologia con drosophila
unique(deg_b$Gene.coincience)#correspondientes a 28 genes

write.table((deg_b), file = " /enhancers_JL/merge_DEMG_enhancers_b.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((deg_b[,6]), file = " /enhancers_JL/merge_DEMG_enhancers_b_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

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
unique(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#8 enhancers presentes definidos por previos autores
unique(dui_b$Gene.coincience)#correspondientes a 8 genes

write.table((dui_b), file = " /enhancers_JL/merge_DUI_enhancers_b.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_b[,6]), file = " /enhancers_JL/merge_DUI_enhancers_b_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
#

# ### now I do a merge with the attack peaks to see how many of them are active and also the peaks are High confidence
#
# ### ATAC-seq peaks
# peaks <- read.delim(" /ATAC/results_JL/ATAC-seq_peaks.txt") #193005
# peak<-peaks[peaks$Infections..intersect. == "I1,I2" ,]
# peak<-peak[,c(1:4)]
# peak$Strand<-"+"
# peak$Score<-0
# peak$Chromosome<-paste("AgamP4",peak$Chromosome,sep = "_")
#
# colnames(peak)<-c( "Chromosome"   ,   "Start"    ,       "End"       ,      "THS.Peak.ID" ,"Strand"      ,    "Score"  )
# write.table((peak), file = " /enhancers_JL/HCS_peaks_all.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
#
#
# dega2<-merge(peak,deg_a,by="THS.Peak.ID")
# dega2_hc<-dega2[,c(7:11)]
# dega2_hc<-dega2_hc[!duplicated(dega2_hc),]
# z<-unique(dega2_hc$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#73
# z<-unique(dega2_hc$Gene.coincience)#41
#
#
# degb2<-merge(peak,deg_b,by="THS.Peak.ID")
# degb2_hc<-degb2[,c(7:11)]
# degb2_hc<-degb2_hc[!duplicated(degb2_hc),]
# z<-unique(degb2_hc$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#16
# z<-unique(degb2_hc$Gene.coincience)#18
#
# duia2<-merge(peak,dui_a,by="THS.Peak.ID")
# duia2_hc<-duia2[,c(7:11)]
# duia2_hc<-duia2_hc[!duplicated(duia2_hc),]
# z<-unique(duia2$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#24
# z<-unique(duia2$Gene.coincience)#13
#
#
# duib2<-merge(peak,dui_b,by="THS.Peak.ID")
# duib2_hc<-duib2[,c(7:11)]
# duib2_hc<-duib2_hc[!duplicated(duib2_hc),]
# z<-unique(duib2_hc$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#5
# z<-unique(duib2_hc$Gene.coincience)#5
#
# ### now I'm going to get the missing genes to see if there are any interesting ones.
# no_hc_dega<-deg_a[which(!(deg_a$Gene.coincience %in% dega2_hc$Gene.coincience)),]
# unique(no_hc_dega$Gene.coincience)## 7 genes se pierden
#
# no_hc_degb<-deg_b[which(!(deg_b$Gene.coincience %in% degb2_hc$Gene.coincience)),]
# unique(no_hc_degb$Gene.coincience)## 10 genes are lost
#
# no_hc_duia<-dui_a[which(!(dui_a$Gene.coincience %in% duia2_hc$Gene.coincience)),]
# unique(no_hc_duia$Gene.coincience)## 2 genes are lost
#
# no_hc_duib<-dui_b[which(!(dui_b$Gene.coincience %in% duib2_hc$Gene.coincience)),]
# unique(no_hc_duib$Gene.coincience)## 3 genes are lost
#
#
#
# #####  I will load the genes present in the gtf and cross them with the DEMG and DUI.
# ## in this way I get the genes with the possibility of having ehnacers
#
# genes<- read.delim(" genomic_data_AgamP4/Genes_AgamP4_release_54.bed", header=FALSE, comment.char="#")
# colnames(genes)<-c("Chr","Start","End","gene","Strand","Score")
#
# deg_gff<-merge(deg_all,genes,by="gene")
# dui_gff<-merge(dui_all,genes,by="gene")
#  write.table(deg_gff[,c(9,10,11,1,12,13)], file = " /enhancers_JL/Regions_genes_DEMG_for_intersect_with_enhancers.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
#  write.table(dui_gff[,c(9,10,11,1,12,13)], file = " /enhancers_JL/Regions_genes_DUI_for_intersect_with_enhancers.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
#
#  #### now I do a bedtools intersect of the gene regions with the enhancers and see which ones are inside the genes.
#
# #
# #  bedtools intersect -a  /enhancers_JL/merge_DEMG_enhancers_a.bed -b  /enhancers_JL/Regions_genes_DEMG_for_intersect_with_enhancers.bed -wa -wb -f 0.51 >  /enhancers_JL/merge_DEMG_enhancers_a_gene_regions.bed
# #  bedtools intersect -a  /enhancers_JL/merge_DUI_enhancers_a.bed -b  /enhancers_JL/Regions_genes_DUI_for_intersect_with_enhancers.bed -wa -wb -f 0.51 >  /enhancers_JL/merge_DUI_enhancers_a_gene_regions.bed
# #
# #  bedtools intersect -a  /enhancers_JL/merge_DEMG_enhancers_b.bed -b  /enhancers_JL/Regions_genes_DEMG_for_intersect_with_enhancers.bed -wa -wb -f 0.51 >  /enhancers_JL/merge_DEMG_enhancers_b_gene_regions.bed
# #  bedtools intersect -a  /enhancers_JL/merge_DUI_enhancers_b.bed -b  /enhancers_JL/Regions_genes_DUI_for_intersect_with_enhancers.bed -wa -wb -f 0.51 >  /enhancers_JL/merge_DUI_enhancers_b_gene_regions.bed
# #
#
#
#  ##  once the intersect is done, I load them here again to see how many enhancers and unique genes appear.
#
#
#  mer_demg_b<- read.delim(" /enhancers_JL/merge_DEMG_enhancers_b_gene_regions.bed", header=F)
#  unique(mer_demg_b$V4)#9 ENHANCERS
#  unique(mer_demg_b$V8)#8 GENES
#
#  mer_demg_a<- read.delim(" /enhancers_JL/merge_DEMG_enhancers_a_gene_regions.bed", header=F)
#  unique(mer_demg_a$V4)#72 ENHANCERS
#  unique(mer_demg_a$V8)#38 GENES
#
#  mer_dui_b<- read.delim(" /enhancers_JL/merge_DUI_enhancers_b_gene_regions.bed", header=F)
#  unique(mer_dui_b$V4)#4 ENHANCERS
#  unique(mer_dui_b$V8)#4 GENES
#
#  mer_dui_a<- read.delim(" /enhancers_JL/merge_DUI_enhancers_a_gene_regions.bed", header=F)
#  unique(mer_dui_a$V4)#21 ENHANCERS
#  unique(mer_dui_a$V8)#12 GENES
#

 ######LOCATING THEM IN THE GENOME

 ####### NOW WE ARE GOING TO SEPARATE THE ENHANCERS THAT ARE INSIDE GENES FROM THOSE THAT ARE NOT AND WITH THE REMAINING ONES WE ARE GOING TO CALCULATE THE DISTANCE WITH THE 5'UTR OF THEIR TARGET GENE.


 deg_a_inter2<-deg_a[which(!(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene. %in% mer_demg_a$V4)),]

 unique(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#100
 unique(mer_demg_a$V4)#72
 unique(deg_a_inter2$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#28


 deg_b_inter2<-deg_b[which(!(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates. %in% mer_demg_b$V4)),]

 unique(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#25
 unique(mer_demg_b$V4)#9
 unique(deg_b_inter2$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#16


 dui_a_inter2<-dui_a[which(!(dui_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene. %in% mer_dui_a$V4)),]

 unique(dui_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#31
 unique(mer_dui_a$V4)#21
 unique(dui_a_inter2$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#10

 dui_b_inter2<-dui_b[which(!(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates. %in% mer_dui_b$V4)),]

 unique(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#8
 unique(mer_dui_b$V4)#4
 unique(dui_b_inter2$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#4

 # ### Once we have the intergenes, we are going to calculate their distance from their target. To do this we are going to subtract the END of each enhancer from the 5'UTR(START) of the gene.

 #
 colnames(deg_a_inter2)[6]<-"gene"
 deg_a_inter3<-merge(genes,deg_a_inter2,by="gene")
 colnames(deg_a_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
 deg_a_inter3<- deg_a_inter3[,c(1:10)]
 deg_a_inter3<-unique( deg_a_inter3)
 deg_a_inter3$distance_target<-deg_a_inter3$Start.gene - ((deg_a_inter3$End.enh + deg_a_inter3$Start.enh)/2)
 deg_a_inter3<-deg_a_inter3[!duplicated(deg_a_inter3),]
 sum(deg_a_inter3$distance_target<0) #11
 sum(deg_a_inter3$distance_target>0) #17

 colnames(deg_b_inter2)[6]<-"gene"
 deg_b_inter3<-merge(genes,deg_b_inter2,by="gene")
 colnames(deg_b_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
 deg_b_inter3<- deg_b_inter3[,c(1:10)]
 deg_b_inter3<-unique( deg_b_inter3)
 deg_b_inter3$distance_target<-deg_b_inter3$Start.gene - ((deg_b_inter3$End.enh + deg_b_inter3$Start.enh)/2)
 deg_b_inter3<-deg_b_inter3[!duplicated(deg_b_inter3),]
 sum(deg_b_inter3$distance_target<0) #12
 sum(deg_b_inter3$distance_target>0) #7

 ### in the case of DEMG B (obtained by drosophila) there are more than the 16 that should be because there are 3 enhancers that are associated to two different genes but there are 16 unique enhancers.


 colnames(dui_a_inter2)[6]<-"gene"
 dui_a_inter3<-merge(genes,dui_a_inter2,by="gene")
 colnames(dui_a_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
 dui_a_inter3<- dui_a_inter3[,c(1:10)]
 dui_a_inter3<-unique( dui_a_inter3)
 dui_a_inter3$distance_target<-dui_a_inter3$Start.gene - ((dui_a_inter3$End.enh + dui_a_inter3$Start.enh)/2)
 dui_a_inter3<-dui_a_inter3[!duplicated(dui_a_inter3),]
 sum(dui_a_inter3$distance_target<0) #1
 sum(dui_a_inter3$distance_target>0) #9

 colnames(dui_b_inter2)[6]<-"gene"
 dui_b_inter3<-merge(genes,dui_b_inter2,by="gene")
 colnames(dui_b_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
 dui_b_inter3<- dui_b_inter3[,c(1:10)]
 dui_b_inter3<-unique( dui_b_inter3)
 dui_b_inter3$distance_target<-dui_b_inter3$Start.gene - ((dui_b_inter3$End.enh + dui_b_inter3$Start.enh)/2)
 dui_b_inter3<-dui_b_inter3[!duplicated(dui_b_inter3),]
 sum(dui_b_inter3$distance_target<0) #1
 sum(dui_b_inter3$distance_target>0) #3



 ##### Once they are located I will study their activity, first I will see which tissue each peak corresponds to, and then I will see in each enhancer where it is active.


 colnames(deg_a_1)[4]<-"Enhancer.ID"
 deg_a_1<-deg_a_1[,c(1:5,10)]
 colnames(deg_b_1)[4]<-"Enhancer.ID"
 deg_b_1<-deg_b_1[,c(1:5,10)]
 ### juntamos los dos tipos de enhancers y los separamos por tejidos

deg_act<-rbind(deg_a_1,deg_b_1)
colnames(peaks)[4]<-"THS.Peak.ID"
deg_act<-merge(deg_act,peaks, by="THS.Peak.ID")
deg_act<-deg_act[,c(1:6,10,11)]
deg_act_sg<-deg_act[(deg_act$Tissue == "S. glands"),]
deg_act_mg<-deg_act[(deg_act$Tissue == "Midguts"),]
deg_act_sg<-unique(deg_act_sg)
deg_act_mg<-unique(deg_act_mg)

deg_act_sg2<-deg_act_sg %>% group_by(Enhancer.ID) %>%tally() #96

deg_act_mg2<-deg_act_mg %>% group_by(Enhancer.ID) %>%tally() #92

deg_act2<-merge(deg_act_mg2,deg_act_sg2,by="Enhancer.ID") ### aqui puedo hacer las comparaciones entre numeros pero hay algunos que solo estan presentes en uno de los tejidos

#### of those appearing in both tissues I will see which tissue has more activity and compare it with their expression.

colnames(deg_act2)<-c("Enhancer.ID", "n.mg","n.sg" )
deg_act2$more.activity<-ifelse(deg_act2$n.mg > deg_act2$n.sg,"MG","SG")
deg_act2$more.activity<-ifelse(deg_act2$n.mg == deg_act2$n.sg,"-",deg_act2$more.activity)

deg_act2<-merge(deg_act2,deg_act,by="Enhancer.ID")
deg_all$gene<- gsub("-R.*","",deg_all$ID)
colnames(deg_act2)[9]<-"gene"

deg_act2<-merge(deg_act2,deg_all,by="gene")
deg_act2<-deg_act2[,c(1:11,18)]
deg_act2<-deg_act2[!duplicated(deg_act2$Enhancer.ID),]

deg_act2$more.expresion<-ifelse(deg_act2$log2FoldChange < 0,"MG","SG")
deg_act2$coincidence<-ifelse(deg_act2$more.activity == deg_act2$more.expresion,"YES","-")# 14 de 63 coinciden

#### I will check if activity and expression coincide for those who are only in one of the 2 tissues.


deg_only_mg<-deg_act_mg2[which(!(deg_act_mg2$Enhancer.ID %in% deg_act2$Enhancer.ID)),]
deg_only_mg<-merge(deg_only_mg,deg_act_mg,by="Enhancer.ID")
deg_all$gene<- gsub("-R.*","",deg_all$ID)
colnames(deg_only_mg)[7]<-"gene"

deg_only_mg<-merge(deg_only_mg,deg_all,by="gene")
deg_only_mg<-deg_only_mg[,c(1:9,17)]
deg_only_mg<-deg_only_mg[!duplicated(deg_only_mg$Enhancer.ID),]### aqui veo que los que tienen mas actividad en MG tienen que tener un logFC<0. de los 29 coinciden 15 y no coinciden 14
sum(deg_only_mg$log2FoldChange > 0)#15 up expressed en SG (no coinciden)
sum(deg_only_mg$log2FoldChange < 0)# 14  up expressed en MG (coinciden)

deg_only_sg<-deg_act_sg2[which(!(deg_act_sg2$Enhancer.ID %in% deg_act2$Enhancer.ID)),]##
deg_only_sg<-merge(deg_only_sg,deg_act_sg,by="Enhancer.ID")
deg_all$gene<- gsub("-R.*","",deg_all$ID)
colnames(deg_only_sg)[7]<-"gene"

deg_only_sg<-merge(deg_only_sg,deg_all,by="gene")
deg_only_sg<-deg_only_sg[,c(1:9,16)]
deg_only_sg<-deg_only_sg[!duplicated(deg_only_sg$Enhancer.ID),] ### aqui veo que los que tienen mas actividad en SG tienen que tener un logFC>0. de los 33 coinciden 19 y no coinciden 14

sum(deg_only_sg$log2FoldChange > 0)#19 up expressed en SG (coinciden)
sum(deg_only_sg$log2FoldChange < 0)# 14 up expressed en MG (no coinciden)


#I do the same for DUI

colnames(dui_a_1)[4]<-"Enhancer.ID"
dui_a_1<-dui_a_1[,c(1:5,10)]
colnames(dui_b_1)[4]<-"Enhancer.ID"
dui_b_1<-dui_b_1[,c(1:5,10)]

dui_act<-rbind(dui_a_1,dui_b_1)
dui_act<-merge(dui_act,peaks, by="THS.Peak.ID")
dui_act<-dui_act[,c(1:6,10,11)]
dui_act_sg<-dui_act[(dui_act$Tissue == "S. glands"),]
dui_act_mg<-dui_act[(dui_act$Tissue == "Midguts"),]
dui_act_sg<-unique(dui_act_sg)
dui_act_mg<-unique(dui_act_mg)

dui_act_sg2<-dui_act_sg %>% group_by(Enhancer.ID) %>%tally() #30

dui_act_mg2<-dui_act_mg %>% group_by(Enhancer.ID) %>%tally() #32

dui_act2<-merge(dui_act_mg2,dui_act_sg2,by="Enhancer.ID") ### aqui puedo hacer las comparaciones entre numeros pero hay algunos que solo estan presentes en uno de los tejidos

#### of those that appear in both tissues I will see which tissue has more activity and compare it with its use.

colnames(dui_act2)<-c("Enhancer.ID", "n.mg","n.sg" )
dui_act2$more.activity<-ifelse(dui_act2$n.mg > dui_act2$n.sg,"MG","SG")
dui_act2$more.activity<-ifelse(dui_act2$n.mg == dui_act2$n.sg,"-",dui_act2$more.activity)

dui_act2<-merge(dui_act2,dui_act,by="Enhancer.ID")
dui_all$gene<- gsub("-R.*","",dui_all$ID)
colnames(dui_act2)[9]<-"gene"

dui_act2<-merge(dui_act2,dui_all,by="gene")
dui_act2<-dui_act2[,c(1:11,18)]
dui_act2<-dui_act2[!duplicated(dui_act2$Enhancer.ID),]

dui_act2$more.expresion<-ifelse(dui_act2$dIF < 0,"MG","SG")
dui_act2$coincidence<-ifelse(dui_act2$more.activity == dui_act2$more.expresion,"YES","-") ## 6 coincidencias de 23

#### I will check if there is a match between activity and use in those that are only in one of the 2 fabrics.

dui_only_mg<-dui_act_mg2[which(!(dui_act_mg2$Enhancer.ID %in% dui_act2$Enhancer.ID)),]
dui_only_mg<-merge(dui_only_mg,dui_act_mg,by="Enhancer.ID")
dui_all$gene<- gsub("-R.*","",dui_all$ID)
colnames(dui_only_mg)[7]<-"gene"

dui_only_mg<-merge(dui_only_mg,dui_all,by="gene")
dui_only_mg<-dui_only_mg[,c(1:9,23)]
dui_only_mg<-dui_only_mg[!duplicated(dui_only_mg$Enhancer.ID),] ### aqui veo que los que tienen mas actividad en MG tienen que tener un dIF<0. de los 9 coinciden 4 y no coinciden 5
sum(dui_only_mg$dIF > 0)#5 up expressed in SG (do not coincide)
sum(dui_only_mg$dIF < 0)# 4 up expressed in MG  (coincide)

dui_only_sg<-dui_act_sg2[which(!(dui_act_sg2$Enhancer.ID %in% dui_act2$Enhancer.ID)),]##
dui_only_sg<-merge(dui_only_sg,dui_act_sg,by="Enhancer.ID")
dui_all$gene<- gsub("-R.*","",dui_all$ID)
colnames(dui_only_sg)[7]<-"gene"

dui_only_sg<-merge(dui_only_sg,dui_all,by="gene")
dui_only_sg<-dui_only_sg[,c(1:9,16)]
dui_only_sg<-dui_only_sg[!duplicated(dui_only_sg$Enhancer.ID),] ### aqui veo que los que tienen mas actividad en SG tienen que tener un logFC>0. de los 7 coinciden 3 y no coinciden 4
sum(dui_only_sg$dIF > 0)# 3 up expressed in SG ( coincide)
sum(dui_only_sg$dIF < 0)# 4  up expressed in MG  (do not coincide)


deg_act_d<-deg_act[,c(2:6)]
deg_act_d<-deg_act_d[!duplicated(deg_act_d$Enhancer.ID),] #125

dui_act_d<-dui_act[,c(2:6)]
dui_act_d<-dui_act_d[!duplicated(dui_act_d$Enhancer.ID),] #39

write.table((deg_act_d), file = " /enhancers_JL/enhancers_activos_DEMG_A_B.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_act_d), file = " /enhancers_JL/enhancers_activos_DUI_A_B.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### now I am going to calculate the coverage of these enhancers with bedtools intersect


 # bedtools intersect -a  /enhancers_JL/enhancers_activos_DEMG_A_B.bed -b  /ATAC/R10-8_d7_nucfree_mapq_10_merged.bam -c >  /enhancers_JL/coverage_enhancers_DEMG_A_B_MG.bed
 # bedtools intersect -a  /enhancers_JL/enhancers_activos_DEMG_A_B.bed -b  /ATAC/R5-1_d14_nucfree_mapq_10_merged.bam -c >  /enhancers_JL/coverage_enhancers_DEMG_A_B_SG.bed
 #
 # bedtools intersect -a  /enhancers_JL/enhancers_activos_DUI_A_B.bed -b  /ATAC/R10-8_d7_nucfree_mapq_10_merged.bam -c >  /enhancers_JL/coverage_enhancers_DUI_A_B_MG.bed
 # bedtools intersect -a  /enhancers_JL/enhancers_activos_DUI_A_B.bed -b  /ATAC/R5-1_d14_nucfree_mapq_10_merged.bam -c >  /enhancers_JL/coverage_enhancers_DUI_A_B_SG.bed

##  I load the tables and counts

cov_demg_mg<- read.delim(" /enhancers_JL/coverage_enhancers_DEMG_A_B_MG.bed", header=F)
colnames(cov_demg_mg)<-c("Chromosome","Start","End","Enhancer.ID","gene","Counts.mg")
cov_demg_sg<- read.delim(" /enhancers_JL/coverage_enhancers_DEMG_A_B_SG.bed", header=F)
colnames(cov_demg_sg)<-c("Chromosome","Start","End","Enhancer.ID","Gene","Counts.sg")

cov_dui_mg<- read.delim(" /enhancers_JL/coverage_enhancers_DUI_A_B_MG.bed", header=F)
colnames(cov_dui_mg)<-c("Chromosome","Start","End","Enhancer.ID","gene","Counts.mg")
cov_dui_sg<- read.delim(" /enhancers_JL/coverage_enhancers_DUI_A_B_SG.bed", header=F)
colnames(cov_dui_sg)<-c("Chromosome","Start","End","Enhancer.ID","Gene","Counts.sg")

H<-cov_dui_sg[!duplicated(cov_dui_sg$Enhancer.ID),]

cov_demg<-merge(cov_demg_mg,cov_demg_sg,by="Enhancer.ID")
cov_demg<-cov_demg[,c(1:6,11)]
cov_demg$more.activity<- ifelse(cov_demg$Counts.mg > cov_demg$Counts.sg,"MG","SG")
cov_demg$more.activity<- ifelse(cov_demg$Counts.mg == cov_demg$Counts.sg,"-",cov_demg$more.activity)



cov_demg<-merge(cov_demg,deg_all,by="gene")
cov_demg<-cov_demg[,c(1:8,15)]
cov_demg<-cov_demg[!duplicated(cov_demg$Enhancer.ID),]

cov_demg$more.expresion<-ifelse(cov_demg$log2FoldChange < 0,"MG","SG")
cov_demg$coincidence<-ifelse(cov_demg$more.activity == cov_demg$more.expresion,"YES","-")
sum(cov_demg$coincidence == "YES")#69 de 125
sum(cov_demg$coincidence == "-")#56 de 125


cov_dui<-merge(cov_dui_mg,cov_dui_sg,by="Enhancer.ID")
cov_dui<-cov_dui[,c(1:6,11)]

cov_dui$more.activity<- ifelse(cov_dui$Counts.mg > cov_dui$Counts.sg,"MG","SG")
cov_dui$more.activity<- ifelse(cov_dui$Counts.mg == cov_dui$Counts.sg,"-",cov_dui$more.activity)

cov_dui<-merge(cov_dui,dui_all,by="gene")
cov_dui<-cov_dui[,c(1:8,15)]
cov_dui<-cov_dui[!duplicated(cov_dui$Enhancer.ID),]

cov_dui$more.expresion<-ifelse(cov_dui$dIF < 0,"MG","SG")
cov_dui$coincidence<-ifelse(cov_dui$more.activity == cov_dui$more.expresion,"YES","-")
sum(cov_dui$coincidence == "YES")#19 de 39
sum(cov_dui$coincidence == "-")#15 de 39
