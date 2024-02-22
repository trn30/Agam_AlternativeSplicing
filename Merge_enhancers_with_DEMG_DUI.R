
### add to the DEMG and DUI the positions of the whole genomic region
deg<- read.delim("/DEMG/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05.txt", header=T) #4933 DEG
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
dui_all<-dui_all[,c(11,12,18,1,19,20,6)]
dui_all$gene<-gsub("-R.*","",dui_all$ID)
unique(dui_all$gene) #144 genes

write.table((deg_all), file = "/Enhancers/DEMG_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_all), file = "/Enhancers/DUI_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)


#### A.GAMBIAE ENHANCERS OBTAINED BY PREVIOUS AUTHORES

enh_a<- read.delim("/enhancers_JL/S11_novel_regulatory_regions_a_JL", header=T)
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


deg_a_1<-rbind(deg_a1,deg_a2) #640
deg_a<-deg_a_1[,c(1:5,10)]
deg_a<-deg_a[!duplicated(deg_a),]
unique(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#88 enhancers present as defined by previous authors
unique(deg_a$Gene.coincience)#corresponding to 44 genes

write.table((deg_a), file = "/Enhancers/merge_DEMG_enhancers_a.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((deg_a[,6]), file = "/Enhancers/merge_DEMG_enhancers_a_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

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
unique(dui_a$Gene.coincience)#corresponding to 17 genes

write.table((dui_a), file = "/Enhancers/merge_DUI_enhancers_a.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_a[,6]), file = "/Enhancers/merge_DUI_enhancers_a_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### ENHANCERS A.GAMBIAE HOMOLOGOUS TO THOSE OF DROSOPHILA

enh_b<- read.delim("/enhancers_JL/S11_novel_regulatory_regions_b_JL", header=T)
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

deg_b_1<-rbind(deg_b1,deg_b2)#250
deg_b<-deg_b_1[,c(1:5,10)]
deg_b<-deg_b[!duplicated(deg_b),]
unique(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#21 enhancers present defined by homology with drosophila
unique(deg_b$Gene.coincience)# corresponding to 22 genes

write.table((deg_b), file = "/Enhancers/merge_DEMG_enhancers_b.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((deg_b[,6]), file = "/Enhancers/merge_DEMG_enhancers_b_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

dui_b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% dui_all$gene),]
b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% dui_all$gene),7]
b1<-dplyr::pull(b1, Annotated.Gene.ID)
dui_b1$Gene.coincience<-b1

dui_b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% dui_all$gene),]
b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% dui_all$gene),8]
b2<-dplyr::pull(b2, Target.by.others..An..gambiae.ortholog.)
dui_b2$Gene.coincience<-b2

dui_b_1<-rbind(dui_b1,dui_b2) ##228
dui_b<-dui_b_1[,c(1:5,10)]
dui_b<-dui_b[!duplicated(dui_b),]
unique(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#13 enhancers present defined by homology with drosophila
unique(dui_b$Gene.coincience)# corresponding to 12 genes

write.table((dui_b), file = "/Enhancers/merge_DUI_enhancers_b.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_b[,6]), file = "/Enhancers/merge_DUI_enhancers_b_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)


# ### now I do a merge with the attack peaks to see how many of them are active and also the peaks are High confidence

### ATAC peaks

peaks <- read.delim("/ATAC/results_JL/ATAC-seq_peaks.txt") #193005
peak<-peaks[peaks$Infections..intersect. == "I1,I2" ,]
peak<-peak[,c(1:4)]
peak$Strand<-"+"
peak$Score<-0
peak$Chromosome<-paste("AgamP4",peak$Chromosome,sep = "_")

colnames(peak)<-c( "Chromosome"   ,   "Start"    ,       "End"       ,      "THS.Peak.ID" ,"Strand"      ,    "Score"  )
write.table((peak), file = "/Enhancers/HCS_peaks_all.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)


# ##### I will load the genes present in the gtf and cross them with the DEMG and DUI.
# ## in this way I get the genes with the possibility of having ehnacers.

genes<- read.delim("/Genomes/genomic_data_AgamP4/Genes_AgamP4_release_54.bed", header=FALSE, comment.char="#")
colnames(genes)<-c("Chr","Start","End","gene","Strand","Score")

deg_gff<-merge(deg_all,genes,by="gene")
dui_gff<-merge(dui_all,genes,by="gene")
write.table(deg_gff[,c(9,10,11,1,12,13)], file = "/Enhancers/Regions_genes_DEMG_for_intersect_with_enhancers.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table(dui_gff[,c(9,10,11,1,12,13)], file = "/Enhancers/Regions_genes_DUI_for_intersect_with_enhancers.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)


#### now I do a bedtools intersect of the gene regions with the enhancers and see which ones are inside the genes.
#### Linux terminal

 # bedtools intersect -a /Enhancers/merge_DEMG_enhancers_a.bed -b /Enhancers/Regions_genes_DEMG_for_intersect_with_enhancers.bed -wa -wb -f 0.51 > /Enhancers/merge_DEMG_enhancers_a_gene_regions.bed
 # bedtools intersect -a /Enhancers/merge_DUI_enhancers_a.bed -b /Enhancers/Regions_genes_DUI_for_intersect_with_enhancers.bed -wa -wb -f 0.51 > /Enhancers/merge_DUI_enhancers_a_gene_regions.bed

 # bedtools intersect -a /Enhancers/merge_DEMG_enhancers_b.bed -b /Enhancers/Regions_genes_DEMG_for_intersect_with_enhancers.bed -wa -wb -f 0.51 > /Enhancers/merge_DEMG_enhancers_b_gene_regions.bed
 # bedtools intersect -a /Enhancers/merge_DUI_enhancers_b.bed -b /Enhancers/Regions_genes_DUI_for_intersect_with_enhancers.bed -wa -wb -f 0.51 > /Enhancers/merge_DUI_enhancers_b_gene_regions.bed


## once the intersect is done I load them here again to see how many enhancers and unique genes appear.

 mer_demg_b<- read.delim("/Enhancers/merge_DEMG_enhancers_b_gene_regions.bed", header=F)
 unique(mer_demg_b$V4)#8 ENHANCERS
 unique(mer_demg_b$V8)#7 GENES

 mer_demg_a<- read.delim("/Enhancers/merge_DEMG_enhancers_a_gene_regions.bed", header=F)
 unique(mer_demg_a$V4)#63 ENHANCERS
 unique(mer_demg_a$V8)#34 GENES

 mer_dui_b<- read.delim("/Enhancers/merge_DUI_enhancers_b_gene_regions.bed", header=F)
 unique(mer_dui_b$V4)#5 ENHANCERS
 unique(mer_dui_b$V8)#5 GENES

 mer_dui_a<- read.delim("/Enhancers/merge_DUI_enhancers_a_gene_regions.bed", header=F)
 unique(mer_dui_a$V4)#20 ENHANCERS
 unique(mer_dui_a$V8)#11 GENES

######LOCALISE THEM IN THE GENOME

####### NOW WE WILL SEPARATE THE ENHANCERS THAT ARE INSIDE GENES FROM THOSE THAT ARE NOT AND WITH THE REMAINING ONES WE WILL CALCULATE THE DISTANCE TO THE 5'UTR OF THEIR TARGET GENE.

deg_a_inter2<-deg_a[which(!(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene. %in% mer_demg_a$V4)),]

unique(deg_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#88
unique(mer_demg_a$V4)#63
unique(deg_a_inter2$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#25


deg_b_inter2<-deg_b[which(!(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates. %in% mer_demg_b$V4)),]

unique(deg_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#21
unique(mer_demg_b$V4)#8
unique(deg_b_inter2$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#13


dui_a_inter2<-dui_a[which(!(dui_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene. %in% mer_dui_a$V4)),]

unique(dui_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#33
unique(mer_dui_a$V4)#20
unique(dui_a_inter2$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#13

dui_b_inter2<-dui_b[which(!(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates. %in% mer_dui_b$V4)),]

unique(dui_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#13
unique(mer_dui_b$V4)#5
unique(dui_b_inter2$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#8

### once we have the intergenes we are going to calculate their distance to their target. To do this we will subtract the END of each enhancer from the 5'UTR(START) of the gene.
 genes<- read.delim("/Genomes/genomic_data_AgamP4/Genes_AgamP4_release_54.bed", header=FALSE, comment.char="#")
 colnames(genes)<-c("Chr","Start","End","gene","Strand","Score")

colnames(deg_a_inter2)[6]<-"gene"
deg_a_inter3<-merge(genes,deg_a_inter2,by="gene")
colnames(deg_a_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
deg_a_inter3<- deg_a_inter3[,c(1:10)]
deg_a_inter3<-unique( deg_a_inter3)
deg_a_inter3$distance_target<-deg_a_inter3$Start.gene - ((deg_a_inter3$End.enh + deg_a_inter3$Start.enh)/2)
deg_a_inter3<-deg_a_inter3[!duplicated(deg_a_inter3),]
sum(deg_a_inter3$distance_target<0) #10
sum(deg_a_inter3$distance_target>0) #15

colnames(deg_b_inter2)[6]<-"gene"
deg_b_inter3<-merge(genes,deg_b_inter2,by="gene")
colnames(deg_b_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
deg_b_inter3<- deg_b_inter3[,c(1:10)]
deg_b_inter3<-unique( deg_b_inter3)
deg_b_inter3$distance_target<-deg_b_inter3$Start.gene - ((deg_b_inter3$End.enh + deg_b_inter3$Start.enh)/2)
deg_b_inter3<-deg_b_inter3[!duplicated(deg_b_inter3),]
sum(deg_b_inter3$distance_target<0) #9
sum(deg_b_inter3$distance_target>0) #6

### in the case of DEMG B (obtained by drosophila) there are more than the 16 that should be because there are 3 enhancers that are associated to two different genes but there are 16 unique enhancers.

colnames(dui_a_inter2)[6]<-"gene"
dui_a_inter3<-merge(genes,dui_a_inter2,by="gene")
colnames(dui_a_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
dui_a_inter3<- dui_a_inter3[,c(1:10)]
dui_a_inter3<-unique( dui_a_inter3)
dui_a_inter3$distance_target<-dui_a_inter3$Start.gene - ((dui_a_inter3$End.enh + dui_a_inter3$Start.enh)/2)
dui_a_inter3<-dui_a_inter3[!duplicated(dui_a_inter3),]
sum(dui_a_inter3$distance_target<0) #4
sum(dui_a_inter3$distance_target>0) #9

colnames(dui_b_inter2)[6]<-"gene"
dui_b_inter3<-merge(genes,dui_b_inter2,by="gene")
colnames(dui_b_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
dui_b_inter3<- dui_b_inter3[,c(1:10)]
dui_b_inter3<-unique( dui_b_inter3)
dui_b_inter3$distance_target<-dui_b_inter3$Start.gene - ((dui_b_inter3$End.enh + dui_b_inter3$Start.enh)/2)
dui_b_inter3<-dui_b_inter3[!duplicated(dui_b_inter3),]
sum(dui_b_inter3$distance_target<0) #4
sum(dui_b_inter3$distance_target>0) #4

##### once they are located I will study their activity, first I will see which tissue each peak corresponds to, and then I will see in each enhancer where it is active.

colnames(deg_a_1)[4]<-"Enhancer.ID"
deg_a_1<-deg_a_1[,c(1:5,10)]
colnames(deg_b_1)[4]<-"Enhancer.ID"
deg_b_1<-deg_b_1[,c(1:5,10)]
### juntamos los dos tipos de enhancers y los separamos por tejidos
peaks <- read.delim("/ATAC/results_JL/ATAC-seq_peaks.txt") #193005
deg_act<-rbind(deg_a_1,deg_b_1)
colnames(peaks)[4]<-"THS.Peak.ID"
deg_act<-merge(deg_act,peaks, by="THS.Peak.ID")
deg_act<-deg_act[,c(1:6,10,11)]
deg_act_sg<-deg_act[(deg_act$Tissue == "S. glands"),]
deg_act_mg<-deg_act[(deg_act$Tissue == "Midguts"),]
deg_act_sg<-unique(deg_act_sg)
deg_act_mg<-unique(deg_act_mg)

deg_act_sg2<-deg_act_sg %>% group_by(Enhancer.ID) %>%tally() #82

deg_act_mg2<-deg_act_mg %>% group_by(Enhancer.ID) %>%tally() #84

deg_act2<-merge(deg_act_mg2,deg_act_sg2,by="Enhancer.ID") ### here I can make comparisons between numbers but there are some that are only present in one of the tissues.
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
deg_act2$coincidence<-ifelse(deg_act2$more.activity == deg_act2$more.expresion,"YES","-")# 16 out of 57 agree

#### I will check if activity and expression coincide for those who are only in one of the 2 tissues.
deg_only_mg<-deg_act_mg2[which(!(deg_act_mg2$Enhancer.ID %in% deg_act2$Enhancer.ID)),]
deg_only_mg<-merge(deg_only_mg,deg_act_mg,by="Enhancer.ID")
deg_all$gene<- gsub("-R.*","",deg_all$ID)
colnames(deg_only_mg)[7]<-"gene"

deg_only_mg<-merge(deg_only_mg,deg_all,by="gene")
deg_only_mg<-deg_only_mg[,c(1:9,17)]
deg_only_mg<-deg_only_mg[!duplicated(deg_only_mg$Enhancer.ID),]### here I see that those with more activity in MG have to have a logFC<0. 15 of the 29 match and 14 do not match.
sum(deg_only_mg$log2FoldChange > 0)# 13 up expressed in SG (do not match)
sum(deg_only_mg$log2FoldChange < 0)# 14  up expressed in MG ( match)

deg_only_sg<-deg_act_sg2[which(!(deg_act_sg2$Enhancer.ID %in% deg_act2$Enhancer.ID)),]##
deg_only_sg<-merge(deg_only_sg,deg_act_sg,by="Enhancer.ID")
deg_all$gene<- gsub("-R.*","",deg_all$ID)
colnames(deg_only_sg)[7]<-"gene"

deg_only_sg<-merge(deg_only_sg,deg_all,by="gene")
deg_only_sg<-deg_only_sg[,c(1:9,16)]
deg_only_sg<-deg_only_sg[!duplicated(deg_only_sg$Enhancer.ID),] ### here I see that those who have more activity in SG must have a logFC>0. 19 of the 33 match and 14 do not match.

sum(deg_only_sg$log2FoldChange > 0)#11 up expressed in SG (match)
sum(deg_only_sg$log2FoldChange < 0)# 14 up expressed in MG (do not match)


# I do the same for DUI

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

dui_act_sg2<-dui_act_sg %>% group_by(Enhancer.ID) %>%tally() #34

dui_act_mg2<-dui_act_mg %>% group_by(Enhancer.ID) %>%tally() #39

dui_act2<-merge(dui_act_mg2,dui_act_sg2,by="Enhancer.ID") ### here I can make comparisons between numbers but there are some that are only present in one of the tissues.
#### of those that appear in both tissues I will see which tissue has more activity and compare it with its use.colnames(dui_act2)<-c("Enhancer.ID", "n.mg","n.sg" )
dui_act2$more.activity<-ifelse(dui_act2$n.mg > dui_act2$n.sg,"MG","SG")
dui_act2$more.activity<-ifelse(dui_act2$n.mg == dui_act2$n.sg,"-",dui_act2$more.activity)

dui_act2<-merge(dui_act2,dui_act,by="Enhancer.ID")
dui_all$gene<- gsub("-R.*","",dui_all$ID)
colnames(dui_act2)[9]<-"gene"

dui_act2<-merge(dui_act2,dui_all,by="gene")
dui_act2<-dui_act2[,c(1:11,18)]
dui_act2<-dui_act2[!duplicated(dui_act2$Enhancer.ID),]

dui_act2$more.expresion<-ifelse(dui_act2$dIF < 0,"MG","SG")
dui_act2$coincidence<-ifelse(dui_act2$more.activity == dui_act2$more.expresion,"YES","-") ## 10 matches out of 27

#### voy a mirar si coinciden actividad y uso en los que solo estan en uno de los 2 tejidos
dui_only_mg<-dui_act_mg2[which(!(dui_act_mg2$Enhancer.ID %in% dui_act2$Enhancer.ID)),]
dui_only_mg<-merge(dui_only_mg,dui_act_mg,by="Enhancer.ID")
dui_all$gene<- gsub("-R.*","",dui_all$ID)
colnames(dui_only_mg)[7]<-"gene"

dui_only_mg<-merge(dui_only_mg,dui_all,by="gene")
dui_only_mg<-dui_only_mg[,c(1:9,16)]
dui_only_mg<-dui_only_mg[!duplicated(dui_only_mg$Enhancer.ID),] ### here I see that those with more activity in MG have to have a dIF<0. 4 of the 9 match and 5 do not match.
sum(dui_only_mg$dIF > 0)#5 up expressed en SG (do not match)
sum(dui_only_mg$dIF < 0)#7  up expressed en MG  (match)

dui_only_sg<-dui_act_sg2[which(!(dui_act_sg2$Enhancer.ID %in% dui_act2$Enhancer.ID)),]##
dui_only_sg<-merge(dui_only_sg,dui_act_sg,by="Enhancer.ID")
dui_all$gene<- gsub("-R.*","",dui_all$ID)
colnames(dui_only_sg)[7]<-"gene"

dui_only_sg<-merge(dui_only_sg,dui_all,by="gene")
dui_only_sg<-dui_only_sg[,c(1:9,16)]
dui_only_sg<-dui_only_sg[!duplicated(dui_only_sg$Enhancer.ID),] ### here I see that those with more activity in SG have to have a logFC>0. 3 of the 7 match and 4 do not match.
sum(dui_only_sg$dIF > 0)# 6 up expressed en SG (match)
sum(dui_only_sg$dIF < 0)# 1  up expressed en MG  (do not match)


deg_act_d<-deg_act[,c(2:6)]
deg_act_d<-deg_act_d[!duplicated(deg_act_d$Enhancer.ID),] #109 corresponding to 61 genes out of the 307 genes that make up the 392 DEMG isoforms.

dui_act_d<-dui_act[,c(2:6)]
dui_act_d<-dui_act_d[!duplicated(dui_act_d$Enhancer.ID),] #46 corresponding to 26 genes out of the 145 genes making up the 247 DUI isoforms


write.table((deg_act_d), file = "/Enhancers/enhancers_activos_DEMG_A_B.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((dui_act_d), file = "/Enhancers/enhancers_activos_DUI_A_B.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### now I am going to calculate the coverage of these enhancers with bedtools intersect
# Linux terminal

# bedtools intersect -a /Enhancers/enhancers_activos_DEMG_A_B.bed -b /ATAC/R10-8_d7_nucfree_mapq_10_merged.bam -c > /Enhancers/coverage_enhancers_DEMG_A_B_MG.bed
# bedtools intersect -a /Enhancers/enhancers_activos_DEMG_A_B.bed -b /ATAC/R5-1_d14_nucfree_mapq_10_merged.bam -c > /Enhancers/coverage_enhancers_DEMG_A_B_SG.bed

# bedtools intersect -a /Enhancers/enhancers_activos_DUI_A_B.bed -b /ATAC/R10-8_d7_nucfree_mapq_10_merged.bam -c > /Enhancers/coverage_enhancers_DUI_A_B_MG.bed
# bedtools intersect -a /Enhancers/enhancers_activos_DUI_A_B.bed -b /ATAC/R5-1_d14_nucfree_mapq_10_merged.bam -c > /Enhancers/coverage_enhancers_DUI_A_B_SG.bed

## I load tables and counts

cov_demg_mg<- read.delim("/Enhancers/coverage_enhancers_DEMG_A_B_MG.bed", header=F)
colnames(cov_demg_mg)<-c("Chromosome","Start","End","Enhancer.ID","gene","Counts.mg")
cov_demg_sg<- read.delim("/Enhancers/coverage_enhancers_DEMG_A_B_SG.bed", header=F)
colnames(cov_demg_sg)<-c("Chromosome","Start","End","Enhancer.ID","Gene","Counts.sg")

cov_dui_mg<- read.delim("/Enhancers/coverage_enhancers_DUI_A_B_MG.bed", header=F)
colnames(cov_dui_mg)<-c("Chromosome","Start","End","Enhancer.ID","gene","Counts.mg")
cov_dui_sg<- read.delim("/Enhancers/coverage_enhancers_DUI_A_B_SG.bed", header=F)
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
sum(cov_demg$coincidence == "YES")#62 of 109
sum(cov_demg$coincidence == "-")#47 of 109


cov_dui<-merge(cov_dui_mg,cov_dui_sg,by="Enhancer.ID")
cov_dui<-cov_dui[,c(1:6,11)]

cov_dui$more.activity<- ifelse(cov_dui$Counts.mg > cov_dui$Counts.sg,"MG","SG")
cov_dui$more.activity<- ifelse(cov_dui$Counts.mg == cov_dui$Counts.sg,"-",cov_dui$more.activity)

cov_dui<-merge(cov_dui,dui_all,by="gene")
cov_dui<-cov_dui[,c(1:8,15)]
cov_dui<-cov_dui[!duplicated(cov_dui$Enhancer.ID),]

cov_dui$more.expresion<-ifelse(cov_dui$dIF < 0,"MG","SG")
cov_dui$coincidence<-ifelse(cov_dui$more.activity == cov_dui$more.expresion,"YES","-")
sum(cov_dui$coincidence == "YES")#27 of 46
sum(cov_dui$coincidence == "-")#19 of 46
