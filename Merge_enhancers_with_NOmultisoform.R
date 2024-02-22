################## which enhancers matches with the NOmultisoformgenes

all_no_multi<- read.delim("/Genomes/genomic_data_AgamP4/Isoforms_NO_multisoform_genes_release54.bed", header=T)
colnames(all_no_multi)<-c("Chr","Start","End","ID","Strand","Score","gene")

### add to the DEMG and DUI the positions of the whole genomic region
region_pro<- read.delim("/Genomes/genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed", header=F)
colnames(region_pro)<-c("Chr","Start","End","ID","Strand","Score")


all<-merge(all_no_multi,region_pro,by="ID")
all<-all[,c(8,9,4,1,11,12)]
all$gene<-gsub("-R.*","",all$ID)

write.table((all), file = "/Enhancers/NOmultisoform_coordinates_all_region.txt",quote = FALSE, sep="\t", row.names = F , col.names = F)


#### A.GAMBIAE ENHANCERS OBTAINED BY PREVIOUS AUTHORES

enh_a<- read.delim("_def/enhancers_JL/S11_novel_regulatory_regions_a_JL", header=T)
enh_a<-enh_a[,c(1:4,7,9:12)]
enh_a$THS.Peak.ID<-ifelse(enh_a$THS.Peak.ID == "",NA,enh_a$THS.Peak.ID)
enh_a<-na.omit(enh_a)#1883 omito las que no tiene pico asociado
library(tidyr)
enh_a<-separate_rows(enh_a,Target.by.others, sep = ",")
enh_a$Target.coincidence <-enh_a$Annotated.Gene.ID == enh_a$Target.by.others
enh_a$Chromosome<-paste("AgamP4",enh_a$Chromosome,sep = "_")
unique(enh_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#706

# now I look at how many coincide in each column of possibilities and take to the gene.coincidence column the ones in the deg group.
all_a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% all$gene),]
a1<-enh_a[which(enh_a$Annotated.Gene.ID %in% all$gene),7]
a1<-dplyr::pull(a1, Annotated.Gene.ID)
all_a1$Gene.coincience<-a1

all_a2<-enh_a[which(enh_a$Target.by.others %in% all$gene),]
a2<-enh_a[which(enh_a$Target.by.others %in% all$gene),8]
a2<-dplyr::pull(a2, Target.by.others)
all_a2$Gene.coincience<-a2


all_a_1<-rbind(all_a1,all_a2) #4545
all_a<-all_a_1[,c(1:5,10)]
all_a<-all_a[!duplicated(all_a),]
unique(all_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#616 enhancers present as defined by previous authorss
unique(all_a$Gene.coincience)#corresponding to 621 genes

write.table((all_a), file = "/Enhancers/merge_NOmultisoforms_enhancers_a.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((all_a[,6]), file = "/Enhancers/merge_NOmultisoforms_enhancers_a_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)


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

all_b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% all$gene),]
b1<-enh_b[which(enh_b$Annotated.Gene.ID %in% all$gene),7]
b1<-dplyr::pull(b1, Annotated.Gene.ID)
all_b1$Gene.coincience<-b1

all_b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% all$gene),]
b2<-enh_b[which(enh_b$Target.by.others..An..gambiae.ortholog. %in% all$gene),8]
b2<-dplyr::pull(b2, Target.by.others..An..gambiae.ortholog. )
all_b2$Gene.coincience<-b2

all_b_1<-rbind(all_b1,all_b2)#4966
all_b<-all_b_1[,c(1:5,10)]
all_b<-all_b[!duplicated(all_b),]
unique(all_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#91 enhancers present defined by homology with drosophila
unique(all_b$Gene.coincience)#corresponding to 644 genes

write.table((all_b), file = "/Enhancers/merge_NOmultisoforms_enhancers_b.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((all_b[,6]), file = "/Enhancers/merge_NOmultisoforms_enhancers_b_VB.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)


# ##### I will load the genes present in the gtf and cross the NO multisoform genes
# ## this way I get the genes with the possibility of having ehnacers

genes<- read.delim("/Genomes/genomic_data_AgamP4/Genes_AgamP4_release_54.bed", header=FALSE, comment.char="#")
colnames(genes)<-c("Chr","Start","End","gene","Strand","Score")

all_gff<-merge(all,genes,by="gene")
write.table(all_gff[,c(8,9,10,1,11,12)], file = "/Enhancers/Regions_genes_NOmultisoforms_for_intersect_with_enhancers.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### now I do a bedtools intersect of the gene regions with the enhancers and see which ones are inside the genes.
# Linux terminal
# bedtools intersect -a /Enhancers/merge_NOmultisoforms_enhancers_a.bed -b /Enhancers/Regions_genes_NOmultisoforms_for_intersect_with_enhancers.bed -wa -wb -f 0.51 > /Enhancers/merge_NOmultisoforms_enhancers_a_gene_regions.bed

# bedtools intersect -a /Enhancers/merge_NOmultisoforms_enhancers_b.bed -b /Enhancers/Regions_genes_NOmultisoforms_for_intersect_with_enhancers.bed -wa -wb -f 0.51 > /Enhancers/merge_NOmultisoforms_enhancers_b_gene_regions.bed

## once the intersect is done I load them here again to see how many enhancers and unique genes appear.

mer_demg_b<- read.delim("/Enhancers/merge_NOmultisoforms_enhancers_b_gene_regions.bed", header=F)
unique(mer_demg_b$V4)#49 ENHANCERS
unique(mer_demg_b$V8)#48 GENES

mer_demg_a<- read.delim("/Enhancers/merge_NOmultisoforms_enhancers_a_gene_regions.bed", header=F)
unique(mer_demg_a$V4)#292 ENHANCERS
unique(mer_demg_a$V8)#209 GENES


######LOCALISE THEM IN THE GENOME

####### NOW WE WILL SEPARATE THE ENHANCERS THAT ARE INSIDE GENES FROM THOSE THAT ARE NOT AND WITH THE REMAINING ONES WE WILL CALCULATE THE DISTANCE TO THE 5'UTR OF THEIR TARGET GENE.


all_a_inter2<-all_a[which(!(all_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene. %in% mer_demg_a$V4)),]

unique(all_a$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#616
unique(mer_demg_a$V4)#292
unique(all_a_inter2$Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.)#324


all_b_inter2<-all_b[which(!(all_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates. %in% mer_demg_b$V4)),]

unique(all_b$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#91
unique(mer_demg_b$V4)#49
unique(all_b_inter2$Enhancer.ID..D..melanogaster.coordinates.database.VT.ID_An..gambiae.coordinates.)#42


### once we have the intergenes we are going to calculate their distance to their target. To do this we will subtract the END of each enhancer from the 5'UTR(START) of the gene.

colnames(all_a_inter2)[6]<-"gene"
all_a_inter3<-merge(genes,all_a_inter2,by="gene")
colnames(all_a_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
all_a_inter3<- all_a_inter3[,c(1:10)]
all_a_inter3<-unique( all_a_inter3)
all_a_inter3$distance_target<-all_a_inter3$Start.gene - ((all_a_inter3$End.enh + all_a_inter3$Start.enh)/2)
all_a_inter3<-all_a_inter3[!duplicated(all_a_inter3),]
sum(all_a_inter3$distance_target<0) #226
sum(all_a_inter3$distance_target>0) #207

colnames(all_b_inter2)[6]<-"gene"
all_b_inter3<-merge(genes,all_b_inter2,by="gene")
colnames(all_b_inter3)<-c("gene", "Chr.enh", "Start.enh","End.enh","Strand","Score", "Chromosome.gene" ,"Start.gene","End.gene", "Enhancer.ID..Previous.author_genomic.coordinates_binding.TFs.annotated.gene.","THS.Peak.ID" )
all_b_inter3<- all_b_inter3[,c(1:10)]
all_b_inter3<-unique( all_b_inter3)
all_b_inter3$distance_target<-all_b_inter3$Start.gene - ((all_b_inter3$End.enh + all_b_inter3$Start.enh)/2)
all_b_inter3<-all_b_inter3[!duplicated(all_b_inter3),]
sum(all_b_inter3$distance_target<0) #210
sum(all_b_inter3$distance_target>0) #208


##### once they are located I will study their activity, first I will see which tissue each peak corresponds to, and then I will see in each enhancer where it is active.


colnames(all_a_1)[4]<-"Enhancer.ID"
all_a_1<-all_a_1[,c(1:5,10)]
colnames(all_b_1)[4]<-"Enhancer.ID"
all_b_1<-all_b_1[,c(1:5,10)]
### juntamos los dos tipos de enhancers y los separamos por tejidos
peaks <- read.delim("/ATAC/results_JL/ATAC-seq_peaks.txt") #193005
all_act<-rbind(all_a_1,all_b_1)
colnames(peaks)[4]<-"THS.Peak.ID"
all_act<-merge(all_act,peaks, by="THS.Peak.ID")
all_act<-all_act[,c(1:6,10,11)]


all_act_d<-all_act[,c(2:6)]
all_act_d<-all_act_d[!duplicated(all_act_d$Enhancer.ID),] #707


write.table((all_act_d), file = "/Enhancers/enhancers_activos_NOmultisoforms_A_B.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
