library(dplyr)
######## Midguts vs. Salivary glands

########PROMOTOR#######
####  I will study the accessibility of ALL DEMGs and then correlate them.

####RNA ####

counts_p<- read.delim("/All_isoforms_filtered_by_cpm_more_than_one_transcript_inf_MG_vs_inf_SG.txt", header=T)
deg<- read.delim("/Isoforms_DEMG_inf_MG_vs_inf_SG_padj_0.05.txt", header=T) #4933 DEG
deg$isoform_ID<-rownames(deg)

counts_p$isoform_ID<-counts_p$ID
deg$Gene_ID<-gsub("-R.*","",deg$isoform_ID)
deg_counts_p<- merge(deg,counts_p,by="isoform_ID")#463

len<- read.delim("/genomic_data_AgamP4/Genes_AgamP4_release_54.bed",header  = F)
colnames(len)<- c("chr","start","End","Gene_ID","strand","score")

deg_f<- merge(len,deg_counts_p,by="Gene_ID")## pierdo genes de t-rna o alguno que no esta anotado
over_deg<- read.delim("/overlapping/DEMG_without_overlapping_genes_ID_definitive.bed",header  = F)
colnames(over_deg)<-"isoform_ID"
deg_f<-deg_f[which((over_deg$isoform_ID %in% deg_f$isoform_ID)),]#337

no_anno<-deg_counts_p[which(!(deg_counts_p$Gene_ID %in% len$Gene_ID)),] #7

### Calculates rpkm and tpm from rna data

rpkm_mg1<-(deg_f$cpm_mg1*1000)/(deg_f$End - deg_f$start)
tpm_mg1<-(rpkm_mg1/sum(rpkm_mg1))*10^6

rpkm_sg1<-deg_f$cpm_sg1*1000/(deg_f$End - deg_f$start)
tpm_sg1<-(rpkm_mg1/sum(rpkm_sg1))*10^6


rpkm_mg3<-deg_f$cpm_mg3*1000/(deg_f$End -deg_f$start)
tpm_mg3<-(rpkm_mg3/sum(rpkm_mg3))*10^6


rpkm_sg3<-deg_f$cpm_sg3*1000/(deg_f$End -deg_f$start)
tpm_sg3<-(rpkm_sg3/sum(rpkm_sg3))*10^6


counts_rna <- as.data.frame(cbind(deg_f,tpm_mg1,tpm_sg1,tpm_mg3,tpm_sg3))
colnames(counts_rna)
counts_rna$Mean_sg_rna<- (counts_rna$tpm_sg1 +counts_rna$tpm_sg3)/2
counts_rna$Mean_mg_rna<- (counts_rna$tpm_mg1+counts_rna$tpm_mg3)/2
counts_rna <- counts_rna[c(1:13,28,29)]
write.table((counts_rna), file = "/Correlation_DEMG/Isoforms_DEMG_MG_vs_SG_0.05_counts_mean.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


#####  separate isoforms by their fold change

up_sg_rna<-subset(counts_rna, log2FoldChange  > 0) ### up in salivary glands 178
up_sg_rna<-up_sg_rna[,c(1:7,13,14)]
up_mg_rna<-subset(counts_rna, log2FoldChange < 0) ### up in midguts 159
up_mg_rna<-up_mg_rna[,c(1:7,13,15)]

##### I extract the accessibility data from the DEMGs isoforms ####

###### ATAC ####
counts_sg1_pr <- read.delim("/ATAC/counts_promoters_all_isoforms/R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed", header=FALSE) #19630058
colnames(counts_sg1_pr)<-c("Chr","Start","End","ID","Strand","Score","Counts")

counts_sg3_pr<- read.delim("/ATAC/counts_promoters_all_isoforms/R5-I3d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed", header=FALSE) #27476048
colnames(counts_sg3_pr)<-c("Chr","Start","End","ID","Strand","Score","Counts")

counts_mg1_pr<- read.delim("/ATAC/counts_promoters_all_isoforms/R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed", header=FALSE) #27476048
colnames(counts_mg1_pr)<-c("Chr","Start","End","ID","Strand","Score","Counts")

counts_mg3_pr <- read.delim("/ATAC/counts_promoters_all_isoforms/R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed", header=FALSE) #21546686
colnames(counts_mg3_pr)<-c("Chr","Start","End","ID","Strand","Score","Counts")

reads_sg1= 19630058
reads_sg3= 27476048
reads_mg1=  27749806
reads_mg3=  21546686

### I will take the co-references of res_sg and res_mg which are the real ones of the peak and not those of the promoter (1Kbp).

rpkm_sg1 <- (10^9*counts_sg1_pr$Counts/(reads_sg1*(counts_sg1_pr$End-counts_sg1_pr$Start)))
tpm_sg1<-(rpkm_sg1/sum(rpkm_sg1))*10^6

rpkm_sg3 <- (10^9*counts_sg3_pr$Counts/(reads_sg3*(counts_sg3_pr$End-counts_sg3_pr$Start)))
tpm_sg3<-(rpkm_sg3/sum(rpkm_sg3))*10^6

rpkm_mg1 <- (10^9*counts_mg1_pr$Counts/(reads_mg1*(counts_mg1_pr$End-counts_mg1_pr$Start)))
tpm_mg1<-(rpkm_mg1/sum(rpkm_mg1))*10^6

rpkm_mg3 <- (10^9*counts_mg3_pr$Counts/(reads_mg3*(counts_mg3_pr$End-counts_mg3_pr$Start)))
tpm_mg3<-(rpkm_mg3/sum(rpkm_mg3))*10^6


counts_sg1_f_pr<-cbind(counts_sg1_pr, rpkm_sg1,tpm_sg1)

counts_sg3_f_pr<-cbind(counts_sg3_pr, rpkm_sg3,tpm_sg3)

counts_sg_final_pr<-merge(counts_sg1_f_pr,counts_sg3_f_pr, by="ID")
counts_sg_final_pr$mean_tpm_sg_atac<-(counts_sg_final_pr$rpkm_sg1+ counts_sg_final_pr$rpkm_sg3)/2
counts_sg_final_pr<-counts_sg_final_pr[,c(1:6,18)]
colnames(counts_sg_final_pr)<-c("isoform_ID", "Chr", "Start","End","Strand","Score","mean_tpm_atac_sg")
counts_sg_final_pr<-counts_sg_final_pr[!duplicated(counts_sg_final_pr), ]
counts_sg_final_pr$Gene_ID<-gsub("-R.*","",counts_sg_final_pr$isoform_ID)
aa_sg <- counts_sg_final_pr %>%
  dplyr::group_by(Gene_ID) %>% dplyr::summarise(Frequency = sum(mean_tpm_atac_sg),NN=dplyr::n())
bb_sg<- rep(aa_sg$Frequency,aa_sg$NN)
counts_sg_final_pr$sum_atac_gene_sg<-bb_sg
counts_sg_final_pr$prop_sg_atac<-counts_sg_final_pr$mean_tpm_atac_sg / counts_sg_final_pr$sum_atac_gene_sg


counts_mg1_f_pr<-cbind(counts_mg1_pr, rpkm_mg1,tpm_mg1)

counts_mg3_f_pr<-cbind(counts_mg3_pr, rpkm_mg3,tpm_mg3)

counts_mg_final_pr<-merge(counts_mg1_f_pr,counts_mg3_f_pr, by="ID")
counts_mg_final_pr$mean_tpm_mg_atac<-(counts_mg_final_pr$rpkm_mg1+ counts_mg_final_pr$rpkm_mg3)/2
counts_mg_final_pr<-counts_mg_final_pr[,c(1:6,18)]
colnames(counts_mg_final_pr)<-c("isoform_ID", "Chr", "Start","End","Strand","Score","mean_tpm_atac_mg")
counts_mg_final_pr<-counts_mg_final_pr[!duplicated(counts_mg_final_pr), ]
counts_mg_final_pr$Gene_ID<-gsub("-R.*","",counts_mg_final_pr$isoform_ID)
aa_mg <- counts_mg_final_pr %>%
  dplyr::group_by(Gene_ID) %>% dplyr::summarise(Frequency = sum(mean_tpm_atac_mg),NN=dplyr::n())
bb_mg<- rep(aa_mg$Frequency,aa_mg$NN)
counts_mg_final_pr$sum_atac_gene_mg<-bb_mg
counts_mg_final_pr$prop_mg_atac<-counts_mg_final_pr$mean_tpm_atac_mg /  counts_mg_final_pr$sum_atac_gene_mg


#### now I will put together the DEMG's with their accessibility data ####

deg_d_sg_pr <- merge(up_sg_rna,counts_sg_final_pr,by="isoform_ID")##178
deg_d_sg_pr<-deg_d_sg_pr[!duplicated(deg_d_sg_pr$isoform_ID), ]

deg_d_mg_pr <- merge(up_mg_rna,counts_mg_final_pr,by="isoform_ID")##158
deg_d_mg_pr<-deg_d_mg_pr[!duplicated(deg_d_mg_pr$isoform_ID), ]


####I study the correlation between these data
cor.test(deg_d_sg_pr$Mean_sg_rna,deg_d_sg_pr$mean_tpm_atac_sg,method="spearman",exact = F)
#0.1205485  p-value = 0.109

cor.test(deg_d_mg_pr$Mean_mg_rna,deg_d_mg_pr$mean_tpm_atac_mg,method="spearman",exact = F)
#0.001201789 ,  p-value  =  0.988

write.table((deg_d_mg_pr), file = "/Correlation_DEMG/Isoforms_up_MG_DEMG_promoter_without_overlap.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
write.table((deg_d_sg_pr), file = "/Correlation_DEMG/Isoforms_up_SG_DEMG_promoter_without_overlap.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


####### Gene Body #######################


###### ATAC ####
counts_sg1_bd <- read.delim("/ATAC/counts_promoters_all_isoforms/R1-I1d14_SG_counts_body_SG_vs_MG_2.bed", header=FALSE) #19630058
colnames(counts_sg1_bd)<-c("Chr","Start","End","ID","Strand","Score","Counts")

counts_sg3_bd<- read.delim("/ATAC/counts_promoters_all_isoforms/R5-I3d14_SG_counts_body_SG_vs_MG_2.bed", header=FALSE) #27476048
colnames(counts_sg3_bd)<-c("Chr","Start","End","ID","Strand","Score","Counts")

counts_mg1_bd<- read.delim("/ATAC/counts_promoters_all_isoforms/R10-I1d7_MG_counts_body_SG_vs_MG_2.bed", header=FALSE) #27476048
colnames(counts_mg1_bd)<-c("Chr","Start","End","ID","Strand","Score","Counts")

counts_mg3_bd <- read.delim("/ATAC/counts_promoters_all_isoforms/R8-I3d7_MG_counts_body_SG_vs_MG_2.bed", header=FALSE) #21546686
colnames(counts_mg3_bd)<-c("Chr","Start","End","ID","Strand","Score","Counts")

reads_sg1= 19630058
reads_sg3= 27476048
reads_mg1=  27749806
reads_mg3=  21546686
### voy a coger las cooredenadas de res_sg y res_mg que son las verdaderas del pico y no las del promotor (1Kbp)

rpkm_sg1 <- (10^9*counts_sg1_bd$Counts/(reads_sg1*(counts_sg1_bd$End-counts_sg1_bd$Start)))
tpm_sg1<-(rpkm_sg1/sum(rpkm_sg1))*10^6

rpkm_sg3 <- (10^9*counts_sg3_bd$Counts/(reads_sg3*(counts_sg3_bd$End-counts_sg3_bd$Start)))
tpm_sg3<-(rpkm_sg3/sum(rpkm_sg3))*10^6

rpkm_mg1 <- (10^9*counts_mg1_bd$Counts/(reads_mg1*(counts_mg1_bd$End-counts_mg1_bd$Start)))
tpm_mg1<-(rpkm_mg1/sum(rpkm_mg1))*10^6

rpkm_mg3 <- (10^9*counts_mg3_bd$Counts/(reads_mg3*(counts_mg3_bd$End-counts_mg3_bd$Start)))
tpm_mg3<-(rpkm_mg3/sum(rpkm_mg3))*10^6


counts_sg1_f_bd<-cbind(counts_sg1_bd, rpkm_sg1,tpm_sg1)

counts_sg3_f_bd<-cbind(counts_sg3_bd, rpkm_sg3,tpm_sg3)

counts_sg_final_bd<-merge(counts_sg1_f_bd,counts_sg3_f_bd, by="ID")
counts_sg_final_bd$mean_tpm_sg_atac<-(counts_sg_final_bd$rpkm_sg1+ counts_sg_final_bd$rpkm_sg3)/2
counts_sg_final_bd<-counts_sg_final_bd[,c(1:6,18)]
colnames(counts_sg_final_bd)<-c("isoform_ID", "Chr", "Start","End","Strand","Score","mean_tpm_atac_sg")
counts_sg_final_bd<-counts_sg_final_bd[!duplicated(counts_sg_final_bd), ]
counts_sg_final_bd$Gene_ID<-gsub("-R.*","",counts_sg_final_bd$isoform_ID)
aa_sg <- counts_sg_final_bd %>%
  dplyr::group_by(Gene_ID) %>% dplyr::summarise(Frequency = sum(mean_tpm_atac_sg),NN=dplyr::n())
bb_sg<- rep(aa_sg$Frequency,aa_sg$NN)
counts_sg_final_bd$sum_atac_gene_sg<-bb_sg
counts_sg_final_bd$prop_sg_atac<-counts_sg_final_bd$mean_tpm_atac_sg / counts_sg_final_bd$sum_atac_gene_sg


counts_mg1_f_bd<-cbind(counts_mg1_bd, rpkm_mg1,tpm_mg1)

counts_mg3_f_bd<-cbind(counts_mg3_bd, rpkm_mg3,tpm_mg3)

counts_mg_final_bd<-merge(counts_mg1_f_bd,counts_mg3_f_bd, by="ID")
counts_mg_final_bd$mean_tpm_mg_atac<-(counts_mg_final_bd$rpkm_mg1+ counts_mg_final_bd$rpkm_mg3)/2
counts_mg_final_bd<-counts_mg_final_bd[,c(1:6,18)]
colnames(counts_mg_final_bd)<-c("isoform_ID", "Chr", "Start","End","Strand","Score","mean_tpm_atac_mg")
counts_mg_final_bd<-counts_mg_final_bd[!duplicated(counts_mg_final_bd), ]
counts_mg_final_bd$Gene_ID<-gsub("-R.*","",counts_mg_final_bd$isoform_ID)
aa_mg <- counts_mg_final_bd %>%
  dplyr::group_by(Gene_ID) %>% dplyr::summarise(Frequency = sum(mean_tpm_atac_mg),NN=dplyr::n())
bb_mg<- rep(aa_mg$Frequency,aa_mg$NN)
counts_mg_final_bd$sum_atac_gene_mg<-bb_mg
counts_mg_final_bd$prop_mg_atac<-counts_mg_final_bd$mean_tpm_atac_mg / counts_mg_final_bd$sum_atac_gene_mg


#### I will now put together the DEMG's with their accessibility data. ####

deg_d_sg_bd <- merge(up_sg_rna,counts_sg_final_bd,by="isoform_ID")
deg_d_sg_bd<-deg_d_sg_bd[!duplicated(deg_d_sg_bd$isoform_ID), ] ## 178

deg_d_mg_bd <- merge(up_mg_rna,counts_mg_final_bd,by="isoform_ID")
deg_d_mg_bd<-deg_d_mg_bd[!duplicated(deg_d_mg_bd$isoform_ID), ] #158


####I study the correlation between these data
cor.test(deg_d_sg_bd$Mean_sg_rna,deg_d_sg_bd$mean_tpm_atac_sg,method="spearman",exact = F)
#0.3029904 ,  p-value  3.94e-05

cor.test(deg_d_mg_bd$Mean_mg_rna,deg_d_mg_bd$mean_tpm_atac_mg,method="spearman",exact = F)
#0.3378486,  p-value   1.416e-05

write.table((deg_d_mg_bd), file = "/Correlation_DEMG/Isoforms_up_MG_DEMG_body_without_overlap.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
write.table((deg_d_sg_bd), file = "/Correlation_DEMG/Isoforms_up_SG_DEMG_body_without_overlap.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
#

 ##################################################


#### SCATTER PLOT

layout(matrix(c(1:4), nrow=2, byrow=FALSE))
layout.show(4)

# gene body
deg_d_mg_bd_c<-deg_d_mg_bd[(deg_d_mg_bd$Mean_mg_rna < 20000),]
p3<-ggplot(deg_d_mg_bd_c,aes(x=Mean_mg_rna,y=mean_tpm_atac_mg))+geom_point()+geom_smooth(method = lm)

deg_d_sg_bd_c<-deg_d_sg_bd[(deg_d_sg_bd$Mean_sg_rna < 6000),]
p4<-ggplot(deg_d_sg_bd_c,aes(x=Mean_sg_rna,y=mean_tpm_atac_sg))+geom_point()+geom_smooth(method = lm)

#promoter

deg_d_mg_pr_c<-deg_d_mg_pr[(deg_d_mg_pr$Mean_mg_rna < 20000),]
p5<-ggplot(deg_d_mg_pr_c,aes(x=Mean_mg_rna,y=mean_tpm_atac_mg))+geom_point(colour="red")+geom_smooth(method = lm)

deg_d_sg_pr_c<-deg_d_sg_pr[(deg_d_sg_pr$Mean_sg_rna < 6000),]
p6<-ggplot(deg_d_sg_pr_c,aes(x=Mean_sg_rna,y=mean_tpm_atac_sg))+geom_point(colour="red")+geom_smooth(method = lm)

library(patchwork)
(p3|p4)/(p5|p6)

#### heatmaps
###SG
all<-cbind(deg_d_sg_bd,deg_d_sg_pr)
all<-all[,c(1,9:15,33)]
colnames(all)<-c("isoform_ID","Mean_sg_rna","chr","start","End","strand","score","mean_sg_atac_bd","mean_sg_atac_pr")


hist(all$mean_sg_atac_bd)
all_clip<-all[(all$Mean_sg_rna < 75000),]
all_clip<-all_clip[(all_clip$mean_sg_atac_pr > 2),]

test_counts_tpm_bd  <- all_clip #### Separate by infection. Sort by mean RNA signal.
colnames(test_counts_tpm_bd)
# Order per mean of rna counts:
test_counts_tpm_bd <- test_counts_tpm_bd[order(test_counts_tpm_bd$Mean_sg_rna),]
test_counts_tpm_bd$order <- 1:dim(test_counts_tpm_bd)[1]
head(test_counts_tpm_bd)
test_counts_tpm2_1 <- test_counts_tpm_bd[,c("isoform_ID","Mean_sg_rna","mean_sg_atac_bd","mean_sg_atac_pr")]
rnames_test_counts_tpm_bd <- test_counts_tpm2_1[,1]
mat_data_test_counts_tpm_bd <- data.matrix(test_counts_tpm2_1[,2:ncol(test_counts_tpm2_1)])
colnames(mat_data_test_counts_tpm_bd) <- c("Mean_sg_rna","mean_sg_atac_bd","mean_sg_atac_pr")
head(test_counts_tpm2_1)
head(mat_data_test_counts_tpm_bd)
my_palette <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)
m_test_counts_tpm_bd <-as.matrix(mat_data_test_counts_tpm_bd)
dim(m_test_counts_tpm_bd) #174


heatmap_test_counts_tpm_bd <- iheatmapr::iheatmap((log2(m_test_counts_tpm_bd[,1:3]+1)), name = "RNA and ATAC TPM mean (log2)",
                                                  colors= my_palette,cluster_rows = "none",
                                                  scale = c("cols"),scale_method="standardize",
                                                  col_labels = T,row_labels = F,
                                                  row_order=test_counts_tpm_bd$order)
heatmap_test_counts_tpm_bd

pdf(file = "/Correlation_DEMG/Correlation_DEMG_Heatmap_SG.pdf")
heatmap_test_counts_tpm_bd
dev.off()


###MG

all<-cbind(deg_d_mg_bd,deg_d_mg_pr)
all<-all[,c(1,9:15,33)]
colnames(all)<-c("isoform_ID","Mean_mg_rna","chr","start","End","strand","score","mean_mg_atac_bd","mean_mg_atac_pr")


hist(all$mean_mg_atac_bd)
hist(all$mean_mg_atac_pr)
hist(all$Mean_mg_rna)

all_clip<-all[(all$mean_mg_atac_bd > 3),]
all_clip<-all_clip[(all_clip$mean_mg_atac_pr > 2),]
all_clip<-all_clip[(all_clip$Mean_mg_rna < 10000),]


test_counts_tpm_bd  <- all_clip ### Separate by infection. Sort by mean RNA signal.
colnames(test_counts_tpm_bd)
# Order per mean of rna counts:
test_counts_tpm_bd <- test_counts_tpm_bd[order(test_counts_tpm_bd$Mean_mg_rna),]
test_counts_tpm_bd$order <- 1:dim(test_counts_tpm_bd)[1]
head(test_counts_tpm_bd)
test_counts_tpm2_1 <- test_counts_tpm_bd[,c("isoform_ID","Mean_mg_rna","mean_mg_atac_bd","mean_mg_atac_pr")]
rnames_test_counts_tpm_bd <- test_counts_tpm2_1[,1]
mat_data_test_counts_tpm_bd <- data.matrix(test_counts_tpm2_1[,2:ncol(test_counts_tpm2_1)])
colnames(mat_data_test_counts_tpm_bd) <- c("Mean_mg_rna","mean_mg_atac_bd","mean_mg_atac_pr")
head(test_counts_tpm2_1)
head(mat_data_test_counts_tpm_bd)
my_palette <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)
m_test_counts_tpm_bd <-as.matrix(mat_data_test_counts_tpm_bd)
dim(m_test_counts_tpm_bd) #134

heatmap_test_counts_tpm_bd2 <- iheatmapr::iheatmap((log2(m_test_counts_tpm_bd[,1:3]+1)), name = "RNA and ATAC TPM mean (log2)",
                                                   colors= my_palette,cluster_rows = "none",
                                                   scale = c("cols"),scale_method="standardize",
                                                   col_labels = T,row_labels = F,
                                                   row_order=test_counts_tpm_bd$order)
heatmap_test_counts_tpm_bd2
pdf(file = "/Correlation_DEMG/Correlation_DEMG_Heatmap_MG_Rerun.pdf")
heatmap_test_counts_tpm_bd2
dev.off()
