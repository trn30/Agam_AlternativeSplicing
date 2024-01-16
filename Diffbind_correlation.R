######## Midguts vs. Salivary glands
####RNA
counts_p<- read.delim(" /Isoform_DEI/Isoforms_DEI_switch_DESeq2/Counts_MG_vs_SG_switchisoform.txt", header=T)
counts_p<- counts_p[,c(6:10,12:15)]
colnames(counts_p)<- c("ID","I1d7", "I1d14", "I3d7", "I3d14","I1d7_lenght", "I1d14_lenght", "I3d7_lenght", "I3d14_lenght")
dei_p<- read.delim(" /Isoform_DEI/Isoforms_DEI_switch_DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T)
dei_p$ID<-rownames(dei_p)
dei_counts_p<- merge(dei_p,counts_p,by="ID")
dei_counts_p<-dei_counts_p[,c(1,3,6,8:15)]


#### Calculates rpkm and tpm from rna data

#EGD-2-I1d14 22100970
#EGD-1-I1d7 28594227
#EGD-7-I3d7 18984545
#EGD-8-I3d14 45308711

reads_mg1<-28594227
rpkm_mg1 <- (10^9*dei_counts_p$I1d7/(reads_mg1*(dei_counts_p$I1d7_lenght)))
tpm_mg1<- (((dei_counts_p$I1d7)*rpkm_mg1)/sum(rpkm_mg1))*10^6

reads_sg1<-22100970
rpkm_sg1 <- (10^9*dei_counts_p$I1d14/(reads_sg1*(dei_counts_p$I1d14_lenght)))
tpm_sg1<- (((dei_counts_p$I1d7)*rpkm_sg1)/sum(rpkm_sg1))*10^6

reads_mg3<-18984545
rpkm_mg3 <- (10^9*dei_counts_p$I3d7/(reads_mg3*(dei_counts_p$I3d7_lenght)))
tpm_mg3<- (((dei_counts_p$I1d7)*rpkm_mg3)/sum(rpkm_mg3))*10^6

reads_sg3<-45308711
rpkm_sg3 <- (10^9*dei_counts_p$I3d14/(reads_sg3*(dei_counts_p$I3d14_lenght)))
tpm_sg3<- (((dei_counts_p$I1d7)*rpkm_sg3)/sum(rpkm_sg3))*10^6




counts_rna <- as.data.frame(cbind(dei_counts_p,tpm_mg1,tpm_sg1,tpm_mg3,tpm_sg3))
colnames(counts_rna)
counts_rna$Mean_TPM_sg_rna<- (counts_rna$tpm_sg1+counts_rna$tpm_sg3)/2
counts_rna$Mean_TPM_mg_rna<- (counts_rna$tpm_mg1+counts_rna$tpm_mg3)/2
counts_rna$Mean_TPM_sg_rna_log<-log2(counts_rna$Mean_TPM_sg+0.1)
counts_rna$Mean_TPM_mg_rna_log<- log2(counts_rna$Mean_TPM_mg+0.1)
#####  separate isoforms by their fold change

up_sg_rna<-subset(counts_rna, log2FoldChange > 0) ### up in salivary glands
up_mg_rna<-subset(counts_rna, log2FoldChange < 0) ### up in midguts

####ATAC

#10-I1d7 27749806
#1 I1d14 19630058
#8 I3d7 21546686
#5 I3d14 27476048


diffbind <- read.delim(" /ATAC/DiffBind/Differential_analysis_output_with_ID.txt", header=T) #19630058
diffbind<-diffbind[,c(1:4,7:10,12)]
diffbind$seqnames<-paste0("AgamP4_",diffbind$seqnames)

diffbind_anno<- read.delim(" /ATAC/DiffBind/Differential_analysis_output_with_ID_annotatePeaks.txt", header=T) #19630058
diffbind_anno$Region_ID<-paste(diffbind_anno$Chr,(diffbind_anno$Start-1),diffbind_anno$End,sep = "_")
diffbind_anno$seqnames<-paste0("AgamP4_",diffbind_anno$seqnames)

diffbind_anno_f<-diffbind_anno[diffbind_anno$Distance.to.TSS > -1000, ]
diffbind_anno_f<-diffbind_anno_f[diffbind_anno_f$Distance.to.TSS < 0, ]


diffbind_f<-merge(diffbind,diffbind_anno_f,by="Region_ID")
diffbind_f<-diffbind_f[,c(1:9,20)]
colnames(diffbind_f)<- c( "Region_ID","Chr","start" , "end", "width","Conc_7d" ,"Conc_14d" , "Fold", "p.value", "ID"  )
diffbind_f<-diffbind_f[!duplicated(diffbind_f$Region_ID), ]


reads_sg= 19630058 + 27476048
reads_mg=  27476048 + 21546686


rpkm_sg_atac<- (10^9*diffbind_f$Conc_14d/(reads_sg*(diffbind_f$width)))
tpm_sg_atac <- (((diffbind_f$width)*rpkm_sg_atac)/sum(rpkm_sg_atac))*10^6

rpkm_mg_atac<- (10^9*diffbind_f$Conc_7d/(reads_mg*(diffbind_f$width)))
tpm_mg_atac <- (((diffbind_f$width)*rpkm_mg_atac)/sum(rpkm_mg_atac))*10^6

##### is already the logarithm


counts_sg_f<-cbind(diffbind_f, rpkm_sg_atac,tpm_sg_atac)


counts_mg_f<-cbind(diffbind_f, rpkm_mg_atac,tpm_mg_atac)
##### The fold change of the ATAC-seq is made at the inverse of the RNA-seq fold change.

up_counts_mg <- merge(counts_mg_f,up_mg_rna,by="ID") # merge the data set by ID
up_counts_mg <- up_counts_mg[,c(1:5,9,12,13,30)]
colnames(up_counts_mg)<-c("ID","Region_ID" ,"Chr","start","end","Foldchange_atac","tpm_mg_atac_log" , "log2FoldChange_rna", "Mean_TPM_mg_rna_log")
up_counts_sg <- merge(counts_sg_f,up_sg_rna,by="ID") # merge the data set by ID
up_counts_sg <- up_counts_sg[,c(1:5,9,12,13,29)]
colnames(up_counts_sg)<-c("ID","Region_ID" ,"Chr","start","end","Foldchange_atac","tpm_mg_atac_log" , "log2FoldChange_rna", "Mean_TPM_mg_rna_log")

up_counts_sg_dif<-up_counts_sg[!duplicated(up_counts_sg$Region_ID), ]
up_counts_sg_dif$Region_ID<-paste("AgamP4",up_counts_sg_dif$Region_ID,sep = "_")
up_counts_mg_dif<-up_counts_mg[!duplicated(up_counts_mg$Region_ID), ]
up_counts_mg_dif$Region_ID<-paste("AgamP4",up_counts_mg_dif$Region_ID,sep = "_")


which(up_counts_mg_dif$Region_ID%in% up_counts_mg$Region_ID)
which(up_counts_mg$Region_ID%in% up_counts_mg_dif$Region_ID)

write.table(up_counts_mg, file = " /Correl_RNA_ATAC/Isoforms_up_MG_comparison_SG_vs_MG.txt",quote = FALSE, sep="\t", row.names = FALSE, col.names = T)
write.table(up_counts_sg, file = " /Correl_RNA_ATAC/Isoforms_up_SG_comparison_SG_vs_MG.txt",quote = FALSE, sep="\t", row.names = FALSE, col.names = T)



hist(up_counts_mg$Mean_TPM_mg_rna_log  )
hist(up_counts_mg$mean_tpm_mg_atac_log)

hist(up_counts_sg$Mean_TPM_sg_rna_log)
hist(up_counts_sg$mean_tpm_sg_atac_log) ### clip <21 at least


up_counts_sg_clip<-unique(up_counts_sg[up_counts_sg$mean_tpm_sg_atac_log < 21, ])


hist(up_counts_sg_clip$mean_tpm_sg_atac_log)


### normality test
shapiro.test(up_counts_mg$mean_tpm_mg_atac_log[1:5000]) #W = 0.97795, p-value = 0.04756
shapiro.test(up_counts_sg$mean_tpm_sg_atac_log[1:5000]) #W = 0.74897, p-value < 2.2e-16
shapiro.test(up_counts_mg$Mean_TPM_mg_rna_log[1:5000]) #W = 0.90725, p-value = 5.148e-07
shapiro.test(up_counts_sg$Mean_TPM_sg_rna_log[1:5000]) #W = 0.94152, p-value = 1.6e-10

#### correlation test
cor.test(up_counts_mg$mean_tpm_mg_atac_log,up_counts_mg$Mean_TPM_mg_rna_log,method="spearman",exact = F) #-0.1349262  p-value = 0.1435


cor.test(up_counts_sg$mean_tpm_sg_atac_log,up_counts_sg$Mean_TPM_sg_rna_log,method="spearman",exact = F) #0.2109188  p-value = 6.99e-05


##### heatmaps
test_counts_tpm_1  <- up_counts_mg #### Separar por infecciones. Ordenar por la señal de RNA media.
colnames(test_counts_tpm_1)
# Order per mean of rna counts:
test_counts_tpm_1 <- test_counts_tpm_1[order(test_counts_tpm_1$mean_tpm_mg_atac_log),]
test_counts_tpm_1$order <- 1:dim(test_counts_tpm_1)[1]
head(test_counts_tpm_1)
test_counts_tpm2_1 <- test_counts_tpm_1[,c("ID","mean_tpm_mg_atac_log","Mean_TPM_mg_rna_log")]
rnames_test_counts_tpm_1 <- test_counts_tpm2_1[,1]
mat_data_test_counts_tpm_1 <- data.matrix(test_counts_tpm2_1[,2:ncol(test_counts_tpm2_1)])
rownames(mat_data_test_counts_tpm_1) <- rnames_test_counts_tpm_1
colnames(mat_data_test_counts_tpm_1) <- c("ATAC_tpm","RNA_tpm")
head(test_counts_tpm2_1)
head(mat_data_test_counts_tpm_1)
my_palette <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)
m_test_counts_tpm_1 <-as.matrix(mat_data_test_counts_tpm_1)
dim(m_test_counts_tpm_1) #367


heatmap_test_counts_tpm_1 <- iheatmapr::iheatmap(((m_test_counts_tpm_1[,1:2])), name = "RNA and ATAC TPM mean (log2)",
                                                 colors= my_palette,cluster_rows = "none",
                                                 scale = c("cols"),scale_method="standardize",
                                                 col_labels = T,row_labels = F,
                                                 row_order=test_counts_tpm_1$order)
heatmap_test_counts_tpm_1



save_iheatmap(heatmap_test_counts_tpm_1," /Correl_RNA_ATAC/heatmap_correl_MG_vs_SG_infeccted_1kb_UP_MG.png")
save_iheatmap(heatmap_test_counts_tpm_1," /Correl_RNA_ATAC/heatmap_correl_MG_vs_SG_infeccted_1kb_UP_MG.pdf")
save_iheatmap(heatmap_test_counts_tpm_1," /Correl_RNA_ATAC/heatmap_correl_MG_vs_SG_infeccted_1kb_UP_MG.html")

test_counts_tpm_1  <- up_counts_sg_clip #### Separar por infecciones. Ordenar por la señal de RNA media.
colnames(test_counts_tpm_1)
# Order per mean of rna counts:
test_counts_tpm_1 <- test_counts_tpm_1[order(test_counts_tpm_1$mean_tpm_sg_atac_log),]
test_counts_tpm_1$order <- 1:dim(test_counts_tpm_1)[1]
head(test_counts_tpm_1)
test_counts_tpm2_1 <- test_counts_tpm_1[,c("ID","mean_tpm_sg_atac_log","Mean_TPM_sg_rna_log")]
rnames_test_counts_tpm_1 <- test_counts_tpm2_1[,1]
mat_data_test_counts_tpm_1 <- data.matrix(test_counts_tpm2_1[,2:ncol(test_counts_tpm2_1)])
rownames(mat_data_test_counts_tpm_1) <- rnames_test_counts_tpm_1
colnames(mat_data_test_counts_tpm_1) <- c("ATAC_tpm","RNA_tpm")
head(test_counts_tpm2_1)
head(mat_data_test_counts_tpm_1)
my_palette <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)
m_test_counts_tpm_1 <-as.matrix(mat_data_test_counts_tpm_1)
dim(m_test_counts_tpm_1) #367


heatmap_test_counts_tpm_1 <- iheatmapr::iheatmap(((m_test_counts_tpm_1[,1:2])), name = "RNA and ATAC TPM mean (log2)",
                                                 colors= my_palette,cluster_rows = "none",
                                                 scale = c("cols"),scale_method="standardize",
                                                 col_labels = T,row_labels = F,
                                                 row_order=test_counts_tpm_1$order)
heatmap_test_counts_tpm_1

save_iheatmap(heatmap_test_counts_tpm_1," /Correl_RNA_ATAC/heatmap_correl_MG_vs_SG_infeccted_1kb_UP_SG.png")
save_iheatmap(heatmap_test_counts_tpm_1," /Correl_RNA_ATAC/heatmap_correl_MG_vs_SG_infeccted_1kb_UP_SG.pdf")
save_iheatmap(heatmap_test_counts_tpm_1," /Correl_RNA_ATAC/heatmap_correl_MG_vs_SG_infeccted_1kb_UP_SG.html")
