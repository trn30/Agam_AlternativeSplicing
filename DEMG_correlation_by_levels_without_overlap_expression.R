
#################################### Correlation testBY LEVELS ######################################
### Expression levels
### Midguts
up_counts_mg<- read.delim("/Correlation_DEMG/Isoforms_up_MG_DEMG_promoter_without_overlap.txt", header=T)
#### levels of expression

cats = split(up_counts_mg, Hmisc::cut2(up_counts_mg$Mean_mg_rna, g=3))
names(cats)
up_mg_levels<-up_counts_mg
#"[   7.87,   293)" "[ 293.28,  1038)" "[1037.77,476207]"


up_mg_levels$rna_categ <- "High"
up_mg_levels$rna_categ[up_counts_mg$Mean_mg_rna < 1038] <- "Medium"
up_mg_levels$rna_categ[up_counts_mg$Mean_mg_rna < 293] <- "Low"
table(up_mg_levels$rna_categ)

#   High    Low Medium
#    51     53     54

high_mg<-subset(up_mg_levels,rna_categ=="High")

medium_mg<- subset(up_mg_levels,rna_categ=="Medium")

low_mg<-subset(up_mg_levels,rna_categ=="Low")

###### I extract the coordinates of the tss to represent in ngsplot

coor<-read.delim(file="/Genomes/genomic_data_AgamP4/Coordinates_only_tss_corrected_all_isoforms.bed",header = F)
colnames(coor)<-c("Chr","Start","End","isoform_ID","Strand","Score")

clas<-merge(up_mg_levels, coor, by="isoform_ID")
clas<-clas[,c(1,9,15,18,19:24)]

high_mg<-subset(clas,rna_categ=="High")
high_mg<-high_mg[,c(6:8,1,9,10,2,3,5)]
medium_mg<- subset(clas,rna_categ=="Medium")
medium_mg<-medium_mg[,c(6:8,1,9,10,2,3,5)]
low_mg<-subset(clas,rna_categ=="Low")
low_mg<-low_mg[,c(6:8,1,9,10,2,3,5)]
colnames(low_mg)

write.table(high_mg, file = "/Correlation_DEMG/Genes_high_expressed_MG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(medium_mg, file = "/Correlation_DEMG/Genes_medium_expressed_MG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(low_mg, file = "/Correlation_DEMG/Genes_low_expressed_MG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

### violin plot ####

p <- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_mg,fill=rna_categ)) +  geom_violin(trim = T) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
pdf(file = "/Correlation_DEMG/Violinplot_levels_expression_MG.pdf")
p
dev.off()
d<- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_mg,fill=rna_categ)) +  geom_boxplot()
pdf(file = "/Correlation_DEMG/Boxplot_levels_expression_MG.pdf")
d
dev.off()


### Correlation test

#### High expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(high_mg$Mean_mg_rna[1:5000]) #

shapiro.test(high_mg$mean_tpm_atac_mg[1:5000]) #

### as they do not follow a normal distribution, we perform the Spearman's correlation

cor.test(high_mg$Mean_mg_rna,high_mg$mean_tpm_atac_mg,method="spearman",exact = F)
#0.02444178  , p-value = 0.8662


#### Medium expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(medium_mg$Mean_mg_rna[1:5000]) #

shapiro.test(medium_mg$mean_tpm_atac_mg[1:5000]) #

### as they do not follow a normal distribution, we perform the Spearman's correlation

cor.test(medium_mg$Mean_mg_rna,medium_mg$mean_tpm_atac_mg,method="spearman",exact = F)
#-0.04939028 p-value = 0.7361


#### Low expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(low_mg$Mean_mg_rna[1:5000]) #

shapiro.test(low_mg$mean_tpm_atac_mg[1:5000]) #


### as they do not follow a normal distribution, we perform the Spearman's correlation

cor.test(low_mg$Mean_mg_rna,low_mg$mean_tpm_atac_mg,method="spearman",exact = F)
# 0.02540946  p-value = 0.8624

######################################################################################################################
#######  SG

up_counts_sg<- read.delim("/Correlation_DEMG/Isoforms_up_SG_DEMG_promoter_without_overlap.txt", header=T)
#### levels of expression

cats = split(up_counts_sg, Hmisc::cut2(up_counts_sg$Mean_sg_rna, g=3))
names(cats)
up_sg_levels<-up_counts_sg
#[1] "[  5.9,  131)" "[131.2,  492)" "[492.2,80428]"


up_sg_levels$rna_categ <- "High"
up_sg_levels$rna_categ[up_counts_sg$Mean_sg_rna < 492] <- "Medium"
up_sg_levels$rna_categ[up_counts_sg$Mean_sg_rna < 131] <- "Low"
table(up_sg_levels$rna_categ)

#   High    Low Medium
#    59     58     61

high_sg<-subset(up_sg_levels,rna_categ=="High")

medium_sg<- subset(up_sg_levels,rna_categ=="Medium")

low_sg<-subset(up_sg_levels,rna_categ=="Low")

###### I extract the coordinates of the TSS to represent in ngsplot

coor<-read.delim(file="/Genomes/genomic_data_AgamP4/Coordinates_only_tss_corrected_all_isoforms.bed",header = F)
colnames(coor)<-c("Chr","Start","End","isoform_ID","Strand","Score")

clas<-merge(up_sg_levels, coor, by="isoform_ID")
clas<-clas[,c(1,9,15,18:24)]

high_sg<-subset(clas,rna_categ=="High")
high_sg<-high_sg[,c(6:8,1,9,10,2,3,5)]
medium_sg<- subset(clas,rna_categ=="Medium")
medium_sg<-medium_sg[,c(6:8,1,9,10,2,3,5)]
low_sg<-subset(clas,rna_categ=="Low")
low_sg<-low_sg[,c(6:8,1,9,10,2,3,5)]
colnames(low_sg)

write.table(high_sg, file = "/Correlation_DEMG/Genes_high_expressed_SG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(medium_sg, file = "/Correlation_DEMG/Genes_medium_expressed_SG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(low_sg, file = "/Correlation_DEMG/Genes_low_expressed_SG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)


###violin plot
p <- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_sg,fill=rna_categ)) +  geom_violin(trim = T) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")

pdf(file = "/Correlation_DEMG/Violinplot_levels_expression_SG.pdf")
p
dev.off()

d<- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_sg,fill=rna_categ)) +  geom_boxplot()
pdf(file = "/Correlation_DEMG/Boxplot_levels_expression_SG.pdf")
d
dev.off()

### Correlation test

#### High expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(high_sg$Mean_sg_rna[1:5000]) #

shapiro.test(high_sg$mean_tpm_atac_sg[1:5000]) #


### as they do not follow a normal distribution, we perform the Spearman's correlation

cor.test(high_sg$Mean_sg_rna,high_sg$mean_tpm_atac_sg,method="spearman",exact = F)
#0.07645099 , p-value = 0.5649


#### Medium expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(medium_sg$Mean_sg_rna[1:5000]) #

shapiro.test(medium_sg$mean_tpm_atac_sg[1:5000]) #

### as they do not follow a normal distribution, we perform the Spearman's correlation

cor.test(medium_sg$Mean_sg_rna,medium_sg$mean_tpm_atac_sg,method="spearman",exact = F)
#0.1363683  p-value = 0.3163


#### Low expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(low_sg$Mean_sg_rna) #

shapiro.test(low_sg$mean_tpm_atac_sg[1:5000]) #

### as they do not follow a normal distribution, we perform the Spearman's correlation

cor.test(low_sg$Mean_sg_rna,low_sg$mean_tpm_atac_sg,method="spearman",exact = F)
# 0.0406546  p-value = 0.7838


############# NGSPLOT
### https://github.com/shenlab-sinai/ngsplot
### a standalone program to visualize enrichment patterns of DNA-interacting proteins at functionally important regions based on next-generation sequencing data.
### Linux terminal
ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_expresed_MG_all.txt -O expression_levels_in_TSS_MG
ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_expresed_SG_all.txt -O expression_levels_in_TSS_SG
