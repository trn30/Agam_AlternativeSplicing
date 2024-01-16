
#################################### CORRELATION BY LEVELS  ######################################
### levels of expression
### MG
up_counts_mg<- read.delim(" /Correl_DEMG/Isoforms_up_MG_DEMG_promoter_without_overlap.txt", header=T)


cats = split(up_counts_mg, Hmisc::cut2(up_counts_mg$Mean_mg_rna, g=3))
names(cats)
up_mg_levels<-up_counts_mg
#"[  7.75,   248)" "[248.38,   934)" "[934.05,469245]"



up_mg_levels$rna_categ <- "High"
up_mg_levels$rna_categ[up_counts_mg$Mean_mg_rna < 934.05] <- "Medium"
up_mg_levels$rna_categ[up_counts_mg$Mean_mg_rna < 248.38] <- "Low"
table(up_mg_levels$rna_categ)

#   High    Low Medium
#   60       61   61

alto_mg<-subset(up_mg_levels,rna_categ=="High")

medio_mg<- subset(up_mg_levels,rna_categ=="Medium")

bajo_mg<-subset(up_mg_levels,rna_categ=="Low")

###### I extract the coordinates of the tss to represent in ngsplot

coor<-read.delim(file=" /genomic_data_AgamP4/Coordinates_only_tss_corrected_all_isoforms.bed",header = F)
colnames(coor)<-c("Chr","Start","End","isoform_ID","Strand","Score")

clas<-merge(up_mg_levels, coor, by="isoform_ID")
clas<-clas[,c(1,10,16,19,20:25)]

alto_mg<-subset(clas,rna_categ=="High")
alto_mg<-alto_mg[,c(6:8,1,9,10,2,3,5)]
medio_mg<- subset(clas,rna_categ=="Medium")
medio_mg<-medio_mg[,c(6:8,1,9,10,2,3,5)]
bajo_mg<-subset(clas,rna_categ=="Low")
bajo_mg<-bajo_mg[,c(6:8,1,9,10,2,3,5)]
colnames(bajo_mg)

write.table(alto_mg, file = " /Correl_DEMG/Genes_high_expressed_MG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(medio_mg, file = " /Correl_DEMG/Genes_medium_expressed_MG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(bajo_mg, file = " /Correl_DEMG/Genes_low_expressed_MG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

### violin plot ####

p <- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_mg,fill=rna_categ)) +  geom_violin(trim = T) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
x11()
p

d<- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_mg,fill=rna_categ)) +  geom_boxplot()
d
### Correlation

####  High expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(alto_mg$Mean_mg_rna[1:5000]) #

shapiro.test(alto_mg$mean_tpm_atac_mg[1:5000]) #

hist(alto_mg$Mean_mg_rna)
hist(alto_mg$mean_tpm_atac_mg)


###  as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(alto_mg$Mean_mg_rna,alto_mg$mean_tpm_atac_mg,method="spearman",exact = F)
#0.008562495 , p-value = 0.9487


#### Half an expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(medio_mg$Mean_mg_rna[1:5000]) #

shapiro.test(medio_mg$mean_tpm_atac_mg[1:5000]) #

hist(medio_mg$Mean_mg_rna)
hist(medio_mg$mean_tpm_atac_mg)


### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(medio_mg$Mean_mg_rna,medio_mg$mean_tpm_atac_mg,method="spearman",exact = F)
#-0.009296285 p-value = 0.9458


#### Low expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(bajo_mg$Mean_mg_rna[1:5000]) #

shapiro.test(bajo_mg$mean_tpm_atac_mg[1:5000]) #

hist(bajo_mg$Mean_mg_rna)
hist(bajo_mg$mean_tpm_atac_mg)

###  as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(bajo_mg$Mean_mg_rna,bajo_mg$mean_tpm_atac_mg,method="spearman",exact = F)
# 0.08786258  p-value = 0.5158

######################################################################################################################
#######  SG
#### levels of expression
up_counts_sg<- read.delim(" /Correl_DEMG/Isoforms_up_SG_DEMG_promoter_without_overlap.txt", header=T)


cats = split(up_counts_sg, Hmisc::cut2(up_counts_sg$Mean_sg_rna, g=3))
names(cats)
up_sg_levels<-up_counts_sg
#[1] "[  5.78,  115)" "[114.65,  429)" "[428.61,78783]"


up_sg_levels$rna_categ <- "High"
up_sg_levels$rna_categ[up_counts_sg$Mean_sg_rna < 428.61] <- "Medium"
up_sg_levels$rna_categ[up_counts_sg$Mean_sg_rna < 114.65] <- "Low"
table(up_sg_levels$rna_categ)

#   High    Low Medium
#    71     72     73

alto_sg<-subset(up_sg_levels,rna_categ=="High")

medio_sg<- subset(up_sg_levels,rna_categ=="Medium")

bajo_sg<-subset(up_sg_levels,rna_categ=="Low")

###### I extract the coordinates of the tss to represent in ngsplot

coor<-read.delim(file=" /genomic_data_AgamP4/Coordinates_only_tss_corrected_all_isoforms.bed",header = F)
colnames(coor)<-c("Chr","Start","End","isoform_ID","Strand","Score")

clas<-merge(up_sg_levels, coor, by="isoform_ID")
clas<-clas[,c(1,10,16,20:25)]

alto_sg<-subset(clas,rna_categ=="High")
alto_sg<-alto_sg[,c(5:7,1,8,9,2:4)]
medio_sg<- subset(clas,rna_categ=="Medium")
medio_sg<-medio_sg[,c(5:7,1,8,9,2:4)]
bajo_sg<-subset(clas,rna_categ=="Low")
bajo_sg<-bajo_sg[,c(5:7,1,8,9,2:4)]
colnames(bajo_sg)

write.table(alto_sg, file = " /Correl_DEMG/Genes_high_expressed_SG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(medio_sg, file = " /Correl_DEMG/Genes_medium_expressed_SG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(bajo_sg, file = " /Correl_DEMG/Genes_low_expressed_SG_DEMG_prom_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

### violin plot
p <- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_sg,fill=rna_categ)) +  geom_violin(trim = T) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
x11()
p
d<- ggplot(clas, aes(x=rna_categ, y=mean_tpm_atac_sg,fill=rna_categ)) +  geom_boxplot()
d
### Correlation

#### High expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(alto_sg$Mean_sg_rna[1:5000]) #

shapiro.test(alto_sg$mean_tpm_atac_sg[1:5000]) #

hist(alto_sg$Mean_sg_rna)
hist(alto_sg$mean_tpm_atac_sg)


### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(alto_sg$Mean_sg_rna,alto_sg$mean_tpm_atac_sg,method="spearman",exact = F)
#0.02418092 , p-value = 0.8425


#### Half an expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(medio_sg$Mean_sg_rna[1:5000]) #

shapiro.test(medio_sg$mean_tpm_atac_sg[1:5000]) #

hist(medio_sg$Mean_sg_rna)
hist(medio_sg$mean_tpm_atac_sg)


### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(medio_sg$Mean_sg_rna,medio_sg$mean_tpm_atac_sg,method="spearman",exact = F)
#-0.08716813  p-value = 0.4797


#### Low expression
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(bajo_sg$Mean_sg_rna[1:5000]) #

shapiro.test(bajo_sg$mean_tpm_atac_sg[1:5000]) #

hist(bajo_sg$Mean_sg_rna)
hist(bajo_sg$mean_tpm_atac_sg)

### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(bajo_sg$Mean_sg_rna,bajo_sg$mean_tpm_atac_sg,method="spearman",exact = F)
# -0.1925588  p-value = 0.1405

###########################################################################################

##########################################################################################
########## levels of accessibility
##### PROMOTER IN MG
up_counts_mg<- read.delim(" /Correl_DEMG/Isoforms_up_MG_DEMG_promoter_without_overlap.txt", header=T)


cats = split(up_counts_mg, Hmisc::cut2(up_counts_mg$mean_tpm_atac_mg , g=3))
names(cats)
up_mg_levels<-up_counts_mg
#[1] "[ 0.036, 6.29)" "[ 6.289,10.17)" "[10.166,24.43]"



up_mg_levels$atac_categ <- "High"
up_mg_levels$atac_categ[up_counts_mg$mean_tpm_atac_mg<10.166] <- "Medium"
up_mg_levels$atac_categ[up_counts_mg$mean_tpm_atac_mg<6.289] <- "Low"
table(up_mg_levels$atac_categ)

#  High  Low  Medium
#   60     62     60


alto_mg<-subset(up_mg_levels,atac_categ=="High")
medio_mg<- subset(up_mg_levels,atac_categ=="Medium")
bajo_mg<-subset(up_mg_levels,atac_categ=="Low")


###### MG
coor<-read.delim(file=" /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed",header = F)
colnames(coor)<-c("Chr","Start","End","isoform_ID","Strand","Score")

clas<-merge(up_mg_levels, coor, by="isoform_ID")
clas<-clas[,c(1,10,16,20:25)]

alto_mg<-subset(clas,atac_categ=="High")
alto_mg<-alto_mg[,c(5:7,1,8,9,2:4)]
medio_mg<- subset(clas,atac_categ=="Medium")
medio_mg<-medio_mg[,c(5:7,1,8,9,2:4)]
bajo_mg<-subset(clas,atac_categ=="Low")
bajo_mg<-bajo_mg[,c(5:7,1,8,9,2:4)]
colnames(bajo_mg)

write.table(alto_mg, file = " /Correl_DEMG/Genes_high_accessible_MG_DEMG_PROM_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(medio_mg, file = " /Correl_DEMG/Genes_medium_accessible_MG_DEMG_PROM_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(bajo_mg, file = " /Correl_DEMG/Genes_low_accessible_MG_DEMG_PROM_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)


### Correlation

#### High accessibility
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(alto_mg$Mean_mg_rna[1:5000]) #

shapiro.test(alto_mg$mean_tpm_atac_mg[1:3000]) #

hist(alto_mg$Mean_mg_rna)
hist(alto_mg$mean_tpm_atac_mg)

### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(alto_mg$mean_tpm_atac_mg,alto_mg$Mean_mg_rna,method="spearman",exact = F)
#-0.1703972 , p-value = 0.193


#### Medium accessibility
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(medio_mg$Mean_mg_rna[1:5000]) #

shapiro.test(medio_mg$mean_tpm_atac_mg[1:5000]) #


hist(medio_mg$Mean_mg_rna)
hist(medio_mg$mean_tpm_atac_mg)

### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(medio_mg$mean_tpm_atac_mg,medio_mg$Mean_mg_rna,method="spearman",exact = F)
#0.09614316   p-value = 0.4649


#### Low  accessibility
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(bajo_mg$Mean_mg_rna[1:5000]) #


shapiro.test(bajo_mg$mean_tpm_atac_mg[1:5000]) #


hist(bajo_mg2$Mean_mg_rna)
hist(bajo_mg$mean_tpm_atac_mg)
bajo_mg2<-bajo_mg[bajo_mg$Mean_mg_rna < 15000,]

### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(bajo_mg$mean_tpm_atac_mg,bajo_mg$Mean_mg_rna,method="spearman",exact = F)
# 0.3808265  p-value = 0.002261

# violin plot
p <- ggplot(clas, aes(x=atac_categ, y=Mean_mg_rna,fill=atac_categ)) +  geom_violin(trim = T) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
x11()
p
d<- ggplot(clas, aes(x=atac_categ, y=Mean_mg_rna,fill=atac_categ)) +  geom_boxplot()
d

######################################################################################################################
#######  PROMOTER IN SG

up_counts_sg<- read.delim(" /Correl_DEMG/Isoforms_up_SG_DEMG_promoter_without_overlap.txt", header=T)

cats = split(up_counts_sg, Hmisc::cut2(up_counts_sg$mean_tpm_atac_sg, g=3))
names(cats)
up_sg_levels<-up_counts_sg
# [1] "[ 0.00, 9.09)" "[ 9.09,23.64)" "[23.64,94.67]"


up_sg_levels$rna_categ <- "High"
up_sg_levels$rna_categ[up_counts_sg$mean_tpm_atac_sg < 23.64] <- "Medium"
up_sg_levels$rna_categ[up_counts_sg$mean_tpm_atac_sg < 9.09] <- "Low"
table(up_sg_levels$rna_categ)

#   High    Low Medium
#   72     73     71



alto_sg<-subset(up_sg_levels,rna_categ=="High")

medio_sg<- subset(up_sg_levels,rna_categ=="Medium")

bajo_sg<-subset(up_sg_levels,rna_categ=="Low")

###### I extract the coordinates of the tss to represent in ngsplot

coor<-read.delim(file=" /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed",header = F)
colnames(coor)<-c("Chr","Start","End","isoform_ID","Strand","Score")

clas<-merge(up_sg_levels, coor, by="isoform_ID")
clas<-clas[,c(1,10,16, 20:25)]

alto_sg<-subset(clas,rna_categ=="High")
alto_sg<-alto_sg[,c(5:7,1,8,9,2:4)]
medio_sg<- subset(clas,rna_categ=="Medium")
medio_sg<-medio_sg[,c(5:7,1,8,9,2:4)]
bajo_sg<-subset(clas,rna_categ=="Low")
bajo_sg<-bajo_sg[,c(5:7,1,8,9,2:4)]
colnames(bajo_sg)

write.table(alto_sg, file = " /Correl_DEMG/Genes_high_accessible_SG_DEMG_PROM_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(medio_sg, file = " /Correl_DEMG/Genes_medium_accessible_SG_DEMG_PROM_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(bajo_sg, file = " /Correl_DEMG/Genes_low_accessible_SG_DEMG_PROM_overlap.bed",quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)


### Correlation

#### High accessibility
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(alto_sg$Mean_sg_rna[1:5000]) #

shapiro.test(alto_sg$mean_tpm_atac_sg[1:5000]) #

hist(alto_sg$Mean_sg_rna)
hist(alto_sg$mean_tpm_atac_sg)


### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(alto_sg$Mean_sg_rna,alto_sg$mean_tpm_atac_sg,method="spearman",exact = F)
#0.1470709 , p-value = 0.2176


#### Medium accessibility
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(medio_sg$Mean_sg_rna[1:5000]) #

shapiro.test(medio_sg$mean_tpm_atac_sg[1:5000]) #

hist(medio_sg$Mean_sg_rna)
hist(medio_sg$mean_tpm_atac_sg)


### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(medio_sg$Mean_sg_rna,medio_sg$mean_tpm_atac_sg,method="spearman",exact = F)
#-0.1843853 p-value =0.1237


#### Low accessibility
#https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
#https://stats.stackexchange.com/questions/252373/pearsons-test-of-correlation-or-spearmans-test?noredirect=1&lq=1
shapiro.test(bajo_sg$Mean_sg_rna[1:5000]) #

shapiro.test(bajo_sg$mean_tpm_atac_sg[1:5000]) #

hist(bajo_sg$Mean_sg_rna)
hist(bajo_sg$mean_tpm_atac_sg)

### as they do not follow a normal distribution, we perform the spearman's correlation

cor.test(bajo_sg$Mean_sg_rna,bajo_sg$mean_tpm_atac_sg,method="spearman",exact = F)
# 0.08380324  p-value = 0.4809

#violinplot
p <- ggplot(clas, aes(x=rna_categ, y=Mean_sg_rna,fill=rna_categ)) +  geom_violin(trim = T) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
x11()
p
d<- ggplot(clas, aes(x=rna_categ, y=Mean_sg_rna,fill=rna_categ)) +  geom_boxplot()
d


############# NGSPLOT

# ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_accessible_MG_PROM.txt -O accessibility_levels_in_gene_body_MG
# ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_accessible_SG_PROM.txt -O accessibility_levels_in_gene_body_SG
#
# ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_expresed_MG_all.txt -O expression_levels_in_TSS_MG
# ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_expresed_SG_all.txt -O expression_levels_in_TSS_SG
