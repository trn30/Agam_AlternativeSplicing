library(IsoformSwitchAnalyzeR)
library(BSgenome.Agam.VectorBase.P412)
load("/isoform/RData/MG_vs_SG_repeticion.RData")


#### all genome isoforms separated by alternative splicing mechanisms
mySwitchList <- dataList_paired
### Step (3) Identify differentially used isoforms
mySwitchList <- isoformSwitchTestDEXSeq(  mySwitchList, reduceToSwitchingGenes=FALSE )
### Step (4) Analyze alternative splicing
mySwitchList <- analyzeAlternativeSplicing( mySwitchList, onlySwitchingGenes=F )
### Step (5) Global splicing analysis
# extractSplicingSummary( mySwitchList, returnResult= T ,  splicingToAnalyze = 'all')
# extractSplicingEnrichment( mySwitchList )
# extractSplicingEnrichmentComparison( mySwitchList )

extractSplicingGenomeWide( mySwitchList ,featureToExtract='all',returnResult=TRUE) ### tengo los número pero no las isoformas

View(mySwitchList$AlternativeSplicingAnalysis)

alt_sp<-mySwitchList$AlternativeSplicingAnalysis
DEMG<- read.delim(" /DESeq2/All_isoforms_more_than_one_transcript_MG_vs_SG.txt", header=T)
colnames(DEMG)[1]<-"isoform_id"
alt_sp_mi<-merge(DEMG,alt_sp, by="isoform_id")

#### I HAVE THE MECHANISMS OF ALL THE ISOFORMS OF THE GENOME, AND I WILL CALCULATE THE DISTRIBUTION OF MECHANISMS FOR EACH CONDITION.
inf_con_mg<- read.delim(" /DESeq2/All_isoforms_filtered_by_cpm_more_than_one_transcript_ctrl_inf_MG.txt", header=T) #848DEG
inf_con_sg<- read.delim(" /DESeq2/All_isoforms_filtered_by_cpm_more_than_one_transcript_ctrl_inf_SG.txt", header=T) #977 DEG
inf_sg_mg<- read.delim(" /DESeq2/All_isoforms_filtered_by_cpm_more_than_one_transcript.txt", header=T) #920 DEG
colnames(inf_con_mg)[1]<-"isoform_id"
colnames(inf_con_sg)[1]<-"isoform_id"
colnames(inf_sg_mg)[1]<-"isoform_id"

#### SEPARATED BY MECHANISM
gw_es<- subset(alt_sp_mi, alt_sp_mi$ES!=0) #971
gw_es<-gw_es[,c(1:4)]

gw_mes<- subset(alt_sp_mi, alt_sp_mi$MES!=0)#427
gw_mes<-gw_mes[,c(1,8:10)]

gw_mee<- subset(alt_sp_mi, alt_sp_mi$MEE!=0)#136
gw_mee<-gw_mee[,c(1,5:7)]

gw_ir<- subset(alt_sp_mi, alt_sp_mi$IR!=0)#95
gw_ir<-gw_ir[,c(1,11:13)]

gw_atts<- subset(alt_sp_mi, alt_sp_mi$ATTS!=0)#348
gw_atts<-gw_atts[,c(1,23:25)]

gw_atss<- subset(alt_sp_mi, alt_sp_mi$ATSS!=0)#1025
gw_atss<-gw_atss[,c(1,20:22)]

gw_a5<- subset(alt_sp_mi, alt_sp_mi$A5!=0)#300
gw_a5<-gw_a5[,c(1,14:16)]

gw_a3<- subset(alt_sp_mi, alt_sp_mi$A3!=0)#286
gw_a3<-gw_a3[,c(1,17:19)]


### multiform genes

df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(286, 300, 136,427,348,1025,971,95)
)
library(scales)
df$delabel<- (percent(df$value/sum(df$value)))
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)
#x11()
a<-pie  + geom_text_repel(aes(label=delabel),
                       position = position_stack(vjust = 0.5),
                       show.legend = FALSE, size=3) + ggtitle("Multisoform genes participation of AS mechanims")+ theme(,plot.title = element_text(size=10), text = element_text(size = 5),legend.key.size = unit(0.5,'cm'))

a
##################################################### Control vs Infeccion MG


##### DEMG

inf_con_mgde<- read.delim(" /DESeq2/Isoforms_DEI_Ctrol_vs_Infec_MG_0.05.txt", header=T)
inf_con_mgde$isoform_id <-row.names(inf_con_mgde)

alt_sp_mgde<-merge(inf_con_mgde,alt_sp_mi, by="isoform_id")

es_mgde<-merge(gw_es,alt_sp_mgde, by= "isoform_id") #7
es_mgde[es_mgde$log2FoldChange > 0 ,1]# INF 1
es_mgde[es_mgde$log2FoldChange < 0 ,1]# CTROL 6

mes_mgde<-merge(gw_mes,alt_sp_mgde, by= "isoform_id") #1
mes_mgde[mes_mgde$log2FoldChange > 0 ,1]# INF 0
mes_mgde[mes_mgde$log2FoldChange < 0 ,1]# CTROL 1

mee_mgde<-merge(gw_mee,alt_sp_mgde, by= "isoform_id") #0
mee_mgde[mee_mgde$log2FoldChange > 0 ,1]# INF 0
mee_mgde[mee_mgde$log2FoldChange < 0 ,1]# CTROL 0

ir_mgde<-merge(gw_ir,alt_sp_mgde, by= "isoform_id") #0
ir_mgde[ir_mgde$log2FoldChange > 0 ,1]# INF 0
ir_mgde[ir_mgde$log2FoldChange < 0 ,1]# CTROL 0

atts_mgde<-merge(gw_atts,alt_sp_mgde, by= "isoform_id") #2
atts_mgde[atts_mgde$log2FoldChange > 0 ,1]# INF 0
atts_mgde[atts_mgde$log2FoldChange < 0 ,1]# CTROL 2

atss_mgde<-merge(gw_atss,alt_sp_mgde, by= "isoform_id") #10
atss_mgde[atss_mgde$log2FoldChange > 0 ,1]# INF 2
atss_mgde[atss_mgde$log2FoldChange < 0 ,1]# CTROL 8

a5_mgde<-merge(gw_a5,alt_sp_mgde, by= "isoform_id") #1
a5_mgde[a5_mgde$log2FoldChange > 0 ,1]# INF 0
a5_mgde[a5_mgde$log2FoldChange < 0 ,1]# CTROL 1

a3_mgde<-merge(gw_a3,alt_sp_mgde, by= "isoform_id") #1
a3_mgde[a3_mgde$log2FoldChange > 0 ,1]# INF 0
a3_mgde[a3_mgde$log2FoldChange < 0 ,1]# CTROL 1

##### DUI

inf_con_mgdui<- read.delim(" /Iso_usage/Isoforms_DUI_Ctrol_vs_Infec_MG_0.05.txt", header=T)

colnames(inf_con_mgdui)[3]<-"isoform_id"

alt_sp_mgdui<-merge(inf_con_mgdui,alt_sp_mi, by="isoform_id")

es_mgdui<-merge(gw_es,alt_sp_mgdui, by= "isoform_id") #1
es_mgdui[es_mgdui$dIF > 0 ,1]# ctrol 1
es_mgdui[es_mgdui$dIF < 0 ,1]# inf 0

mes_mgdui<-merge(gw_mes,alt_sp_mgdui, by= "isoform_id") #4
mes_mgdui[mes_mgdui$dIF > 0 ,1]# ctrol 2
mes_mgdui[mes_mgdui$dIF < 0 ,1]# inf 2

mee_mgdui<-merge(gw_mee,alt_sp_mgdui, by= "isoform_id") #0
mee_mgdui[mee_mgdui$dIF > 0 ,1]# ctrol 0
mee_mgdui[mee_mgdui$dIF < 0 ,1]# inf 0

ir_mgdui<-merge(gw_ir,alt_sp_mgdui, by= "isoform_id") #0
ir_mgdui[ir_mgdui$dIF > 0 ,1]# ctrol 0
ir_mgdui[ir_mgdui$dIF < 0 ,1]# inf 0

atts_mgdui<-merge(gw_atts,alt_sp_mgdui, by= "isoform_id") #2
atts_mgdui[atts_mgdui$dIF > 0 ,1]# ctrol 1
atts_mgdui[atts_mgdui$dIF < 0 ,1]# inf 1

atss_mgdui<-merge(gw_atss,alt_sp_mgdui, by= "isoform_id") #3
atss_mgdui[atss_mgdui$dIF > 0 ,1]# ctrol 2
atss_mgdui[atss_mgdui$dIF < 0 ,1]# inf 1

a5_mgdui<-merge(gw_a5,alt_sp_mgdui, by= "isoform_id") #5
a5_mgdui[a5_mgdui$dIF > 0 ,1]# ctrol 3
a5_mgdui[a5_mgdui$dIF < 0 ,1]# inf 2

a3_mgdui<-merge(gw_a3,alt_sp_mgdui, by= "isoform_id") #1
a3_mgdui[a3_mgdui$dIF > 0 ,1]# ctrol 0
a3_mgdui[a3_mgdui$dIF < 0 ,1]# inf 1

######## pie chart
### DEMG
df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(1, 1, 0,1,2,10,7,0)
)
df$delabel<- (percent(df$value/sum(df$value)))
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")


pie <- bp + coord_polar("y", start=0)
#x11()
b<-pie + geom_text_repel(aes(label=delabel),
                      position = position_stack(vjust = 0.5),
                      show.legend = FALSE, size=3) + ggtitle("Participation of AS mechanims in Ctrl vs. Inf MG (DEMG)")+ theme(legend.position = 'none',plot.title = element_text(size=10), text = element_text(size = 1))
b
### DUI

df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(1,5,0,4,2,3,1,0)
)
df$delabel<- (percent(df$value/sum(df$value)))
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")


pie <- bp + coord_polar("y", start=0)
#x11()
c<-pie + geom_text_repel(aes(label=delabel),
                      position = position_stack(vjust = 0.5),
                      show.legend = FALSE, size=3)  + ggtitle("Participation of AS mechanims in Ctrl vs. Inf MG (DUI)")+ theme(legend.position = 'none',plot.title = element_text(size=10), text = element_text(size = 1))


c
#################################################################### Control vs Infeccion SG


##### DEMG

inf_con_sgde<- read.delim(" /DESeq2/Isoforms_DEI_Ctrol_vs_Infec_SG_0.05.txt", header=T) #56DEG
inf_con_sgde$isoform_id <-row.names(inf_con_sgde)


alt_sp_sgde<-merge(inf_con_sgde,alt_sp_mi, by="isoform_id")

es_sgde<-merge(gw_es,alt_sp_sgde, by= "isoform_id") #15
es_sgde[es_sgde$log2FoldChange > 0 ,1]# INF 2
es_sgde[es_sgde$log2FoldChange < 0 ,1]# CTROL 13

mes_sgde<-merge(gw_mes,alt_sp_sgde, by= "isoform_id") #10
mes_sgde[mes_sgde$log2FoldChange > 0 ,1]# INF 4
mes_sgde[mes_sgde$log2FoldChange < 0 ,1]# CTROL 6

mee_sgde<-merge(gw_mee,alt_sp_sgde, by= "isoform_id") #3
mee_sgde[mee_sgde$log2FoldChange > 0 ,1]# INF 0
mee_sgde[mee_sgde$log2FoldChange < 0 ,1]# CTROL 3

ir_sgde<-merge(gw_ir,alt_sp_sgde, by= "isoform_id") #2
ir_sgde[ir_sgde$log2FoldChange > 0 ,1]# INF 1
ir_sgde[ir_sgde$log2FoldChange < 0 ,1]# CTROL 1

atts_sgde<-merge(gw_atts,alt_sp_sgde, by= "isoform_id") #5
atts_sgde[atts_sgde$log2FoldChange > 0 ,1]# INF 3
atts_sgde[atts_sgde$log2FoldChange < 0 ,1]# CTROL 2

atss_sgde<-merge(gw_atss,alt_sp_sgde, by= "isoform_id") #32
atss_sgde[atss_sgde$log2FoldChange > 0 ,1]# INF 4
atss_sgde[atss_sgde$log2FoldChange < 0 ,1]# CTROL 28

a5_sgde<-merge(gw_a5,alt_sp_sgde, by= "isoform_id") #6
a5_sgde[a5_sgde$log2FoldChange > 0 ,1]# INF 0
a5_sgde[a5_sgde$log2FoldChange < 0 ,1]# CTROL 6

a3_sgde<-merge(gw_a3,alt_sp_sgde, by= "isoform_id") #3
a3_sgde[a3_sgde$log2FoldChange > 0 ,1]# INF 1
a3_sgde[a3_sgde$log2FoldChange < 0 ,1]# CTROL 2

##### DUI


inf_con_sgdui<- read.delim(" /Iso_usage/Isoforms_DUI_Ctrol_vs_Infec_SG_0.05.txt", header=T) #35

colnames(inf_con_sgdui)[3]<-"isoform_id"

alt_sp_sgdui<-merge(inf_con_sgdui,alt_sp_mi, by="isoform_id")

es_sgdui<-merge(gw_es,alt_sp_sgdui, by= "isoform_id") #9
es_sgdui[es_sgdui$dIF > 0 ,1]# ctrol 6
es_sgdui[es_sgdui$dIF < 0 ,1]# inf 3

mes_sgdui<-merge(gw_mes,alt_sp_sgdui, by= "isoform_id") #2
mes_sgdui[mes_sgdui$dIF > 0 ,1]# ctrol 1
mes_sgdui[mes_sgdui$dIF < 0 ,1]# inf 1

mee_sgdui<-merge(gw_mee,alt_sp_sgdui, by= "isoform_id") #2
mee_sgdui[mee_sgdui$dIF > 0 ,1]# ctrol 1
mee_sgdui[mee_sgdui$dIF < 0 ,1]# inf 1

ir_sgdui<-merge(gw_ir,alt_sp_sgdui, by= "isoform_id") #5
ir_sgdui[ir_sgdui$dIF > 0 ,1]# ctrol 2
ir_sgdui[ir_sgdui$dIF < 0 ,1]# inf 3

atts_sgdui<-merge(gw_atts,alt_sp_sgdui, by= "isoform_id") #5
atts_sgdui[atts_sgdui$dIF > 0 ,1]# ctrol 3
atts_sgdui[atts_sgdui$dIF < 0 ,1]# inf 2

atss_sgdui<-merge(gw_atss,alt_sp_sgdui, by= "isoform_id") #11
atss_sgdui[atss_sgdui$dIF > 0 ,1]# ctrol 4
atss_sgdui[atss_sgdui$dIF < 0 ,1]# inf 7

a5_sgdui<-merge(gw_a5,alt_sp_sgdui, by= "isoform_id") #5
a5_sgdui[a5_sgdui$dIF > 0 ,1]# ctrol 1
a5_sgdui[a5_sgdui$dIF < 0 ,1]# inf 4

a3_sgdui<-merge(gw_a3,alt_sp_sgdui, by= "isoform_id") #4
a3_sgdui[a3_sgdui$dIF > 0 ,1]# ctrol 2
a3_sgdui[a3_sgdui$dIF < 0 ,1]# inf 2

######## pie chart
### DEMG
df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(3,6,3,10,5,32,15,2)
)
df$delabel<- (percent(df$value/sum(df$value)))
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")


pie <- bp + coord_polar("y", start=0)
#x11()
d<-pie + geom_text_repel(aes(label=delabel),
                      position = position_stack(vjust = 0.5),
                      show.legend = FALSE, size=3) + ggtitle("Participation of AS mechanims in Ctrl vs. Inf SG (DEMG)")+ theme(legend.position = 'none',plot.title = element_text(size=10), text = element_text(size = 1))

d

### DUI

df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(4,5,2,2,5,11,9,5)
)
df$delabel<- (percent(df$value/sum(df$value)))
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")


pie <- bp + coord_polar("y", start=0)
#x11()
e<-pie + geom_text_repel(aes(label=delabel),
                      position = position_stack(vjust = 0.5),
                      show.legend = FALSE, size=3)  + ggtitle("Participation of AS mechanims in Ctrl vs. Inf SG (DUI)")+ theme(legend.position = 'none',plot.title = element_text(size=10), text = element_text(size = 1))

e


############################################ Infeccion Mg vs Infeccion SG

#####DEMG
sg_sgmg_de<- read.delim(" /DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T) #463
sg_sgmg_de$isoform_id <-row.names(sg_sgmg_de)

alt_sp_sgmg_de<-merge(alt_sp_mi,sg_sgmg_de, by="isoform_id")

es_sgmg_de<-merge(gw_es,alt_sp_sgmg_de, by= "isoform_id") #170
es_sgmg_de[es_sgmg_de$log2FoldChange > 0 ,1]# sg 102
es_sgmg_de[es_sgmg_de$log2FoldChange < 0 ,1]# mg 68

mes_sgmg_de<-merge(gw_mes,alt_sp_sgmg_de, by= "isoform_id") #73
mes_sgmg_de[mes_sgmg_de$log2FoldChange > 0 ,1]# sg 35
mes_sgmg_de[mes_sgmg_de$log2FoldChange < 0 ,1]# mg 38

mee_sgmg_de<-merge(gw_mee,alt_sp_sgmg_de, by= "isoform_id") #17
mee_sgmg_de[mee_sgmg_de$log2FoldChange > 0 ,1]# sg 6
mee_sgmg_de[mee_sgmg_de$log2FoldChange < 0 ,1]# mg 11

ir_sgmg_de<-merge(gw_ir,alt_sp_sgmg_de, by= "isoform_id") #16
ir_sgmg_de[ir_sgmg_de$log2FoldChange > 0 ,1]# sg 6
ir_sgmg_de[ir_sgmg_de$log2FoldChange < 0 ,1]# mg 10

atts_sgmg_de<-merge(gw_atts,alt_sp_sgmg_de, by= "isoform_id") #51
atts_sgmg_de[atts_sgmg_de$log2FoldChange > 0 ,1]# sg 29
atts_sgmg_de[atts_sgmg_de$log2FoldChange < 0 ,1]# mg 22

atss_sgmg_de<-merge(gw_atss,alt_sp_sgmg_de, by= "isoform_id") #194
atss_sgmg_de[atss_sgmg_de$log2FoldChange > 0 ,1]# sg 105
atss_sgmg_de[atss_sgmg_de$log2FoldChange < 0 ,1]# mg 89

a5_sgmg_de<-merge(gw_a5,alt_sp_sgmg_de, by= "isoform_id") #44
a5_sgmg_de[a5_sgmg_de$log2FoldChange > 0 ,1]# sg 22
a5_sgmg_de[a5_sgmg_de$log2FoldChange < 0 ,1]# mg 22

a3_sgmg_de<-merge(gw_a3,alt_sp_sgmg_de, by= "isoform_id") #34
a3_sgmg_de[a3_sgmg_de$log2FoldChange > 0 ,1]# sg 9
a3_sgmg_de[a3_sgmg_de$log2FoldChange < 0 ,1]# mg 25


#####DUI
inf_sgmg_dui<- read.delim(" /Iso_usage/Isoforms_DUI_MG_vs_SG_0.05_counts.txt", header=T) #848DEG

colnames(inf_sgmg_dui)[1]<-"isoform_id"


alt_sp_sgmg_dui<-merge(alt_sp_mi,inf_sgmg_dui, by="isoform_id")


es_sgmg_dui<-merge(gw_es,alt_sp_sgmg_dui, by= "isoform_id") #90
es_sgmg_dui[es_sgmg_dui$dIF > 0 ,1]# sg 46
es_sgmg_dui[es_sgmg_dui$dIF < 0 ,1]# mg 44

mes_sgmg_dui<-merge(gw_mes,alt_sp_sgmg_dui, by= "isoform_id") #51
mes_sgmg_dui[mes_sgmg_dui$dIF > 0 ,1]# sg 27
mes_sgmg_dui[mes_sgmg_dui$dIF < 0 ,1]# mg 24

mee_sgmg_dui<-merge(gw_mee,alt_sp_sgmg_dui, by= "isoform_id") #13
mee_sgmg_dui[mee_sgmg_dui$dIF > 0 ,1]# sg 5
mee_sgmg_dui[mee_sgmg_dui$dIF < 0 ,1]# mg 8

ir_sgmg_dui<-merge(gw_ir,alt_sp_sgmg_dui, by= "isoform_id") #6
ir_sgmg_dui[ir_sgmg_dui$dIF > 0 ,1]# sg 2
ir_sgmg_dui[ir_sgmg_dui$dIF < 0 ,1]# mg 4

atts_sgmg_dui<-merge(gw_atts,alt_sp_sgmg_dui, by= "isoform_id") #18
atts_sgmg_dui[atts_sgmg_dui$dIF > 0 ,1]# sg 11
atts_sgmg_dui[atts_sgmg_dui$dIF < 0 ,1]# mg 7

atss_sgmg_dui<-merge(gw_atss,alt_sp_sgmg_dui, by= "isoform_id") #94
atss_sgmg_dui[atss_sgmg_dui$dIF > 0 ,1]# sg 43
atss_sgmg_dui[atss_sgmg_dui$dIF < 0 ,1]# mg 51

a5_sgmg_dui<-merge(gw_a5,alt_sp_sgmg_dui, by= "isoform_id") #16
a5_sgmg_dui[a5_sgmg_dui$dIF > 0 ,1]# sg 10
a5_sgmg_dui[a5_sgmg_dui$dIF < 0 ,1]# mg 6

a3_sgmg_dui<-merge(gw_a3,alt_sp_sgmg_dui, by= "isoform_id") #13
a3_sgmg_dui[a3_sgmg_dui$dIF > 0 ,1]# sg 2
a3_sgmg_dui[a3_sgmg_dui$dIF < 0 ,1]# mg 2

write.table((es_sgmg_dui), file = " /Isoforms_DUI_MG_vs_SG_ES_ID.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
write.table((atss_sgmg_dui), file = " /Isoforms_DUI_MG_vs_SG_ATSSS_ID.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)

write.table((atss_sgmg_de), file = " /Isoforms_DEMG_MG_vs_SG_ATSS_ID.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)
write.table((es_sgmg_de), file = " /Isoforms_DEMG_MG_vs_SG_ES_ID.txt",quote = FALSE, sep="\t", row.names = T , col.names = T)


######## pie chart
### DEMG
df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(34,44,17,73,51,194,170,16)
)
df$delabel<- (percent(df$value/sum(df$value)))
#df$delabel<- paste(percent(df$value/sum(df$value)),df$Mechanism,sep=" ")

head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")


pie <- bp + coord_polar("y", start=0)
#x11()
f<-pie + geom_text_repel(aes(label=delabel),
                      position = position_stack(vjust = 0.5),
                      show.legend = FALSE, size=3) + ggtitle("Participation of AS mechanims in Inf MG vs. Inf SG (DEMG)") + theme(legend.position = 'none',plot.title = element_text(size=10), text = element_text(size = 1))

f

### DUI

df <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  value = c(13,16,13,51,18,94,90,6)
)
df$delabel<- (percent(df$value/sum(df$value)))
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=Mechanism))+
  geom_bar(width = 1, stat = "identity")


pie <- bp + coord_polar("y", start=0)

g<-pie + geom_text_repel(aes(label=delabel),
                      position = position_stack(vjust = 0.5),
                      show.legend = F, size=3)  + ggtitle("Participation of AS mechanims in Inf MG vs. Inf SG (DUI)") + theme(legend.position = 'none',plot.title = element_text(size=10), text = element_text(size = 1))

g
library(tableExtra)
h<- rectGrob(gp=gpar(fill="grey90"))

library(gridExtra)
grid.arrange(b,c,d,e,f,g,a,h,ncol=2)
pdf( file = " /figures/Piechart_composition_AS_mechanisms.pdf", onefile = FALSE, height=6, width = 9)
grid.arrange(b,c,d,e,f,g,a,h,ncol=2)
dev.off()

########################################################################################
#    accumulated hitogram
#########################################################################################
all<- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(286, 300, 136,427,348,1025,971,95)
)
all$Comparison<-"Multisoform Genes"
all$percent<-(all$Frequency/sum(all$Frequency))*100

demg_mg <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(1, 1, 0,1,2,10,7,0)
)
demg_mg$Comparison<-"Ctrl vs Inf MG"
demg_mg$percent<-(demg_mg$Frequency/sum(demg_mg$Frequency))*100

dui_mg<- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(1,5,0,4,2,3,1,0))

  dui_mg$Comparison<-"Ctrl vs Inf MG"
  dui_mg$percent<-(dui_mg$Frequency/sum(dui_mg$Frequency))*100

  demg_sg<- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(3,6,3,10,5,32,15,2)
)
  demg_sg$Comparison<-"Ctrl vs Inf SG"
  demg_sg$percent<-(demg_sg$Frequency/sum(demg_sg$Frequency))*100

dui_sg <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(4,5,2,2,5,11,9,5)
)
dui_sg$Comparison<-"Ctrl vs Inf SG"
dui_sg$percent<-(dui_sg$Frequency/sum(dui_sg$Frequency))*100

demg_mgsg <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(34,44,17,73,51,194,170,16)
)
demg_mgsg$Comparison<-"Inf MG vs Inf SG"
demg_mgsg$percent<-(demg_mgsg$Frequency/sum(demg_mgsg$Frequency))*100

dui_mgsg <- data.frame(
  Mechanism = c("A3", "A5", "MEE","MES", "ATTS", "ATSS","ES","IR"),
  Frequency = c(13,16,13,51,18,94,90,6)
)
dui_mgsg$Comparison<-"Inf MG vs Inf SG"
dui_mgsg$percent<-(dui_mgsg$Frequency/sum(dui_mgsg$Frequency))*100

datoscompletodui<-rbind(dui_sg,dui_mg,dui_mgsg,all)
datoscompletodemg<-rbind(demg_sg,demg_mg,demg_mgsg,all)

ggplot(datoscompletodemg,aes(x=Comparison,y=Frequency,fill= Mechanism))+geom_bar(stat = "identity", position = position_fill())
ggplot(datoscompletodui,aes(x=Comparison,y=Frequency,fill= Mechanism))+geom_bar(stat = "identity", position = position_fill())
