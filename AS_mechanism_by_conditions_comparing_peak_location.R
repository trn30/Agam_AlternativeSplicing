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

##### MG vs SG
#####DUI
inf_sgmg_dui<- read.delim(" /Iso_usage/Isoforms_DUI_MG_vs_SG_0.05_counts.txt", header=T) #848DEG

colnames(inf_sgmg_dui)[1]<-"isoform_id"


alt_sp_sgmg_dui<-merge(alt_sp_mi,inf_sgmg_dui, by="isoform_id")


es_sgmg_dui<-merge(gw_es,alt_sp_sgmg_dui, by= "isoform_id") #90

mes_sgmg_dui<-merge(gw_mes,alt_sp_sgmg_dui, by= "isoform_id") #51

mee_sgmg_dui<-merge(gw_mee,alt_sp_sgmg_dui, by= "isoform_id") #13

ir_sgmg_dui<-merge(gw_ir,alt_sp_sgmg_dui, by= "isoform_id") #6

atts_sgmg_dui<-merge(gw_atts,alt_sp_sgmg_dui, by= "isoform_id") #18

atss_sgmg_dui<-merge(gw_atss,alt_sp_sgmg_dui, by= "isoform_id") #94

a5_sgmg_dui<-merge(gw_a5,alt_sp_sgmg_dui, by= "isoform_id") #16

a3_sgmg_dui<-merge(gw_a3,alt_sp_sgmg_dui, by= "isoform_id") #13

##### I am going to associate each mechanism with a zone

es_sgmg_dui$Region<- "Body"
es_sgmg_dui<-es_sgmg_dui[,c(1:4,33,65)]
###### now with the peaks I will see to which zones these peaks are noted compared to the zone expected by the mechanism
peaks <- read.delim(" /ATAC/results_JL/ATAC-seq_peaks.txt") #193005
peak<-peaks[peaks$Infections..intersect. == "I1,I2" ,]

peak$Strand<-"+"
peak$Score<-0
peak$Chromosome<-paste("AgamP4",peak$Chromosome,sep = "_")
peak$isoform_id<-peak$Nearest.Promoter.ID..HOMER.

# es
es_sgmg_dui$Region<- "Body"
es_sgmg_dui<-es_sgmg_dui[,c(1:4,33,65)]
es_sgmg_dui_2<-merge(es_sgmg_dui,peak,by="isoform_id")
es_sgmg_dui_2<-es_sgmg_dui_2[,c(1:10,12,18)]

es_sgmg_dui_3<-es_sgmg_dui_2 %>% group_by(Annotation.final) %>%tally() #93

#atss
atss_sgmg_dui$Region<- "Promoter"
atss_sgmg_dui<-atss_sgmg_dui[,c(1:4,33,65)]
atss_sgmg_dui_2<-merge(atss_sgmg_dui,peak,by="isoform_id")
atss_sgmg_dui_2<-atss_sgmg_dui_2[,c(1:10,12,18)]

atss_sgmg_dui_3<-atss_sgmg_dui_2 %>% group_by(Annotation.final) %>%tally() #93


#atts
atts_sgmg_dui$Region<- "Downstream"
atts_sgmg_dui<-atts_sgmg_dui[,c(1:4,33,65)]
atts_sgmg_dui_2<-merge(atts_sgmg_dui,peak,by="isoform_id")
atts_sgmg_dui_2<-atts_sgmg_dui_2[,c(1:10,12,18)]

atts_sgmg_dui_3<-atts_sgmg_dui_2 %>% group_by(Annotation.final) %>%tally() #93


#####DEMG
inf_sgmg_de<- read.delim(" /DESeq2/Isoforms_DEI_MG_vs_SG_0.05.txt", header=T) #848DEG
inf_sgmg_de$isoform_id<-rownames(inf_sgmg_de)


alt_sp_sgmg_de<-merge(alt_sp_mi,inf_sgmg_de, by="isoform_id")

es_sgmg_de<-merge(gw_es,alt_sp_sgmg_de, by= "isoform_id") #170

mes_sgmg_de<-merge(gw_mes,alt_sp_sgmg_de, by= "isoform_id") #73

mee_sgmg_de<-merge(gw_mee,alt_sp_sgmg_de, by= "isoform_id") #17

ir_sgmg_de<-merge(gw_ir,alt_sp_sgmg_de, by= "isoform_id") #16

atts_sgmg_de<-merge(gw_atts,alt_sp_sgmg_de, by= "isoform_id") #51

atss_sgmg_de<-merge(gw_atss,alt_sp_sgmg_de, by= "isoform_id") #194

a5_sgmg_de<-merge(gw_a5,alt_sp_sgmg_de, by= "isoform_id") #44

a3_sgmg_de<-merge(gw_a3,alt_sp_sgmg_de, by= "isoform_id") #34


# es
es_sgmg_de$Region<- "Body"
es_sgmg_de<-es_sgmg_de[,c(1:4,30,35)]
es_sgmg_de_2<-merge(es_sgmg_de,peak,by="isoform_id")
es_sgmg_de_2<-es_sgmg_de_2[,c(1:10,12,18)]

es_sgmg_de_3<-es_sgmg_de_2 %>% group_by(Annotation.final) %>%tally() #93

#atss
atss_sgmg_de$Region<- "Promoter"
atss_sgmg_de<-atss_sgmg_de[,c(1:4,30,35)]
atss_sgmg_de_2<-merge(atss_sgmg_de,peak,by="isoform_id")
atss_sgmg_de_2<-atss_sgmg_de_2[,c(1:10,12,18)]

atss_sgmg_de_3<-atss_sgmg_de_2 %>% group_by(Annotation.final) %>%tally() #93


#atts
atts_sgmg_de$Region<- "Downstream"
atts_sgmg_de<-atts_sgmg_de[,c(1:4,30,35)]
atts_sgmg_de_2<-merge(atts_sgmg_de,peak,by="isoform_id")
atts_sgmg_de_2<-atts_sgmg_de_2[,c(1:10,12,18)]

atts_sgmg_de_3<-atts_sgmg_de_2 %>% group_by(Annotation.final) %>%tally() #93

#    cumulative histogram

#########################################################################################
es_sgmg_dui_3$Type_mechanism<-"Exon skipping"
atss_sgmg_dui_3$Type_mechanism<-"Alternative start site"
atts_sgmg_dui_3$Type_mechanism<-"Alternative termination site"

es_sgmg_de_3$Type_mechanism<-"Exon skipping"
atss_sgmg_de_3$Type_mechanism<-"Alternative start site"
atts_sgmg_de_3$Type_mechanism<-"Alternative termination site"

datoscompletodui<-rbind(es_sgmg_dui_3,atss_sgmg_dui_3,atts_sgmg_dui_3)
datoscompletodemg<-rbind(es_sgmg_de_3,atss_sgmg_de_3,atts_sgmg_de_3)

ggplot(datoscompletodemg,aes(x=Type_mechanism,y=n,fill= Annotation.final))+geom_bar(stat = "identity", position = position_fill())
ggplot(datoscompletodui,aes(x=Type_mechanism,y=n,fill= Annotation.final))+geom_bar(stat = "identity", position = position_fill())
