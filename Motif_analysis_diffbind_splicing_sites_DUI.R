###### MOTIF ANALYSIS
### promoter + gene body

# http://homer.ucsd.edu/homer/motif/
# http://homer.ucsd.edu/homer/ngs/peakMotifs.html
# http://homer.ucsd.edu/homer/motif/practicalTips.html
# findMotifsGenome.pl peaks.txt mm8r peakAnalysis -size 200 -len 8

### we will analyse the diffbind peaks that coincide with DUI isoforms located with splicing sites (98 in total = 95 + 3).

dui_eve<-read.delim(" /Correl_by_mechanisms/DUI_Isoforms_events_whitout_overlapping.bed",header = T)

diffb_MG<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed", header=F)
colnames(diffb_MG)<-c( "Chr_region","Start_region","End_event","Region_ID","Strand_region","Score_region","Chr_event","Start_event","End_region","event_ID","strand","score")
diffb_MG_def<-merge(diffb_MG,dui_eve,by="event_ID")
unique(diffb_MG_def$Region_ID)#5

diffb_SG<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed", header=F)
colnames(diffb_SG)<-c( "Chr_region","Start_region","End_event","Region_ID","Strand_region","Score_region","Chr_event","Start_event","End_region","event_ID","strand","score")
diffb_SG_def<-merge(diffb_SG,dui_eve,by="event_ID")
unique(diffb_SG_def$Region_ID)#86
### junto los dos tejidos
diffb_def_all<-rbind(diffb_MG_def,diffb_SG_def)
### For the motif analysis I will keep the regions of the diffbind peaks coinciding with splicing sites and DUI.

diffb_def_all<-diffb_def_all[,c(2,3,10,5,11:12)]
diffb_def_all2<-unique(diffb_def_all)

write.table((diffb_def_all2), file = " /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

#### findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]

## findMotifsGenome.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE -p 20 -len 8,10,12 -size -100,100 -mset insects -dumpFasta -preparse -N 98
#######  now I would do the same for the peaks that meet the above conditions AND are exon skipping events.
diffb_ES_all<- read.delim(" /ATAC/DiffBind/Diffbind_peaks_coincident_with_splicing_sites_ES_and_DUI.txt", header=T)
### For the motif analysis I will keep the regions of the diffbind peaks coinciding with splicing sites ES of the DUIs.
diffb_ES_all<-diffb_ES_all[,c(6,11,12,1,7,8)]
diffb_ES_all<-unique(diffb_ES_all)
write.table((diffb_ES_all), file = " /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_exon_skipping_motif_analysis.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)


##  findMotifsGenome.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_exon_skipping_motif_analysis.bed  genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE_ES -p 20 -len 8,10,12 -size -100,100 -mset insects -dumpFasta -preparse -N 22



##########################################################################
# ANOTATION OF REASONS concerning the first analysis (all splicing events)

#### KNOWN MOTIFS

### MOTIF 1 CATCMCTA

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/knownResults/known1.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_KNWON_annotatePeaks_motif1_CATCMCTA.txt # 1e-3

### MOTIF 2 GGYCATAAAW

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/knownResults/known2.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_KNOWN_annotatePeaks_motif2_GGYCATAAAW.txt # 1e-2

#### HOMER MOTIFS

### MOTIF 1 CACCWATT

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif1.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif1_CACCWATT.txt # 1e-34

### MOTIF 2 CTTGCAAGGT

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif2.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif2_CTTGCAAGGT.txt # 1e-31

### MOTIF 3 ACARAGAAAT

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif3.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif3_ACARAGAAAT.txt # 1e-28

### MOTIF 4 ACACAAAG

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif4.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif4_ACACAAAG.txt # 1e-27

### MOTIF 5 GAGTCTTTAT

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif5.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif5_GAGTCTTTAT.txt # 1e-27

### MOTIF 6 GAAAACAGTT

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif6.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif6_GAAAACAGTT.txt # 1e-27
### MOTIF 7 AAACAGTAAAAG

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif7.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif7_AAACAGTAAAAG.txt # 1e-27
### MOTIF 8 CAAGAAGA

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa  -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif8.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif8_CAAGAAGA.txt # 1e-27
### MOTIF 9 GAMCGAGH

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif9.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif9_GAMCGAGH.txt # 1e-27

### MOTIF 10 TTATAAYT

annotatePeaks.pl  /ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/homerResults/motif10.motif >  /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif10_TTATAAYT.txt # 1e-27

###########################################################################################

###### now I will see with which DUI isoforms these motifs are annotated.

#### KNOWN MOTIFS

### MOTIF 1 CATCMCTA
known_1<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_KNOWN_annotatePeaks_motif1_CATCMCTA.txt", header=T)
colnames(known_1)[1]<- "Peak_ID"
known_1<-known_1[,c(1:6,11,22)]
known_1$Nearest.PromoterID<-gsub("ID=","",known_1$Nearest.PromoterID)
known_1$Nearest.PromoterID<-gsub(";.*","",known_1$Nearest.PromoterID)
known_1$Nearest.PromoterID<-gsub("-E.*","",known_1$Nearest.PromoterID)
known_1$Nearest.PromoterID<-gsub("exon_","",known_1$Nearest.PromoterID)
colnames(known_1)[7]<- "ID"

known_12<-known_1
known_12[known_12==""]<-NA
known_1d<-na.omit(known_12)#14

known_1dui<-merge(known_1d,dui,by="ID")
unique(known_1dui$ID)#10
unique(known_1dui$Peak_ID)#10 peaks with the motif

### MOTIF 2 GGYCATAAAW
known_2<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_KNOWN_annotatePeaks_motif2_GGYCATAAAW.txt", header=T)
colnames(known_2)[1]<- "Peak_ID"
known_2<-known_2[,c(1:6,11,22)]
known_2$Nearest.PromoterID<-gsub("ID=","",known_2$Nearest.PromoterID)
known_2$Nearest.PromoterID<-gsub(";.*","",known_2$Nearest.PromoterID)
known_2$Nearest.PromoterID<-gsub("-E.*","",known_2$Nearest.PromoterID)
known_2$Nearest.PromoterID<-gsub("exon_","",known_2$Nearest.PromoterID)
colnames(known_2)[7]<- "ID"

known_22<-known_2
known_22[known_22==""]<-NA
known_2d<-na.omit(known_22)#12 peaks with the motif

known_2dui<-merge(known_2d,dui,by="ID")
unique(known_2dui$ID)#10
unique(known_2dui$Peak_ID)#10

#### HOMER MOTIFS

### MOTIF 1 CACCWATT
homer_1<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif1_CACCWATT.txt", header=T)
colnames(homer_1)[1]<- "Peak_ID"
homer_1<-homer_1[,c(1:6,11,22)]
homer_1$Nearest.PromoterID<-gsub("ID=","",homer_1$Nearest.PromoterID)
homer_1$Nearest.PromoterID<-gsub(";.*","",homer_1$Nearest.PromoterID)
homer_1$Nearest.PromoterID<-gsub("-E.*","",homer_1$Nearest.PromoterID)
homer_1$Nearest.PromoterID<-gsub("exon_","",homer_1$Nearest.PromoterID)
colnames(homer_1)[7]<- "ID"

homer_12<-homer_1
homer_12[homer_12==""]<-NA
homer_1d<-na.omit(homer_12)#35 peaks with the motif

homer_1dui<-merge(homer_1d,dui,by="ID")
unique(homer_1dui$ID)#22
unique(homer_1dui$Peak_ID)#25

### MOTIF 2 CTTGCAAGGT
homer_2<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif2_CTTGCAAGGT.txt", header=T)
colnames(homer_2)[1]<- "Peak_ID"
homer_2<-homer_2[,c(1:6,11,22)]
homer_2$Nearest.PromoterID<-gsub("ID=","",homer_2$Nearest.PromoterID)
homer_2$Nearest.PromoterID<-gsub(";.*","",homer_2$Nearest.PromoterID)
homer_2$Nearest.PromoterID<-gsub("-E.*","",homer_2$Nearest.PromoterID)
homer_2$Nearest.PromoterID<-gsub("exon_","",homer_2$Nearest.PromoterID)
colnames(homer_2)[7]<- "ID"

homer_22<-homer_2
homer_22[homer_22==""]<-NA
homer_2d<-na.omit(homer_22)#30 peaks with the motif

homer_2dui<-merge(homer_2d,dui,by="ID")
unique(homer_2dui$ID)#21
unique(homer_2dui$Peak_ID)#25

### MOTIF 3 ACARAGAAAT

  homer_3<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif3_ACARAGAAAT.txt", header=T)
colnames(homer_3)[1]<- "Peak_ID"
homer_3<-homer_3[,c(1:6,11,22)]
homer_3$Nearest.PromoterID<-gsub("ID=","",homer_3$Nearest.PromoterID)
homer_3$Nearest.PromoterID<-gsub(";.*","",homer_3$Nearest.PromoterID)
homer_3$Nearest.PromoterID<-gsub("-E.*","",homer_3$Nearest.PromoterID)
homer_3$Nearest.PromoterID<-gsub("exon_","",homer_3$Nearest.PromoterID)
colnames(homer_3)[7]<- "ID"

homer_32<-homer_3
homer_32[homer_32==""]<-NA
homer_3d<-na.omit(homer_32)#36 peaks with the motif

homer_3dui<-merge(homer_3d,dui,by="ID")
unique(homer_3dui$ID)#26
unique(homer_3dui$Peak_ID)#30


### MOTIF 4 ACACAAAG
homer_4<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif4_ACACAAAG.txt", header=T)
colnames(homer_4)[1]<- "Peak_ID"
homer_4<-homer_4[,c(1:6,11,22)]
homer_4$Nearest.PromoterID<-gsub("ID=","",homer_4$Nearest.PromoterID)
homer_4$Nearest.PromoterID<-gsub(";.*","",homer_4$Nearest.PromoterID)
homer_4$Nearest.PromoterID<-gsub("-E.*","",homer_4$Nearest.PromoterID)
homer_4$Nearest.PromoterID<-gsub("exon_","",homer_4$Nearest.PromoterID)
colnames(homer_4)[7]<- "ID"

homer_42<-homer_4
homer_42[homer_42==""]<-NA
homer_4d<-na.omit(homer_42)#39 peaks with the motif

homer_4dui<-merge(homer_4d,dui,by="ID")
unique(homer_4dui$ID)#26
unique(homer_4dui$Peak_ID)#31


### MOTIF 5 GAGTCTTTAT
homer_5<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif5_GAGTCTTTAT.txt", header=T)
colnames(homer_5)[1]<- "Peak_ID"
homer_5<-homer_5[,c(1:6,11,22)]
homer_5$Nearest.PromoterID<-gsub("ID=","",homer_5$Nearest.PromoterID)
homer_5$Nearest.PromoterID<-gsub(";.*","",homer_5$Nearest.PromoterID)
homer_5$Nearest.PromoterID<-gsub("-E.*","",homer_5$Nearest.PromoterID)
homer_5$Nearest.PromoterID<-gsub("exon_","",homer_5$Nearest.PromoterID)
colnames(homer_5)[7]<- "ID"

homer_52<-homer_5
homer_52[homer_52==""]<-NA
homer_5d<-na.omit(homer_52)#28 peaks with the motif

homer_5dui<-merge(homer_5d,dui,by="ID")
unique(homer_5dui$ID)#20
unique(homer_5dui$Peak_ID)#24

### MOTIF 6 GAAAACAGTT
homer_6<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif6_GAAAACAGTT.txt", header=T)
colnames(homer_6)[1]<- "Peak_ID"
homer_6<-homer_6[,c(1:6,11,22)]
homer_6$Nearest.PromoterID<-gsub("ID=","",homer_6$Nearest.PromoterID)
homer_6$Nearest.PromoterID<-gsub(";.*","",homer_6$Nearest.PromoterID)
homer_6$Nearest.PromoterID<-gsub("-E.*","",homer_6$Nearest.PromoterID)
homer_6$Nearest.PromoterID<-gsub("exon_","",homer_6$Nearest.PromoterID)
colnames(homer_6)[7]<- "ID"

homer_62<-homer_6
homer_62[homer_62==""]<-NA
homer_6d<-na.omit(homer_62)#31 peaks with the motif

homer_6dui<-merge(homer_6d,dui,by="ID")
unique(homer_6dui$ID)#21
unique(homer_6dui$Peak_ID)#23

### MOTIF 7 AAACAGTAAAAG
homer_7<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif7_AAACAGTAAAAG.txt", header=T)
colnames(homer_7)[1]<- "Peak_ID"
homer_7<-homer_7[,c(1:6,11,22)]
homer_7$Nearest.PromoterID<-gsub("ID=","",homer_7$Nearest.PromoterID)
homer_7$Nearest.PromoterID<-gsub(";.*","",homer_7$Nearest.PromoterID)
homer_7$Nearest.PromoterID<-gsub("-E.*","",homer_7$Nearest.PromoterID)
homer_7$Nearest.PromoterID<-gsub("exon_","",homer_7$Nearest.PromoterID)
colnames(homer_7)[7]<- "ID"

homer_72<-homer_7
homer_72[homer_72==""]<-NA
homer_7d<-na.omit(homer_72)#31 peaks with the motif

homer_7dui<-merge(homer_7d,dui,by="ID")
unique(homer_7dui$ID)#20
unique(homer_7dui$Peak_ID)#22

### MOTIF 8 CAAGAAGA
homer_8<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif8_CAAGAAGA.txt", header=T)
colnames(homer_8)[1]<- "Peak_ID"
homer_8<-homer_8[,c(1:6,11,22)]
homer_8$Nearest.PromoterID<-gsub("ID=","",homer_8$Nearest.PromoterID)
homer_8$Nearest.PromoterID<-gsub(";.*","",homer_8$Nearest.PromoterID)
homer_8$Nearest.PromoterID<-gsub("-E.*","",homer_8$Nearest.PromoterID)
homer_8$Nearest.PromoterID<-gsub("exon_","",homer_8$Nearest.PromoterID)
colnames(homer_8)[7]<- "ID"

homer_82<-homer_8
homer_82[homer_82==""]<-NA
homer_8d<-na.omit(homer_82)#34 peaks with the motif

homer_8dui<-merge(homer_8d,dui,by="ID")
unique(homer_8dui$ID)#25
unique(homer_8dui$Peak_ID)#27

### MOTIF 9 GAMCGAGH
homer_9<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif9_GAMCGAGH.txt", header=T)
colnames(homer_9)[1]<- "Peak_ID"
homer_9<-homer_9[,c(1:6,11,22)]
homer_9$Nearest.PromoterID<-gsub("ID=","",homer_9$Nearest.PromoterID)
homer_9$Nearest.PromoterID<-gsub(";.*","",homer_9$Nearest.PromoterID)
homer_9$Nearest.PromoterID<-gsub("-E.*","",homer_9$Nearest.PromoterID)
homer_9$Nearest.PromoterID<-gsub("exon_","",homer_9$Nearest.PromoterID)
colnames(homer_9)[7]<- "ID"


homer_92<-homer_9
homer_92[homer_92==""]<-NA
homer_9d<-na.omit(homer_92)#35 peaks with the motif

homer_9dui<-merge(homer_9d,dui,by="ID")
unique(homer_9dui$ID)#24
unique(homer_9dui$Peak_ID)#28


### MOTIF 10 TTATAAYT
homer_10<- read.delim(" /Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results_DEFINITIVE/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif10_TTATAAYT.txt", header=T)
colnames(homer_10)[1]<- "Peak_ID"
homer_10<-homer_10[,c(1:6,11,22)]
homer_10$Nearest.PromoterID<-gsub("ID=","",homer_10$Nearest.PromoterID)
homer_10$Nearest.PromoterID<-gsub(";.*","",homer_10$Nearest.PromoterID)
homer_10$Nearest.PromoterID<-gsub("-E.*","",homer_10$Nearest.PromoterID)
homer_10$Nearest.PromoterID<-gsub("exon_","",homer_10$Nearest.PromoterID)
colnames(homer_10)[7]<- "ID"

homer_102<-homer_10
homer_102[homer_102==""]<-NA
homer_10d<-na.omit(homer_102)#30 peaks with the motif

homer_10dui<-merge(homer_10d,dui,by="ID")
unique(homer_10dui$ID)#22
unique(homer_10dui$Peak_ID)#25
