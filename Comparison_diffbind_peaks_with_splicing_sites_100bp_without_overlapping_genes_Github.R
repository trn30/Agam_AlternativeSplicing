####### I will perform the chromatin study with DUi isoforms at the splicing sites.
#### 1. with all the splicign events I am going to keep only those occurring in the DUI isoforms.
all_eve<- read.delim(" /Correl_by_mechanisms/Isoforms_Genome_wide_all_events_all_info.txt", header=T)
all_eve$gene<-gsub("-R.*","",all_eve$isoform_id)
colnames(all_eve)[1]<- "ID"
dui<- read.delim(" /Iso_usage/Isoforms_DUI_MG_vs_SG_0.05 (1).txt", header=T) #211 dui
unique(dui$gene) #116 genes
sum(dui$dIF<0)##112 en MG
sum(dui$dIF>0)## 99 en SG

## I eliminate overlapping genes
over<-read.delim(" /overlapping/DUI_without_overlapping_genes_ID_definitive.bed", header=F) #180
p<-unique(over$V1)#180 genes
dui_ov<-dui[which(( dui$ID %in%  over$V1)),]
z<-unique(dui_ov$gene)#97 genes
dui_eve<-merge(all_eve,dui_ov,by="ID")
dui_eve<-dui_eve[,c(1,3:7,9,10,15,16:18,19)]

write.table((dui_eve[,c(4,2,3,7,6,5)]), file = " /Correl_by_mechanisms/DUI_Isoforms_events_whitout_overlapping.bed",quote = FALSE, sep="\t", row.names =F , col.names = T)
## I will check which isoforms do not have localised elements.
dui_ov[which(!(dui_ov$ID %in% dui_eve$ID)),1] #15 isoformas

### Now what I'm going to do is to add about 100 base pairs to the ends of the splicing sites to increase the area a bit and take the peripheral areas. I do it with bedtools slop.
#### I do this in the linux terminal
#bedtools slop -i  /Correl_by_mechanisms/DUI_Isoforms_events_whitout_overlapping.bed -g  /genomic_data_AgamP4/chrom.sizes_release_54 -b 100 >  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed

### with this I will make a bedtools intersect of the diffbind peaks

# bedtools intersect -a  /ATAC/DiffBind/Diffbind_peaks_JL.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wa   >  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping.bed
# bedtools intersect -a  /ATAC/DiffBind/Diffbind_peaks_JL.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wb -wa >  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_b.bed
#

### I now load the splicing sites that match a diffbind peak.
diff_spli<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_b.bed", header=F)
colnames(diff_spli)[10]<-"event_ID"
diff_def<-merge(diff_spli,dui_eve,by="event_ID")
unique(dui_eve$ID)#165
unique(diff_def$ID)#98
unique(diff_def$gene.x)#55

#### I will do the same by separating the difbind peaks into MG and SG.
# I separate the peaks in MG and SG
diff_log<-read.delim(" /ATAC/DiffBind/Differential_analysis_output_with_ID.txt",header = T)
diff_mg<-diff_log[(diff_log$Fold>0),12]
diff_sg<-diff_log[(diff_log$Fold<0),12]
### now I make the merge with the bed we are interested in to make the intersect and save it.
diff_pe<-read.delim(" /ATAC/DiffBind/Diffbind_peaks_JL.bed",header = F)
pe_mg<-diff_pe[(diff_pe$V4 %in% diff_mg),]
pe_sg<-diff_pe[(diff_pe$V4 %in% diff_sg),]

write.table((pe_mg), file = " /ATAC/DiffBind/Diffbind_peaks_JL_MG.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((pe_sg), file = " /ATAC/DiffBind/Diffbind_peaks_JL_SG.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

### with this I will make a bedtools intersect of the diffbind peaks.

# bedtools intersect -a  /ATAC/DiffBind/Diffbind_peaks_JL_MG.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wa   >  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG.bed
# bedtools intersect -a  /ATAC/DiffBind/Diffbind_peaks_JL_SG.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wa   >  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG.bed

### to see how many genes have differential peaks in splicing sites I use the -wb option.
# bedtools intersect -a  /ATAC/DiffBind/Diffbind_peaks_JL_MG.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wb -wa >  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed
# bedtools intersect -a  /ATAC/DiffBind/Diffbind_peaks_JL_SG.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wb -wa >  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed

###### IMPORTANT
diffb_MG<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed", header=F)
colnames(diffb_MG)[10]<-"event_ID"
diffb_MG_def<-merge(diffb_MG,dui_eve,by="event_ID")
unique(dui_eve$ID)#165
unique(diffb_MG_def$ID)#5

diffb_SG<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed", header=F)
colnames(diffb_SG)[10]<-"event_ID"
diffb_SG_def<-merge(diffb_SG,dui_eve,by="event_ID")
unique(dui_eve$ID)#165
unique(diffb_SG_def$ID)#95


#### 95 of SG + 3 of MG = 98 OF 180 dui HAVE DIFFBIND PEAK as 2 of the 5 peaks in MG are of the same isoform as those in SG.

##########################################################################################################

### another check we are going to do is with MACS2 peaks and the bdgdiff and bdgcmp functions.

##### with bdgdiff we are going to obtain the doferential peaks and for them we are going to obtain the bdg files with the callpeaks function (in the JL files I can't find it).

# macs2 callpeak -t  /ATAC/R10-8_d7_nucfree_mapq_10_merged.bam -f BAMPE -g 273109044 --keep-dup all -n  /ATAC/MACS2_0323/R10-8_d7_macs2_call_peaks_for_bdgdiff --trackline -B -q 0.01 --nomodel --call-summits
#
# macs2 callpeak -t  /ATAC/R5-1_d14_nucfree_mapq_10_merged.bam -f BAMPE -g 273109044 --keep-dup all -n  /ATAC/MACS2_0323/R5-1_d14_macs2_call_peaks_for_bdgdiff --trackline -B -q 0.01 --nomodel --call-summits

### I now use bdgdiff (https://github.com/macs3-project/MACS/wiki/Call-differential-binding-events)

#macs2 bdgdiff --t1  /ATAC/MACS2_0323/R10-8_d7_macs2_call_peaks_for_bdgdiff_treat_pileup.bdg  --c1  /ATAC/MACS2_0323/R10-8_d7_macs2_call_peaks_for_bdgdiff_control_lambda.bdg --t2  /ATAC/MACS2_0323/R5-1_d14_macs2_call_peaks_for_bdgdiff_treat_pileup.bdg  --c2  /ATAC/MACS2_0323/R5-1_d14_macs2_call_peaks_for_bdgdiff_control_lambda.bdg --d1 60994 --d2 64309 --outdir  /ATAC/MACS2_0323/out_macs2_bdgdiff --o-prefix Inf_7d_vs_d14


#### I use bdgcmp to have run track with the difference between the two macs2 tracks. (https://github.com/macs3-project/MACS/wiki/Build-Signal-Track)
# macs2 bdgcmp -t   /ATAC/MACS2_0323/R10-8_d7_macs2_call_peaks_for_bdgdiff_treat_pileup.bdg  -c  /ATAC/MACS2_0323/R5-1_d14_macs2_call_peaks_for_bdgdiff_treat_pileup.bdg -o  /ATAC/MACS2_0323/Inf_7d_vs_d14_substract_TRACK.bdg -m subtract -p 0.00001


### once I have the differential peaks (cond1 for MG and cond 2 for SG) I will see how many of them merge with the splice sites.


# bedtools intersect -a  /ATAC/MACS2_0323/out_macs2_bdgdiff/Inf_7d_vs_d14_c3.0_cond1.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed  -wa   >  /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_MG.bed
# bedtools intersect -a  /ATAC/MACS2_0323/out_macs2_bdgdiff/Inf_7d_vs_d14_c3.0_cond2.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed  -wa   >  /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_SG.bed

### to see how many genes have differential peaks in splicing sites I use the -wb option.
# bedtools intersect -a  /ATAC/MACS2_0323/out_macs2_bdgdiff/Inf_7d_vs_d14_c3.0_cond1.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed  -wb   >  /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed
# bedtools intersect -a  /ATAC/MACS2_0323/out_macs2_bdgdiff/Inf_7d_vs_d14_c3.0_cond2.bed -b  /Correl_by_mechanisms/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed  -wb  >  /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed


### I now load the splicing sites that match a differential macs2 peak.
madif_MG<- read.delim(" /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed", header=F)
colnames(madif_MG)[9]<-"event_ID"
madif_MG_def<-merge(madif_MG,dui_eve,by="event_ID")
unique(dui_eve$ID)#165
unique(madif_MG_def$ID)#1

madif_SG<- read.delim(" /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed", header=F)
colnames(madif_SG)[9]<-"event_ID"
madif_SG_def<-merge(madif_SG,dui_eve,by="event_ID")
unique(dui_eve$ID)#165
unique(madif_SG_def$ID)#64

#### 65 out of 180 DUI HAVE MACS2 DIFFERENTIAL PEAK

#### To compare the two methods what we are going to do is to make an intersect between the peaks of both methods. The intersect is done with 25% of the peaks must coincide.

# bedtools intersect -a  /ATAC/MACS2_0323/out_macs2_bdgdiff/Inf_7d_vs_d14_c3.0_cond1.bed -b  /ATAC/DiffBind/Diffbind_peaks_JL_MG.bed -f 0.25 -wa   >  /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_MG.bed
# bedtools intersect -a  /ATAC/MACS2_0323/out_macs2_bdgdiff/Inf_7d_vs_d14_c3.0_cond2.bed -b  /ATAC/DiffBind/Diffbind_peaks_JL_SG.bed -f 0.25 -wa   >  /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_SG.bed

coin_MG<- read.delim(" /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_MG.bed", header=F) #681
coin_SG<- read.delim(" /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_SG.bed", header=F) #9934

### now I do the same for the differential peaks of both methods which also coincide with a splicing site.
#
# bedtools intersect -a  /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_MG.bed -b  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG.bed -f 0.25 -wa -wb   >  /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_in_splicing_sites_whitout_overlapping_MG.bed
# bedtools intersect -a  /ATAC/DiffBind/merge_macsdiff_peaks_with_DUI_splicing_sites_whitout_overlapping_SG.bed -b  /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG.bed -f 0.25 -wa  -wb  >  /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_in_splicing_sites_whitout_overlapping_SG.bed

coin_ss_MG<- read.delim(" /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_in_splicing_sites_whitout_overlapping_MG.bed", header=F) #1
coin_ss_SG<- read.delim(" /ATAC/DiffBind/coincident_peaks_macsdiff_diffbind_in_splicing_sites_whitout_overlapping_SG.bed", header=F) #170
coin_ss_SG<-unique(coin_ss_SG)#51


##########################################################################################################################################################################################
##########################################################################################################################################################################################

#### I will now see if the prediction created for Es is correct.
#### the prediction is that at the exon skipping splicing sites there will be some with accessibility proportional to the use of that isoform due to the presence of spliceosomal proteins or the opposite will happen due to the binding of repressor proteins.
#### 1. I will take DUI isoforms that have diffbind peaks at splicing sites.
#### 2. I'm going to see where these MG or SG isoforms are most used.
#### 3. I'm going to look at the peak diffbind where it is more accessible.


diffb_MG<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed", header=F)
diffb_MG<-diffb_MG[,c(1:4,7:12)]
colnames(diffb_MG)<-c("Chr_peak","Start_peak","End_peak","Region_ID","Chr_event","Start_event","End_event","event_ID","Score","Strand")
diffb_SG<- read.delim(" /ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed", header=F)
diffb_SG<-diffb_SG[,c(1:4,7:12)]
colnames(diffb_SG)<-c("Chr_peak","Start_peak","End_peak","Region_ID","Chr_event","Start_event","End_event","event_ID","Score","Strand")

### to see where they are most used I combine the data sets with the DUI data sets.
diffb_MG_d<-merge(dui_eve,diffb_MG,by="event_ID")#11
diffb_SG_d<-merge(dui_eve,diffb_SG,by="event_ID")#213

diffb_MG_d<-unique(diffb_MG_d)#7
diffb_MG_d<-diffb_MG_d[,c(1:7,9,10,15:17,19,20)]
diffb_MG_d$more.used<-ifelse(diffb_MG_d$dIF > 0 , "SG","MG")

diffb_SG_d<-unique(diffb_SG_d)#141
diffb_SG_d<-diffb_SG_d[,c(1:7,9,10,15:17,19,20)]
diffb_SG_d$more.used<-ifelse(diffb_SG_d$dIF > 0 , "SG","MG")

### I add the information of the fold change of the diffbind peaks.
diff_log<-read.delim(" /ATAC/DiffBind/Differential_analysis_output_with_ID.txt",header = T)

diffb_MG_d2<-merge(diffb_MG_d,diff_log,by="Region_ID")
diffb_MG_d2<-diffb_MG_d2[,c(1:15,24)]
diffb_MG_d2$more.access<-ifelse(diffb_MG_d2$Fold > 0 , "MG","SG")
diffb_MG_d2$coincidence<-ifelse(diffb_MG_d2$more.access == diffb_MG_d2$more.used, "YES","NO")
sum(diffb_MG_d2$coincidence == "YES")#1
sum(diffb_MG_d2$coincidence == "NO")#6

diffb_SG_d2<-merge(diffb_SG_d,diff_log,by="Region_ID")
diffb_SG_d2<-diffb_SG_d2[,c(1:15,24)]
diffb_SG_d2$more.access<-ifelse(diffb_SG_d2$Fold > 0 , "MG","SG")
diffb_SG_d2$coincidence<-ifelse(diffb_SG_d2$more.access == diffb_SG_d2$more.used, "YES","NO")
sum(diffb_SG_d2$coincidence == "YES")#68
sum(diffb_SG_d2$coincidence == "NO")#73

### filter for the ES
diffb_MG_d2$type_event<-gsub("_.*","",diffb_MG_d2$event_ID)
diff_MG_ES<-diffb_MG_d2[diffb_MG_d2$type_event == "ES",] #4
sum(diff_MG_ES$coincidence == "YES")#0
sum(diff_MG_ES$coincidence == "NO")#4

diffb_SG_d2$type_event<-gsub("_.*","",diffb_SG_d2$event_ID)
diff_SG_ES<-diffb_SG_d2[diffb_SG_d2$type_event == "ES",] #20
sum(diff_SG_ES$coincidence == "YES")#7
sum(diff_SG_ES$coincidence == "NO")#13
### I will now join the MG and SG matches to make the motif analysis with all peaks, removing duplicates manually.
diff_ES_all<-rbind(diff_SG_ES,diff_MG_ES)
write.table((diff_ES_all), file = " /ATAC/DiffBind/Diffbind_peaks_coincident_with_splicing_sites_ES_and_DUI.txt",quote = FALSE, sep="\t", row.names = F , col.names = T)

##### I do the same for ATSS events where I expect + expression, + accessibility.

diff_MG_ATSS<-diffb_MG_d2[diffb_MG_d2$type_event == "ATSS",] #3
sum(diff_MG_ATSS$coincidence == "YES")#1
sum(diff_MG_ATSS$coincidence == "NO")#2

diff_SG_ATSS<-diffb_SG_d2[diffb_SG_d2$type_event == "ATSS",] #83
sum(diff_SG_ATSS$coincidence == "YES")#40
sum(diff_SG_ATSS$coincidence == "NO")#43



##### I do the same for IR events where I expect + accessibility == a + intron retention.

diff_MG_IR<-diffb_MG_d2[diffb_MG_d2$type_event == "IR",] #0
sum(diff_MG_IR$coincidence == "YES")#
sum(diff_MG_IR$coincidence == "NO")#

diff_SG_IR<-diffb_SG_d2[diffb_SG_d2$type_event == "IR",] #1
sum(diff_SG_IR$coincidence == "YES")
sum(diff_SG_IR$coincidence == "NO")#1
