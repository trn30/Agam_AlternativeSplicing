####### I will perform the chromatin study with DUi isoforms at the splicing sites.
#### 1. with all the splicign events I am going to keep only those occurring in the DUI isoforms.
all_eve<- read.delim("/Correl_by_mechanisms/Isoforms_Genome_wide_all_events_all_info.txt", header=T)
all_eve$gene<-gsub("-R.*","",all_eve$isoform_id)
colnames(all_eve)[1]<- "ID"
dui<- read.delim("/DUI/Isoforms_DUI_inf_Mg_vs_inf_SG_padj_0.05.txt", header=T) #211 dui
colnames(dui)[3]<-"ID"

sum(dui$dIF<0)## 123 en MG
sum(dui$dIF>0)## 124 en SG

## I eliminate overlapping genes
over<-read.delim("/overlapping/DUI_without_overlapping_genes_ID_definitive.bed", header=F) #180
p<-unique(over$V1)#200 genes
dui_ov<-dui[which(( dui$ID %in%  over$V1)),]
z<-unique(dui_ov$gene)#144 genes
dui_eve<-merge(all_eve,dui_ov,by="ID")
dui_eve<-dui_eve[,c(1,3:7,9,10,15,16:18,19)]

write.table((dui_eve[,c(4,2,3,7,6,5)]), file = "/Motif_analysis/DUI_Isoforms_events_without_overlapping.bed",quote = FALSE, sep="\t", row.names =F , col.names = F)
## I will check which isoforms do not have localised elements.
dui_ov[which(!(dui_ov$ID %in% dui_eve$ID)),3] #17 isoformas

### Now what I'm going to do is to add about 100 base pairs to the ends of the splicing sites to increase the area a bit and take the peripheral areas. I do it with bedtools slop.
#### I do this in the linux terminal
#bedtools slop -i /Motif_analysis/DUI_Isoforms_events_without_overlapping.bed -g /Genomes/genomic_data_AgamP4/chrom.sizes_release_54 -b 100 > /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed

### with this I will make a bedtools intersect of the diffbind peaks
#linux terminal
# bedtools intersect -a /ATAC/DiffBind/Diffbind_peaks_JL.bed -b /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wa   > /Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_without_overlapping.bed
# bedtools intersect -a /ATAC/DiffBind/Diffbind_peaks_JL.bed -b /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wb -wa > /Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_without_overlapping_b.bed


### I now load the splicing sites that match a diffbind peak.
diff_spli<- read.delim("/ATAC/DiffBind/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_b.bed", header=F)
colnames(diff_spli)[10]<-"event_ID"
diff_def<-merge(diff_spli,dui_eve,by="event_ID")
unique(dui_eve$ID)#183
unique(diff_def$ID)#56
unique(diff_def$gene)#37

#### I will do the same by separating the difbind peaks into MG and SG.
# I separate the peaks in MG and SG
diff_log<-read.delim("/ATAC/DiffBind/Differential_analysis_output_with_ID.txt",header = T)
diff_mg<-diff_log[(diff_log$Fold>0),12]
diff_sg<-diff_log[(diff_log$Fold<0),12]
### ahora hago el merge con el bed que nos interesa para hacer el intersect y lo guardo
diff_pe<-read.delim("/ATAC/DiffBind/Diffbind_peaks_JL.bed",header = F)
pe_mg<-diff_pe[(diff_pe$V4 %in% diff_mg),]
pe_sg<-diff_pe[(diff_pe$V4 %in% diff_sg),]

write.table((pe_mg), file = "/Motif_analysis/Diffbind_peaks_JL_MG.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)
write.table((pe_sg), file = "/Motif_analysis/Diffbind_peaks_JL_SG.bed",quote = FALSE, sep="\t", row.names = F , col.names = F)

### with this I will make a bedtools intersect of the diffbind peaks.
# linux terminal
# bedtools intersect -a /Motif_analysis/Diffbind_peaks_JL_MG.bed -b /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wa   > /Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG.bed
# bedtools intersect -a /Motif_analysis/Diffbind_peaks_JL_SG.bed -b /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wa   > /Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG.bed

### to see how many genes have differential peaks in splicing sites I use the -wb option.
# bedtools intersect -a /Motif_analysis/Diffbind_peaks_JL_MG.bed -b /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wb -wa > /Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed
# bedtools intersect -a /Motif_analysis/Diffbind_peaks_JL_SG.bed -b /Motif_analysis/DUI_Isoforms_events_slop_100bp_whitout_overlapping.bed -wb -wa > /Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed

###### IMPORTANT

diffb_MG<- read.delim("/Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_MG_b.bed", header=F)
colnames(diffb_MG)[10]<-"event_ID"
diffb_MG_def<-merge(diffb_MG,dui_eve,by="event_ID")
unique(dui_eve$ID)#183
unique(diffb_MG_def$ID)#5

diffb_SG<- read.delim("/Motif_analysis/merge_diffbind_peaks_with_DUI_splicing_sites_whitout_overlapping_SG_b.bed", header=F)
colnames(diffb_SG)[10]<-"event_ID"
diffb_SG_def<-merge(diffb_SG,dui_eve,by="event_ID")
unique(dui_eve$ID)#183
unique(diffb_SG_def$ID)#94


#### 94 de SG + 4 de MG = 98 DE 247 DUI HAVE DIFFBIND PEAK as 1 of the 5 peaks in MG are of the same isoform as those in SG.
