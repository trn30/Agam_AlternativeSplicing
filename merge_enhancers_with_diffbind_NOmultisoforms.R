####################### I will try to see how many enhacers match diffbind peaks and NOmultisoforms.

### now with bedtools I'm going to see how many diffbind peaks match with demg and DUI enhancers.
# Linux Terminal

# bedtools intersect -a /Enhancers/merge_NOmultisoforms_enhancers_a.bed -b /ATAC/DiffBind/Diffbind_peaks_JL.bed -wa -wb -f 0.51 > /Enhancers/merge_NOmultisoforms_enhancers_a_diffbind_peaks.bed
#
# bedtools intersect -a /Enhancers/merge_NOmultisoforms_enhancers_b.bed -b /ATAC/DiffBind/Diffbind_peaks_JL.bed -wa -wb -f 0.51 > /Enhancers/merge_NOmultisoforms_enhancers_b_diffbind_peaks.bed


mer_demg_b<- read.delim("/Enhancers/merge_NOmultisoforms_enhancers_b_diffbind_peaks.bed", header=F)
unique(mer_demg_b$V4)#4 ENHANCERS
unique(mer_demg_b$V6)#57 GENES

mer_demg_a<- read.delim("/Enhancers/merge_NOmultisoforms_enhancers_a_diffbind_peaks.bed", header=F)
p<- unique(mer_demg_a$V4)#121 ENHANCERS
unique(mer_demg_a$V6)#169 GENES
