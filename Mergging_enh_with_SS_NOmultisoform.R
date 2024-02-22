all_eve<- read.delim("/Correl_by_mechanisms/Isoforms_Genome_wide_all_events_all_info.txt", header=T)

all_eve_d<- all_eve[,c(5,3,4,9)]
all_eve_d$strand<-"+"
all_eve_d$score<-0

write.table((all_eve_d), file = "/Enhancers/Isoforms_Genome_wide_all_events_all_info.bed" ,quote = FALSE, sep="\t", row.names = F , col.names = F)

#### now I do a bedtools intersect of the SS with the enhancers and see which ones are inside the MO multisoform genes
### Linux terminal
bedtools intersect -wa -wb -f 0.25 -a /Enhancers/Isoforms_Genome_wide_all_events_all_info.bed -b /Enhancers/merge_NOmultisoforms_enhancers_b.bed > /Enhancers/merge_NOmultisoforms_enhancers_b_SS.bed

bedtools intersect -wa -wb -f 0.25 -a /Enhancers/Isoforms_Genome_wide_all_events_all_info.bed -b /Enhancers/merge_NOmultisoforms_enhancers_a.bed > /Enhancers/merge_NOmultisoforms_enhancers_a_SS.bed

a<- read.delim("/Enhancers/merge_NOmultisoforms_enhancers_a_SS.bed", header=F)
b<- read.delim("/Enhancers/merge_NOmultisoforms_enhancers_b_SS.bed", header=F)

a2<-unique(a$V10) # 21
b2<-unique(b$V10) # 2


def<- c(a2,b2)
unique(def)#23
