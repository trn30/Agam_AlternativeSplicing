#### Load all isoforms from the genome
all_iso<- read.delim("/Genomes/genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed", header=F)
colnames(all_iso)[4]<-"ID"

all_iso$gene<-gsub("-R.","",all_iso$ID)

#### extract the genes which have more than 1 isoform
a<-gsub("AGAP[0-9]{6}-RA","FALSE",all_iso$ID) ### make -RA isoforms FALSE
pos<-which(all_iso$ID == a) # we take out RB's positions and subtract 1 from his position to get his RA as well.
pos_ra<- pos-1

pos_def<-c(pos,pos_ra)
pos_def<-sort(pos_def)# I order from lowest to highest
pos_def<-unique(pos_def)# I remove duplicates
all_iso_def<-all_iso[pos_def,]#3505
## these are genes with more than one isoform.
cc<-unique(all_iso$gene)
pp<-unique(all_iso_def$gene)
write.table(all_iso_def, file = "/Genomes/genomic_data_AgamP4/Isoforms_of_multisoform_genes_release54.bed",quote = FALSE, sep="\t", row.names = T , col.names = T)


all_no_multi<-all_iso[which(!(all_iso$gene %in% all_iso_def$gene)),]
write.table(all_no_multi, file = "/Genomes/genomic_data_AgamP4/Isoforms_NO_multisoform_genes_release54.bed",quote = FALSE, sep="\t", row.names = T , col.names = T)
