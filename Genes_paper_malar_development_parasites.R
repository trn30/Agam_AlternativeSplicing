#GENES PAPER MALAR

dui<- read.delim("Isoforms_DUI_MG_vs_SG_0.05 (1).txt", header=T) #211 dui
deg<- read.delim("Isoforms_DEI_MG_vs_SG_0.05.txt", header=T) #4933 DEG
deg$isoform_ID<-rownames(deg)
deg$gene<-gsub("-R.*","",deg$isoform_ID)
malar<- read.delim("Genes_paper_malar_2013.txt", header=F) #4933 DEG


malar[which(malar$V2 %in% dui$gene),] #1
malar[which(malar$V2 %in% deg$gene),] #10
