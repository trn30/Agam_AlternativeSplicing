####################################################################################################################################################################################################
################################################################################### ALTERNATIVE SPLICING ###########################################################################################
####################################################################################################################################################################################################

### mapping of samples and counts with salmon are  made following  Expresion diferencial.odt and Alternative_splicing.odt
#1- analysis of the quality of reads with Fastqc
#2- trimming in case of adapters with bbmap(specifically with bbduk)

bduk.sh -Xmx1g in= 329_1.fastq in=329_2.fastq out1=329t_1.fastq out2=329t_2.fastq qtrim=l trimq=10

##NOT NECESSARY IF THEY DON'T HAVE ADAPTERS
#3- In case the GC content is high, the sortmerna programme will be used.
#I give you the database of eukaryotes and bacteria as well as the short and long rRNAs of the species in question (always from the silva rrna database).

sortmerna --ref  /rRNA_databases/silva-bac-16s-id90.fasta --ref  /rRNA_databases/silva-bac-23s-id98.fasta --ref  /rRNA_databases/silva-euk-28s-id98.fasta --ref  /rRNA_databases/silva-euk-18s-id95.fasta --ref  /rRNA_databases/Agambiae_silva_lsu_r132_28_170420_nogap.fasta --ref  /rRNA_databases/Agambiae_silva_ssu_r138_28_170420_nogap.fasta --ref  /rRNA_databases/Pfalciparum_silva_LSU_290818_nogap.fasta --ref  /rRNA_databases/Pfalciparum_silva_SSU_290818_nogap.fasta -reads  /subset/329_1.fq -reads  /subset/329_2.fq -fastx -paired_in -v -num_alignments 1 -other  /non-rRNA/non_rRNA_329

# after doing the sortmerna I separate the files again into 1 and 2

reformat.sh in= /non-rRNA/non_rRNA_327.fastq.fastq out1=327s_1.fq out2=327s_2.fq t=$cores

## 4- The output of sortmetna is Mapped with SALMON for the 3 comparisons and also with STAR to use IGV
## The code is in the script Alternative_splicing_STAR_SALMON_QC.sh
#5- After the mapping we do a QC of the bam files with quialimap, also we use MultiQC to summarise all the quality controls performed
## The code is in the script Alternative_splicing_STAR_SALMON_QC.sh

#6- extraction of multiform genes from the anopheles gambiae genome (vectobase release 54)
## R script extract_isoforms_from_multisoform_genes.R

#7- Quantification of differential expression (DEseq2)
##### obtain Differential expressed isoforms of multisoform genes (DEMG)
### results of the DEMG in the comparions of Ctrol vs. Inf MG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_ctrol.R
### results of the DEMG in the comparions of Ctrol vs. Inf SG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_ctrol.R
### results of the DEMG in the comparions of Inf MG vs. Inf SG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_Inf.R


#8- Quantification of Usage (IsoformSwitchAnalyzeR)
### results of the DUI in the comparions of Ctrol vs. Inf MG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_ctrol.R
### results of the DUI in the comparions of Ctrol vs. Inf SG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_ctrol.R
### results of the DUI in the comparions of Inf MG vs. Inf SG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_Inf.R

#9- Comparison of our DEMG and DUI isoforms with genes present in the paper https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-12-216.
### results of the matches in the comparions of Ctrol vs. Inf MG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_ctrol.R
### results of the matches in the comparions of Ctrol vs. Inf SG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_ctrol.R
### results of the matches in the comparions of Inf MG vs. Inf SG are in the R script Repeticion_DESEQ_Isoformswitch_Inf_vs_Inf.R

#10- Obtain volcano plot DEMG vs DUI
# The results are in the R script script_volcanos_percent_DUI_vs_DEMG.R

#11- Obtain use of events in each comparison
# The results are in the R script AS_mechanism_by_conditions.R

#12- Correlation between accessibility and expression in comparison of Inf MG vs. Inf SG

##### adaptation of bam file to the Vectorbase relesa 54

for f in $( ls | egrep -E 'out.bam$' ); do  samtools sort -@ 20 "${f}" -o "${f}_sort.bam"; done

for f in $( ls | egrep -E '_sort.bam$' ); do samtools view -H "${f}"  | sed -e 's/2L/AgamP4_2L/g' -e 's/2R/AgamP4_2R/g' -e 's/3L/AgamP4_3L/g' -e 's/3R/AgamP4_3R/g' -e 's/Mt/AgamP4_Mt/g' -e 's/UNKN/AgamP4_UNKN/g' -e 's/X/AgamP4_X/g' -e  's/Y_unplaced/AgamP4_Y_unplaced/g' | samtools reheader - "${f}"  > "${f}_correct_chr.bam" ;done

for f in $( ls | egrep -E '_chr.bam$' ); do  samtools index "${f}" ;done

for f in $( ls | egrep -E '_merged.bam$' ); do samtools flagstat "${f}" > "${f}_flagstat.txt" ; done

#### count the number of reads

samtools flagstat R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam > R5-I3d14_nucfree_mapq_10.bam.cram.bam_chr.bam_flagstat.txt

samtools flagstat R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam > R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam_flagstat.txt

samtools flagstat R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam > R1-I1d14_nucfree_mapq_10.bam.cram.bam_chr.bam_flagstat.txt

samtools flagstat R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam > R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam_flagstat.txt

## Rna-seq del NAR21

samtools flagstat EGD-1-I1d7_S1_L001_R1_001_mapped.bam > EGD-1-I1d7_S1_L001_R1_001_mapped.bam_flagstat.txt

samtools flagstat EGD-2-I1d14_S5_L001_R1_001_mapped.bam > EGD-2-I1d14_S5_L001_R1_001_mapped.bam_flagstat.txt

samtools flagstat EGD-7-I3d7_S3_L001_R1_001_mapped.bam > EGD-7-I3d7_S3_L001_R1_001_mapped.bam_flagstat.txt

samtools flagstat EGD-8-I3d14_S4_L001_R1_001_mapped.bam > EGD-8-I3d14_S4_L001_R1_001_mapped.bam_flagstat.txt

### NORMALIZATION TO REPRESENT THE SAMPLES IN IGV
### Midguts inf vs. control
#RNA
samtools view -b -h -s 0.3650883 Midgut_Infected_7d_M2Aligned.out.bam_sort.bam_correct_chr.bam > Midgut_Infected_7d_M2Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam
samtools view -b -h -s 0.313299 Midgut_Infected_7d_M3Aligned.out.bam_sort.bam_correct_chr.bam > Midgut_Infected_7d_M3Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam
samtools view -b -h -s 0.8415303 Midgut_Non-Infected_7d_M3Aligned.out.bam_sort.bam_correct_chr.bam > Midgut_Non-Infected_7d_M3Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam
samtools view -b -h -s 1 Midgut_Non-Infected_7d_M2Aligned.out.bam_sort.bam_correct_chr.bam > Midgut_Non-Infected_7d_M2Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam

####salivary Glands inf vs. control
###RNA
samtools view -b -h -s 1 Salivary_glands_Infected_14d_M2Aligned.out.bam_sort.bam_correct_chr.bam > Salivary_glands_Infected_14d_M2Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam
samtools view -b -h -s 0.9707 Salivary_glands_Infected_14d_M3Aligned.out.bam_sort.bam_correct_chr.bam > Salivary_glands_Infected_14d_M3Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam
samtools view -b -h -s 0.7407 Salivary_glands_Non-Infected_14d_M1Aligned.out.bam_sort.bam_correct_chr.bam > Salivary_glands_Non-Infected_14d_M1Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam
samtools view -b -h -s 0.8357 Salivary_glands_Non-Infected_14d_M2Aligned.out.bam_sort.bam_correct_chr.bam > Salivary_glands_Non-Infected_14d_M2Aligned.out.bam_sort.bam_correct_chr.bam_norm.bam


####Midguts vs. Salivary Glands
###RNA
samtools view -b -h -s 0.8589  EGD-2-I1d14_S5_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam > EGD-2-I1d14_S5_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam_norm.bam
samtools view -b -h -s 0.6639 EGD-1-I1d7_S1_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam > EGD-1-I1d7_S1_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam_norm.bam
samtools view -b -h -s 1 EGD-7-I3d7_S3_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam > EGD-7-I3d7_S3_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam_norm.bam
samtools view -b -h -s 0.419  EGD-8-I3d14_S4_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam > EGD-8-I3d14_S4_L001_R1_001_mapped.bam_correct_chr.bam_sorted.bam_norm.bam

###ATAC
samtools view -b -h -s 0.6841  R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam > R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam_norm.bam
samtools view -b -h -s 0.9671 R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam > R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam_norm.bam
samtools view -b -h -s 0.8810 R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam > R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam_norm.bam
samtools view -b -h -s 0.6909 R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  > R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam_norm.bam

### To quantify the Accessibility in the region of interest we use Bedtools
### The code is in the script files_correlation_DEMG_DUI.sh

### To perform the correlation we eliminate the isoforms DEMG corresponding to genes that overlap with other genes
### The information to check overlapping zones is in check_overlaping_zones.R
### first we eliminate overlapping gene with other genes and calculate the correlation of DEMG group without overlapping genes
# The results are in the R script correlation_DEMG_without_overlap.R

### Also we divide the DEMG isoforms into levels according to their expression and see if they correlate with accessibility. We expect that the higher the expression, the higher the accessibility.
## to represent the levels of expression in the TSS we use NGSPLOT
# The results are in the R script DEMG_correlation_by_levels_without_overlap_expression.R


#13- See of there is a relation between DiffBind peaks and our DUI isoforms
### We searched for diffbind peaks that matched splicing sites of our DUI isoforms.
###The process performed can be found in : Comparison_diffbind_peaks_with_splicing_sites_100bp_without_overlapping_genes
### With the diffbind peaks coincident with splicing sites in the DUI we will do a motif analysis and we annotate the resulting motifs
# The results are in the R script Motif_analysis_diffbind_splicing_sites_DUI.R


#14- Study the relation between our DEMG and DUI isoforms and enhancers
#### We use enhancers known by other authors and by our group taking drosophila orthologs from the paper https://academic.oup.com/nargab/article/3/1/lqaa113/6104947#246280024 (NAR2021) and we look how many of them are present and active in our groups.
#### We look into how many enhancers are annotated to our DEMG and DUI groups
#The results are in the R script merge_enhancers_with_DEMG_DUI.R

#### How many enhancers matches with DiffBind peaks an also with our isoform groups
#The results are in the R script merge_enhancers_with_Diffbind_DEMG_DUI.R

#### Which enhancers overlap with the Splicing Sites of our isoforms DEMG and DUI
#The results are in the R script mergging_enh_with_SS.R

# Also we do the same process with the NO multisoform genes, to compare the prevalence of coincidences

#### We look into how many enhancers are annotated to NO multisoform genes
#The results are in the R script merge_enhancers_with_DEMG_DUI.R

#### How many enhancers matches with DiffBind peaks an also with NO multisoform genes
#The results are in the R script merge_enhancers_with_Diffbind_DEMG_DUI.R

#### Which enhancers overlap with the Splicing Sites of NO multisoform genes
#The results are in the R script mergging_enh_with_SS.R
