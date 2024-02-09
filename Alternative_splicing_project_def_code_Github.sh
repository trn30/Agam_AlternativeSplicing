### mapping of samples and counts with salmon are  made following  Expresion diferencial.odt and Alternative_splicing.odt
#1- analysis of the quality of reads with Fastqc
#2- trimming in case of adapters with bbmap(specifically with bbduk)

bduk.sh -Xmx1g in= $PWD/329_1.fastq in=$PWD/329_2.fastq out1=329t_1.fastq out2=329t_2.fastq qtrim=l trimq=10

##NOT NECESSARY IF THEY DON'T HAVE ADAPTERS
#3- In case the GC content is high, the sortmerna programme will be used.
#I give you the database of eukaryotes and bacteria as well as the short and long rRNAs of the species in question (always from the silva rrna database).


sortmerna --ref  /rRNA_databases/silva-bac-16s-id90.fasta --ref  /rRNA_databases/silva-bac-23s-id98.fasta --ref  /rRNA_databases/silva-euk-28s-id98.fasta --ref  /rRNA_databases/silva-euk-18s-id95.fasta --ref  /rRNA_databases/Agambiae_silva_lsu_r132_28_170420_nogap.fasta --ref  /rRNA_databases/Agambiae_silva_ssu_r138_28_170420_nogap.fasta --ref  /rRNA_databases/Pfalciparum_silva_LSU_290818_nogap.fasta --ref  /rRNA_databases/Pfalciparum_silva_SSU_290818_nogap.fasta -reads  /subset/329_1.fq -reads  /subset/329_2.fq -fastx -paired_in -v -num_alignments 1 -other  /non-rRNA/non_rRNA_329


# after doing the sortmerna I separate the files again into 1 and 2

reformat.sh in= /non-rRNA/non_rRNA_327.fastq.fastq out1=327s_1.fq out2=327s_2.fq t=$cores

## 4- Mapping with SALMON
# mapped mode
# In this mode, the programme creates its own alignments and then performs the count. To use it, the first thing we have to do is to create an index of the genome we are interested in.

salmon index -t transcripts.fa -i transcripts_index  -k 31

# The number of k-mers should be adapted according to the size of the reads to make it more sensitive. 31 is fine for reads of 75 bp or more.
salmon quant -i transcripts_index_AgamP4_1.3.0 -l A -1  /RNA-seq/raw_data/2567-VGC-0011_Salivary_glands_Infected_14d_S3_1.fastq.gz -2  /RNA-seq/raw_data/2567-VGC-0011_Salivary_glands_Infected_14d_S3_2.fastq.gz  --threads 4 -o  /RNA-seq/out_salmon/2567-VGC-0011_Salivary_glands_Infected_14d_S3  --validateMappings --go 4 --mismatchSeedSkip 5  # OPTIMISED


#5- deeptools for PCA representation

##plotcorrelation- how similar the genes or samples are to each other

plotCorrelation -in readCounts.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers  -o heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab

##plotPCA- is a way to reduce the number of variables, leaving the ones that have the most influence on the results, so that you keep as much information as possible but with fewer variables.

plotPCA -in readCounts.npz -o PCA_readCounts.png -T "PCA of read counts"



#6- Quantification of differential expression (DEseq2)
##### obtain Differential expressed isoforms of multisoform genes (DEMG)
### the R script is  DESeq2_analysis_MG_vs_SG.R Ctrl_vd_Infec_MG_analysis Ctrl_vd_Infec_SG_analysis

#7- Quantification of Usage (IsoformSwitchAnalyzeR)
### obtain isoforms differentially used
## the code is in  Isoform_switch_MG_vs_SG.R and  Isoform_switch_Ctrol_vs_Infec.R


### Obtain volcano plot DEMG vs DUI
# script_volcanos_percent_DUI_vs_DEG.R

###  obtain use of events by group
# AS_mechanism_by_conditions.R

### comparison of our DEMG and DUI isoforms with genes present in the paper https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-12-216.
###Genes_paper_malar_development_parasites.R

###### correlation ATAC vs. RNA

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

### to quantifay the accessibiliy in the region of interest we use Bedtools
########### ISOFORMAS #############
##### accessibility in the promoter ######
### I DO READ COUNTS ON THE ENTIRE PROMOTER

bedtools intersect -c -a Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  > R5-I3d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  > R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  > R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam > R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed


##### accessibility in the isoform body ######

### I DO READ COUNTS ON THE ENTIRE BODY

bedtools intersect -c -a Isoforms_AgamP4_release_54.bed -b R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  > R5-I3d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a Isoforms_AgamP4_release_54.bed -b R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  > R1-I1d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a Isoforms_AgamP4_release_54.bed -b R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  > R8-I3d7_MG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a Isoforms_AgamP4_release_54.bed -b R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam > R10-I1d7_MG_counts_body_SG_vs_MG_2.bed


### the information to check overlapping zones is in check_overlaping_zones.R
#### first we eliminate overlapping gene with isoforms Promoters and calculate the correlation of DEMG group whitout overlapping genes
### THE COUNTS OF THE ENRICHMENT ARE IN files_correlations_DEG_DUI.sh or above this information
#### inrformation about the correlation is in : correlation_DEMG_without_overlap

#### divide the expression and accessibility by levels to perform the correlation and make Ngsplot
# DEMG_correlation_by_levels_without_overlap.R


 ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_accesible_SG.txt -O expression_genebody_SG
 ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_expression_SG.txt -O accesibility_tss_SG

#### We searched for diffbind peaks that matched splicing sites of our DUI isoforms. The process performed can be found in : Comparison_diffbind_peaks_with_splicing_sites_100bp_without_overlapping_genes.R
#### With the diffbind peaks coincident with splicing sites in the DUI we will do a motif analysis. We do the same with these peaks that are also located in splicing events of the exon skipping type.
# preparation and execution Motif_analysis_diffbind_splicing_sites_DUI.R
### annotation

#### KNOWN MOTIFS

### MOTIF 1 CATCMCTA

annotatePeaks.pl Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff VectorBase-54_AgambiaePEST.gff -m /peaks_with_DUI_ans_splicing_sites_without_overlapping_results/knownResults/known1.motif > /peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_KNWON_ aks_motif1_CATCMCTA.txt # 1e-3

### MOTIF 2 GGYCATAAAW

annotatePeaks.pl Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff VectorBase-54_AgambiaePEST.gff -m /peaks_with_DUI_ans_splicing_sites_without_overlapping_results/knownResults/known2.motif > /peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_KNOWN_annotatePeaks_motif2_GGYCATAAAW.txt # 1e-2

#### HOMER MOTIFS

### MOTIF 1 CACCWATT

annotatePeaks.pl Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif1.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif1_CACCWATT.txt # 1e-34

### MOTIF 2 CTTGCAAGGT

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif2.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif2_CTTGCAAGGT.txt # 1e-31

### MOTIF 3 ACARAGAAAT

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif3.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif3_ACARAGAAAT.txt # 1e-28

### MOTIF 4 ACACAAAG

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif4.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif4_ACACAAAG.txt # 1e-27

### MOTIF 5 GAGTCTTTAT

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif5.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif5_GAGTCTTTAT.txt # 1e-27

### MOTIF 6 GAAAACAGTT

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif6.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif6_GAAAACAGTT.txt # 1e-27
### MOTIF 7 AAACAGTAAAAG

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif7.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif7_AAACAGTAAAAG.txt # 1e-27
### MOTIF 8 CAAGAAGA

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa  -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif8.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif8_CAAGAAGA.txt # 1e-27
### MOTIF 9 GAMCGAGH

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif9.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif9_GAMCGAGH.txt # 1e-27

### MOTIF 10 TTATAAYT

annotatePeaks.pl    Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed   Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff  VectorBase-54_AgambiaePEST.gff -m   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif10.motif >   peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif10_TTATAAYT.txt # 1e-27

### to study if any type of event is related to specific accessibility of a region, we use the script AS_mechanism_by_conditions_comparing_peak_location.R
## finally we extrac enhancers known by other authors and by our group taking drosophila orthologs and we look how many of them are present and active in our groups.
# merge_enhancers_with_DEMG_DUI.R
