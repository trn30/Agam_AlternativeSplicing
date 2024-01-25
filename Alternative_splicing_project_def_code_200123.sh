### mapping of samples and counts with salmos are  made following  /data/amd_HDD_4TB/Barbara/Alternative_splicing/Expresion diferencial.odt and Alternative_splicing.odt
#1- analysis of the quality of reads with Fastqc
#2- trimming in case of adapters with bbmap(specifically with bbduk)

bduk.sh -Xmx1g in= $PWD/329_1.fastq in=$PWD/329_2.fastq out1=329t_1.fastq out2=329t_2.fastq qtrim=l trimq=10

##NOT NECESSARY IF THEY DON'T HAVE ADAPTERS
#3- In case the GC content is high, the sortmerna programme will be used.
#I give you the database of eukaryotes and bacteria as well as the short and long rRNAs of the species in question (always from the silva rrna database).


sortmerna --ref /home/barbara/Programs/sortmerna/rRNA_databases/silva-bac-16s-id90.fasta --ref /home/barbara/Programs/sortmerna/rRNA_databases/silva-bac-23s-id98.fasta --ref /home/barbara/Programs/sortmerna/rRNA_databases/silva-euk-28s-id98.fasta --ref /home/barbara/Programs/sortmerna/rRNA_databases/silva-euk-18s-id95.fasta --ref /home/barbara/Programs/sortmerna/rRNA_databases/Agambiae_silva_lsu_r132_28_170420_nogap.fasta --ref /home/barbara/Programs/sortmerna/rRNA_databases/Agambiae_silva_ssu_r138_28_170420_nogap.fasta --ref /home/barbara/sortmerna/rRNA_databases/Pfalciparum_silva_LSU_290818_nogap.fasta --ref /home/barbara/sortmerna/rRNA_databases/Pfalciparum_silva_SSU_290818_nogap.fasta -reads /home/barbara/prueba/subset/329_1.fq -reads /home/barbara/prueba/subset/329_2.fq -fastx -paired_in -v -num_alignments 1 -other /home/barbara/prueba/non-rRNA/non_rRNA_329


# after doing the sortmerna I separate the files again into 1 and 2

reformat.sh in=/home/barbara/prueba/non-rRNA/non_rRNA_327.fastq.fastq out1=327s_1.fq out2=327s_2.fq t=$cores

## 4- Mapping with SALMON
# mapped mode
# In this mode, the programme creates its own alignments and then performs the count. To use it, the first thing we have to do is to create an index of the genome we are interested in.

salmon index -t transcripts.fa -i transcripts_index  -k 31

# The number of k-mers should be adapted according to the size of the reads to make it more sensitive. 31 is fine for reads of 75 bp or more.
salmon quant -i transcripts_index_AgamP4_1.3.0 -l A -1 /media/i7_station/i7_HDD_6TB/FMP_2020_Barbara/Training/RNA-seq/raw_data/2567-VGC-0011_Salivary_glands_Infected_14d_S3_1.fastq.gz -2 /media/i7_station/i7_HDD_6TB/FMP_2020_Barbara/Training/RNA-seq/raw_data/2567-VGC-0011_Salivary_glands_Infected_14d_S3_2.fastq.gz  --threads 4 -o /media/i7_station/i7_HDD_6TB/FMP_2020_Barbara/Training/RNA-seq/out_salmon/2567-VGC-0011_Salivary_glands_Infected_14d_S3  --validateMappings --go 4 --mismatchSeedSkip 5  # OPTIMISED


#5- deeptools for PCA representation

##plotcorrelation- how similar the genes or samples are to each other

    plotCorrelation -in readCounts.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers  -o heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab

##plotPCA- is a way to reduce the number of variables, leaving the ones that have the most influence on the results, so that you keep as much information as possible but with fewer variables.

plotPCA -in readCounts.npz -o PCA_readCounts.png -T "PCA of read counts"



#6- Quantification of differential expression (DEseq2)
##### obtain Differential expressed isoforms of multisoform genes (DEMG)
### the R script is /data/amd_HDD_4TB/Barbara/Alternative_splicing/DESeq2/DESeq2_analysis_MG_vs_SG_070422.R Ctrl_vd_Infec_MG_analysis_210722 Ctrl_vd_Infec_SG_analysis_210722

#7- Quantification of Usage (IsoformSwitchAnalyzeR)
### obtain isoforms differentially used
## the code is in /data/amd_HDD_4TB/Barbara/Alternative_splicing/Iso_usage/Isoform_switch_MG_vs_SG.R and /data/amd_HDD_4TB/Barbara/Alternative_splicing/Iso_usage/Isoform_switch_Ctrol_vs_Infec.R


### Obtain volcano plot DEMG vs DUI
# /data/amd_HDD_4TB/Barbara/Alternative_splicing/script_volcanos_percent_DUI_vs_DEG_def_210722.R

###  obtain use of events by roup ad if the changes are statistically relevants
# /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/AS_mechanism_by_conditions.R AS_mechanism_by_conditions_CHI_test.R

### comparison of our DEMG and DUI isoforms with genes present in the paper https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-12-216.
### /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Genes_paper_malar_development_parasites.R

###### correlation ATAC vs. RNA

##### adaptation of bam file to the Vectorbase relesa 54
/out_star_tfm

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
###midguts inf vs. control
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

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R5-I3d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed


##### accessibility in the isoform body ######

### I DO READ COUNTS ON THE ENTIRE BODY

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R5-I3d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R1-I1d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R8-I3d7_MG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam > /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/counts_promoters_all_isoforms/R10-I1d7_MG_counts_body_SG_vs_MG_2.bed


### the information to check overlapping zones is in /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/overlapping/check_overlaping_zones.R
#### first we eliminate overlapping gene with isoforms Promoters and calculate the correlation of DEMG group whitout overlapping genes
# /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/overlapping/correlation_DEMG_without_overlap_definitive_120423.R


#### divide the espression and accessibility by levels to perform the correlation and make ngsplot
### THE COUNTS OF THE ENRICHMENT ARE IN files_correlations_DEG_DUI.sh or above this information
# THE CODE RELATED WITH THIS IS CALLED /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/overlapping/DEMG_correlation_by_levels_without_overlap.r

 /home/ipbln/egdiaz/bin/ngsplot/bin/ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_accesible_SG.txt -O expression_genebody_SG
 /home/ipbln/egdiaz/bin/ngsplot/bin/ngs.plot.r -G AgamP4 -R bed -C config_ngsplot_expression_SG.txt -O accesibility_tss_SG

#### We searched for diffbind peaks that matched splicing sites of our DUI isoforms. The process performed can be found in : /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Comparison_diffbind_peaks_with_splicing_sites_100bp_without_overlapping_genes.R
#### With the diffbind peaks coincident with splicing sites in the DUI we will do a motif analysis. We do the same with these peaks that are also located in splicing events of the exon skipping type.
# preparation and execution /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/Motif_analysis_diffbind_splicing_sites_DUI_DEF_130423.R
### annotation

#### KNOWN MOTIFS

### MOTIF 1 CATCMCTA

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/knownResults/known1.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_KNWON_ aks_motif1_CATCMCTA.txt # 1e-3

### MOTIF 2 GGYCATAAAW

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/knownResults/known2.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_KNOWN_annotatePeaks_motif2_GGYCATAAAW.txt # 1e-2

#### HOMER MOTIFS

### MOTIF 1 CACCWATT

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif1.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif1_CACCWATT.txt # 1e-34

### MOTIF 2 CTTGCAAGGT

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif2.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif2_CTTGCAAGGT.txt # 1e-31

### MOTIF 3 ACARAGAAAT

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif3.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif3_ACARAGAAAT.txt # 1e-28

### MOTIF 4 ACACAAAG

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif4.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif4_ACACAAAG.txt # 1e-27

### MOTIF 5 GAGTCTTTAT

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif5.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif5_GAGTCTTTAT.txt # 1e-27

### MOTIF 6 GAAAACAGTT

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif6.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif6_GAAAACAGTT.txt # 1e-27
### MOTIF 7 AAACAGTAAAAG

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif7.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif7_AAACAGTAAAAG.txt # 1e-27
### MOTIF 8 CAAGAAGA

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa  -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif8.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif8_CAAGAAGA.txt # 1e-27
### MOTIF 9 GAMCGAGH

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif9.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif9_GAMCGAGH.txt # 1e-27

### MOTIF 10 TTATAAYT

annotatePeaks.pl /data/amd_HDD_4TB/Barbara/Alternative_splicing/ATAC/DiffBind/Diffbind_peaks_with_DUI_and_splicing_sites_whitout_overlapping_motif_analysis.bed  /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -gff /data/amd_HDD_4TB/Barbara/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff -m /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/homerResults/motif10.motif > /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/Motif_analysis/peaks_with_DUI_ans_splicing_sites_without_overlapping_results/annotation/peaks_with_DUI_ans_splicing_sites_HOMER_annotatePeaks_motif10_TTATAAYT.txt # 1e-27

### to study if any type of event is related to specific accessibility of a region, we use the script /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/AS_mechanism_by_conditions_comparing_peak_location.R
## finally we extrac enhancers known by other authos and by our group taking drosophila orthologs and we look how many of them are present and active in our groups.
# /data/amd_HDD_4TB/Barbara/Alternative_splicing_def/enhancers_JL/merge_enhancers_with_DEMG_DUI_def_310123.R
