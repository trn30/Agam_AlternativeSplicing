


########### GENES #############
##### Accessibility at the promoter ######

### I DO READ COUNTS THROUGHOUT THE PROMOTER.

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_genes/R5-I3d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_genes/R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_genes/R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_promoters_genes/R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed



##### Accessibility at the PEAK LOCATED IN THE PROMOTER ######
### I INTERSECT THE HCS PEAKS WITH THE PROMOTERS AND THEN PERFORM THE COUNTS.

bedtools intersect -a  /ATAC/HCS_salivary_glands.bed  -b  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb.bed -f 0.51 -wa >  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_genes_sg_0.51_2_peaks.bed

bedtools intersect -a  /ATAC/HCS_midguts.bed -b  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb.bed -f 0.51 -wa >  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_genes_mg_0.51_2_peaks.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_genes_mg_0.51_2_peaks.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2_peaks_gen.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_genes_sg_0.51_2_peaks.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2_peaks_gen.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_genes_mg_0.51_2_peaks.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2_peaks_gen.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_genes_sg_0.51_2_peaks.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R5-I314_SG_counts_1kb_upstream_SG_vs_MG_2_peaks_gen.bed

##### Accessibility in the gene body ######

### I DO READ COUNTS ON THE WHOLE BODY


bedtools intersect -c -a  /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_genes/R5-I3d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_genes/R1-I1d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_genes/R8-I3d7_MG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_promoters_genes/R10-I1d7_MG_counts_body_SG_vs_MG_2.bed


##### Accessibility at PEAK LOCATED ON THE GENE BODY  ######
### I INTERSECT THE HCS PEAKS WITH THE BODY AND THEN PERFORM THE COUNTS.

bedtools intersect -a  /ATAC/HCS_salivary_glands.bed -b  /genomic_data_AgamP4/Genes_AgamP4_release_54.bed  -f 0.51  -wa >  /ATAC/counts_in_peaks/reads_HCS_gene_body_sg_0.51_2.bed

bedtools intersect -a  /ATAC/HCS_midguts.bed  -b  /genomic_data_AgamP4/Genes_AgamP4_release_54.bed -f 0.51  -wa >  /ATAC/counts_in_peaks/reads_HCS_gene_body_mg_0.51_2.bed


bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_gene_body_mg_0.51_2.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R10-I1d7_MG_counts_body_SG_vs_MG_2_peaks_gen.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_gene_body_sg_0.51_2.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R1-I1d14_SG_counts_body_SG_vs_MG_2_peaks_gen.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_gene_body_mg_0.51_2.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R8-I3d7_MG_counts_body_SG_vs_MG_2_peaks_gen.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_gene_body_sg_0.51_2.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R5-I314_SG_counts_body_SG_vs_MG_2_peaks_gen.bed


############################################################################################################################################################################################################################################################################

########### ISOFORMS #############
##### Accessibility at the promoter ######
### I DO READ COUNTS THROUGHOUT THE PROMOTER.

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_all_isoforms/R5-I3d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_all_isoforms/R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_all_isoforms/R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_promoters_all_isoforms/R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2.bed



##### Accessibility at the PEAK LOCATED IN THE PROMOTER ######
### I INTERSECT THE HCS PEAKS WITH THE PROMOTERS AND THEN PERFORM THE COUNTS.


bedtools intersect -a  /ATAC/HCS_salivary_glands.bed  -b  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -f 0.51 -wa >  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_iso_sg_0.51_2_peaks.bed

bedtools intersect -a  /ATAC/HCS_midguts.bed -b  /genomic_data_AgamP4/Promoters_AgamP4_release_54_1kb_all_isoforms.bed -f 0.51 -wa >  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_iso_mg_0.51_2_peaks.bed


bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_iso_mg_0.51_2_peaks.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R10-I1d7_MG_counts_1kb_upstream_SG_vs_MG_2_peaks_iso.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_iso_sg_0.51_2_peaks.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R1-I1d14_SG_counts_1kb_upstream_SG_vs_MG_2_peaks_iso.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_iso_mg_0.51_2_peaks.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R8-I3d7_MG_counts_1kb_upstream_SG_vs_MG_2_peaks_iso.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_regions_1kb_upstream_iso_sg_0.51_2_peaks.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R5-I314_SG_counts_1kb_upstream_SG_vs_MG_2_peaks_iso.bed


##### Accessibility in the ISOFORM body ######

### I DO READ COUNTS ON THE WHOLE BODY

bedtools intersect -c -a  /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_all_isoforms/R5-I3d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_all_isoforms/R1-I1d14_SG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_promoters_all_isoforms/R8-I3d7_MG_counts_body_SG_vs_MG_2.bed

bedtools intersect -c -a  /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_promoters_all_isoforms/R10-I1d7_MG_counts_body_SG_vs_MG_2.bed



##### Accessibility at the PEAK LOCATED ON THE ISOFORM BODY  ######
### I INTERSECT THE HCS PEAKS WITH THE BODY AND THEN PERFORM THE COUNTS.

bedtools intersect -a  /ATAC/HCS_salivary_glands.bed -b  /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed  -f 0.51 -wa >  /ATAC/counts_in_peaks/reads_HCS_iso_body_sg_0.51_2_peaks.bed

bedtools intersect -a  /ATAC/HCS_midguts.bed  -b  /genomic_data_AgamP4/Isoforms_AgamP4_release_54.bed -f 0.51  -wa >  /ATAC/counts_in_peaks/reads_HCS_iso_body_mg_0.51_2_peaks.bed


bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_iso_body_mg_0.51_2_peaks.bed -b  /ATAC/R10-I1d7_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R10-I1d7_MG_counts_body_SG_vs_MG_2_peaks_iso.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_iso_body_sg_0.51_2_peaks.bed -b  /ATAC/R1-I114d_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R1-I1d14_SG_counts_body_SG_vs_MG_2_peaks_iso.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_iso_body_mg_0.51_2_peaks.bed -b  /ATAC/R8-I3d7_nucfree_mapq_10.bam.cram.bam_chr.bam  >  /ATAC/counts_in_peaks/R8-I3d7_MG_counts_body_SG_vs_MG_2_peaks_iso.bed

bedtools intersect -c -a  /ATAC/counts_in_peaks/reads_HCS_iso_body_sg_0.51_2_peaks.bed -b  /ATAC/R5-I314d_nucfree_mapq_10.bam.cram.bam_chr.bam >  /ATAC/counts_in_peaks/R5-I314_SG_counts_body_SG_vs_MG_2_peaks_iso.bed
