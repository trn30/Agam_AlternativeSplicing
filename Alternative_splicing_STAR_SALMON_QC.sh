#### 	MG vs SG


#################################################          STAR            #######################################################################################################


cd /Alternative_splicing/raw_data/paired ## same as /A_gambiae_P_falciparum_infected_RNA-seq_NARGB_2021/RNA_Agam_0119/out_sortmerna

dsrc d -t10 EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq

dsrc d -t10 EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_2.fq.fastq.dsrc.fastq


STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /P_falciparum_ATAC-seq_NAR_2018/ATAC_seq/igenomes/Anopheles-gambiae_AgamP4_STAR_idx --genomeFastaFiles /P_falciparum_ATAC-seq_NAR_2018/ATAC_seq/igenomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa --sjdbGTFfile /P_falciparum_ATAC-seq_NAR_2018/ATAC_seq/igenomes/VectorBase-54_AgambiaePEST.gff.gtf --sjdbOverhang 150 --genomeSAindexNbases 13



### Mapping all files:
     cd /Alternative_splicing/raw_data/paired ## same as /A_gambiae_P_falciparum_infected_RNA-seq_NARGB_2021/RNA_Agam_0119/out_sortmerna
        # sortmerna with default databases + custom databases: (esta es la mejor opción basada en los fastqc)
        # for f in $(ls | egrep -E '_no_rRNA_def_cust_dbs.fastq$' ); do unmerge-paired-reads.sh "${f}" "${f}_1.fq" "${f}_2.fq"; done &&

  for f in $(ls | egrep -E '_1.fastq.dsrc.fastq$' );do STAR --genomeDir /P_falciparum_ATAC-seq_NAR_2018/ATAC_seq/igenomes/Anopheles-gambiae_AgamP4_STAR_idx --sjdbGTFfile /P_falciparum_ATAC-seq_NAR_2018/ATAC_seq/igenomes/VectorBase-54_AgambiaePEST.gff.gtf --readFilesIn "${f}" "$(sed 's/_1.fq.fastq.dsrc.fastq/_2.fq.fastq.dsrc.fastq/g' <<<"${f}")" --outSAMtype BAM Unsorted --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --alignIntronMin 1 --alignIntronMax 249416 --outFileNamePrefix /Alternative_splicing/out_STAR_2024/"${f}_mapped.bam" --runThreadN 20 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3; done



## QC:
        cd /Alternative_splicing/out_STAR_2024/ &&
            for f in $(ls | egrep -E '.bam$'); do samtools sort -@ 20 -o "${f}_sorted.bam" "${f}"; done &&
        for f in $(ls | egrep -E '_sorted.bam$'); do samtools index -@ 20 "${f}"; done &&
        for f in $(ls | egrep -E '_sorted.bam$'); do qualimap bamqc -bam "${f}" -nt 20 -c; done &&
        for f in $( ls | egrep -E '_sorted.bam$' ); do samtools stats -@ 20 "${f}" | grep ^SN | cut -f 2- && echo "${f}"; done &&
        for f in $( ls | egrep -E '_sorted.bam$' ); do samtools flagstat -@ 20 "${f}" && echo "${f}"; done

###########################################################           salmon       ###########################################################################################

# mapped mode
# In this mode, the programme creates its own alignments and then performs the count. To use it, the first thing we have to do is to create an index of the genome we are interested in.

salmon index -t /Genomes/genomic_data_AgamP4/Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.12.fa -i transcripts_index_24  -k 31

# The number of k-mers should be adapted according to the size of the reads to make it more sensitive. 31 is fine for reads of 75 bp or more.
cd /Alternative_splicing/raw_data/paired/

for f in $(ls | egrep -E '_1.fq.fastq.dsrc.fastq$' );do salmon quant -i /Alternative_splicing/transcripts_index_24 -l A -1  "${f}" -2 "$(sed 's/_1.fq.fastq.dsrc.fastq/_2.fq.fastq.dsrc.fastq/g' <<<"${f}")" --threads 20 -o  /Alternative_splicing/out_salmon_2024/"${f}"  --validateMappings --go 4 --mismatchSeedSkip 5 ; done  # OPTIMISED


###########################################################           QUALIMAP       ###########################################################################################


### SINGLE

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-1-I1d7_Single -a proportional -bam /Alternative_splicing/out_star_tfm/EGD-1-I1d7_S1_L001_R1_001_mapped.bam_correct_chr.bam -p strand-specific-forward -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf  --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-2-I1d14_Single -a proportional -bam /Alternative_splicing/out_star_tfm/EGD-2-I1d14_S5_L001_R1_001_mapped.bam_correct_chr.bam -p strand-specific-forward -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf  --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-7-I3d7_Single -a proportional -bam /Alternative_splicing/out_star_tfm/EGD-7-I3d7_S3_L001_R1_001_mapped.bam_correct_chr.bam -p strand-specific-forward -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf  --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-8-I3d14_Single -a proportional -bam /Alternative_splicing/out_star_tfm/EGD-8-I3d14_S4_L001_R1_001_mapped.bam_correct_chr.bam -p strand-specific-forward -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf  --java-mem-size=16G


### PAIRED
qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-1-I1d7_Paired -a proportional -bam /Alternative_splicing/out_STAR_2024/EGD-1-I1d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq_mapped.bamAligned.out.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf  --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-2-I1d14_Paired  -a proportional -bam /Alternative_splicing/out_STAR_2024/EGD-2-I1d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastqmapped.bamAligned.out.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-7-I3d7_Paired  -a proportional -bam /Alternative_splicing/out_STAR_2024/EGD-7-I3d7_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq_mapped.bamAligned.out.bam -p strand-specific-forward -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/EGD-8-I3d14_Paired  -a proportional -bam /Alternative_splicing/out_STAR_2024/EGD-8-I3d14_merged_reads.trimmed.fastq_def_cust_dbs_no_rRNA.fastq_1.fq.fastq.dsrc.fastq_mapped.bamAligned.out.bam -p strand-specific-forward -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G






###########################################################           MULTIQC       ###########################################################################################



multiqc  /Alternative_splicing_def/qualimap/Paired /Alternative_splicing_def/FastQC/Paired  /Alternative_splicing/out_STAR_2024/Paired  /Alternative_splicing/out_salmon_2024/Paired  -o /Alternative_splicing_def/MultiQC/Paired






###############################################################################################################################################################################################################
###############################################################################################################################################################################################################


#### Ctrol vs Inf SG

cd /Alternative_splicing/raw_data/sortmerna/Salivary_glands

 for f in $(ls | egrep -E '.dsrc$' ); do dsrc d -t20 "${f}" "${f}.fastq"   ; done



###########################################################           salmon       ###########################################################################################

# mapped mode
# In this mode, the programme creates its own alignments and then performs the count. To use it, the first thing we have to do is to create an index of the genome we are interested in.


# The number of k-mers should be adapted according to the size of the reads to make it more sensitive. 31 is fine for reads of 75 bp or more.
cd /Alternative_splicing/raw_data/sortmerna/Salivary_glands

 for f in $(ls | egrep -E '_1.fastq.dsrc.fastq$' );do salmon quant -i /Alternative_splicing/transcripts_index_24 -l A -1  "${f}" -2 "$(sed 's/_1.fastq.dsrc.fastq/_2.fastq.dsrc.fastq/g' <<<"${f}")" --threads 20 -o  /Alternative_splicing/out_salmon_2024/Salivary_glands/"${f}"  --validateMappings --go 4 --mismatchSeedSkip 5 ; done  # OPTIMISED





###########################################################           QUALIMAP       ###########################################################################################


qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Salivary_glands/Salivary_glands_Infected_14d_M3 -a proportional -bam /Alternative_splicing/out_star_tfm/Salivary_glands_bam/Salivary_glands_Infected_14d_M3Aligned.out.bam_sort.bam_correct_chr.bam -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf  --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Salivary_glands/Salivary_glands_Infected_14d_M2 -a proportional -bam /Alternative_splicing/out_star_tfm/Salivary_glands_bam/Salivary_glands_Infected_14d_M2Aligned.out.bam_sort.bam_correct_chr.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Salivary_glands/Salivary_glands_Non-Infected_14d_M1 -a proportional -bam /Alternative_splicing/out_star_tfm/Salivary_glands_bam/Salivary_glands_Non-Infected_14d_M1Aligned.out.bam_sort.bam_correct_chr.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Salivary_glands/Salivary_glands_Non-Infected_14d_M2 -a proportional -bam /Alternative_splicing/out_star_tfm/Salivary_glands_bam/Salivary_glands_Non-Infected_14d_M2Aligned.out.bam_sort.bam_correct_chr.bam -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G


###########################################################           MULTIQC       ###########################################################################################



multiqc  /Alternative_splicing_def/qualimap/Salivary_glands /Alternative_splicing_def/FastQC/Salivary_glands /Alternative_splicing/out_star_tfm/Salivary_glands_bam  /Alternative_splicing/out_salmon_2024/Salivary_glands -o /Alternative_splicing_def/MultiQC/Salivary_glands







###############################################################################################################################################################################################################
###############################################################################################################################################################################################################




#### Ctrol vs Inf MG

cd /Alternative_splicing/raw_data/sortmerna/Midguts

 for f in $(ls | egrep -E '.dsrc$' ); do dsrc d -t18 "${f}" "${f}.fastq"   ; done
 for f in $(ls | egrep -E '.fastq$' ); do reformat.sh in= "${f}" out1="${f}_1.fastq" out2="${f}_2.fastq" t=20  ; done




###########################################################           salmon       ###########################################################################################

# mapped mode
# In this mode, the programme creates its own alignments and then performs the count. To use it, the first thing we have to do is to create an index of the genome we are interested in.


# The number of k-mers should be adapted according to the size of the reads to make it more sensitive. 31 is fine for reads of 75 bp or more.
cd /Alternative_splicing/raw_data/sortmerna/Midguts

 for f in $(ls | egrep -E '_1.fastq.dsrc.fastq$' );do salmon quant -i /Alternative_splicing/transcripts_index_24 -l A -1  "${f}" -2 "$(sed 's/_1.fastq.dsrc.fastq/_2.fastq.dsrc.fastq/g' <<<"${f}")" --threads 20 -o  /Alternative_splicing/out_salmon_2024/Midguts/"${f}"  --validateMappings --go 4 --mismatchSeedSkip 5 ; done  # OPTIMISED






###########################################################           QUALIMAP       ###########################################################################################


qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Midguts/Midgut_Infected_7d_M2 -a proportional -bam /Alternative_splicing/out_star_tfm/Midguts_bam/Midgut_Infected_7d_M2Aligned.out.bam_sort.bam_correct_chr.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Midguts/Midgut_Infected_7d_M3 -a proportional -bam /Alternative_splicing/out_star_tfm/Midguts_bam/Midgut_Infected_7d_M3Aligned.out.bam_sort.bam_correct_chr.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Midguts/Midgut_Non-Infected_7d_M2 -a proportional -bam /Alternative_splicing/out_star_tfm/Midguts_bam/Midgut_Non-Infected_7d_M2Aligned.out.bam_sort.bam_correct_chr.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G

qualimap rnaseq -outdir /Alternative_splicing_def/qualimap/Midguts/Midgut_Non-Infected_7d_M3  -a proportional -bam /Alternative_splicing/out_star_tfm/Midguts_bam/Midgut_Non-Infected_7d_M3Aligned.out.bam_sort.bam_correct_chr.bam  -gtf /Genomes/genomic_data_AgamP4/VectorBase-54_AgambiaePEST.gff.gtf -pe --java-mem-size=16G



###########################################################           MULTIQC       ###########################################################################################


multiqc  /Alternative_splicing_def/qualimap/Midguts /Alternative_splicing_def/FastQC/Midguts /Alternative_splicing/out_star_tfm/Midguts_bam  /Alternative_splicing/out_salmon_2024/Midguts -o /Alternative_splicing_def/MultiQC/Midguts


#################################################################################################################################################


### repetition of the PCA and the PCOAs are in PCOA_AS_ctrl_vs_inf_050224.R
