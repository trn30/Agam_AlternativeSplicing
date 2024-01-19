# Agam_AlternativeSplicing
Analysis pipeline used in the article "Alternative splicing and its regulation in the malaria vector Anopheles gambiae"(https://doi.org/10.1101/2023.07.18.549290)
## Authors
Bárbara Díaz-Terenti, Elena Gómez-Díaz*
###### Instituto de Parasitología y Biomedicina López-Neyra, Consejo Superior de Investigaciones Científicas (IPBLN, CSIC), Granada, Spain

## Abstract
 
Alternative splicing (AS) is a highly conserved mechanism that allows to expand the coding capacity of the genome, by modifying the way multiple isoforms are expressed or used to generate different phenotypes. Despite its importance in physiology and disease, genome-wide studies of AS are lacking in most insects, including mosquitoes. Even for model organisms, chromatin associated processes involved in the regulation AS are poorly known. In this study, we investigated AS in the mosquito Anopheles gambiae in the context of tissue-specific gene expression and mosquito responses to a Plasmodium falciparum infection, as well as the relationship between patterns of differential isoform expression and usage with chromatin accessibility changes.  For this, we combined RNA-seq and ATAC-seq data from A. gambiae midguts and salivary glands, and from infected and non-infected midguts. We report differences between tissues in the expression of 456 isoforms and in the use of 211 isoforms. Secondly, we find a clear and significant association between chromatin accessibility states and tissue-specific patterns of AS. The analysis of differential accessible regions located at splicing sites permitted the identification of several motifs resembling the binding sites of Drosophila transcription factors. Finally, the genome-wide analysis of tissue-dependent enhancer activity revealed that approximately 20% of A. gambiae transcriptional enhancers annotate to a differentially expressed or used isoform and that their activation status is linked to AS differences between tissues. This research elucidates the role of AS in mosquito vector gene expression and identifies regulatory regions potentially involved in AS regulation, which could be important in the development of novel strategies for vector control.


## Code Explanation
  The full pipeline of the analysis is in the file **Alternative_splicing_project_def_code_Github.sh**, this file will redirect you to other scripts also present in this repository. An outline summary of the contents of the scripts is:

 - AS_mechanism_by_conditions.R
   
  See how many type of events apear in the DEMG and DUI isoforms of the _Inf MG vs. Inf SG_
 + AS_mechanism_by_conditions_CHI_test.R


 * AS_mechanism_by_conditions_comparing_peak_location.R


 - check_overlaping_zones.R


 - Comparison_diffbind_peaks_with_splicing_sites_100bp_without_overlapping_genes_Github.R


 - correlation_DEMG_without_overlap_definitive_120423.R


 - Ctrl_vd_Infec_MG_analysis_210722_Github.R


 - Ctrl_vd_Infec_SG_analysis_210722_Github.R


 - DEMG_correlation_by_levels_without_overlap.r


 - DESeq2_analysis_MG_vs_SG_070422.R


 - Diffbind_correlation.R


 - files_correlations_DEG_DU_Github.sh


 - Genes_paper_malar_development_parasites.R


 - isoform_switch_Ctrl_vs_Infec.R


 - Isoform_switch_MG_vs_SG.R


 - merge_enhancers_with_DEMG_DUI_def_310123_Github.R


 - Motif_analysis_diffbind_splicing_sites_DUI_DEF_130423_Github.R


 - script_volcanos_percent_DUI_vs_DEG_def_210722.R













