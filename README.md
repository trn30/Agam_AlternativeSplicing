# Agam_AlternativeSplicing
Analysis pipeline used in the article "Alternative splicing and its regulation in the malaria vector Anopheles gambiae"(https://doi.org/10.1101/2023.07.18.549290)
## Authors
Bárbara Díaz-Terenti, Elena Gómez-Díaz*
###### Instituto de Parasitología y Biomedicina López-Neyra, Consejo Superior de Investigaciones Científicas (IPBLN, CSIC), Granada, Spain

## Abstract
 
Alternative splicing (AS) is a highly conserved mechanism that allows to expand the coding capacity of the genome, by modifying the way multiple isoforms are expressed or used to generate different phenotypes. Despite its importance in physiology and disease, genome-wide studies of AS are lacking in most insects, including mosquitoes. Even for model organisms, chromatin associated processes involved in the regulation AS are poorly known. In this study, we investigated AS in the mosquito Anopheles gambiae in the context of tissue-specific gene expression and mosquito responses to a Plasmodium falciparum infection, as well as the relationship between patterns of differential isoform expression and usage with chromatin accessibility changes.  For this, we combined RNA-seq and ATAC-seq data from A. gambiae midguts and salivary glands, and from infected and non-infected midguts. We report differences between tissues in the expression of 456 isoforms and in the use of 211 isoforms. Secondly, we find a clear and significant association between chromatin accessibility states and tissue-specific patterns of AS. The analysis of differential accessible regions located at splicing sites permitted the identification of several motifs resembling the binding sites of Drosophila transcription factors. Finally, the genome-wide analysis of tissue-dependent enhancer activity revealed that approximately 20% of A. gambiae transcriptional enhancers annotate to a differentially expressed or used isoform and that their activation status is linked to AS differences between tissues. This research elucidates the role of AS in mosquito vector gene expression and identifies regulatory regions potentially involved in AS regulation, which could be important in the development of novel strategies for vector control.


## Code Explanation
 The full pipeline of the analysis is in the file **Alternative_splicing_project_def_code_Github.sh**, this file will redirect you to other scripts also present in this repository. An outline summary of the contents of the scripts is:

 - **DESeq2_analysis_MG_vs_SG.R**

   Obtain Differential expressed isoforms of multisoform genes (DEMG) in the _Inf MG vs. Inf SG_ comparison

 - **Ctrl_vd_Infec_MG_analysis.R**

   Obtain Differential expressed isoforms of multisoform genes (DEMG) in the _Ctrl vs. Inf MG_ comparison

 - **Ctrl_vd_Infec_SG_analysis.R**

   Obtain Differential expressed isoforms of multisoform genes (DEMG) in the _Ctrl MG vs. Inf SG_ comparison

 - **Isoform_switch_MG_vs_SG.R**

   Quantification of the usage with IsoformSwitchAnalyzeR, obtaing isoforms differentially used (DUI) in the _Inf MG vs. Inf SG_ comparison

 - **Isoform_switch_Ctrl_vs_Infec.R**

   Quantification of the usage with IsoformSwitchAnalyzeR, obtaing isoforms differentially used (DUI) in the _Ctrl vs. Inf MG_ and _Ctrl vs. Inf SG_ comparisons
   
 - **Script_volcanos_percent_DUI_vs_DEG.R**

   Obtain  the volcano plot of DEMG vs DUI in the _Inf MG vs. Inf SG_ comparison
   
 - **AS_mechanism_by_conditions.R**
   
   See how many type of events apear in the DEMG and DUI isoforms of the _Inf MG vs. Inf SG_

 - **Genes_paper_malar_development_parasites.R**

   Comparison of our DEMG and DUI isoforms with genes present in the paper https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-12-216

 + **Check_overlaping_zones.R**

   How to check for overlapping areas by removing overlapping genes with promoters of our isoforms and calculate the correlation of the DEMG cluster without these overlapping genes.

 - **Files_correlations_DEG_DUI.sh**

   Enrichment counts at promoters and in the body of isoforms, these counts are necessary for the study of the correlation between expression and accessibility
   
 - **Correlation_DEMG_without_overlap.R**

   To see if there is a correlation between the expression and accessibility of the DEMG isoforms in the comparison _Inf MG vs. Inf SG_ 

 - **Comparison_diffbind_peaks_with_splicing_sites_100bp_without_overlapping_genes.R**

   We searched for diffbind peaks that matched splicing sites of our DUI isoforms

 - **Motif_analysis_diffbind_splicing_sites_DUI.R**

   With the diffbind peaks coincident with splicing sites in the DUI we will do a motif analysis. We do the same with these peaks that are also located in splicing events of the exon skipping type

 * **AS_mechanism_by_conditions_comparing_peak_location.R**

   To study if any type of event is related to specific accessibility of a region, we use the script AS_mechanism_by_conditions_comparing_peak_location.R

 - **Merge_enhancers_with_DEMG_DUI.R**

   How we extrac enhancers known by other authors and by our group taking _Drosophila_ orthologs and we look how many of them are present and active in our groups of DEMG and DUI






 















