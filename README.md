# Code and data to reproduce results from Saurabh's PhD  thesis

<table>
  <tr>
   <td>Step No
   </td>
   <td>Step
   </td>
   <td>Details
   </td>
   <td>Input
   </td>
   <td>Method
   </td>
   <td>Program file
   </td>
   <td>Output
   </td>
   <td>Output file
   </td>
   <td>Notes
   </td>
  </tr>
  <tr>
   <td>1
   </td>
   <td>Phylogenies curation
   </td>
   <td>GTDB phylogeny
   </td>
   <td>
   </td>
   <td>Manual download
<p>
Pruning by custom script
   </td>
   <td>/Codes/tree_pruning_GTDB.Rmd
   </td>
   <td>GTDB bacteria phylogeny, release 80
<p>
Pruned GTDB phylogeny
   </td>
   <td>/Input_data/phylogenies/GTDB_bac120_r80.tree
<p>
/Input_data/phylogenies/GTDB_bac120_r80_pruned0.01.newick
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Order level phylogenies
   </td>
   <td>Pruned GTDB phylogeny
   </td>
   <td>Manual extraction
   </td>
   <td>NA
   </td>
   <td>Order-level phylogenies
   </td>
   <td>Example [/Input_data/phylogenies/Bacteroidale_GTDB_14k_tree.newick]
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>2
   </td>
   <td>GC data collection
   </td>
   <td>Metadata from GTDB dataset
   </td>
   <td>
   </td>
   <td>Manual download
   </td>
   <td>NA
   </td>
   <td>Metadata files
   </td>
   <td>/Input_data/phylogenies/GTDB_bac_metadata_r80.tsv
<p>
/Input_data/phylogenies/GTDB_bac_metadata_r80_pruned0.01.tsv
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>3
   </td>
   <td>Phylogenetic model fitting
   </td>
   <td>BM
   </td>
   <td>Order-level phylogenies
<p>
GC data [/Input_data/levolution_files_inference_actual_data/Bacteroidale_GTDB_14k_GC.txt]
   </td>
   <td>Geiger
   </td>
   <td>/Codes/fitting_Brownian_and_OU_to_order_level_clades.Rmd
   </td>
   <td>Maximum likelihood and best-fit parameters
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>OU
   </td>
   <td>Order-level phylogenies
<p>
GC data [/Input_data/levolution_files_inference_actual_data/Bacteroidale_GTDB_14k_GC.txt]
   </td>
   <td>Geiger
   </td>
   <td>/Codes/fitting_Brownian_and_OU_to_order_level_clades.Rmd
   </td>
   <td>Maximum likelihood and best-fit parameters
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Levy jumps
   </td>
   <td>Order-level phylogenies
<p>
GC data [/Input_data/levolution_files_inference_actual_data/Bacteroidale_GTDB_14k_GC.txt]
<p>
Config file [/Input_data/levolution_config_files/levolution_Bacteroidale_GTDB_14k_GC_config1_alpha0.1.param]
   </td>
   <td>Levolution
   </td>
   <td>Direct in unix shell according to levolution website
   </td>
   <td>Best-fit parameter values
<p>
Log files
<p>
Parameter used for model inference
<p>
Phylogenies with branch specific posterior probabilities based on actual data
   </td>
   <td>/Input_data/levolution_files_inference_actual_data/Bacteroidale_GTDB_14k_GC_config1_alpha0.1_oneLineSummary.txt
<p>
/Input_data/levolution_files_inference_actual_data/levolution_Bacteroidale_GTDB_14k_GC_config1_alpha0.1.log
<p>
/Input_data/levolution_files_inference_actual_data/levolution_Bacteroidale_GTDB_14k_GC_config1_alpha0.1.param
<p>
/Input_data/levolution_files_inference_actual_data/Bacteroidale_GTDB_14k_GC_config1_alpha0.1.post
   </td>
   <td>These are just examples for Bacteroidale, for one value of alpha. Same command can be run for multiple values of alpha and multiple order-level clades.
<p>
Branch length in posterior probability phylogeny represent the posterior probability of occurrence of jump(s)
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Summary of model fitting
   </td>
   <td>Model fitting outputs
   </td>
   <td>Custom script
   </td>
   <td>/Codes/summary_of_levolution_Levy_jumps_best_fit_parameters.Rmd
   </td>
   <td>File with best fit model parameters
   </td>
   <td>/Input_data/levolution_files_inference_actual_data/levolution_best_fit_parameters_summary.txt
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>4
   </td>
   <td>Jump identification
   </td>
   <td>Simulations to identify thresholds
   </td>
   <td>Order level phylogenies
<p>
Best fit parameter values
   </td>
   <td>Custom scripts
   </td>
   <td>/Codes/GC_Levy_jumps_simulations_GTDB.Rmd
   </td>
   <td>Simulated GC data
<p>
Simulated jump details
   </td>
   <td>/Input_data/levolution_files_inference_simulated_data/Bacteroidale_GTDB_14k_GC_config1_alpha0.25_geigerSim1.traits
   </td>
   <td>5 independent simulations for each order level clade
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Jump inference on simulated data
   </td>
   <td>Order-level phylogenies
<p>
Simulated GC data
   </td>
   <td>Levolution
   </td>
   <td>Direct in unix shell according to levolution website
   </td>
   <td>Posterior probabilities of jumps in simulated data
   </td>
   <td>/Input_data/levolution_files_inference_simulated_data/Bacteroidale_GTDB_14k_GC_config1_alpha0.25_geigerSim1.post
   </td>
   <td>For 5 independent sets of simulated data for each order-level clade
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Precision recall calculations
   </td>
   <td>Posterior probabilities of jumps in simulated data
<p>
Actually simulated jump details
   </td>
   <td>Custom script
   </td>
   <td>/Codes/precision_recall_curves_simulated_jumps.Rmd
   </td>
   <td>Precision recall curves
<p>
Precision recall summaries
   </td>
   <td>/Input_data/levolution_files_inference_simulated_data/Bacteroidale_GTDB_14k_GC_config1_alpha0.25_geigerSim1.details
<p>
/Input_data/levolution_files_inference_simulated_data/PRC_all_facets_Jan2021.png
   </td>
   <td>Pooled data from 5 independent simulations for each order-level clade
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Identification of thresholds
   </td>
   <td>Precision-recall curves
   </td>
   <td>Manual
   </td>
   <td>
   </td>
   <td>Posterior probability thresholds for each order-level clade
   </td>
   <td>
   </td>
   <td>For each order-level clade
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Branch identification and  visual mapping
   </td>
   <td>Phylogenies with branch specific posterior probabilities based on actual data
<p>
Posterior probability thresholds for each order-level clade
   </td>
   <td>Custom script
   </td>
   <td>/Codes/levolution_jump_mapping_phylogeny.Rmd
   </td>
   <td>Images of jumps mapped on the phylogenies
   </td>
   <td>/levolution_files_inference_actual_data/levolution_Bacteroidale_GTDB_14k_GC_config1_alpha0.25_pp0.95_unlabelled_jumpNum_magma_scale.png
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Jump magnitude analysis
   </td>
   <td>Phylogenies with branch specific posterior probabilities based on actual data
<p>
Posterior probability thresholds for each order-level clade
<p>
GC data
   </td>
   <td>Custom script
   </td>
   <td>/Codes/levolution_jump_magnitudes_and_direction.Rmd
   </td>
   <td>GC content distribution
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>5
   </td>
   <td>Ecological correlates
   </td>
   <td>Lifestyle and habitat data
   </td>
   <td>Collected from literature
   </td>
   <td>Manually
   </td>
   <td>NA
   </td>
   <td>Data on lifestyle and habitats of focal and sister clades
   </td>
   <td>Table 6.2
<p>
AS1, oxygen dependence sheet
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Data on oxygen dependence
   </td>
   <td>Downloaded from Madin et al (2020)
   </td>
   <td>Custom script
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>/Input_data/trait_data/
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Analysis of changes in lifestyle/habitat/oxygen dependence
   </td>
   <td>Data on lifestyle and habitats of focal and sister clades
   </td>
   <td>Manually
   </td>
   <td>NA
   </td>
   <td>Changes in lifestyle/habitat/oxygen dependence
   </td>
   <td>Table 6.2
<p>
AS1, oxygen dependence sheet
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>6
   </td>
   <td>dN/dS analysis using the method of Gueguen and Duret (2018)
   </td>
   <td>Estimation of dN/dS
   </td>
   <td>Summary details of datasets and parameters used [/Input_data/GD18_datasets_summary2.tab].
<p>
Bppml config template file [/Codes/ml_nonhom_gammaHet.bpp]
<p>
MapNH config template file [/Codes/map_dNdS_per_branch_template.bpp]
   </td>
   <td>Custom script
<p>
(Dependencies: Python, Bpp suite, Standalone BLAST, PRANK, Gblocks)
   </td>
   <td>Master pipeline notebook:
<p>
/Codes/dNdS_estimation_pipeline2a.Rmd
   </td>
   <td>Best fit model parameters
<p>
Counts of dN and dS
   </td>
   <td>Example, model fit parameters:
<p>
Example, dN, dS counts:
<p>
/Results/mapNH/Paracoccus/Paracoccus_HEG_orthologs40_set1.fa.counts_dN.dnd
<p>
/Results/mapNH/Paracoccus/Paracoccus_HEG_orthologs40_set1.fa.counts_dS.dnd
   </td>
   <td>Performs full pipeline for dNdS estimation by GD18 method.
<p>
Change the name of the clade for each analysis in the master pipeline notebook.
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>Meta-analysis of dNdS estimated using GD18 method
   </td>
   <td>dN and dS counts from mapNH
   </td>
   <td>Custom script
   </td>
   <td>/Codes/dNdS_GC3_metaanalysis_revised.R
   </td>
   <td>Compiled data for all clades analyzed
   </td>
   <td>/Results/dNdS/focal_clades_dNdS_gc3_metadata_orthologs40_set1_revised.tab
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>7
   </td>
   <td>Simulations of codon evolution models, analysis and comparison with YN00
   </td>
   <td>
   </td>
   <td>Bppml config files 
<p>
[/Results/mapNH/Buchnera/Buchnera7a_simulations_new/Buchnera_seq_simul_eOmegaSmaller_eGCLess.bpp etc.]
<p>
Buchnera simulation phylogenies
<p>
[/Results/mapNH/Buchnera/Buchnera7a_simulations_new/Buchnera7a_simul_intree_eBranchesEq.newick etc.]
<p>
Bppml config files [/Results/mapNH/Buchnera/Buchnera7a_simulations_new/ml_nonhom_Buchnera7a_simul.bpp]
<p>
MapNH config file [/Results/mapNH/Buchnera/Buchnera7a_simulations_new/map_dNdS_per_branch_Buchnera7a_simul.bpp]
<p>
PAML config file
<p>
[/Results/YN00/Buchnera/Buchnera_simul_yn00_new/Buchnera_yn00_unEq_weights_simul_template.ctl]
   </td>
   <td>Custom script (Dependencies: bpp suite, PAML)
   </td>
   <td>/Codes/bppml_behavior_simulations_revised.Rmd
   </td>
   <td>Comparison figure of bppml, YN00, and simulated omega/dNdS values
   </td>
   <td>/Results/dNdS/Buchnera/simulation_results_dN_dS_omega_newer.svg
   </td>
   <td>Buchnera folders in dNdS, mapNH, and YN00 directories are compressed.
   </td>
  </tr>
  <tr>
   <td>8
   </td>
   <td>Genome features analysis, including dN/dS estimation by YN00, clustering analysis
   </td>
   <td>
   </td>
   <td>Phylogenies with branch specific posterior probabilities based on actual data
<p>
Metadata from GTDB dataset
   </td>
   <td>Custom script
<p>
(Dependencies: PAML)
   </td>
   <td>/Codes/levolution_jumps_pairwise_dNdS_comparisons.Rmd
   </td>
   <td>Comparison of dN/dS using YN00
<p>
Genome feature changes
<p>
Clustering by genome features
   </td>
   <td>/Results/dNdS/jump_analysis/all_features_changes_table.tab
<p>
/Results/dNdS/jump_analysis/all_features_comparison_table.tab
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>9
   </td>
   <td>Analysis of repair enzymes gain/loss
   </td>
   <td>
   </td>
   <td>HMM profiles of enzymes
<p>
Phylogenies with branch specific posterior probabilities based on actual data
<p>
Metadata from GTDB dataset
   </td>
   <td>Manually downloaded from Uniprot
<p>
Custom scripts
<p>
(Dependencies: HMMER3)
   </td>
   <td>/Codes/repair_enzyme_analysis_GTDB_focal_clades.Rmd
<p>
 
<p>
/Codes/repair_enzyme_detection_GTDB_focal_clades.Rmd
   </td>
   <td>Unique gain, loss events in downward and upward jumps
   </td>
   <td>/Results/dNdS/jump_analysis/levolution_downwardJumps_deviant_enz.tab
<p>
/Results/dNdS/jump_analysis/levolution_upwardJumps_deviant_enz.tab
<p>
/Results/dNdS/jump_analysis/levolution_enzyme_gain_loss.xls
   </td>
   <td>
   </td>
  </tr>
</table>

