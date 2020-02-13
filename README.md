# Multi-Trait Genome-wide Association Studies Methods Reveal Loci Associated with Maize Inflorescence and Leaf Architecture
This github contains the files and scripts as well as additional results for the manuscript:
Multi-Trait Genome-wide Association Studies Methods Reveal Loci Associated with Maize Inflorescence and Leaf Architecture

Folder: BLUPS
  "BLUP Model Trait Pipeline Summaries.csv" -> Full model trait details 
  
 "GLMMER_BLUPs.R" -> R script for pipeline to evalute and take raw public data avalible at panzea.org to BLUPS
  "PanzeaDataMaizeBLUPS.csv" -> The resulting trait BLUPs from GLMMER_BLUPs.R
  "Reference GBS taxanames.csv" -> Taxa name reference file for GBS taxa names
  "trait_names.csv" -> Trait name csv file required for GLMMER_BLUPs.R

Folder: Simulation_Study
  "simplePHENOTYPES_logfiles" -> Folder containing log files from simple phenotypes
  "Detriming False Positive rate of Simulation.R" -> R script for determining false and true postive for multiple GAPIT GWAS  output files, requires simplePHENOTYPE output
  "Function to detrimine false and true postive.R" -> Function used in "Detriming False Positive rate of Simulation.R" and "getting false postive rate for gemma output.R"
  "GAPIT.Kin.VanRadenEntireGennome55KSNPFILTEREDIMPUTED.csv" -> A VanRaden Kinship matrix using the all ten chromosomes from the 55k legacy SNP chip data from panzea.org
  "LOCOVanRadenKinship55kLinkImputeFiltered.Rdata" -> 10 VanRaden Kinship matrix from the 55k legacy SNP chip data from panzea.org, where each kinship corresponds to a single chromosome.
  "SNP55kv2_282_GenoFiles_GapitInput.RData"-> 55k SNP chip numeric myGM and myGD files required for GAPIT input 
  "getting false postive rate for gemma output.R" -> R script for detriming false and true postive for multiple GEMMA GWAS output files, require simplePHENOTYPE output
  "run_gemma.R" -> R script for running GEMMA for simplePHENOTYPE output
  "running_gwas.R" -> R script for running GAPIT for simplePHENOTYPE output 
  
Files listed as GAPIT.MLM... correspond to the MLM GWAS model results from GAPIT for each trait BLUP and each PC
"pheno.assoc.txt.zip" -> Multivariate MLM GWAS results from GEMMA
"pheno.log.txt" -> GEMMA log file
"Manhattan plots and QQ plots.zip" ->  zip folder containing the QQ plots and Manhattan plots from each GWAS 
