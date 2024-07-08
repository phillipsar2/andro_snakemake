# Description of various scripts used in big bluestem analyses

## Variant calling and filtering
`genoDPfilter.R` -- filter SNPs based on genotype depth and missingness   
`bwa_index.sh` -- index genome for use in bwa-mem2  
`vcf2ADmatrix.R` -- extract allele depth from VCF and convert to input format for EBG  
`gltable2mpgl.R` -- convert a table of PL values from GATK VariantsToTable to mpgl format for EBG  


## Population genetics
`bigblue_fst_rho.R` -- Calculate Fst and Rho
`bigblue_Fis.R` -- Plot individual inbreeding coefficients estimated from ngsF
`bigblue_STRUCTURE_plots.R` -- Create input file for STRUCTURE and plot output
`bigblue_thetas.R` -- Calculate per-site values of thetas estimated by ANGSD
`singleread_pca.R` -- Plot single read PCA
`kinship_matrix.R` -- Calculate and plot kinship matrix

## Common garden
`GDD_Calculations.R` -- Daymetr calcualtions of growing degree days at 10 degrees C
`androclimatena.R`-- Extracing ClimateNA data
`heritability.R` -- Estimate narrow-sense heritability for each trait with MegaLMM
`traits_pca.R` -- Run PCA on standarized traits
`gmatrix.R` -- Estimate G matrix with MegaLMM
`linear_mod_all.R` -- linear mixed models (Model 1 and 2)
`transfer_distance.R` -- Estimate effect of transfer distance on growth
