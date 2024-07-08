# The consequences of polyploidy in adaptation of a dominant prairie grass

This repo details the bioinformatic methods for {doi}. 

The majority of the methods are encapsulated in a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow, but scripts related to analysis of the common garden data and some population genetics calculations are in stand-alone R scripts.

## Project organization
<pre>
├── README.md <- The top-level README.md for Andropogon or polyploidy enthusiasts that want to replicate these analyses.   
├── rules    
|   ├── angsd.smk   
|   ├── calling.smk
|   ├── determine_ploidy.smk
|   ├── filtering.smk
|   ├── mapping.smk
|   ├──    pop_struc.smk
|   ├──    pop_struc_highcov.smk   
|   └── procesam.smk   
├─environment.yml   
├── scripts   
│   ├──allelebalance_filter.sh   
│   ├─allelebalance_filter.R   
│   ├── process_trees.R <- Script for tree figures.   
│   └── poa_pca_pi.R <- Script for PCA and pi figures.    
├── data   
│   ├── raw <- The original WGS data dump.   
│   ├── interm  <- Intermediate data that has been transformed.   
│   ├── processed <- The final datasets.
│   ├── vcf <- The unfiltered vcfs.   
│   ├── gene <- Genes downloaded from NCBI for phylogenetic analyses.    
│   └── genome <- The reference genome.   
├── reports <- Generated analyses as HTML, PDF, or .txt.    
├── Snakefile   
├── config.py   
├── submit.json   
└── submit.sh   
</pre>


## Overview of the process
For more details, see the publication.

### 1. Alignment
- Short read data was verified as A. gerardi. Any data that had mixed identiy as a different species or identified entirely as another species was not alinged.
- reads from different lanes were merged into a single file (JGI individuals)
- pair-end reads were aligned with bwa-mem2/2.0

## 2. Post-alignment processing
- BAMs were sorted, read groups were added, and duplicates were marked with GATK and samtools
- Low-coverage (< 5X) BAMs from the same genotype were merged together (PanAnd & UCD) with samtools merge and read groups were re-added with gatk. Bams that were merged together were moved to the folder data/interm/raw-unmerged
- Low quality individuals (< 0.5X) were removed from the low-coverage data set and moved into the data/interm/low-coverage folder
- High coverage bams (from JGI) from genotypes that didn't have seperately sequenced low coverage data (PanAnd & UCD) were subsampled down to 2-4X coverage (samtools view -s) and placed in data/final_bams/lowcov
- All remaining low-coverage bams (PanAnd & UCD) that were not discarded were moved to `data/final_bams/lowcov`
- All remaining high-coverage bams (JGI) that were not discarded were moved to `data/final_bams/highcov`

## 3. Determine ploidy
- nQuire was utilized to determine ploidy for all high-coverage Andropogon individuals. \
I attempted to cluster the normalized maximum log-likelihood with mclust5 in R but it identified 8 and 9 clusters. \
Instead, I plotted the normalized values and colored the points but known/unknown ploidy to identify groups. \
Scripts are local in `~/Andropogon/nQuire`. Ploidy was identified for 16 of 17 unknown individuals.
- If genome size could not be determined or inferred for an individual, the bam was moved to the folder data/interm/unknown-ploidy

## 4. Identify good sites in SNP filtering
- Genotypes are called with bcftools mpileup & bcftools call for 2 sets of individuals: all high coverage (JGI) and all low covera$
- SNPs are extracted using GATK. VCFs are converted to table format with GATK to examine quality distributions and set hard filter$
    -Low coverage SNPs (unfiltered): 463,763,393
    - High coverage SNPs (unfiltered): 480,857,741
- SNPs are hard filtered using GATK.
    - Lowcov: QUAL >= 30, MQ >= 30, biallelic
    - Highcov: QUAL >= 30, MQ >= 30, biallelic
- SNPs are filtered for genotype depth and missingness using a custom script (filter written by Mitra Melon)
    - Lowcov: min = 1, max = `qpois(p = 0.99)`, < 20% missing data
    - Highcov: min = 8, max = `qpois(p = 0.99)`, < 20% missing data
- Final filtered SNP counts:
    - Lowcov: 11,707,655
    - Highcov: 102,242,020

## 5. Single-read genotypes (low-coverage only)
- Single-read genotypes were extracted using ANGSD for all 11,707,655 positions
- 50k SNPs were randomly selected using shuf

## 6. PCA
### High cov SNPs:
- ANGSD was utilized to estimate GLs and then run a PCA (PCangsd) as all individuals are 6x
- Genotype likelihoods were estimated directly from the BAMs using the Samtools method (`-GL 1`) assuming the reference allele is the major allele.
- Reads with multiple best hits (`-uniqueOnly 1`) and flags above 255 (`-remove_bads 1`) were removed and only proper pairs were included.
- Sites with a mapping quality below 30 (`-minMapQ 30`), a minimum base quality score of 30 (`-minQ 30`), \
and a total sequencing depth below 8 and above 70 (`-doCounts 1 -setMinDepthInd 8 -setMaxDepthInd 70`) were dropped. 
- Genotype likelihoods were output as a beagle file (`-doGlf 2`).
- 30K sites were randomly grabbed from the beagle file using shuf.
- The PCA was ran with PCangsd using default settings. The PCA was plotted in R with ggplot2.

### Low cov SNPs:
- The 50k SNPs with single-read genotypes were used to run a PCA in ANGSD. Single-read gentoypes are ploidy-neutral.
    - PCA was also run with 20k random SNPs.
    - I should run multiple iterations of the 50k SNPs PCA
- Results of the PCA were plotted with a custom script in R.
- Seperate PCAs were ran for the CG and All Andropogon

## 7. Kinship matrix
- A single read was randomly sampled (`-doIBS`) for each sites with < 20% missing data. Then, 50k (CG) or 100k (All Andro) random SNPs were \
selected for the kinship matrix.
- The single-read genotypes were used to estimate the kinshp matrix in a custom script (`kinship_matrix.R`)
- To estimate the diagonal of the kinship matrix (aka inbreeding), only sites with genotype depth > 2 were included. A single read was sampled\
twice at each site and both draws were used to calculate inbreeding (the probability of sharing an allele IBD with yourself). After drawing \
the reads twice, I randomly selected the same set of 100k sites from both runs. These 100k snps were used to calculate inbreeding as the \
probability of sharing an allele IBD with yourself.

## 8. STRUCTURE
### Common garden individuals:
- 50k random SNPs with < 20% missing data were used for STRUCTURE

### All individuals
- 100k random SNPs with < 20% missing data used for STRUCTURE

## 9. Inbreeding coefficient
- Estimated only in the 6x genotypes
- Bams were subsampled to less than 1.5X coverage, if they weren't already below that threshold, (`samtools -view -b -s`) in order to reduce\
 bias caused by coverage. Coverage percentages to subsample to are stored in `./6x_bams_to_subset.tsv`.
- Genotype likelihoods were estimated and SNPs were called in ANGSD using the above SNP filtering criteria.
- Inbreeding coefficient was estimated with `ngsF` and plotted in R.

## 10. Aneuploidy detection
- lowcov BAM files were converted to bed files using `bedtools bamtobed`
- average coverage was calculated in 1 Mbp windows across the genome for each genotype using `bedtools coverage`
-- Coverage was not calculated for the following genotypes due to data quality:`AN20N023`,`AN20NSCR000363`,`AN20NSCR000409`

## 11. Thetas
- Estiamted only in the 6x genotypes utilizing the subsampled bams.
- The SAF, SFS, and thetas were estimated for (1) each population with more than 3 genotypes and (2) the West and East\
genetic group. Populations with more than 3 genotypes were randomly subsampled to 3. The West and East groups were run will their full sample sizes.

# ---- Other ----


#6 Estimating genotype probabilities
##Low coverage SNPs:
- The vcfs are split by ploidy using bcftools view
- Using a custom R script that utilizes vcfR (scripts/vcf2ADmatrix.R), I extracted the matrices for total read count and alternate read count
- I ran EBG independently on each ploidy to calculate phred-scaled genotype probabilites (1000 interactions; error = 0.01)
- The phred-scaled genotype likelihood EBG output was converted to a matrix and the files for each ploidy were combined using a custom script (scripts/ebg2glmatrix_2ploidy.R)

#8 PCA
##High cov SNPs:
- ANGSD was utilized to estimate GLs and then run a PCA (PCangsd) as all individuals are 6x
- Genotype likelihoods were estimated directly from the BAMs using the Samtools method (-GL 1) assuming the reference allele is the major allele (-doMajorMinor 4). 
- Reads with multiple best hits (-uniqueOnly 1) and flags above 255 ( -remove_bads 1) were removed and only proper pairs were included (-only_proper_pairs 1).
- Sites with a mapping quality below 30 (-minMapQ 30), a minimum base quality score of 30 (-minQ 30), \
and a total sequencing depth below 8 and above 70 (-doCounts 1 -setMinDepthInd 8 -setMaxDepthInd 70) were dropped. The maximum number of missing genotypes allowed at a site was 8 (-minInd 8).
- Genotype likelihoods were output as a beagle file (-doGlf 2).
- 30K sites were randomly grabbed from the beagle file using shuf.
- The PCA was ran with PCangsd using default settings. The PCA was plotted in R with ggplot2.

####
# Analyses that were unsuccessful
* ENTROPY
* Treemix
