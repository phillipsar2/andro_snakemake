# The consequences of polyploidy in adaptation of a dominant prairie grass

This repo details the bioinformatic methods for {doi}. 

The majority of the methods are encapsulated in a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow, but scripts related to analysis of the common garden data and some population genetics calculations are in stand-alone R scripts.

## Project organization
<pre>
├── README.md <- The top-level README.md for Andropogon or polyploidy enthusiasts that want to replicate these analyses.   
├── rules
|   ├── calling.smk
|   ├── determine_ploidy.smk
|   ├── filtering.smk
|   ├── mapping.smk
|   ├── pop_struc.smk
|   ├── pop_struc_highcov.smk
|   └── procesam.smk
├─  environment.yml 
├─  scripts
│   ├── README.md
│   ├── filtering <- custom scripts for variant filtering
│   ├── common_garden <- scripts related to linear mixed models & phenotype data
│   └── popgen <- analysis and plotting scripts for pop gen analyses
├── data <- each subfolder has a similar structure as angsd
│   ├── angsd <- intermediate and final files for analyses run with ANGSD
│   |	├── 6x_subsample <- files related to only 6x samples
│   |   ├── highcov <- files related to all high-coverage samples
│   |   └── lowcov <- files related to all low-coverage samples
│   ├── final_bams <- merged, deduplicated bam files with known ploidy
│   ├── genome <- reference genome
│   ├── interm  <- intermediate files in read mapping and SNP calling
│   ├── nQuire
│   └── structure
├── reports <- Generated analyses as HTML, PDF, or .txt.
├── qc <- quality check output for raw data
├── Snakefile
├── config.py
├── submit.json
└── submit.sh
</pre>


## Overview of the process
For more details, see the publication.

### 1. Alignment
- Short read data was verified as A. gerardi. Any data that had mixed identiy as a different species or identified entirely as another species was not alinged.
- Reads from different lanes were merged into a single file (JGI individuals)
- Pair-end reads were aligned with bwa-mem2/2.0

### 2. Post-alignment processing
- BAMs were sorted, read groups were added, and duplicates were marked with GATK and samtools
- Low-coverage (< 5X) BAMs from the same genotype were merged together (PanAnd & UCD) with samtools merge and read groups were re-added with gatk. Bams that were merged together were moved to the folder `data/interm/raw-unmerged`
- Low quality individuals (< 0.5X) were removed from the low-coverage data set and moved into the `data/interm/low-coverage` folder
- High coverage bams (from JGI) from genotypes that didn't have seperately sequenced low coverage data (PanAnd & UCD) were subsampled down to 2-4X coverage (`samtools view -s`) and placed in `data/final_bams/lowcov`
- All remaining low-coverage bams (PanAnd & UCD) that were not discarded were moved to `data/final_bams/lowcov`
- All remaining high-coverage bams (JGI) that were not discarded were moved to `data/final_bams/highcov`

### 3. Determine ploidy
- nQuire was utilized to determine ploidy for all high-coverage Andropogon individuals. \
I attempted to cluster the normalized maximum log-likelihood with mclust5 in R but it identified 8 and 9 clusters. \
Instead, I plotted the normalized values and colored the points but known/unknown ploidy to identify groups. Ploidy was identified for 16 of 17 unknown individuals.
- If genome size could not be determined or inferred for an individual, the bam was moved to the folder `data/interm/unknown-ploidy`

### 4. Identify good sites in SNP filtering
- Genotypes are called with bcftools mpileup & bcftools call for 2 sets of individuals: all high coverage (JGI) and all low coverage
- SNPs are extracted using GATK. VCFs are converted to table format with GATK to examine quality distributions and set hard filters.
    - Low coverage SNPs (unfiltered): 463,763,393
    - High coverage SNPs (unfiltered): 480,857,741
- SNPs are hard filtered using GATK.
    - Lowcov: QUAL >= 30, MQ >= 30, biallelic
    - Highcov: QUAL >= 30, MQ >= 30, biallelic
- SNPs are filtered for genotype depth and missingness using `scripts/filtering/genoDPfilter.R`
    - Lowcov: min = 1, max = `qpois(p = 0.99)`, < 20% missing data
    - Highcov: min = 8, max = `qpois(p = 0.99)`, < 20% missing data
- Final filtered SNP counts:
    - Lowcov (all): 11,707,655
    - Highcov: 102,242,020

### 5. Single-read genotypes (low-coverage only)
- Single-read genotypes were extracted using ANGSD for all 11,707,655 positions
- 100k SNPs were randomly selected using shuf

### 6. PCA
#### High cov SNPs:
- ANGSD was utilized to estimate GLs and then run a PCA (PCangsd) as all individuals are 6x
- Genotype likelihoods were estimated directly from the BAMs using the Samtools method (`-GL 1`) assuming the reference allele is the major allele.
- Reads with multiple best hits (`-uniqueOnly 1`) and flags above 255 (`-remove_bads 1`) were removed and only proper pairs were included.
- Sites with a mapping quality below 30 (`-minMapQ 30`), a minimum base quality score of 30 (`-minQ 30`), \
and a total sequencing depth below 8 and above 70 (`-doCounts 1 -setMinDepthInd 8 -setMaxDepthInd 70`) were dropped. 
- Genotype likelihoods were output as a beagle file (`-doGlf 2`).
- 30K sites were randomly grabbed from the beagle file using shuf.
- The PCA was ran with PCangsd using default settings. 

#### Low cov SNPs:
- The 100k SNPs with single-read genotypes were used to run a PCA in ANGSD.
    - PCA was also run with 50k random SNPs.
- Results of the PCA were plotted with `scripts/popgen/singleread_pca.R`

### 7. Kinship matrix
- A single read was randomly sampled (`-doIBS`) for each sites with < 20% missing data. Then, 100k random SNPs were \
selected for the kinship matrix with `shuf`.
- The single-read genotypes were used to estimate the kinshp matrix in a custom script (`scripts/popgen/kinship_matrix.R`)

### 8. STRUCTURE
- 100k random SNPs with < 20% missing data used for STRUCTURE

### 9. Inbreeding coefficient (Fis)
- Estimated only in the 6x genotypes with high coverage data using `ngsF` from ngstools.
- Genotype likelihoods were estimated and SNPs were called in ANGSD using the above SNP filtering criteria.
- Approximate values were generated independtly 3 times and final values were run 3 times using the approximate values as priors.

### 10. Thetas
- Estiamted only in the 6x genotypes utilizing the subsampled bams.
- The SAF, SFS, and thetas were estimated for (1) each population with more than 3 genotypes and (2) the West and East\
genetic group. Populations with more than 3 genotypes were randomly subsampled to 3. The West and East groups were run with their full sample sizes.


## Data availability
The short-read WGS data is available on NCBI Sequence Read Archive (SRA) under BioProject PRJNA1109389.
Supplementary data including genotype metadata, raw phenotype data, and flow cytometry data is available on [Dryad](https://doi.org/10.5061/dryad.gxd2547v1).
The *A. gerardi* reference genome is available on Phytozome under genome IDs 784 and 783. 
