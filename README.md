# Population genetics of  Andropogon gerardii

A snakemake workflow for all of the things looking at origins of polyploidy.

#1 Alignment
- Short read data was verified as A. gerardi. Any data that had mixed identiy as a different species or identified entirely as another species was not alinged.
- reads from different lanes were merged into a single file (JGI individuals)
- pair-end reads were aligned with bwa-mem2/2.0

#2 Post-alignment processing
- BAMs were sorted, read groups were added, and duplicates were marked with GATK and samtools
- Low-coverage (< 5X) BAMs from the same genotype were merged together (PanAnd & UCD) with samtools merge and read groups were re-added with gatk. Bams that were merged together were moved to the folder data/interm/raw-unmerged
- Low quality individuals (< 0.5X) were removed from the low-coverage data set and moved into the data/interm/low-coverage folder
- High coverage bams (from JGI) from genotypes that didn't have seperately sequenced low coverage data (PanAnd & UCD) were subsampled down to 2-4X coverage (samtools view -s) and placed in data/final_bams/lowcov
- All remaining low-coverage bams (PanAnd & UCD) that were not discarded were moved to data/final_bams/lowcov
- All remaining high-coverage bams (JGI) that were not discarded were moved to data/final_bams/highcov

#3 Determine ploidy
- nQuire was utilized to determine ploidy for all high-coverage Andropogon individuals. \
I attempted to cluster the normalized maximum log-likelihood with mclust5 in R but it identified 8 and 9 clusters. \
Instead, I plotted the normalized values and colored the points but known/unknown ploidy to identify groups. \
Scripts are local in ~/Andropogon/nQuire. Ploidy was identified for 16 of 17 unknown individuals.
- If genome size could not be determined or inferred for an individual, the bam was moved to the folder data/interm/unknown-ploidy

#4 Identify good sites in SNP filtering
- Genotypes are called with bcftools mpileup & bcftools call for 2 sets of individuals: all high coverage (JGI) and all low covera$
- SNPs are extracted using GATK. VCFs are converted to table format with GATK to examine quality distributions and set hard filter$
        -- Low coverage SNPs (unfiltered): 463,763,393
        -- High coverage SNPs (unfiltered): 480,857,741
- SNPs are hard filtered using GATK.
        -- Lowcov: QUAL >= 30, MQ >= 30, biallelic
        -- Highcov: QUAL >= 30, MQ >= 30, biallelic
- SNPs are filtered for genotype depth and missingness using a custom script (filter written by Mitra Melon)
        -- Lowcov: min = 1, max = qpois(p = 0.99), < 20% missing data
        -- Highcov: min = 8, max = qpois(p = 0.99), < 20% missing data
- Final filtered SNP counts:
        -- Lowcov: 11,707,655
        -- Highcov: 102,242,020

#5 Single-read genotypes (low-coverage only)
- Single-read genotypes were extracted using ANGSD for all 11,707,655 positions
- 50k SNPs were randomly selected using shuf

#6 PCA
High cov SNPs:
- ANGSD was utilized to estimate GLs and then run a PCA (PCangsd) as all individuals are 6x
- Genotype likelihoods were estimated directly from the BAMs using the Samtools method (-GL 1) assuming the reference allele is th$
- Reads with multiple best hits (-uniqueOnly 1) and flags above 255 ( -remove_bads 1) were removed and only proper pairs were incl$
- Sites with a mapping quality below 30 (-minMapQ 30), a minimum base quality score of 30 (-minQ 30), \
and a total sequencing depth below 8 and above 70 (-doCounts 1 -setMinDepthInd 8 -setMaxDepthInd 70) were dropped. The maximum num$
- Genotype likelihoods were output as a beagle file (-doGlf 2).
- 30K sites were randomly grabbed from the beagle file using shuf.
- The PCA was ran with PCangsd using default settings. The PCA was plotted in R with ggplot2.

Low cov SNPs:
- The 50k SNPs with single-read genotypes were used to run a PCA in ANGSD. Single-read gentoypes are ploidy-neutral.
	-- PCA was also run with 20k random SNPs.
	-- I should run multiple iterations of the 50k SNPs PCA
- Results of the PCA were plotted with a custom script in R.

#7 Kinship matrix
- Kinship matrix requires sites without any missing data. Sites without any missing data for the common garden individuals \
were extracted with ANGSD. The input was all of the sites with 20% or less missing data and these sites were further subsetted by \
requiring sites have data for all 79 individuals (-minInd 79). A single read was randomly selected for each site (`-doIBS`).
- The single-read genotypes were used to estimate the kinshp matrix in a custom script (`kinship_matrix.R`)
- Number of SNPs used to generate kinship matrix: 28,633

# --------


#6 Estimating genotype probabilities
Low coverage SNPs:
- The vcfs are split by ploidy using bcftools view
- Using a custom R script that utilizes vcfR (scripts/vcf2ADmatrix.R), I extracted the matrices for total read count and alternate read count
- I ran EBG independently on each ploidy to calculate phred-scaled genotype probabilites (1000 interactions; error = 0.01)
- The phred-scaled genotype likelihood EBG output was converted to a matrix and the files for each ploidy were combined using a custom script (scripts/ebg2glmatrix_2ploidy.R)


#7 ENTROPY
Low cov SNPs:
- 10k random SNPs were selected and the subsequent GL matrix was converted to mpgl format using a custom script (scripts/glmat2mpgl.R)
- ENTROPY was run with default burn-in and MCMC steps. Every 5th step was recorded after burn-in (-t 5) and ENTROPY was run for k=c(2:14)


#8 PCA
High cov SNPs:
- ANGSD was utilized to estimate GLs and then run a PCA (PCangsd) as all individuals are 6x
- Genotype likelihoods were estimated directly from the BAMs using the Samtools method (-GL 1) assuming the reference allele is the major allele (-doMajorMinor 4). 
- Reads with multiple best hits (-uniqueOnly 1) and flags above 255 ( -remove_bads 1) were removed and only proper pairs were included (-only_proper_pairs 1).
- Sites with a mapping quality below 30 (-minMapQ 30), a minimum base quality score of 30 (-minQ 30), \
and a total sequencing depth below 8 and above 70 (-doCounts 1 -setMinDepthInd 8 -setMaxDepthInd 70) were dropped. The maximum number of missing genotypes allowed at a site was 8 (-minInd 8).
- Genotype likelihoods were output as a beagle file (-doGlf 2).
- 30K sites were randomly grabbed from the beagle file using shuf.
- The PCA was ran with PCangsd using default settings. The PCA was plotted in R with ggplot2.

Low cov SNPs:
- The GL table was filtered for X% missing data. This left a total of X SNPs.
- 20K SNPs were randomly selected using 'shuf'.
- Individual allele frequencies were calculated in R using a custom script. 
- Allele frequencies were standardized seperately for each ploidy using a Patterson's standardization.\
 This standardization method properly accounts for the difference in variation between the two cytotypes.
- A PCA was ran in R using 'prcomp' on the standardized allele frequencies. Results were plotted with ggplot2.
