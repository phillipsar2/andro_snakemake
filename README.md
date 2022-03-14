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

#4 SNP calling
- Genotypes are called with bcftools mpileup & bcftools call for 2 sets of individuals: all high coverage (JGI) and all low coverage (subsampled JGI, PanAnd, & UCD)

#5 SNP filtering
- SNPs are extracted using GATK. VCFs are converted to table format with GATK to examine quality distributions and set hard filtering cutoffs.
	-- Low coverage SNPs (unfiltered): 463,763,393
        -- High coverage SNPs (unfiltered): 480,857,741
- SNPs are hard filtered using GATK.
	-- Lowcov: QUAL >= 30, MQ >= 30, biallelic 
	-- Highcov:
- SNPs are filtered for genotype depth and missingness using a custom script (filter written by Mitra Melon)
	-- Lowcov: min = 1, max = qpois(p = 0.99), < 20% missing data
	-- Highcov:
- Final filtered SNP counts:
	-- Lowcov: 11,707,655
	-- Highcov:

#6 Estimating genotype probabilities
- Using a custom R script that utilizes vcfR, I extract the matrices for total read count and alternate read count
- Then I use EBG to calculate phred-scaled genotype probabilites
- The phred-scaled genotype probs are converted to mpgl format using a custom script

