# Population genetics of  Andropogon gerardii

A snakemake workflow for all of the things looking at origins of polyploidy.

#1 Alignment
- Short read data was verified as A. gerardi. Any data that had mixed identiy as a different species or identified entirely as another species was not alinged.
- reads from different lanes were merged into a single file (JGI individuals)
- pair-end reads were aligned with bwa-mem2/2.0

#2 Post-alignment processing
- BAMs were sorted, read groups were added, and duplicates were marked with GATK and samtools
- Low-coverage (< 5X) BAMs from the same genotype were merged together (PanAnd) with samtools merge and read groups were re-added with gatk. Bams that were merged together were moved to the folder data/interm/raw-unmerged
- Low quality individuals (< 0.5X) were removed from the low-coverage data set and moved into the data/interm/low-coverage folder
- If genome size could not be determined or inferred for an individual, the bam was moved to the folder data/interm/unknown-ploidy
- High coverage bams (from JGI) were subsampled down to 2-4X coverage (samtools view -s) and placed in data/final_bams/lowcov
- All remaining low-coverage bams (PanAnd) that were not discarded were moved to data/final_bams/lowcov
- All remaining high-coverage bams (JGI) that were not discarded were moved to data/final_bams/highcov

#3 Determine ploidy
- nQuire was utilized to determine ploidy for all high-coverage Andropogon individuals. \
I attempted to cluster the normalized maximum log-likelihood with mclust5 in R but it identified 8 and 9 clusters. \
Instead, I plotted the normalized values and colored the points but known/unknown ploidy to identify groups. \
Scripts are local in ~/Andropogon/nQuire. Ploidy was identified for 16 of 17 unknown individuals. 

# SNP calling
- Genotypes are called with bcftools mpileup for 2 sets of individuals: all high coverage (JGI) and all low coverage (subsampled + PanAnd)
- Using a custom R script that utilizes vcfR, I extract the matrices for total read count and alternate read count
- Then I use EBG to calculate phred-scaled genotype probabilites
- The phred-scaled genotype probs are converted to mpgl format using a custom script
