# Population genetics of  Andropogon gerardii

A snakemake workflow for all of the things looking at origins of polyploidy.

#1 Alignment
- Short read data was verified as A. gerardi. Any data that had mixed identiy as a different species or identified entirely as another species was not alinged.
- reads from different lanes were merged into a single file (JGI individuals)
- pair-end reads were aligned with bwa-mem2/2.0

#2 Post-alignment processing
- BAMs were sorted, read groups were added, and duplicates were marked with GATK and samtools
- Low-coverage (< 5X) BAMs from the same genotype were merged together (PanAnd) with samtools merge and read groups were re-added with gatk
- Low quality individuals (< 0.5X) were removed from the low-coverage data set and moved into the data/interm/low-coverage folder

