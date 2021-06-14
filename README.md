# Population genetics of  Andropogon gerardii

A snakemake workflow for all of the things looking at origins of polyploidy.

#1 Alignment
- Short read data was verified as A. gerardi. Any data that had mixed identiy as a different species or identified entirely as another species was not alinged.
- reads from different lanes were merged into a single file (JGI individuals)
- pair-end reads were aligned with bwa-mem2/2.0
