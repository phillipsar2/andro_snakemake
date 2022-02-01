# Bed_coverage
#bedcov_in = "data/sorted_bam/{induv}.merge.bam"
#bedcov_out = "data/coverage/{induv}.cov"


# Genome
ref = "data/genome/JGIgenome/v1/Andropogon_gerardii_var_Kellogg_1272_HAP1_V1_release/Andropogon_gerardii_var_Kellogg_1272/sequences/Andropogon_gerardii_var_Kellogg_1272.mainGenome.fasta"
contig_list = "data/genome/JGIgenome/v1/Andropogon_gerardii_var_Kellogg_1272_HAP1_V1_release/Andropogon_gerardii_var_Kellogg_1272/sequences/contig.list"

# Merge low-coverage bams
#bam_file = "bams_to_merge.tsv"

# Sets of finalized genotypes with their associated ploidy
highcov ="highcov_w_ploidy.csv"
sample_map = "highcov_subsample_map.txt"

lowcov = "allgenos_w_ploidy.csv"


# --- old stuff ---

# HaplotypeCaller
#haplo_in = "data/sorted_bam/{sample}.merge.bam"
haplo_in = "data/interm/mark_dups/{sample}.dedup.bam"


# Joint genotyping
#joint_out = "data/raw/vcf/andro1278.raw.vcf"
joint_out = "data/raw/vcf/test.raw.vcf"


# Filter galliart vcf
vcf = "data/galliart_data/spipsforPCA.vcf"

