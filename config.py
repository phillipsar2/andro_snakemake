# Bed_coverage
#bedcov_in = "data/sorted_bam/{induv}.merge.bam"
#bedcov_out = "data/coverage/{induv}.cov"


# Prepare the genome
ref = "data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta"

# HaplotypeCaller
#haplo_in = "data/sorted_bam/{sample}.merge.bam"
haplo_in = "data/interm/mark_dups/{sample}.dedup.bam"


# Joint genotyping
#joint_out = "data/raw/vcf/andro1278.raw.vcf"
joint_out = "data/raw/vcf/test.raw.vcf"


# Filter galliart vcf
vcf = "data/galliart_data/spipsforPCA.vcf"

