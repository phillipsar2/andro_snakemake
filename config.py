# Bed_coverage
#bedcov_in = "data/sorted_bam/{induv}.merge.bam"
#bedcov_out = "data/coverage/{induv}.cov"


# Add read groups
addrg_in = "data/mergensort/{sample}.merge.sorted.bam"

# MarkDuplicates
#mark_in = "data/sorted_bam/{sample}.merge.bam"
##mark_in = "/home/aphillip/Andro_snakemake/data/sorted_reads/{sample}.sorted.bam"
mark_in = "data/interm/addrg/{sample}.rg.bam"


# Prepare the genome
ref = "data/genome/ANDRO_contigs1278.fasta"
ref_dict = "data/genome/ANDRO_contigs1278.dict"
ref_fai = "data/genome/ANDRO_contigs1278.fasta.fai"
#ref = "/home/aphillip/Andropogoneae/Andro_sequences/ANDRO.asmAll.contigs.fasta"
#ref_dict = "/home/aphillip/Andropogoneae/Andro_sequences/ANDRO.asmAll.contigs.dict"
#ref_fai = "/home/aphillip/Andropogoneae/Andro_sequences/ANDRO.asmAll.contigs.fasta.fai"


# HaplotypeCaller
#haplo_in = "data/sorted_bam/{sample}.merge.bam"
haplo_in = "data/interm/mark_dups/{sample}.dedup.bam"


# Joint genotyping
#joint_out = "data/raw/vcf/andro1278.raw.vcf"
joint_out = "data/raw/vcf/test.raw.vcf"
