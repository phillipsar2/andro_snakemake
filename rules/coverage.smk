
#SAMPLE = glob_wildcards("data/sorted_bam/{sample}.merge.bam").sample

# Calculate coverage of each bam file
rule bed_coverage:
    input:
        "data/sorted_bam/{sample}.merge.bam",
    output:
        "data/coverage/{sample}.cov",
    shell:
        "samtools depth {input} | "
        """awk '{{sum+=$3;sumsq+=$3*$3}} END {{ print sum/592152502,",",sqrt(sumsq/592152502 - (sum/592152502)**2)}}' > {output}"""
