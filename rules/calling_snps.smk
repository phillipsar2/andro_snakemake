# Calling SNPs using GATK (probably maybe)
# Follows GATK best practices

# keep invariant sites http://seqanswers.com/forums/showthread.php?t=43006
## include --emitRefConfidence GVCF in haplotype HaplotypeCaller

rule haplotype_caller:
    input:
        ref = config.ref,
        bam = "data/sorted_bam/{sample}.merge.bam"
    output:
        gvcf = "data/gvcfs/{sample}.g.vcf",
        index = "data/gvcfs/{sample}.g.vcf.idx"
    run:
        # -G specifies annotations to apply variant calls
        # --emitRefConfidence keeps invariant sites
        shell("gatk HaplotypeCaller \
        -I {input.bam} \
        -O {output.gvcf} \
        -R {input.ref} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --emitRefConfidence GVCF")

# Combine gvcfs for joint variant calling
# Use trick from https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
rule combine_gvcfs:
    input:
        expand("data/gvcfs/{sample}.g.vcf", sample=SAMPLE)
    output:
        "data/combined_database/chr{chr}/vcfheader.vcf"
    # Join multiple files with -V.join trick - takes all the inputs and inputs them after a -V in the shell command
    params:
        files = lambda wildcards, input: " -V ".join(input),
        dir = "data/combined_database/chr{chr}",
        region = "{chr}"
    run:
        shell("gatk GenomicsDBImport \
        -V {params.files} \
        --genomicsdb-workspace-path {params.dir} \
        --intervals {params.region}")
