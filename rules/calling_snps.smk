# Calling SNPs using GATK (probably maybe)
# Follows GATK best practices


# Prep reference genome for GATK by indexing and creating a dictionary
rule genome_index:
    input:
        config.ref
    output:
        dict = config.ref_dict,
        fai = config.ref_fai
    run:
        shell("gatk CreateSequenceDictionary \
        -R {input} \
        -O {output.dict}")
        shell("samtools faidx {input}")

# keep invariant sites http://seqanswers.com/forums/showthread.php?t=43006
## include --emitRefConfidence GVCF in haplotype HaplotypeCaller
rule haplotype_caller:
    input:
        ref = config.ref,
        bam = config.haplo_in
    output:
        gvcf = "data/gvcf/{sample}.g.vcf",
        #index = "data/gvcf/{sample}.g.vcf.idx"
    run:
        # -G specifies annotations to apply variant calls
        # --emitRefConfidence keeps invariant sites
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output.gvcf} \
        --reference {input.ref} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -ERC GVCF")

#gatk HaplotypeCaller --input ~/Andro_snakemake/data/interm/INCQ_PCRfree_1-McKain_252-WiDiv_43_CCGCGGTT_Andropogon_gerardii_I842_L3.dedup.bam
# --output ~/Andro_snakemake/data/interm/INCQ_PCRfree_1-McKain_252-WiDiv_43_CCGCGGTT_Andropogon_gerardii_I842_L3.g.vcf
# --reference data/genome/ANDRO_contigs1278.fasta
# -G StandardAnnotation
# -G AS_StandardAnnotation
# -ERC GVCF

# Combine gvcfs for joint variant calling
# Use trick from https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
##### NOT RECOGNIZING THE DIRECTORY
#rule combine_gvcfs:
#    input:
#        expand("data/gvcf/{sample}.g.vcf", sample=SAMPLE),
#    output:
#        config.dir_out
    # Join multiple files with -V.join trick - takes all the inputs and inputs them after a -V in the shell command
#    params:
#        files = lambda wildcards, input: " -V ".join(input),
        ##dir = "data/combined_database/gvcf",
#        int = "data/genome/list1278"
    # Must specify workspace directory with --genomicsdb-workspace-path (must be empty)
    # Have to specify intervals to work with - probably use contigs?
#    run:
#        shell("gatk GenomicsDBImport \
#        -V {params.files} \
#        --genomicsdb-workspace-path {output} \
#        --batch-size 100 \
#        --intervals {params.int}")

# Joint-calling SNPs - can be rerun at any point whe samples are added to the cohort
rule joint_genotyping:
    input:
        expand("data/gvcf/{sample}.g.vcf", sample=SAMPLE),
        #dir = directory("data/gvcf/combined_dataset"),
        #ref = config.ref
    output:
        config.joint_out
    params:
        files = lambda wildcards, input: " -V ".join(input),
        int = "data/genome/list1278",
        db = "gendb://data/gvcf/combined"
    run:
        shell("gatk GenomicsDBImport \
        -V {params.files} \
        --genomicsdb-workspace-path data/gvcf/combined \
        --batch-size 100 \
        --intervals {params.int}")
        shell("gatk GenotypeGVCFs \
        -R config.ref \
        -V {params.db} \
        -new-qual \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --includeNonVariantSites\
        -O {output}")
