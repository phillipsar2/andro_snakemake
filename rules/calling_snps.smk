# Calling SNPs using GATK (probably maybe)
# Follows GATK best practices


# Samtools mpileup 
#rule samtools_mpileup:
#    input:
#        bam = "/group/jrigrp6/andropogon_sequence/data/interm/mark_dups/{sample}.dedup.bam.gz", 
#        ref = "/group/jrigrp6/andropogon_sequence/data/genome/ANDRO_contigs1278.fasta",
#    output:
#        temp("data/mpileup/{sample}.mpileup"),
##        "data/bcfs/{sample}.bcf",
#    run:
#        shell("zcat {input.bam} | bcftools mpileup -d100000 -a AD -Ou -Bf {input.ref} - -o {output}")
##        bcftools call -mv -Ob -o {output} - ")

#rule bcftools_call:
#    input:
#        "data/mpileup/{sample}.mpileup",
#    output:
#        "data/bcfs/{sample}.vcf",
#    run:
#        shell("bcftools call -mv -Ov {input} > {output}")

# bgzip vcfs
rule bgzip_vcfs:
    input:
        "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz"
#        "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf"
    output:
#        "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz",
        "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz.tbi"
    params:
        sample = "{sample}" 
    run:
#        shell("gunzip {input}")
#        shell("bgzip /group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{params.sample}.g.vcf")
#        shell("bgzip {input}")
#        shell("tabix -p vcf /group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{params.sample}.g.vcf.gz")
        shell("tabix -p vcf {input}")

# index bgzipped vcfs
#rule bcftools_index: 
#    input: 
#        "data/bcfs/{sample}.new.vcf.gz",
#    output:
#        index = "data/bcfs/{sample}.new.vcf.gz.csi",
#    run:
#        shell("bcftools index -c {input}")

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
        bam = "/group/jrigrp10/andropogon_shortreads/oldgenome_bams/{sample}.dedup.bam"
    output:
#        outdir = "data/gvcf/{sample}.g.vcf",
        outdir = "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz"
    params:
        regions = "data/genome/contigs.list",
        sample = "{sample}"
    run:
        # -G specifies annotations to apply variant calls
        # --emitRefConfidence keeps invariant sites
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output /group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{params.sample}.g.vcf \
        --reference {input.ref} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -L {params.regions} \
        -ploidy 6 \
        -ERC BP_RESOLUTION")
        shell("bgzip /group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{params.sample}.g.vcf")

#rule combine_gvcfs:
#    input:
#        ref = config.ref,
#        gvcf = expand("/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz", sample=SAMPLE)
#    output:
#        "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/combined.g.vcf.gz"
#    params:
#        files = lambda wildcards, input: " --variant ".join(input.gvcf)
#    run:
#        shell("gatk CombineGVCFs \
#        -R {input.ref} \
#        --variant {params.files} \
#        -O {output} \
#        -RF MappingQualityReadFilter --minimum-mapping-quality 20 \
#        ")

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
#        directory("data/interm/combined_database/{count}"),
    # Join multiple files with -V.join trick - takes all the inputs and inputs them after a -V in the shell command
#    params:
#        files = lambda wildcards, input: " -V ".join(input),
        ##dir = "data/combined_database/gvcf",
#        int = "data/data/processed/scattered_intervals/{count}-scattered.intervals",
#        tmp = "/scratch/aphillip/genomicsdbimport",
#        map = "data/processed/sample_map",
    # Must specify workspace directory with --genomicsdb-workspace-path (must be empty)
    # Have to specify intervals to work with - probably use contigs?
#    run:
#        shell("mkdir -p {params.tmp}")
#        shell("gatk GenomicsDBImport \
#        -V {params.files} \
#        --genomicsdb-workspace-path {output} \
#        --batch-size 50 \
#        --reader-threads 8 \
#        --sample-name-map {params.map} \
#        --intervals {params.int} \
#        --tmp-dir {params.tmp}")
#        shell("rm -rf {params.tmp}")

# Joint-calling SNPs - can be rerun at any point whe samples are added to the cohort
#rule joint_genotyping:
#    input:
        #expand("data/gvcf/{sample}.g.vcf", sample=SAMPLE),
#        dir = directory("data/interm/combined_database"),
#        ref = config.ref,
#    output:
#        "data/raw/vcf/andro.raw.snps.indels.vcf",
#    params:
        #files = lambda wildcards, input: " -V ".join(input),
        #int = "data/genome/list1278",
#        db = "gendb://data/interm/combined_database",
#    run:
        #shell("gatk GenomicsDBImport \
        #-V {params.db} \
        #--genomicsdb-workspace-path data/gvcf/combined \
        #--batch-size 100 \
        #--intervals {params.int}")
#        shell("gatk GenotypeGVCFs \
#        -R {input.ref} \
#        -V {params.db} \
#        -new-qual \
#        -G StandardAnnotation \
#        -G AS_StandardAnnotation \
#       --include-non-variant-sites \
#        -O {output}")
