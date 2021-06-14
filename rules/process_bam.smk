# Proccessing merged bam files

#rule samtools_resort:
#    input:
#        "data/sorted_bam/{sample}.merge.bam"
#    output:
#        temp("data/mergensort/{sample}.merge.sorted.bam")
#    params:
#        "-m 16G"
#    run:
#        shell("samtools sort {input} > {output}")

#rule samtools_bamtosam:
#    input:
#        "data/mergensort/{sample}.merge.sorted.bam"
#    output:
#        temp("data/mergensort/{sample}.merge.sorted.sam")
#    run:
#        shell("samtools view -h -o {output} {input}")


#rule samtools_samtobam:
#    input:
#        "data/mergensort/{sample}.merge.sorted.sam"
#    output:
#        "data/mergensort/{sample}.NEW.bam"
#    run:
#        shell("samtools view -bS {input} > {output}")

rule add_rg:
    input:
#        "data/mergensort/{sample}.NEW.bam"
        "data/sorted_bam/{sample}.sorted.bam"
    output:
        bam = temp(touch("data/interm/addrg/{sample}.rg.bam"))
        #index = temp(touch("data/interm/addrg/{sample}.rg.bai"))
    params:
        tmp = "/scratch/aphillip/addrg/{sample}",
        sample = "{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk AddOrReplaceReadGroups \
        -I={input} \
        -O={output.bam} \
        -RGID=4 \
        -RGLB=lib1 \
        -RGPL=illumina \
        -RGPU=unit1 \
        -RGSM={params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")


rule mark_dups:
    input:
        config.mark_in
    output:
        # Make bam and index temporary
#        bam = "data/interm/mark_dups/{sample}.dedup.bam",
        bam = "/group/jrigrp10/andropogon_shortreads/oldgenome_bams/{sample}.dedup.bam",
        index = "/group/jrigrp10/andropogon_shortreads/oldgenome_bams/{sample}.dedup.bai",
#        index = "data/interm/mark_dups/{sample}.dedup.bai",
        metrics = "qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/scratch/aphillip/mark_dups/{sample}"
    run:
        # Create a scratch directory
        shell("mkdir -p {params.tmp}")
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={params.tmp}")
        # Remove scratch directory
        shell("rm -rf {params.tmp}")


# Quality metrics with qualimap

rule bamqc:
    input:
        "data/interm/mark_dups/{sample}.dedup.bam",
    output:
        "reports/bamqc/{sample}_stats/qualimapReport.html"
    params:
        dir = "reports/bamqc/{sample}_stats"
#    threads: 8
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 8 \
        -nr 100000 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=64G")


# Combine qualimap results with 'python magic' maybe

#rule multibamqc:
#    input:
#        all = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES)
#    output:
#        "reports/multisampleBamQcReport.html"
#    params:
#        outdir = "reports",
#        infile = "models/bamqc_list.txt"
#    run:
#        shell("find reports/bamqc -mindepth 1 -maxdepth 1 -type d | grep SamC > models/ALL.bamqclist.txt")
#        import pandas as pd
#        data = pd.read_csv("models/ALL.bamqclist.txt", sep = " ", header = None, names = ['filename'])
#        data['sample'] = data['filename'].str.split('.').str[0].str.split('/').str[2].str.split('_stats').str[0]
#        data = data.sort_values('sample', axis = 0, ascending = True)
#        data = data[['sample','filename']]
#        data.to_csv(r'models/bamqc_list.txt', header = None, index = None, sep = ' ', mode = 'a')
#        shell("qualimap multi-bamqc \
#        --data {params.infile} \
#        --paint-chromosome-limits \
#        -outdir {params.outdir} \
#        -outformat html \
#        --java-mem-size=40G")

# Realign indels

