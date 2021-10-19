#  https://github.com/KamilSJaron/smudgeplot


# Calculate kmer frequences from the fasta files and then store in a histogram

rule extract_kmers:
    input:
        r1 = "/group/jrigrp10/andropogon_shortreads/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii.merge.R1.fastq.gz",
        r2 = "/group/jrigrp10/andropogon_shortreads/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii.merge.R2.fastq.gz"
    output:
        hist = "data/kmer/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii.hist"
    params:
        geno = "INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii",
        list = "INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii.FILES",
        tmp = "/scratch/aphillip/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii",
        kmcdb = "data/kmer/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii"
    rule:
        # kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
        # tmp = temp directory
        # kmcdb = output file name prefix for KMC database
        # shell("module load KMC")
        shell("mkdir {params.tmp}")
        shell("ls {params.geno}*.fastq.gz > {params.list}")
        shell("kmc -k21 -t{threads} -m64 -ci1 -cs10000 {params.list} {params.kmcdb} {params.tmp}")
        shell("kmc_tools transform {params.kmcdb} histogram {output.hist} -cx10000")
        shell("rm -r {params.tmp}")


# extract genomic kmers within a coverage threshold then generate smudgeplot

rule smudgeplot:
    input:
        "data/kmer/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii.hist"
    output:
        "data/smudgeplot/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii_smudgeplot.png"
    params:
        kmcdb = "data/kmer/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii"
        prefix = "data/smudgeplot/INER_PCRfree_Loretta_Web-2_16_ACGCACCT_Andropogon_gerardii"
    shell:
        """
        L=$(smudgeplot.py cutoff {input} L)
        U=$(smudgeplot.py cutoff {input} U)
        echo $L $U # these need to be sane values
        # L should be like 20 - 200
        # U should be like 500 - 300
        kmc_tools transform {params.kmcdb} -ci"$L" -cx"$U" reduce {params.kmcdb}_L"$L"_U"$U"
        smudge_pairs {params.kmcdb}_L"$L"_U"$U" {params.kmcdb}_L"$L"_U"$U"_coverages.tsv {params.kmcdb}_L"$L"_U"$U"_pairs.tsv > {params.kmcdb}_familysizes.tsv
        smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv -o {params.prefix}
        """

