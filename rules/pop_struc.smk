# Evaluating population structure with ENTROPY

# (1) Generate ENTROPY input files
rule ent_in:
    input:
        "data/ebg/lowcov/10k_lowcov-GL.txt"
    output:
        "data/entropy/lowcov/10k_lowcov.qk2inds.txt",
        "data/entropy/lowcov/10k_lowcov_pntest_meangl.txt",
        "data/entropy/lowcov/10k_lowcov.mpgl"
    params:
        pre = "data/entropy/lowcov/10k_lowcov"
    shell:
        """
        Rscript scripts/glmat2mpgl.R --gl {input} --out {params.pre}
        """

# (2) Run ENTROPY
# not sure if I can actually activate a conda env in a snakemake rule?? might have to run it sepearately outside snakemake
# run with default burn-in of 1000 (-b 1000) and 1000 MCMC steps (-l 1000); store every 5 step after burn-in (-t 5)
# recommends running three times
#rule entropy:
    input:
        mpgl = "data/entropy/lowcov/10k_lowcov.mpgl",
        ploidy = "data/entropy/lowcov/10k_lowcov.ploidy_inds.txt",
        init_admix = "data/entropy/lowcov/10k_lowcov.qk{k}inds.txt",
    output:
        "data/entropy/lowcov/10k_lowcov.{k}.{chain}.hdf5"
    params:
        k = {k}
    shell:
        """
        module load bio3 
        source activate entropy-2.0
        entropy -i {input.mpgl} \
        -n {input.ploidy} \
        -k {params.k} \
        -q {input.init_admix} \
        -t 5 \ 
        -o {output}
        """
