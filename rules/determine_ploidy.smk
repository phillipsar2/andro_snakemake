#  https://github.com/KamilSJaron/smudgeplot

## nQuire approach ##

rule calc_coverage:
    input:
        "data/interm/mark_dups/{unknown}.dedup.bam"
    output:
        f = "data/nQuire/{unknown}.bin",
        d = "data/nQuire/{unknown}_denoised.bin",
        hist = "data/nQuire/{unknown}_denoised.hist"
    params:
        "data/nQuire/{unknown}"
    shell:
        """
        ~/toolsfordayz/nQuire/nQuire create -q 30 -b {input} -o {params}
        ~/toolsfordayz/nQuire/nQuire denoise {output.f} -o {params}_denoised
        ~/toolsfordayz/nQuire/nQuire histo {output.d} > {output.hist}
        """

# > denoised histograms should be manually inspected before continuing. Inspect for % sites kept and pattern of hist.

rule run_GMM:
    input:
        expand("data/nQuire/{unknown}_denoised.bin", unknown = UNKNOWN)
    output:
        "data/nQuire/nquire_results_denoised.txt"
    threads:
        8
    shell:
        """
        ~/toolsfordayz/nQuire/nQuire lrdmodel -t {threads} {input} > {output}
        """
        
