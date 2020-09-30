#!/bin/bash

snakemake --jobs 100 --use-conda \
--rerun-incomplete \
--latency-wait 120 \
--cluster-config submit.json \
--cluster "sbatch --mem {cluster.mem} -J {cluster.name} --time {cluster.time} -p {cluster.p}"

