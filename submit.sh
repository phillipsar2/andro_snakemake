#!/bin/bash

snakemake --jobs 120 --use-conda \
--rerun-incomplete \
--latency-wait 120 \
--cluster-config submit.json \
--cluster "sbatch --mem {cluster.mem} --time {cluster.time} -p {cluster.p}"

