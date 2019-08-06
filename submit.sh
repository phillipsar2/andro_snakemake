#!/bin/bash

snakemake --jobs 7 --use-conda \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -p {cluster.p} --time {cluster.time} -n {cluster.nCPUs}"

