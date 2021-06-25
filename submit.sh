#!/bin/bash

snakemake --jobs 50 --use-conda \
--rerun-incomplete \
--latency-wait 120 \
--cluster-config submit.json \
--cluster "sbatch --mem {cluster.mem} -J {cluster.name} --time {cluster.time} -p {cluster.p} -o {cluster.o} --cpus-per-task {cluster.cpus-per-task}"

