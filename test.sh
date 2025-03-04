#!/bin/bash

CONFIG='config/cluster'
CORES=36

### export the environment 
conda env export > environment.yml

## build the workflow from the most current snakefile
cp Snakefile workflow.smk

# run it
snakemake --profile ${CONFIG} --use-conda --cores ${CORES} --rerun-incomplete --latency-wait 90 --verbose -s workflow.smk 