from datetime import datetime
import pandas as pd
import yaml
from pathlib import Path
import re
import os
import glob
import sys
from utils import utils

BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"

# big picture variables
OUTPUT = config['output_path']
print("\nOUTPUT PATH:")
print(OUTPUT)

# load in fastq path
input_path = os.path.abspath(config['inputs'])
input_df = pd.read_csv(input_path, comment="#")
samples = input_df['sample_id'].to_list()

# get input names 
input_names = utils.get_input_names(input_df, OUTPUT)

print("\nINPUT FILES:")
[print(x) for x in samples]

# timestamp
now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
test_file = OUTPUT + "test/" + now + ".txt"


# biolegend reads
hash_dir = config['hashing_read_dir']
hash_samples = [os.path.basename(x).replace('.fastq.gz', '') for x in glob.glob(hash_dir + "*.fastq.gz")]


################ RULE FILES ################
include: "rules/reference.smk"
include: "rules/demultiplex.smk"
include: "rules/hashing.smk"
include: "rules/core.smk"
include: "rules/anndata.smk"
include: "rules/isoquant.smk"
include: "rules/rna_velocity.smk"

rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/transcripts.fa',
        OUTPUT + 'references/annotations.gtf',
        OUTPUT + 'references/geneTable.csv',
        expand(OUTPUT + "fastq/{sid}.raw.fastq.gz", sid=samples),
        expand(OUTPUT + "demultiplex/{sid}.done", sid=samples),
        expand(OUTPUT + "reports/fastqc/{sid}.report.html", sid=samples),
        expand(OUTPUT + "mapping/{sid}.bam.bai", sid=samples),
        expand(OUTPUT + "mapping/{sid}.tagged.bam.bai", sid=samples),
        expand(OUTPUT + "reports/bamstats/{sid}.bamstats", sid=samples),
        expand(OUTPUT + "individual_counts/{sid}.counts.txt", sid=samples),
        OUTPUT + 'reports/seqkit_stats/raw_report.txt',
        OUTPUT + 'reports/seqkit_stats/demultiplexed_report.txt',
        OUTPUT + 'merged/merged.bam.bai',
        OUTPUT + 'merged/merged.stats',
        OUTPUT + 'merged/merged.bamstats',
        OUTPUT + 'merged/merged.counts.txt',
        OUTPUT + 'scanpy/raw.anndata.h5ad',
        OUTPUT + "reports/anndata/library_summary.txt",
    

rule isoforms:
    input:
        expand(OUTPUT + "isoquant/{sid}.done", sid=samples),
        OUTPUT + "isoquant/annotations.db",
        OUTPUT + "isoquant_prepared/gene_counts.csv",
        OUTPUT + "isoquant_prepared/transcript_counts.csv",
        OUTPUT + "isoquant_prepared/isoforms.csv",


rule rna_velocity:
    input:
        OUTPUT + 'velocyto/merged.tagged.bam',
        OUTPUT + 'velocyto/run_velocyto.done',


rule hashing:
    input:
        expand(OUTPUT + "hash_fastq/{sid}.fastq.gz", sid=hash_samples),
        OUTPUT + "hash_reference/hash.fasta",
        OUTPUT + "hash_reference/barcode_whitelist.txt",
        OUTPUT + "hash_reference/barcode_translation.txt",
        OUTPUT + "reports/nanoplexer/hashing_report.txt", # demutliplex with nanoplexer
        OUTPUT + "hash_reference/detected_barcodes.fasta",
        OUTPUT + "hash_reference/barcode_mapping.csv",
        OUTPUT + "hash_reference/translated_barcodes.fasta",
        OUTPUT + "hash_reference/translated_barcodes.fasta.amb",
        expand(OUTPUT + "hash_alignment/{sid}.bam", sid=tag_list),
        expand(OUTPUT + "reports/hash_count/{sid}.csv", sid=tag_list),
        OUTPUT + "reports/hash_map/hashmap.csv",

rule test:
    output:
        touch(test_file),