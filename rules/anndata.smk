rule make_anndata:
    input:
        counts=OUTPUT + 'merged/merged.counts.txt',
        gtf=OUTPUT + "references/annotations.gtf",
    conda:
        "scanpy"
    output:
        OUTPUT + "scanpy/raw.anndata.h5ad",
    shell:
        """python scripts/make_anndata.py {input.counts} {input.gtf} {output}"""


rule anndata_report:
    input:
        OUTPUT + "scanpy/raw.anndata.h5ad",
    output:
        OUTPUT + "reports/anndata/library_summary.txt",
    conda:
        "scanpy"
    shell:
        """python scripts/library_summary.py {input} > {output}"""


rule process_anndata:
    input:
        raw=OUTPUT + "scanpy/raw.anndata.h5ad",
        v5=OUTPUT + 'v5_tagged/v5_result.table.csv'
    output:
        OUTPUT + "scanpy/processed.anndata.h5ad",
    conda:
        "scanpy"
    params:
        annotations="/home/cstansbu/git_repositories/ONT-single-cell/config/gene_annotations/"
    shell:
        """python scripts/process_anndata.py {input.raw} {output} {params.annotations} {input.v5}"""
        
        
rule cluster_anndata:
    input:
        OUTPUT + "scanpy/processed.anndata.h5ad",
    output:
        OUTPUT + "scanpy/clustered.anndata.h5ad",
    conda:
        "scanpy"
    shell:
        """python scripts/cluster_cells.py {input} {output}"""