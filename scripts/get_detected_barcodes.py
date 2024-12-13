import scanpy as sc
import os
import sys
import pandas as pd

def anndata_to_fasta(fpath, outpath):
    """
    Converts anndata file to fasta format.

    Args:
        fpath: Path to the anndata file.
        outpath: Path to write the fasta file.
    """

    adata = sc.read_h5ad(fpath)

    # Overwrite existing file if it exists
    if os.path.exists(outpath):
        os.remove(outpath)

    with open(outpath, "a") as outfile:
        for bc in adata.obs_names:
            print(f">{bc}", file=outfile)
            # Assuming 'bc' itself is the sequence you want to write
            print(bc, file=outfile) 

if __name__ == "__main__":    
    anndata_file_path = sys.argv[1]
    fasta_output_path = sys.argv[2]

    anndata_to_fasta(anndata_file_path, fasta_output_path)