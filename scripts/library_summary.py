import scanpy as sc
import numpy as np
import os
import sys
from datetime import datetime

def anndata_info(adata, anndata_path):
  """
  Prints basic information about an AnnData object.

  Args:
    adata: An AnnData object.
    anndata_path: Path to the AnnData object file.
  """

  # --- Calculate metrics ---

  # Total reads
  total_reads = np.sum(adata.X)

  # Unique genes
  unique_genes = adata.shape[1]

  # Reads per gene
  reads_per_gene = np.sum(adata.X, axis=0)
  avg_reads_per_gene = total_reads / unique_genes
  median_reads_per_gene = np.median(reads_per_gene)
  min_reads_per_gene = np.min(reads_per_gene)
  max_reads_per_gene = np.max(reads_per_gene)

  # Total cells
  total_cells = adata.shape[0]

  # Reads per cell
  reads_per_cell = np.sum(adata.X, axis=1)
  avg_reads_per_cell = total_reads / total_cells
  median_reads_per_cell = np.median(reads_per_cell)
  min_reads_per_cell = np.min(reads_per_cell)
  max_reads_per_cell = np.max(reads_per_cell)

  # Genes per cell
  genes_per_cell = np.sum(adata.X > 0, axis=1)
  avg_genes_per_cell = np.mean(genes_per_cell)
  median_genes_per_cell = np.median(genes_per_cell)
  min_genes_per_cell = np.min(genes_per_cell)
  max_genes_per_cell = np.max(genes_per_cell)

  # Sparsity
  sparsity = np.sum(adata.X == 0) / adata.X.size * 100

  # Highly expressed genes
  top_n = 10  # Change this to the desired number of top genes
  top_genes_indices = np.argpartition(np.sum(adata.X, axis=0), -top_n)[-top_n:]
  top_genes = adata.var_names[top_genes_indices]
  top_genes_total_reads = reads_per_gene[top_genes_indices]

  # --- Print information ---

  print("-" * 30)
  print("AnnData Object Information")
  print("-" * 30)

  print("\n--- General Information ---")
  print(f"File path:         {anndata_path}")
  print(f"Date and Time:     {datetime.now()}")
  print(f"Total reads:       {total_reads:.3f}")
  print(f"Unique genes:      {unique_genes}")
  print(f"Total cells:       {total_cells}")
  print(f"Sparsity:          {sparsity:.2f}%") 

  print("\n--- Reads per Gene ---")
  print(f"Avg. reads/gene:   {avg_reads_per_gene:.3f}")
  print(f"Median reads/gene: {median_reads_per_gene:.3f}")
  print(f"Min reads/gene:    {min_reads_per_gene:.3f}")
  print(f"Max reads/gene:    {max_reads_per_gene:.3f}")

  print("\n--- Reads per Cell ---")
  print(f"Avg. reads/cell:   {avg_reads_per_cell:.3f}")
  print(f"Median reads/cell: {median_reads_per_cell:.3f}")
  print(f"Min reads/cell:    {min_reads_per_cell:.3f}")
  print(f"Max reads/cell:    {max_reads_per_cell:.3f}")

  print("\n--- Genes per Cell ---")
  print(f"Avg. genes/cell:   {avg_genes_per_cell:.3f}")
  print(f"Median genes/cell: {median_genes_per_cell:.3f}")
  print(f"Min genes/cell:    {min_genes_per_cell}")
  print(f"Max genes/cell:    {max_genes_per_cell}")

  print("\n--- Highly Expressed Genes ---")
  for i, gene in enumerate(top_genes):
    print(f"{gene}: {top_genes_total_reads[i]:.3f} reads")

  print("-" * 30)


if __name__ == "__main__":
  anndata_path = sys.argv[1]

  adata = sc.read_h5ad(anndata_path)
  adata.var_names = adata.var['gene_name'].values
  sc.logging.print_memory_usage()
  anndata_info(adata, anndata_path) 