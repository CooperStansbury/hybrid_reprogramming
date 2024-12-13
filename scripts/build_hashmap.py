import pandas as pd
import numpy as np
import os
import sys


if __name__ == "__main__":
    output_path = sys.argv[1]
    file_list = sys.argv[2:]

    df = []

    for fpath in file_list:
        tag = os.path.basename(fpath).replace(".csv", "")
        tmp = pd.read_csv(fpath)
        tmp.columns = ['barcode', 'reads', 'mapq']
        tmp['tag'] = tag

        total_barcodes = len(tmp)
        mapped_barcodes = len(tmp[tmp['reads'] > 0])

        # Nicer print statements with more information and percentages
        print(f"\n{'='*20} Processing File: {tag} {'='*20}")
        print(f"File path: {fpath}")
        print(f"Total reads: {tmp['reads'].sum()}")
        print(f"Total barcodes: {total_barcodes}")
        print(f"Barcodes with reads mapped: {mapped_barcodes} ({(mapped_barcodes/total_barcodes)*100:.2f}%)")
        print(f"Mean reads per barcode: {tmp['reads'].mean():.2f}")
        print(f"Median reads per barcode: {tmp['reads'].median():.2f}")

        df.append(tmp)

    df = pd.concat(df)

    pdf = df.copy()
    pdf = pdf[pdf['barcode'].notna()]
    pdf = pdf[pdf['mapq'] > 0]

    pdf = pd.pivot_table(
        pdf, 
        index="barcode",
        columns="tag",
        values="reads",
    )

    pdf['assigned_tag'] = pdf.idxmax(axis=1)
    pdf = pdf.reset_index()

    # Calculate and print percentages for tag assignments
    print("\n===== Tag Assignment Summary =====")
    tag_counts = pdf['assigned_tag'].value_counts()
    total_assigned = tag_counts.sum()
    for tag, count in tag_counts.items():
        percentage = (count / total_assigned) * 100
        print(f"{tag}: {count} ({percentage:.2f}%)")

    pdf.to_csv(output_path, index=False)