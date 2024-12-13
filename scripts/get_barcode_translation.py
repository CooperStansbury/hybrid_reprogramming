import os
import sys
import pandas as pd
import pysam


def main():
    if len(sys.argv) != 5:
        print("Usage: script_name detected_barcodes translation_path mapping_csv translation_fasta")
        sys.exit(1)

    detected_barcodes = sys.argv[1]
    translation_path = sys.argv[2]
    mapping_csv = sys.argv[3]
    translation_fasta = sys.argv[4]

    # Load detected barcodes
    print("\nLoading detected barcodes...")
    detected_fasta = pysam.FastaFile(detected_barcodes)
    detected = list(detected_fasta.references)
    print(f"Number of detected barcodes: {len(detected)}")

    # Load translation path
    print("\nLoading translation path...")
    df = pd.read_csv(translation_path, sep='\t', header=None, names=['cell_id', 'feature_id'])
    print(f"Shape of translation dataframe: {df.shape}")

    # Build the mapping
    print("\nBuilding mapping...")
    mapping = df[df['cell_id'].isin(detected)].reset_index(drop=True)
    print(f"Shape of mapping dataframe: {mapping.shape}")
    mapping.to_csv(mapping_csv, index=False)
    print(f"Mapping saved to {mapping_csv}")

    # Write the translation FASTA
    print("\nWriting translation FASTA...")
    if os.path.exists(translation_fasta):
        os.remove(translation_fasta)
        print(f"Existing file '{translation_fasta}' has been removed.")

    with open(translation_fasta, "a") as outfile:
        for _, row in mapping.iterrows():
            cell_id = row['cell_id']
            feature_id = row['feature_id']
            print(f">{cell_id}", file=outfile)
            print(feature_id, file=outfile)

    print(f"\nTranslation FASTA saved to {translation_fasta}")
    print("Done!")

if __name__ == "__main__":
    main()