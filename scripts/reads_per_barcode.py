import pysam
import pandas as pd
import sys
from collections import defaultdict

def count_reads_and_avg_mapq_per_contig(sam_file_path):
    """
    Counts the number of reads and calculates average mapping quality per contig in a SAM file.

    Args:
        sam_file_path: The path to the SAM file.

    Returns:
        A pandas DataFrame with columns 'Contig', 'Read Count', and 'Average Mapping Quality'.
    """
    samfile = pysam.AlignmentFile(sam_file_path, "r")

    # Use defaultdict for efficient counting and summing
    data_per_contig = defaultdict(lambda: {'count': 0, 'mapq_sum': 0})
    for read in samfile:
        contig = read.reference_name
        data_per_contig[contig]['count'] += 1
        data_per_contig[contig]['mapq_sum'] += read.mapping_quality

    samfile.close()

    # Create DataFrame and calculate average
    df = pd.DataFrame.from_dict(data_per_contig, orient='index').reset_index()
    df.columns = ['Contig', 'Read Count', 'Mapping Quality Sum']
    df['Average Mapping Quality'] = df['Mapping Quality Sum'] / df['Read Count']
    df.drop(columns=['Mapping Quality Sum'], inplace=True)

    return df

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <sam_file_path> <output_csv_path>")
        sys.exit(1)

    sam_file_path = sys.argv[1]
    output_csv_path = sys.argv[2]

    results_df = count_reads_and_avg_mapq_per_contig(sam_file_path)
    results_df.to_csv(output_csv_path, index=False)
    print(f"Results saved to {output_csv_path}")
    

# import pysam
# import pandas as pd
# import sys
# from collections import Counter

# def count_reads_per_contig(sam_file_path):
#     """
#     Counts the number of reads aligned to each contig in a SAM file.

#     Args:
#         sam_file_path: The path to the SAM file.

#     Returns:
#         A pandas DataFrame with columns 'Contig' and 'Read Count'.
#     """
#     samfile = pysam.AlignmentFile(sam_file_path, "r")

#     # Use Counter for efficient counting
#     reads_per_contig = Counter(read.reference_name for read in samfile)

#     samfile.close()

#     return pd.DataFrame.from_dict(reads_per_contig, orient='index', columns=['Read Count']).reset_index().rename(columns={'index':'Contig'})

# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python script_name.py <sam_file_path> <output_csv_path>")
#         sys.exit(1)

#     sam_file_path = sys.argv[1]
#     output_csv_path = sys.argv[2]

#     results_df = count_reads_per_contig(sam_file_path)
#     results_df.to_csv(output_csv_path, index=False)
#     print(f"Results saved to {output_csv_path}")