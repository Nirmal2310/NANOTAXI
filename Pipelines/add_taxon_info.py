import pandas as pd
import sys
import argparse

def add_taxa_data(count_file, taxa_file):
    try:
        taxa_data = pd.read_csv(taxa_file, sep = "\t")

        count_column_name = ['REF_ID', 'COUNT']

        count_data = pd.read_csv(count_file, sep = "\t", header=None, names=count_column_name)

        final_data = pd.merge(count_data, taxa_data, on='REF_ID')

        final_data = final_data.groupby('Species', as_index = False)["COUNT"].sum()

        final_data = pd.merge(final_data, taxa_data, on='Species')

        final_data.to_csv(sys.stdout, sep = "\t", index=False)
    
    except Exception as e:
        print(f"Error occured: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add Taxon Data to Counts File.")
    parser.add_argument("-c", "--counts", type=str, required=True, help="Input Counts File.")
    parser.add_argument("-t", "--taxa", type=str, required=True, help="Input GSR DB Taxonomy file.")
    args = parser.parse_args()
    
    add_taxa_data(args.counts, args.taxa)