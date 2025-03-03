import pysam
import pandas as pd
import sys
import argparse

def extract_alignment_data(bam_file, taxa_file, identity, coverage):
    try:
        header_data = []
        alignment_data = []
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                read_id = read.query_name
                if read.is_unmapped:
                    continue
                cigar = read.cigartuples
                aligned_seq = len(read.query_alignment_sequence)
                ref_id = read.reference_name
                if cigar:
                    match = 0
                    for i in range(len(cigar)):
                        if cigar[i][0] == 7:
                            match += cigar[i][1]
                    perc_match = (100 * match / aligned_seq)
                    alignment_data.append({
                            "READ_ID": read_id,
                            "REF_ID": ref_id,
                            "MATCH": match,
                            "ALIGNED_LENGTH": aligned_seq,
                            "PERCENTAGE_IDENTITY": perc_match
                    })

        bam = pysam.AlignmentFile(bam_file, "rb")

        header = bam.header.to_dict ()

        for i in range(len(header['SQ'])):
            header_data.append({
                "REF_ID": header.get('SQ')[i]['SN'],
                "REF_LENGTH": header.get('SQ')[i]['LN']
                })

        alignment_dataframe = pd.DataFrame(alignment_data)

        header_dataframe = pd.DataFrame(header_data)

        final_data = pd.merge(header_dataframe, alignment_dataframe, on='REF_ID')

        final_data = final_data.assign(PERCENTAGE_COVERAGE = 100*final_data['ALIGNED_LENGTH']/final_data['REF_LENGTH'])

        identity = int(identity)

        coverage = int(coverage)

        final_data = final_data.query('PERCENTAGE_IDENTITY >= @identity and PERCENTAGE_COVERAGE >= @coverage').groupby('REF_ID').size().reset_index(name='COUNT')

        taxa_data = pd.read_csv(taxa_file, sep = "\t")

        taxa_data = taxa_data.rename(columns={'FeatureID': 'REF_ID'})

        final_data = pd.merge(final_data, taxa_data, on='REF_ID')

        final_data = final_data.groupby('Species', as_index = False)["COUNT"].sum()

        final_data.to_csv(sys.stdout, sep = "\t", index=False)
    
    except Exception as e:
        print(f"Error occured: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Alignment Information from BAM file.")
    parser.add_argument("-b", "--bam", type=str, required=True, help="Input sorted BAM file mapped to GSR DB.")
    parser.add_argument("-t", "--taxa", type=str, required=True, help="Input GSR DB Taxonomy file.")
    parser.add_argument("-i", "--identity", type=str, required=True, help="Percentage Identity Cutoff.")
    parser.add_argument("-c", "--coverage", type=str, required=True, help="Percentage Coverage Cutoff.")
    args = parser.parse_args()
    
    extract_alignment_data(args.bam, args.taxa, args.identity, args.coverage)