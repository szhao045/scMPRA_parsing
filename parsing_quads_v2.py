import gzip
import argparse
import pandas as pd
from fuzzysearch import find_near_matches
'''
This code is largely copied from our scTRIP parsing code, credit: A Ramu. 
'''

def check_bulk_plasmid_lib(read1):
    '''
    General information: spike-in library checks are done with 2by150 reads, so
    Read1 Alone is enough for giving all the information.
    Sample Read is stored at ./sample_read.dna in the same folder 
    Input: read-in line from read1. 
    Output: parsed promBC and rBC
    '''
    # Global Variables for extracting the barcodes
    # They are all 8 pbs, but I can always increase them anytime
    BEFORE_PBC = 'AATCTAGA'
    AFTER_PBC = 'GTCGAGAT'
    BEFORE_RBC = 'AAGTTATG'
    AFTER_RBC = 'GCTTTAAG'
    # Initiate a dictionary to hold results
    pop = {}
    # Initiate zero values for counting 
    pop['wrong'] = 0
    # Find the before pBC, and the after pBC, and get the promoter BC, and return length. 
    pBC_left_lim = find_near_matches(BEFORE_PBC, read1, max_l_dist = 1)
    pBC_right_lim = find_near_matches(AFTER_PBC, read1, max_l_dist = 1)
    rBC_left_lim = find_near_matches(BEFORE_RBC, read1, max_l_dist = 1)
    rBC_right_lim = find_near_matches(AFTER_RBC, read1, max_l_dist = 1)
    if len(pBC_left_lim) !=0 and len(pBC_right_lim) != 0 :
        if pBC_left_lim[0].end < pBC_right_lim[0].start:
            putative_pBC = read1[pBC_left_lim[0].end:pBC_right_lim[0].start]
            if len(rBC_left_lim) !=0 and len(rBC_right_lim)!=0:
                if rBC_left_lim[0].end < rBC_right_lim[0].start:
                    putative_rBC = read1[rBC_left_lim[0].end:rBC_right_lim[0].start]
                else:
                    putative_rBC = []
            else:
                putative_rBC = []
        else: 
            putative_pBC = []
    else:
        putative_pBC = []
    if len(putative_pBC) == 12:
        if len(putative_rBC) == 25:
            # Get the counts
            pop['bc_pair'] = putative_pBC + '\t' + putative_rBC
        else:
            pop['wrong'] += 1
    else:
        pop['wrong'] += 1
    # check
    return pop

def extract_10xBCs(line2):
    if len(line2) != 0:
        cellBC = line2[0:16]
        umi = line2[16:28]
        pop = cellBC + '\t' + umi
    else:
        pop = 'N'*16 + '\t' + 'N'*12
    return pop

def parse_fastq(r1_file, r2_file):
    '''
    Function to parse fastq file from bulk RNA-seq. 
    This code should be able to be expanded to run with the later parsing
    So it's important to make it modular. 
    '''
    # Initiate a dict for holding the parsed barcodes, key is promBC + rBC
    # value is the number of reads
    total_reads = 0
    wrong_reads = 0
    line_num = 0
    promBCrBC = {}
    with gzip.open(r1_file, 'rt') as r1:
        with gzip.open(r2_file, 'rt') as r2:
        # Call bulk plasmid parsing 
        # TODO: make it general purpose to call
            for line1,line2 in zip(r1,r2):    
                line1 = line1.rstrip("\n")
                line2 = line2.rstrip("\n")
                line_num += 1
                if line_num % 4 == 2:
                    total_reads += 1
                    # We first deal with line1 which comes from R2 of the sequencing file
                    bc_pair_info = check_bulk_plasmid_lib(line1)
                    # Get the parcode pair
                    # We then deal with line 2 which comes from R1 of the sequencer.
                    sc_info = extract_10xBCs(line2)
                    # bad cell
                    bad_trio = 'N'*16 + '\t' + 'N'*12
                    if 'bc_pair' in bc_pair_info:
                        if sc_info != bad_trio:
                            bc_pair = bc_pair_info['bc_pair']
                            bc_full = sc_info + '\t' + bc_pair
                            if bc_full not in promBCrBC:
                                promBCrBC[bc_full] = 0
                            promBCrBC[bc_full] += 1
                    # Get the stats
                    wrong_reads += bc_pair_info['wrong']
    print(f'The total number of reads is {total_reads}')
    print(f'The number of wrong is {wrong_reads}')
    df = pd.DataFrame(list(promBCrBC.items()),columns = ['bc_info','counts']) 
    df[['cellBC','umi','pBC','rBC']] = df['bc_info'].str.split('\t', expand=True)
    df = df.drop(columns=['bc_info'])
    return df

def parse_arguments():
    parser = argparse.ArgumentParser(description='Save a csv file for the pBC, rBC, and read counts ')
    parser.add_argument("--R1",
                        help = "Read1 of the FASTQ pair")
    parser.add_argument("--R2",
                        help = "Read2 of the FASTQ pair")
    parser.add_argument("--name",
                    help = "name for saving the parsed barcodes")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    pop_df = parse_fastq(args.R1, args.R2)
    pop_df.to_csv(args.name + '.tsv', index = None, sep = '\t')

if __name__ == "__main__":
    main()