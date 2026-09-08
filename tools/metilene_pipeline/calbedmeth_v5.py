import numpy as np
from collections import defaultdict
import argparse
import logging
import os
import pandas as pd
import csv

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def bed2poslist(bed_file):
    """
    Read bed file and convert it to a dictionary.
    """
    chr_pos_list = []
    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, *bed_value = line.strip().split()
            for pos in range(int(start), int(end)+1):
                chr_pos_list.append(chrom + ":" + str(pos))
    return set(chr_pos_list)

def pos_met2subset(met_file, chr_pos_list, tmp_outf):
    """
    Extract methylation values for positions in chr_pos_list to subset pos methylation file.
    """
    with open(met_file, 'r') as f, open(tmp_outf, 'w') as outf:
        for line in f:
            if not line.startswith('chrom'):
                chrom, pos, *met_values = line.strip().split()
                if chrom + ":" + pos in chr_pos_list:
                    outf.write(line)



def main():
    parser = argparse.ArgumentParser(description='Calculate methylation values in BED regions.')
    parser.add_argument('-b', '--bed', required=True, help='BED file')
    parser.add_argument('-m', '--met', required=True, help='Methylation file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()
    logging.info("Reading BED file and get positions")
    chr_pos_list = bed2poslist(args.bed)
    logging.info("Extract methylation values for positions in chr_pos_list to subset pos methylation file.")
    pos_met2subset(args.met, chr_pos_list, args.output+'.tmp')
    
    fileBed = pd.read_csv(args.output+'.tmp',header=None,index_col=[0,1],sep="\t")
    ohandle = open(args.output,"w")
    out = csv.writer(ohandle,delimiter='\t')
    for rec in csv.reader(open(args.bed,"r"),delimiter="\t"):
        value = fileBed.loc[(rec[0],),:].loc[int(rec[1]):int(rec[2])].mean()
        out.writerow(rec+value.tolist())
    os.remove(args.output+'.tmp') 
    logging.info("Done!")
    
if __name__ == '__main__':
    main()
