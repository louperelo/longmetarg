#!/usr/bin/env python

import csv
import sys
import argparse
import pandas as pd
from ete3 import NCBITaxa

parser = argparse.ArgumentParser()
parser.add_argument("--blast_out", metavar="<BLAST>", dest="blast", help="enter the path to the blast output file. By default '/out_dir/result_blast.tsv' will be used",
                    type=str, default="result_tax.tsv")

# print help message for user
parser.print_help()

# get command line arguments
args = parser.parse_args()

# read files from path
blast_path = args.blast

#adapted from https://stackoverflow.com/questions/43867631/how-can-i-get-taxonomic-rank-names-from-taxid

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

if __name__ == '__main__':
    lst = pd.read_csv(blast_path, sep='\t', header=None).iloc[:,1].fillna(2).to_list()
    taxids = [int(i) for i in lst]
    #possible ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    desired_ranks = ['species', 'class']
    
    results = list()
    for taxid in taxids:
        results.append(list())
        results[-1].append(str(taxid))
        ranks = get_desired_ranks(taxid, desired_ranks)
        for key, rank in ranks.items():
            if rank != '<not present>':
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else:
                results[-1].append(rank)

    with open('taxon.tsv', 'w') as f:
        sys.stdout = f
        #generate the header
        header = ['taxid']
        header.extend(desired_ranks)
        print('\t'.join(header))
        #print the results
        for result in results:
            print('\t'.join(result))
    