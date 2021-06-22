#!/usr/bin/env python

import pandas as pd
from scipy import stats
import numpy as np
#import seaborn as sns
#import matplotlib.pyplot as plt
import math
from Bio import SeqIO
import io 
import re 
import pysam
from functools import reduce 
import argparse
import os 

parser = argparse.ArgumentParser()
parser.add_argument("--bam_file", metavar="<BAM>", dest="bam", help="enter the path to the alignment.bam file. By default 'aln_F4.bam' will be used",
                    type=str, default="aln_F4.bam")
parser.add_argument("--reads_fasta", metavar="<FASTA>", dest="fasta", help="enter the path to the original fasta file being analysed. By default 'reads.fasta' will be used",
                    type=str, default="reads.fasta")
parser.add_argument("--ident", metavar="<IDENT>", dest="ident", help="enter the int value for minimum identity. By default 80 will be used",
                    type=int, default= 80)
parser.add_argument("--cov_length", metavar="<COV>", dest="cov", help="enter the int value for minimum coverage length. By default 95 will be used",
                    type=int, default= 95)
parser.add_argument("--folder_out", metavar="<OUT>", dest="out", help="enter name for output files. By default 'arg_results' will be used",
                    type=str, default="../out_dir/")
parser.add_argument("--aro_idx", metavar="<IDX>", dest="idx", help="enter the path to the aro_index.csv file. By default 'aro_index.tsv' will be used",
                    type=str, default="aro_index.tsv")


# print help message for user
parser.print_help()

# get command line arguments
args = parser.parse_args()

# read files from path
bam = args.bam
fasta = args.fasta
ident = args.ident
covlen = args.cov
folder = args.out
idx = args.idx

#read list of cigar tuples and get number of matches (0), insertions (1) or deletions (2)
#auxiliary function in parse_bam()
def read_cigar(lof_tup, idnum):
    x = 0
    for t in lof_tup:
        if(t[0]==idnum):
            x += t[1]
    return x

#Joins information from BAM file in pandas dataframe
#query sequence: query_name, query_length
#reference sequence: reference_name (gives one string, is split into ARO, ID, gene name and NCBI reference id), reference_start, reference_length
#alignment: query_alignment_length, number of mismatches and gaps (tag 'NM)
#calculates sequence identity % (identity(A,B)=100*(identical nucleotides / min(length(A),length(B)))), with identical nucleotides = query_alignment_length - NM
#calculates cover length % (query_alignment_length*100 / reference_length)
pd.options.mode.chained_assignment = None
def parse_bam(bam_path):
    aln_file = pysam.AlignmentFile(bam_path, "rb")
    lst = []
    # loop over alignments, get values per contig and store in list of lists (lst)
    for index, aln in enumerate(aln_file.fetch(until_eof = True)): #index = int(0 ... n), aln = all information on read
        substr = [aln.query_name, aln.query_length, aln.query_alignment_length, aln.get_tag('NM'), aln.reference_length, aln.reference_start, aln.cigartuples]
        #divide information in reference_name
        string = str(aln.reference_name)
        start=[]
        stop=[]
        for i, c in enumerate(string):
            if ((c==':')):
                start.append(i+1)
            elif (c=='|'):
                stop.append(i)
            else:
                continue
        stop.append(len(string)) 

        for i in range(0, len(start)):
            #substr = []
            substr.append(string[start[i]:stop[i]])
        lst.append(substr)
    #print(lst[0:10])
    df = pd.DataFrame(lst, columns=('contig_name', 'contig_length', 'aln_length', 'aln_nm', 'ref_length', 'ref_start', 'c_tuples', 'ref_ARO', 'ref_ID', 'ref_genename', 'ref_NCBI'))
    #get number of matches from cigar tuples
    df['matches'] = df['c_tuples'].apply(lambda x: read_cigar(x, 0))
    df['insertions'] = df['c_tuples'].apply(lambda x: read_cigar(x, 1))
    df['deletions'] = df['c_tuples'].apply(lambda x: read_cigar(x, 2))
    #infer contig_length in repetitions of same contig_name (otherwise the value is 0)
    for i in range(1, df.shape[0]-1):
        if (df['contig_name'].iloc[i+1]==df['contig_name'].iloc[i]):
            df['contig_length'].iloc[i+1] = df['contig_length'].iloc[i]
    
    #calculate coverage length
    df['cov_length'] = df['aln_length']*100/df['ref_length']

    #Sequence identity is the amount of characters which match exactly between two different sequences.
    #identity(A,B)=100% (num identical nucleotides / min(length(A),length(B)))
    df['cov_identity'] = 100*df['matches']/(df.loc[:,['aln_length','ref_length']].min(axis=1))
    return df

#Filter df for highest identity and coverlength rates
def filter_best(df, ident, cov_l):
    return df[(df['cov_identity']>=ident) & (df['cov_length']>=cov_l)]

#Filter assembly fasta for contigs of interest (data) and save to out_name.fasta
#for taxonomic analysis
def arg_contigs(data, fasta, out_name):
    #filter contigs with antibiotic resistance genes
    arg_contigs = data['contig_name'].drop_duplicates().to_list()
    # filter contig sequence information from original fasta file
    #filter fasta for contigs with antibiotic resistance genes (arg) for taxonomic analysis
    fasta_sequences = SeqIO.parse(open(fasta),'fasta')
    with open(out_name, 'w') as out_file:
        for fasta in fasta_sequences:
            #name, sequence = fasta.id, fasta.seq.tostring() #tostring() should be replaced by str(fasta.seq), but is not working on my computer
            name, sequence = fasta.id, str(fasta.seq) 
            for c in arg_contigs:
                if (name==c):
                    out_file.write('>'+ name + '\n' + sequence + '\n')

#check for and eliminate less significant (lower cover identity) overlaps
#generate list of index numbers of non-overlapping hits from df sorted by coverage identity (highest first)
#in case of overlaps, keep the hit with the highest coverage identity
def overlaps(df_in):
    df = df_in.reset_index()
    #list of contig_names
    reads = df['contig_name'].unique()
    #list of indices to keep
    keep = []
    #check overlaps for one contig_name at a time
    for read in reads:
        #create dataframe for each contig_name, sorted by cov_identity, highest value first
        readdf = df[df['contig_name']==read].sort_values(by='cov_identity', ascending=False)
        #list of indices to keep for each read
        k=[]
        #iterate over each enty for one read
        for i in range(0, readdf.shape[0]-1):
            #append first entry of sorted readdf (highest cov_identity) to list of indices to keep for this contig_name      
            k.append(readdf['index'].iloc[0])
            #list for indices of contigs not overlapping with first entry
            lst=[]
            #compare first entry with all other entries
            for j in range (i+1, readdf.shape[0]):
                #get start s and end e position of two resistance gene hits 
                s1, e1 = readdf['ref_start'].iloc[i], readdf['ref_start'].iloc[i] + readdf['ref_length'].iloc[i]
                s2, e2 = readdf['ref_start'].iloc[j], readdf['ref_start'].iloc[j] + readdf['ref_length'].iloc[j]               
                #if there is no overlap, add the entry index to lst
                if (e1<s2 or e2<s1):
                    lst.append(readdf['index'].iloc[j])                                         
            #update readdf, only keep entries with index in lst
            readdf = readdf[readdf['index'].isin(lst)]
            #if updated readdf only contains one entry, add index to k and pass on to next read
            if (readdf.shape[0]==1):
                k.append(readdf['index'].iloc[0])
                break
            #if updated readdf is empty, pass on to next read
            if(readdf.shape[0]==0):
                break
        #append indices for each read to lst keep    
        keep.append(k)
    #flatten list of lists (keep)
    keep = reduce(lambda x,y: x+y,keep)
    return(df[df['index'].isin(keep)])

if __name__ == "__main__":

    #extract data of interest from bam file, filter best hits and eliminate overlaps
    result_df = overlaps(filter_best(parse_bam(bam), ident, covlen))

    #add corresponding drug class from CARD aro_index.tsv to result_df
    rgdrug_dict = pd.read_csv(idx, sep='\t').set_index('ARO Name').to_dict()['Drug Class']
    result_df['drug_class'] = result_df['ref_genename'].map(rgdrug_dict)
    #save result_df as tsv
    result_df.to_csv("argHitsDf.tsv", sep='\t')
    
    #save reads/contigs of hits in result_df in 'result.fasta' for further analysis with PlasFlow or Blast/Diamond
    arg_contigs(result_df, fasta, "argHits.fasta") 


