#!/usr/bin/env python

import sys
import pandas as pd
import os.path
from pathlib import Path


#get the arguments
in_file = sys.argv[1]
out_file = sys.argv[2]
chrom = sys.argv[3] 
method = sys.argv[4]

#transform strings of file and path to Path
my_file = Path(in_file)

#check if the vcf exists
if my_file.is_file():
    # file exists
    print("input IBD file: ", in_file)
else:
    print("input IBD file: ", in_file , " doesn't exist")
    quit()

#check chromosome number
if (int(chrom) > 0 and int(chrom) <= 22):
    #chromosome is a valid autosome
    print("chromosome number: ", chrom )
else: 
    print("chromosome number: ", chrom , " is invalid, please launch with autosomes")
    quit()




ibd_segments = pd.read_csv(my_file, sep='\t', lineterminator='\n')

if method=="HapIBD":
    ibd_segments.columns = ['SAMPLE1','HAP_INDEX1','SAMPLE2','HAP_INDEX2','CHROM','START','END','GEN_LENGTH']
    #create new columns merging the ids in order
    iibd_segments['SAMPLE1'] = ibd_segments['SAMPLE1'].astype(str)
    ibd_segments['SAMPLE2'] = ibd_segments['SAMPLE2'].astype(str)
    ibd_segments['id1_id2'] = ibd_segments['SAMPLE1'] + ':' + ibd_segments['SAMPLE2']
    #create the IBD summ per pair of individuals
    pair_len = ibd_segments.groupby('id1_id2',as_index=False).agg({'id1_id2':'first' ,'GEN_LENGTH':'sum'})
    pair_num = ibd_segments['id1_id2'].value_counts()
    df_pair_num = pd.DataFrame({'id1_id2':pair_num.index , 'segment_count':pair_num.values})
    pair_IBD = pd.merge(pair_len, df_pair_num, on="id1_id2")

if method=="PhaseIBD":
    ibd_segments.columns = ['chromosome','id1','id2','id1_haplotype','id2_haplotype','start','end','start_cm','end_cm','start_bp','end_bp','size_cM','size_bp']
    #create new columns merging the ids in order
    ibd_segments['id1'] = ibd_segments['id1'].astype(str)
    ibd_segments['id2'] = ibd_segments['id2'].astype(str)
    ibd_segments['id1_id2'] = ibd_segments['id1'] + ':' + ibd_segments['id2']
    #create the IBD summ per pair of individuals
    pair_len = ibd_segments.groupby('id1_id2',as_index=False).agg({'id1_id2':'first' ,'size_cM':'sum'})
    pair_num = ibd_segments['id1_id2'].value_counts()
    df_pair_num = pd.DataFrame({'id1_id2':pair_num.index , 'segment_count':pair_num.values})
    pair_IBD = pd.merge(pair_len, df_pair_num, on="id1_id2")


#write output file
pair_IBD.to_csv(out_file, sep='\t', index=False)













