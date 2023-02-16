#!/usr/bin/env python

import phasedibd as ibd
import sys
import pandas as pd
import os.path
from pathlib import Path

print("TEST")

#get the arguments
vcf_file = sys.argv[1] 
chrom = int(sys.argv[2]) 
out_path = sys.argv[3]
if sys.argv[4].endswith('.map'):
    genetic_map = sys.argv[4]
else:
    genetic_map="NA"
sample_id=sys.argv[5]



print("TEST2")

#transform strings of file and path to Path
my_file = Path(vcf_file)
my_path = Path(out_path)
my_map = Path(genetic_map)
my_ids = Path(sample_id)

#check if the vcf exists
if my_file.is_file():
    # file exists
    print("input vcf file: ", vcf_file)
else:
    print("input vcf file: ", vcf_file , " doesn't exist")
    quit()

#check chromosome number
if (chrom > 0 and chrom <= 22):
    #chromosome is a valid autosome
    print("chromosome number: ", chrom )
else: 
    print("chromosome number: ", chrom , " is invalid, please launch with autosomes")
    quit()


#check if out_path exists
if my_path.is_dir():
    # directory exists
    print("output directory: ", out_path )
else: 
    print("output directory: ", out_path , " doesn't exist")
    quit()


#check wheter there is a genetic map specified or not
if my_map.is_file():
    # file exists
    print("genetic map file: ", genetic_map)
    haplotypes = ibd.VcfHaplotypeAlignment(vcf_file, genetic_map)
else:
    print("genetic map file not specified")
    haplotypes = ibd.VcfHaplotypeAlignment(vcf_file)

#check if id file exists
if my_ids.is_file():
    # file exists
    print("ID corresponding table file: ", sample_id)
else:
    print("ID corresponding table file: ", sample_id , " doesn't exist")
    quit()



#instanciate a tpbwt objet and launch its method
tpbwt = ibd.TPBWTAnalysis()
ibd_results = tpbwt.compute_ibd(haplotypes, chromosome=str(chrom), L_f=2.0 , L_m=200 )

#print number of segment
#print("number of IBD segments: ", ibd_results.shape[0])


#add columns for size of IBD segment
ibd_results['size_cM'] = ibd_results['end_cm'] - ibd_results['start_cm']
ibd_results['size_bp'] = ibd_results['end_bp'] - ibd_results['start_bp']

#copy the dataframe
#ibd_results_original=ibd_results.copy()


#convert back id1 and id2 back to CaG ids
correspondance_id = pd.read_csv(sample_id, sep='\t', lineterminator='\n')
print(len(correspondance_id.index))


#change the id1 and id2 for the sample ids in the initial VCF
ibd_results['id1'] = ibd_results['id1'].map(correspondance_id.set_index('index')['VCF_ID'])
ibd_results['id2'] = ibd_results['id2'].map(correspondance_id.set_index('index')['VCF_ID'])


print("output file")
#write output file
out_file=out_path+"/chr"+str(chrom)+".ibd"
ibd_results.to_csv(out_file, sep='\t', index=False)

#out_file=out_path+"/chr"+str(chrom)+"original.ibd"
#ibd_results_original.to_csv(out_file, sep='\t', index=False)


exit()

