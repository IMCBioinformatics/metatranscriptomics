# Host Contamination calcualtor
# This script takes the three seqkit files and generates read counts of low quality and clean reads and host contamination 

### Importing required packages

import pandas as pd
import os
import logging


### Reading in input files and parameters

fwd = snakemake.params[0]
rev = snakemake.params[1]

raw = pd.read_csv(snakemake.input[0],sep="\s+").astype(str)

prinseq = pd.read_csv(snakemake.input[1],sep="\s+").astype(str)

bmtagger = pd.read_csv(snakemake.input[2],sep="\s+").astype(str)

### Creating a function to assign read types
def assign_read_type(x):
    if '_R1_' in x or '_R1.' in x or '_1.' in x or '_1_' in x:
        return "1"
    # Check reverse patterns
    elif '_R2_' in x or '_R2.' in x or '_2.' in x or '_2_' in x:
        return "2"
    else:
        return "Unknown"


### Adding new columns for read types and sample IDs to existing dataframes

##### RAW READS        
raw['Type'] = raw.iloc[:,0].apply(assign_read_type)
raw['Sample']=raw.iloc[:, 0].str.split(fwd,expand=True).loc[:,0].str.split(rev,expand=True).loc[:,0]


##### READS after prinseq
prinseq['Type'] = prinseq.iloc[:,0].apply(assign_read_type)
prinseq['Sample']=prinseq.iloc[:,0].str.split("_filtered_",expand=True,).iloc[:,0]


##### READS after bmtagger
bmtagger['Type'] = bmtagger.iloc[:,0].apply(assign_read_type)
bmtagger['Sample']=bmtagger.iloc[:,0].str.split("_bmtagged_",expand=True,).iloc[:,0]


### Multiindexing and changing row lables to sample and type
raw=raw.set_index(['Sample',"Type"])

prinseq=prinseq.set_index(['Sample',"Type"])

bmtagger=bmtagger.set_index(['Sample',"Type"])


### Joining all three dataframes
merged=bmtagger.join(prinseq,lsuffix="_bmtagged", rsuffix="_prinseq").join(raw)


### Removing row lables        
merged=merged.reset_index()

### Restructuring the merged table and calculation of host contamination, clean reads, and low quality reads
merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')] = merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')].applymap(lambda x: str(x).replace(',', '')).astype(float)

merged.loc[:,'Host_contamination']=(merged.loc[:,'num_seqs_prinseq']-merged.loc[:,'num_seqs_bmtagged'])

merged.loc[:,'Clean_reads']=merged.loc[:,'num_seqs_bmtagged']

merged.loc[:,'Low_quality_reads']=merged.loc[:,'num_seqs']-merged.loc[:,'num_seqs_prinseq']

### Saving the final table in csv format        
merged.loc[:,['Sample','Type','Low_quality_reads','Clean_reads','Host_contamination','num_seqs','num_seqs_prinseq','num_seqs_bmtagged']].to_csv(snakemake.output[0],index=False)
