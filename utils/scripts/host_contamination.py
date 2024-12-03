# Host Contamination calcualtor
# This script takes the three seqkit files and generates read counts of  clean reads and host contamination 
# The script can sometimes fail at merge step, when raw R1 and R2 suffix are different from bmtagger and prinseq files. 


import pandas as pd
import os
import logging
import re

logging.basicConfig(level=logging.DEBUG)

logging.info("The raw file R1 suffix:"+ str(snakemake.params[0]))
logging.info("The raw file R2 suffix:"+ str(snakemake.params[1]))

raw=pd.read_csv(snakemake.input[0],sep="\s+").astype(str)
# After reading raw
logging.debug(f"Raw columns: {raw.columns.tolist()}")
logging.debug(f"Raw sample data:\n{raw.head()}")


prinseq=pd.read_csv(snakemake.input[1],sep="\s+").astype(str)
# After reading prinseq
logging.debug(f"Prinseq columns: {prinseq.columns.tolist()}")
logging.debug(f"Prinseq sample data:\n{prinseq.head()}")


bmtagger=pd.read_csv(snakemake.input[2],sep="\s+").astype(str)
# After reading bmtagger
logging.debug(f"Bmtagger columns: {bmtagger.columns.tolist()}")
logging.debug(f"Bmtagger sample data:\n{bmtagger.head()}")


#raw['Sample']=raw.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split(snakemake.params[0],expand=True).loc[:,0].str.split(snakemake.params[1],expand=True).loc[:,0]

##########How to create the Sample name
def extract_sample(filepath, delimiter='_'):
    # Extract the basename (e.g., "Library_negative_batch1_2.fastq.gz")
    basename = os.path.basename(filepath)

    # Remove multiple extensions by applying splitext twice
    # First split removes '.gz', second removes '.fastq'
    name_without_ext = os.path.splitext(os.path.splitext(basename)[0])[0]

    # Split from the right once to remove the trailing number
    sample = name_without_ext.rsplit(delimiter, 1)[0]

    return sample


# Apply the function to create the 'Sample' column
raw['Sample'] = raw['file'].apply(extract_sample)




#def find_str(x):
# for p in snakemake.params:
#  if p in x:
#   return p.split('_')[1].replace('R','')
# return x


###########How to create the type name
def find_str(x):
    """
    Extracts the Type identifier (e.g., '1' or '2') from the filename.
    Assumes filenames are in the format: SampleName_1.fastq.gz or SampleName_2.fastq.gz
    """
    basename = os.path.basename(x)
    pattern = re.compile(r'_(\d)\.fastq(?:\.gz)?$')
    match = pattern.search(basename)
    if match:
        return match.group(1)  # Returns '1' or '2'
    else:
        logging.warning(f"Regex did not match for file: {basename}")
        return None  # or a default value if preferred



raw['Type']=raw.iloc[:,0].apply(find_str)
logging.debug(f"raw columns: {raw.columns.tolist()}")
logging.debug(f"raw sample data:\n{raw.head()}")



prinseq['Sample']=prinseq.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_filtered",expand=True,).iloc[:,0]

prinseq['Type']=prinseq.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_filtered_",expand=True,).iloc[:,1].str.split(".fastq",expand=True,).iloc[:,0]

bmtagger['Sample']=bmtagger.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_bmtagged",expand=True,).iloc[:,0].str.split(".fastq",expand=True,).iloc[:,0]

bmtagger['Type']=bmtagger.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_bmtagged_",expand=True,).iloc[:,1].str.split(".fastq",expand=True,).iloc[:,0]


raw=raw.set_index(['Sample','Type'])
logging.debug(f"raw columns: {raw.columns.tolist()}")
logging.debug(f"raw sample data:\n{raw.head()}")

prinseq=prinseq.set_index(['Sample','Type'])
logging.debug(f"prinseq columns: {prinseq.columns.tolist()}")
logging.debug(f"prinseq sample data:\n{prinseq.head()}")

bmtagger=bmtagger.set_index(['Sample','Type'])
logging.debug(f"bmtagger columns: {bmtagger.columns.tolist()}")
logging.debug(f"bmtagger sample data:\n{bmtagger.head()}")

logging.warning("Make sure these three indices match in raw, prinseq and bmtagger files because they are used in the merge step:")

logging.warning(f'Raw (Sample,Type): {raw.index[0]}')
logging.warning(f'Prinseq (Sample,Type): {prinseq.index[0]}')
logging.warning(f'Bmtagger (Sample,Type): {bmtagger.index[0]}')





if raw.isnull().values.any():
    logging.debug("There are null values in the raw DataFrame.")
else:
    logging.debug("No null values found in the raw DataFrame.")



if prinseq.isnull().values.any():
    logging.debug("There are null values in the prinseq DataFrame.")
else:
    logging.debug("No null values found in the prinseq DataFrame.")



if bmtagger.isnull().values.any():
    logging.debug("There are null values in the bmtagger DataFrame.")
else:
    logging.debug("No null values found in the bmtagger DataFrame.")



# Check the index of each DataFrame
print("bmtagger Index:", bmtagger.index)
print("prinseq Index:", prinseq.index)
print("raw Index:", raw.index)


merged=bmtagger.join(prinseq,lsuffix="_bmtagged", rsuffix="_prinseq").join(raw)

logging.debug(f"merged columns: {merged	.columns.tolist()}")
logging.debug(f"merged sample data:\n{merged.head()}")


merged=merged.reset_index()


print("=== Head of Columns 0 to 6 ===")
print(merged.iloc[:, 0:7].head(7))
print("\n")



print("=== Head of Columns 7 to 13 ===")
print(merged.iloc[:, 7:14].head(7))
print("\n")



print("=== Head of Columns 14 to 19 ===")
print(merged.iloc[:, 14:20].head(7))
print("\n")



print("=== Head of Columns 20 to 25 ===")
print(merged.iloc[:, 20:25].head(7))
print("\n")



merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')] = merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')].applymap(lambda x: str(x).replace(',', '')).astype(float)


merged.loc[:,'Host_contamination']=(merged.loc[:,'num_seqs_prinseq']-merged.loc[:,'num_seqs_bmtagged'])



print("=== Head of Columns 0 to 6 ===")
print(merged.iloc[:, 0:7].head(7))
print("\n")



print("=== Head of Columns 7 to 13 ===")
print(merged.iloc[:, 7:14].head(7))
print("\n")



print("=== Head of Columns 14 to 19 ===")
print(merged.iloc[:, 14:20].head(7))
print("\n")



print("=== Head of Columns 20 to 26 ===")
print(merged.iloc[:, 20:28].head(7))
print("\n")


merged.loc[:,'Clean_reads']=merged.loc[:,'num_seqs_bmtagged']


merged.loc[:,'Low_quality_reads']=merged.loc[:,'num_seqs']-merged.loc[:,'num_seqs_prinseq']


print("=== Head of Columns 0 to 6 ===")
print(merged.iloc[:, 0:7].head(7))
print("\n")



print("=== Head of Columns 7 to 13 ===")
print(merged.iloc[:, 7:14].head(7))
print("\n")



print("=== Head of Columns 14 to 19 ===")
print(merged.iloc[:, 14:20].head(7))
print("\n")



print("=== Head of Columns 20 to 29 ===")
print(merged.iloc[:, 20:29].head(7))
print("\n")



merged.loc[:,['Sample','Type','Low_quality_reads','Clean_reads','Host_contamination','num_seqs','num_seqs_prinseq','num_seqs_bmtagged']].to_csv(snakemake.output[0],index=False)
