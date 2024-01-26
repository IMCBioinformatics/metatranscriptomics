#!/usr/bin/python
import os
import sys
import re
import csv
import argparse
import json
import glob

### Sample command
# python sgb_to_gtdb_profile_abundances.py --metaphlan_SGB_profile_infile output_with_controls/metaphlan/2016_8622_S18_profile.txt --sgb_to_gtdb_tsv_file mpa_vOct22_CHOCOPhlAnSGB_202212_SGB2GTDB.tsv --output_dir GTDB_profiles

parser = argparse.ArgumentParser()

metaphlan_SGB_profile_infile = None
sgb_to_gtdb_tsv_file = None
output_dir = None

parser.add_argument('--metaphlan_SGB_profile_infile', action='store', dest='metaphlan_SGB_profile_infile',
                    help='input metaphlan SGB profile file as input. (i.e. sample_name_profile.txt)')
parser.add_argument('--sgb_to_gtdb_tsv_file', action='store', dest='sgb_to_gtdb_tsv_file',
                    help='input metaphlan SGB to GTDB taxonomy file as input. (i.e. mpa_vOct22_CHOCOPhlAnSGB_202212_SGB2GTDB.tsv)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

metaphlan_SGB_profile_infile = results.metaphlan_SGB_profile_infile
sgb_to_gtdb_tsv_file = results.sgb_to_gtdb_tsv_file
output_dir = results.output_dir

if(metaphlan_SGB_profile_infile == None):
	print('\n')
	print('error: please use the --metaphlan_SGB_profile_infile option to specify the input metaphlan SGB profile file directory as input')
	print('metaphlan_SGB_profile_infile =' + ' ' + str(metaphlan_SGB_profile_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(sgb_to_gtdb_tsv_file == None):
    print('\n')
    print('error: please use the --sgb_to_gtdb_tsv_file option to specify the input metaphlan SGB to GTDB taxonomy file as input')
    print('sgb_to_gtdb_tsv_file =' + ' ' + str(sgb_to_gtdb_tsv_file))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

# Create the output directory if it does not exist.
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# Sort input list using natural alphanumeric sort.
def alphanumeric_sort(list):

    # Sort the given iterable in the way that humans expect.
    convert = lambda text: int(text) if text.isdigit() else text
    
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    
    return sorted(list, key = alphanum_key)

#def merge_metaphlan_tables(metaphlan_profile_dir, merged_metaphlan_table_file):

# The SGB ID to GTDB taxonomy dictionary.
sgb_to_gtdb_dict = {}

i = 0
with open(sgb_to_gtdb_tsv_file, "r") as sgb_to_gtdb_tsv_input_file:
    csv_reader = csv.reader(sgb_to_gtdb_tsv_input_file, delimiter='\t')
    for row in csv_reader:

#        print(row)
        if(i != 0):
            (sgb_id, gtdb_lineage) = row
            sgb_to_gtdb_dict[sgb_id] = gtdb_lineage
        i += 1

#print(sgb_to_gtdb_dict)
#sys.exit()

 # {'0': 'clade_name', '1': 'clade_taxid', '2': 'relative_abundance', '3': 'coverage', '4': 'estimated_number_of_reads_from_the_clade'}

# The list of clade names.
clade_names_list = []

# The metaphlan profile file dictionary.
metaphlan_profile_dict = {}

# The metaphlan profile input file basename.
basename = os.path.basename(metaphlan_SGB_profile_infile)

# Filename without the extension ".txt".
filename = os.path.splitext(basename)[0]

#        print(filename)

# Open the metaphlan profile input file.
metaphlan_profile_filehandle = open(metaphlan_SGB_profile_infile, "r")

# The comment lines count to use as an index for the comments list containing all the comment lines.
comment_lines_count = 0

# the comments list to store all the comment lines.
comments_list = []

# The header row dictionary to get the column names and the location index positioning in the row.
header_dict = {}

# The row count number to count the number of rows. More specifically, when zero we need to find the header comment line and process the row to obtain all the column names for the row dictionary row_dict.
row_count = 0

# Iterate through the metaphlan profile file for sample_name.
for line in metaphlan_profile_filehandle.readlines():

    # Search for the comment lines first so that we can figure out which is the header line and where the data starts.
    if(re.search("^#", line)):
    
        comment_lines_count += 1
#                print(line)
        
        # Get the comment line and remove the newline.
        comment_line = line.replace("\n", "")
        
        # Append the comment line to the comments list.
        comments_list.append(comment_line)
        
    else: # Data starts at this line number. We want to get the header line from the comment_line list.
    
        # Get the row and remove the newline.
        row = line.replace("\n", "")
        if(row_count == 0):
        
            header = comments_list[comment_lines_count-1]
            #print(header.replace("#", "").split("\t"))
            column_list = header.replace("#", "").split("\t")
            column_count = 0
            
            # Order of the SGB header.
            # {'0': 'clade_name', '1': 'clade_taxid', '2': 'relative_abundance', '3': 'coverage', '4': 'estimated_number_of_reads_from_the_clade'}
            for column_name in column_list:
                header_dict[str(column_count)] = str(column_name)
                column_count += 1
     
#                print(header_dict)
        column_index = 0
#                print(row)
        
        # For each row entry make a row dictionary with keys for each column based on the header index.
        row_dict = {}
        for column_value in row.split("\t"):
            if(column_index < len(header_dict)):
                column_name = header_dict[str(column_index)]
                row_dict[str(column_name)] = column_value

            column_index += 1

        # Get the clade name from the row dictionary.
        clade_name = row_dict["clade_name"]
        
        # If the clade name is "UNCLASSIFIED" and contains the SGB ID "t__SGB"
        if((clade_name == "UNCLASSIFIED") or (re.search("t__SGB", clade_name))):
        
            # Append the clade names to a list to sort.
            clade_names_list.append(clade_name)
            
            # Initialize the metaphlan_profile_dict and keys and assign rows for each file based on clade_name.
            if(not(str(clade_name) in metaphlan_profile_dict)):
            
                metaphlan_profile_dict[str(clade_name)] = {}
                
                # Storing the row_dict for each row at the clade_name, filename keypair of the metaphlan_profile_dict dictionary for the first entry after initalization.
                metaphlan_profile_dict[str(clade_name)] = row_dict
            else:
                # Storing the row_dict for each row at the clade_name, filename keypair of the metaphlan_profile_dict dictionary
                metaphlan_profile_dict[str(clade_name)] = row_dict
                
        row_count += 1

#print(metaphlan_profile_dict)
#sys.exit()
clade_names_set_sorted = alphanumeric_sort(set(clade_names_list))
del clade_names_list

csv_writer_file_handle = open(os.path.join(output_dir, filename + ".txt"), "w+")
csv_writer = csv.writer(csv_writer_file_handle, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)

# Write the metaphlan database comment entry.
csv_writer.writerow([comments_list[0]])

# Write the metaphlan profile header comment entry.
csv_writer.writerow(["#clade_name", "estimated_number_of_reads_from_the_clade"])

metaphlan_GTDB_dict = {}
for clade_name in clade_names_set_sorted:
    if(clade_name == "UNCLASSIFIED"):
#        print(metaphlan_profile_dict[clade_name])
        if("UNCLASSIFIED" in metaphlan_profile_dict):
        
            clade_dict = metaphlan_profile_dict[clade_name]
            
            metaphlan_GTDB_dict[clade_name] = int(clade_dict["estimated_number_of_reads_from_the_clade"])
            
    elif(re.search("t__SGB", clade_name)):
        # k__Bacteria|p__Actinobacteria|c__Actinomycetia|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Gardnerella|s__Gardnerella_vaginalis|t__SGB17302
        clade_dict = metaphlan_profile_dict[clade_name]
#        print(clade_dict)
        
        sgb_id = clade_name.split("|")[-1].replace("t__", "")
        gtdb_tax_lineage = sgb_to_gtdb_dict[sgb_id]
#        print(gtdb_tax_lineage)

        gtdb_tax_level_list = gtdb_tax_lineage.split(";")
#        print(gtdb_tax_level_list)
        
        # The GTDB taxonomy lineage list.
        gtdb_taxonomy_lineage_list = []
        
        # Iterate over the GTDB taxonomy levels starting from domain "d__" to species "s__"
        # Add the abundance to the current abundance value for this taxonomic level.
        for gtdb_tax_level in gtdb_tax_level_list:
        
            # Add the GTDB taxonomy level to the lineage list.
            gtdb_taxonomy_lineage_list.append(gtdb_tax_level)
            
            # Concatenate this part of the lineage using the ";" delimiter.
            gtdb_taxonomy = ";".join(gtdb_taxonomy_lineage_list)
            
            # If the gtdb_taxonomy string is not a key in the metaphlan_GTDB_dict dictionary then set to zero 0.
            if(not(gtdb_taxonomy in metaphlan_GTDB_dict)):
                metaphlan_GTDB_dict[gtdb_taxonomy] = 0
                
            # If the gtdb_taxonomy string is a key in the metaphlan_GTDB_dict dictionary then add to the previous count.
            if(gtdb_taxonomy in metaphlan_GTDB_dict):
                metaphlan_GTDB_dict[gtdb_taxonomy] += int(clade_dict["estimated_number_of_reads_from_the_clade"])

# Delete the metaphlan profile dictionary because we do not need it anymore.
del metaphlan_profile_dict
                
# Get a list the the GTDB clade names so we can sort the list.
gtdb_clade_names_list = metaphlan_GTDB_dict.keys()
gtdb_clade_names_set_sorted = alphanumeric_sort(set(gtdb_clade_names_list))
del gtdb_clade_names_list

# Write all the new GTDB taxonomy abundance counts for the GTDB metaphlan profile.
for gtdb_clade_name in gtdb_clade_names_set_sorted:
    csv_writer.writerow([gtdb_clade_name,str(metaphlan_GTDB_dict[gtdb_clade_name])])

# Closing csv_writer_file_handle.
csv_writer_file_handle.close()

# Delete the clade names sorted set because we are finished.
del clade_names_set_sorted

# Delete the GTDB clade names sorted set because we are finished.
del gtdb_clade_names_set_sorted

# Delete the metaphlan GTDB profile dictionary because we are finished.
del metaphlan_GTDB_dict
