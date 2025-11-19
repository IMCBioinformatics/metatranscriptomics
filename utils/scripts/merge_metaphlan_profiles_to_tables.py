#!/usr/bin/python
import os
import sys
import re
import csv
import argparse
import json
import glob

# Sample command.
# python merge_metaphlan_profiles_to_tables.py --metaphlan_SGB_profile_dir output_with_controls/metaphlan --metaphlan_GTDB_profile_dir output_with_controls/metaphlan/GTDB --output_dir output_with_controls

parser = argparse.ArgumentParser()

metaphlan_SGB_profile_dir = None
metaphlan_GTDB_profile_dir = None
output_dir = None

parser.add_argument('--metaphlan_SGB_profile_dir', action='store', dest='metaphlan_SGB_profile_dir',
                    help='input metaphlan SGB profile file directory as input. (i.e. output/metaphlan)')
parser.add_argument('--metaphlan_GTDB_profile_dir', action='store', dest='metaphlan_GTDB_profile_dir',
                    help='input metaphlan GTDB profile file directory as input. (i.e. output/metaphlan/GTDB)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

metaphlan_SGB_profile_dir = results.metaphlan_SGB_profile_dir
metaphlan_GTDB_profile_dir = results.metaphlan_GTDB_profile_dir
output_dir = results.output_dir

if(metaphlan_SGB_profile_dir == None):
	print('\n')
	print('error: please use the --metaphlan_SGB_profile_dir option to specify the input metaphlan SGB profile file directory as input')
	print('metaphlan_SGB_profile_dir =' + ' ' + str(metaphlan_SGB_profile_dir))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(metaphlan_GTDB_profile_dir == None):
    print('\n')
    print('error: please use the --metaphlan_GTDB_profile_dir option to specify the input metaphlan GTDB profile file directory as input as input')
    print('metaphlan_GTDB_profile_dir =' + ' ' + str(metaphlan_GTDB_profile_dir))
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

def merge_metaphlan_tables(metaphlan_profile_dir, merged_metaphlan_table_file):
    # {'0': 'clade_name', '1': 'clade_taxid', '2': 'relative_abundance', '3': 'coverage', '4': 'estimated_number_of_reads_from_the_clade'}

    # Dictionary used to organize all the samples into a table. Dictionary where the first key is the lineage (clade_name) and the second key is the sample name. If NCBI/SGB taxonomy metaphlan the value for each key1,key2 pair is the row from that sample with the clade_name, clade_taxid, relative_abundance, coverage, and estimated_number_of_reads_from_the_clade. If GTDB metaphlan
    metaphlan_table_dict = {}

    # The list of clade names that will be used in the merging step.
    clade_names_list = []

    # The list of sample names that will be used in the merging step.
    sample_names_list = []

    # The clade taxid dictionary.
    clade_taxid_dict = {}

    # Iterate over the list of the metaphlan profile files from the metaphlan file directory.
    #metaphlan_profile_list = os.listdir(metaphlan_profile_dir)
    #for metaphlan_profile_filename in metaphlan_profile_list:
    for metaphlan_profile_infile in glob.glob(os.path.join(metaphlan_profile_dir,"*_profile.txt")):
    #    print(metaphlan_profile_infile)
    #    sys.exit()

        # The metaphlan profile input file basename.
        basename = os.path.basename(metaphlan_profile_infile)

        # Filename without the extension ".txt".
        filename = os.path.splitext(basename)[0]

        # Append the sample names to a list to sort.
        sample_names_list.append(filename)

#        print(filename)

        # Open the metaphlan profile input file.
        metaphlan_profile_filehandle = open(metaphlan_profile_infile, "r")

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
#                print(row_dict)

                # Get the clade name from the row dictionary.
                clade_name = row_dict["clade_name"]

                # Append the clade names to a list to sort.
                clade_names_list.append(clade_name)

                # If the clade_taxid column exists in the metaphlan profile file. This is the NCBI taxonomy in the SGB metaphlan profiles.
                if("clade_taxid" in row_dict):
                    clade_taxid = row_dict["clade_taxid"]
                    clade_taxid_dict[str(clade_name)] = clade_taxid

                # Initialize the metaphlan_table_dict and keys and assign rows for each file based on clade_name.
                if(not(str(clade_name) in metaphlan_table_dict)):
                    metaphlan_table_dict[str(clade_name)] = {}
                if(not(str(filename) in metaphlan_table_dict[str(clade_name)])):
                    # Need to initialize the value at the clade_name, filename keypair as a dictionary.
                    metaphlan_table_dict[str(clade_name)][str(filename)] = {}

                    # Storing the row_dict for each row at the clade_name, filename keypair of the metaphlan_table_dict dictionary for the first entry after initalization.
                    metaphlan_table_dict[str(clade_name)][str(filename)] = row_dict
                else:
                    # Storing the row_dict for each row at the clade_name, filename keypair of the metaphlan_table_dict dictionary
                    metaphlan_table_dict[str(clade_name)][str(filename)] = row_dict

                row_count += 1

    clade_names_set_sorted = alphanumeric_sort(set(clade_names_list))
    del clade_names_list

    sample_names_set_sorted = alphanumeric_sort(set(sample_names_list))
    del sample_names_list

    csv_writer_file_handle = open(os.path.join(output_dir, merged_metaphlan_table_file), "w+")
    csv_writer = csv.writer(csv_writer_file_handle, delimiter='\t', quoting=csv.QUOTE_NONE)
    if(clade_name in clade_taxid_dict):
        csv_writer.writerow(["clade_name", "clade_taxid"] + list(sample_names_set_sorted))
    else:
        csv_writer.writerow(["clade_name"] + list(sample_names_set_sorted))

    # Iterate over the clade names sorted set.
    for clade_name in clade_names_set_sorted:
#        print(clade_name)

        metaphlan_table_row = []
        metaphlan_table_row.append(clade_name)

        # If the clade_taxid exists add the clade_taxid. Not added for the GTDB output.
        if(clade_name in clade_taxid_dict):
            clade_taxid = clade_taxid_dict[clade_name]
            metaphlan_table_row.append(clade_taxid)

        # Iterate over the sample names sorted set.
        for sample_name in sample_names_set_sorted:
#            print(sample_name)

            if(clade_name in metaphlan_table_dict):
                if(sample_name in metaphlan_table_dict[clade_name]):

                    # Get the metaphlan_profile_entry using the clade name and sample name as keys.
                    metaphlan_profile_entry = metaphlan_table_dict[clade_name][sample_name]

#                    print(metaphlan_profile_entry)

                    # The SGB metaphlan profile file header order.
                    # 'clade_name', 'clade_taxid', 'relative_abundance', 'coverage', 'estimated_number_of_reads_from_the_clade'

                    # Construct the metaphlan abundance tables. If "estimated_number_of_reads_from_the_clade" in metaphlan_profile_entry use abundance_count.
                    if(("estimated_number_of_reads_from_the_clade" in metaphlan_profile_entry) and ((merged_metaphlan_table_file == "metaphlan/results/merged_abundance_table_SGB.txt") or (merged_metaphlan_table_file == "metaphlan/results/merged_abundance_table_GTDB.txt"))):
#                    if(("estimated_number_of_reads_from_the_clade" in metaphlan_profile_entry) and ((merged_metaphlan_table_file == "merged_abundance_table_SGB.txt"))):
                        abundance_count = metaphlan_profile_entry["estimated_number_of_reads_from_the_clade"]
#                        print(abundance_count)
                        metaphlan_table_row.append(abundance_count)
                    # Construct the metaphlan relative abundance tables.
                    elif((merged_metaphlan_table_file == "metaphlan/results/merged_abundance_table_SGB_relab.txt") or (merged_metaphlan_table_file == "metaphlan/results/merged_abundance_table_GTDB_relab.txt")): # Otherwise "estimated_number_of_reads_from_the_clade" not in metaphlan_profile_entry use "relative_abundance" instead.
                        relative_abundance = metaphlan_profile_entry["relative_abundance"]
#                        print(relative_abundance)
                        metaphlan_table_row.append(relative_abundance)

                else: # Otherwise clade_name not in sample_name assign to zero count 0.
                    metaphlan_table_row.append("0")
    #        sample_count += 1

        # Write the entire row with all the metaphlan information.
        csv_writer.writerow(metaphlan_table_row)
        #print("\t".join(metaphlan_table_row))
        #sys.exit()

    csv_writer_file_handle.close()

    # Delete the clade names sorted set because we are finished.
    del clade_names_set_sorted

    # Delete the clade names sorted set because we are finished.
    del sample_names_set_sorted

    # Delete the clade names sorted set because we are finished.
    del metaphlan_table_dict

### END FUNCTION ###

### MAIN ###

## Merge the metaphlan profile tables.

# Merge the metaphlan SGB taxonomy profiles into a large abundance table with all samples.
merge_metaphlan_tables(metaphlan_SGB_profile_dir, "metaphlan/results/merged_abundance_table_SGB.txt")

# Merge the metaphlan GTDB taxonomy profiles into a large abundance table with all samples.
merge_metaphlan_tables(metaphlan_GTDB_profile_dir, "metaphlan/results/merged_abundance_table_GTDB.txt")

# Merge the metaphlan SGB taxonomy profiles into a large relative abundance table with all samples.
#merge_metaphlan_tables(metaphlan_SGB_profile_dir, "metaphlan/results/merged_abundance_table_SGB_relab.txt")

# Merge the metaphlan GTDB taxonomy profiles into a large relative abundance table with all samples.
#merge_metaphlan_tables(metaphlan_GTDB_profile_dir, "metaphlan/results/merged_abundance_table_GTDB_relab.txt")
