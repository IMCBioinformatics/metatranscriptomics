#!/usr/bin/env python3

import argparse
import os
import sys
import re
import pandas as pd
from itertools import takewhile

def merge(aaastrIn, ostm):

    merged_tables = pd.DataFrame()

    for f in aaastrIn:

        # --- Step 1: find the header line starting with #clade_name ---
        header_line = None
        mpaVersion = None
        with open(f) as fin:
            for line in fin:
                # capture MetaPhlAn DB version if present
                if line.startswith('#mpa_v'):
                    mpaVersion = line.strip()
                if line.startswith("#clade_name"):
                    header_line = line.strip()
                    break

        if header_line is None:
            raise RuntimeError(f"No header starting with #clade_name found in {f}")

        # --- Step 2: remove the leading '#' and extract all column names ---
        names = header_line.lstrip("#").strip().split("\t")

        # --- Step 3: read only required two columns ---
        required_cols = ["clade_name", "relative_abundance"]
        for col in required_cols:
            if col not in names:
                raise RuntimeError(f"Column {col} not found in file {f}")

        idx_clade = names.index("clade_name")
        idx_relab = names.index("relative_abundance")
        usecols = [idx_clade, idx_relab]

        # Read data: skip all lines starting with '#', use our column names
        iIn = pd.read_csv(
            f,
            sep="\t",
            comment="#",      # ignore all header lines
            header=None,      # no header in data block
            names=names,      # assign the names we parsed
            usecols=usecols,  # only clade_name + relative_abundance
        ).fillna("")

        # Set clade_name as index
        iIn = iIn.set_index("clade_name")

        # Choose sample name
        sample_name = os.path.splitext(os.path.basename(f))[0]
        this_sample = iIn.iloc[:, 0].rename(sample_name).to_frame()

        # Debug line (optional)
        print(f"Adding sample column {sample_name} from file {f} (shape {this_sample.shape})",
              file=sys.stderr)

        if merged_tables.empty:
            merged_tables = this_sample
        else:
            merged_tables = merged_tables.join(this_sample, how="outer")

    # Write output (no version, unless you decide to aggregate mpaVersion values)
    merged_tables.fillna("0").reset_index().to_csv(
        ostm, index=False, sep="\t"
    )


argp = argparse.ArgumentParser( prog = "merge_metaphlan_tables_relab.py",
    description = """Performs a table join on one or more metaphlan output files.""")
argp.add_argument( "aistms",    metavar = "input.txt", nargs = "+",
    help = "One or more tab-delimited text tables to join" )
argp.add_argument( '-o',    metavar = "output.txt", nargs = 1,
    help = "Name of output file in which joined tables are saved" )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

argp.usage = argp.format_usage()[7:]+"\n\n\tPlease make sure to supply file paths to the files to combine. If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n\n\t\tpython merge_metaphlan_tables_relab.py Table1.txt Table2.txt Table3.txt > output.txt\n\n\tA wildcard to indicate all .txt files that start with Table can be used as follows:\n\n\t\tpython merge_metaphlan_tables_relab.py Table*.txt > output.txt"


def main( ):
    args = argp.parse_args( )
    if args.o is None:
        merge(args.aistms, sys.stdout)
    else:
        with open(args.o[0], 'w') as fout:
            merge(args.aistms, fout)

if __name__ == '__main__':
    main()
