#!/usr/bin/env bash

source ../venv/bin/activate

# Download GFF files containing gene definitions. The script expects
# access to Google Drive. If this is not the case refer to GFF file
# in the publication - Additional data file 4.

python download_gff_files.py


# Once the GFF file are downloaded, do pairwise alignments. This step
# needs additional annotation documents downloaded from Google Sheets.
# Without access refer to supplementary tables. Finaly, reference 
# assemblies are also needed, update the path to their location.


python analyze_array_palindrome.py
python permutation_test.py
