#!/usr/bin/env python3
"""
Parse alignment file (BAM file) and collect "signature" for every read.
Variant positions read from external file. JSON format is expected with
the following stucture:

[
    {"pos": 739,  "source": "var_call_hu2", "alt1": "A", "alt2": "G", "type":"SNV", "found" : "."},
    {"pos": 816,  "source": "ref_assembly", "alt1": "T", "alt2": "TGGCGGAGGTGGAGGTGG", "type":"INS", "found" : "-"},
]

pos: integer, position relative to single gene reference.
source: informative, why this position was added to the lsit
alt1: reference allele
alt2: variant allele
type: [SNV, INS]
found: string set to initial value "." for SNVs and "-" for INS.
       filled out by the script.

"""

import pysam
import json
import argparse
import copy
import sys
from typing import List, Dict

# Constants for CIGAR operations
CIGAR_MATCH = 0
CIGAR_INSERTION = 1
CIGAR_DELETION = 2
CIGAR_SKIP = 3
CIGAR_SOFT_CLIP = 4
CIGAR_HARD_CLIP = 5


def parse_variants( variants_file: str) -> (List[Dict], List[int]):
    """Parse the variants JSON file and return variants and their positions."""
    variant_positions = []
    variants_json = []
    with open(variants_file, 'r') as file:
        variants_json = json.load(file)
        for variant in variants_json:
            variant['pos'] = variant['pos'] - 1 #convert to 0 based
            variant_positions.append(variant['pos'])
            del variant
        variant_positions.sort()

    return variants_json, variant_positions


def process_read(read, local_variants, variant_positions, report_all):
    """Process a single read and update local_variants with alleles."""
    ref_pos = read.reference_start
    query_pos = 0

    if read.is_unmapped:
        return None

    if len(variant_positions) < 1:
        return None

    if (not report_all) and ref_pos > variant_positions[0]:
        return None

    for cigar_op, length in read.cigartuples:
        if cigar_op == CIGAR_MATCH:
            for _ in range(length):
                if ref_pos in variant_positions:
                    nucleotide = read.query_sequence[query_pos]
                    for x in local_variants:
                        if x["pos"] == ref_pos and x["type"] == "SNV":
                            x['found'] = nucleotide

                ref_pos += 1
                query_pos += 1
        elif cigar_op == CIGAR_INSERTION:
            if ref_pos in variant_positions:
                for x in local_variants:
                    if x['pos'] == ref_pos and x['type'] == "INS" and length == len(x['alt2']):
                        x['found'] = "I"
            query_pos += length

        elif cigar_op == CIGAR_DELETION:
            ref_pos += length
        elif cigar_op == CIGAR_SKIP:
            ref_pos += length
        elif cigar_op == CIGAR_SOFT_CLIP:
            query_pos += length
        elif cigar_op == CIGAR_HARD_CLIP:
            continue

    return "".join([x['found'] for x in local_variants])


def main():
    parser = argparse.ArgumentParser(description="Collect uniq combinations of SNPs in a bam file")
    parser.add_argument("-i", "--input", help="input bam file", required=True)
    parser.add_argument("-v", "--variants", help="list of positions where SNPs are located", required=True)
    parser.add_argument("-o", "--json_out", help="output json")
    parser.add_argument("-t", "--tab_out", help="tab output with all sequenes")
    parser.add_argument("-f", "--filter_FP", help="if provided will write a filtered SAM file to stdout (filtered for all reads that contain given signature)")
    parser.add_argument("-r", "--report_all_positions", action='store_true', default=False)

    args = parser.parse_args()

    variants, variant_positions = parse_variants(args.variants)
    signature_counts = {}

    bam_file = args.input

    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        with (open(args.tab_out,"w") if args.tab_out is not None else sys.stderr) as TAB_FILE:

            if args.filter_FP != None:
                print(bamfile.header, end="")

            for read in bamfile:

                quals = read.query_qualities
                quality_string = ""
                if quals is not None:
                    quality_string = ''.join(chr(q + 33) for q in quals)

                local_variants = copy.deepcopy(variants)
                signature = process_read(read, local_variants, variant_positions, args.report_all_positions)

                if signature is None:
                    continue

                if args.filter_FP == signature:
                    print(read.to_string())

                print(f"{read.query_name}\t{read.query_sequence}\t{signature}\t{quality_string}", file=TAB_FILE)

                signature_counts[signature] = signature_counts.get(signature,0) + 1

    if args.json_out:
        with open(args.json_out, 'w') as outfile:
            json.dump(signature_counts, outfile)
    else:
        print(signature_counts, file=sys.stderr)

if __name__ == "__main__":
    main()
