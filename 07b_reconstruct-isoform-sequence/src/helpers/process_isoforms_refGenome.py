#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict, Counter
import traceback
import pysam
import re
import json

def get_jacard_similarity(intervals_A, intervals_B):
    if intervals_A[0][0] > intervals_B[0][0]:
        if intervals_A[0][0] < intervals_B[0][1]:
            intervals_B[0] = (intervals_A[0][0], intervals_B[0][1])

    if intervals_B[0][0] > intervals_A[0][0]:
        if intervals_B[0][0] < intervals_A[0][1]:
            intervals_A[0] = (intervals_B[0][0], intervals_A[0][1])

    if intervals_A[-1][1] < intervals_B[-1][1]:
        if intervals_A[-1][1] > intervals_B[-1][0]:
            intervals_B[-1] = (intervals_B[-1][0], intervals_A[-1][1])

    if intervals_B[-1][1] < intervals_A[-1][1]:
        if intervals_B[-1][1] > intervals_A[-1][0]:
            intervals_A[-1] = (intervals_A[-1][0], intervals_B[-1][1])

    def segments_to_positions(segments):
        positions = set()
        for start, end in segments:
            positions.update(range(start, end))
        return positions

    pos_A = segments_to_positions(intervals_A)
    pos_B = segments_to_positions(intervals_B)

    intersection = len(pos_A & pos_B)
    union = len(pos_A | pos_B)

    if union == 0:
        return 0.0
    
    return intersection / union

def sam_to_interval(alignment):
    """Convert a SAM alignment to a list of [start, end] intervals."""
    intervals = []
    if alignment.is_unmapped:
        return intervals

    # Parse CIGAR to get spliced intervals
    ref_pos = alignment.reference_start
    exon_length = 0
    for op, length in alignment.cigartuples or []:
        if op == 0: #Match/mismatch
            exon_length += length
        elif op == 2:
            exon_length += length
        elif op == 3:
            intervals.append([ref_pos, ref_pos + exon_length])
            ref_pos += length + exon_length
            exon_length = 0

    intervals.append([ref_pos, ref_pos + exon_length])
    return intervals

def intervals_to_splice_pattern(intervals):
    """Convert list of intervals to splice pattern (remove first and last coordinate)"""
    if not intervals:
        return tuple()

    modified_intervals = []

    if len(intervals) == 1:
        #Single interval - no structure to compare
        modified_intervals = [(intervals[0][0], intervals[0][1])]
        return modified_intervals

    for i, interval in enumerate(intervals):
        if i == 0:
            modified_intervals.append((None, interval[1]))
        elif i == len(intervals) -1:
            modified_intervals.append((interval[0], None))
        else:
            modified_intervals.append(tuple(interval))

    return tuple(modified_intervals)

def read_signature_file(signature_file):
    """read signatrues identified for each individual read"""

    signatures = {}

    try:
        with open(signature_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                try:
                    parts = line.split('\t')
                    if len(parts) < 3:
                        print(f"Warning: Line {line_num} has fewer than 3 columns, skipping: {line}")
                        continue
                    
                    key = parts[0]
                    value = parts[2]
                    signatures[key] = value

                except Exception as e:
                    if line.startswith("#"):
                        pass
                    else:
                        raise Exception(e)

        print(f"Loaded {len(signatures)} signatures from {signature_file}")
        return signatures

    except FileNotFoundError:
        print(f"Error: Signature file {signature_file} not found!")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading signature file: {e}")
        sys.exit(1)

def get_substring(pattern, string):
    result = ""
    match = re.search(pattern, string)
    if match:
        extracted = match.group(0)
        result = extracted
    return result

def read_isoforms(isoform_file):
    """Read isoforms identified by StringTie"""

    isoforms = []

    try:
        with open(isoform_file, 'r') as f:

            transcript_id = ""
            isoform = ""
            intervals = []
            for line in f:
                line = line.strip()

                try:
                    parts = line.split('\t')
                    if parts[2] == "transcript":
                        print("trans")
                        if intervals != []:
                            
                            if isoform != "NA":
                                splice_pattern = intervals_to_splice_pattern(intervals)
                                isoforms.append({
                                    "id": transcript_id, 
                                    "isoform": isoform, 
                                    "intervals": intervals, 
                                    "splice_pattern": splice_pattern})

                            intervals = []
                            transcript_id = ""
                            isoform = ""
                        isoform = get_substring(r"[A-Z]+\d{0,1}_isoform_\d+", parts[8])
                        if isoform == "":
                            isoform = "NA"
                        transcript_id = get_substring(r"TCONS_\d+", parts[8])
                        if transcript_id == "":
                            raise Exception("Error: did not match transcript identifier")


                    if parts[2] == "exon": 
                        # first coordinate is -1 to get BED like 0 based first interval
                        intervals.append((int(parts[3]) -1, int(parts[4])))
                except Exception as e:
                    if line.startswith("#"):
                        pass
                    else:
                        raise Exception(e)

            
            if isoform != "NA":
                splice_pattern = intervals_to_splice_pattern(intervals)
                isoforms.append({
                    "id": transcript_id, 
                    "isoform": isoform, 
                    "intervals": intervals, 
                    "splice_pattern": splice_pattern})

    except FileNotFoundError:
        print(f"Error: Isoforms file {isoform_file} not found!")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading isoforms file: {e}")
        sys.exit(1)

    print(isoforms)
    return(isoforms)

def compare_intervals(intervals_A, intervals_B):
    """ Compare two intervals for exact match """
    
    if len(intervals_A) != len(intervals_B):
        return False

    for index, element in enumerate(intervals_A):
        if intervals_A[index][0] != intervals_B[index][0]:
            return False
        if intervals_A[index][1] != intervals_B[index][1]:
            return False
    
    return True


def process_read_alignments(bam_file, signatures, isoforms, log_file, jacard):
    """Process read alignments and match to transcript intervals."""

    transcript_to_reads = defaultdict(list)
    read_ids = set()
    matched_reads = 0
    total_reads = 0
    discarded_reads = {}
    blacklist = []

    reads_to_isoform = defaultdict(str)
    

    # print(f"\t{'\t'.join([x['isoform'] for x in isoforms])}")
    # for isoform_A in isoforms:
    #     print(f"{isoform_A['isoform']}", end = "\t")
    #     for isoform_B in isoforms:
    #         jacard = get_jacard_similarity(isoform_A['intervals'], isoform_B['intervals'])
    #         print(f"{jacard}", end= "\t")

    #     print()



    with pysam.AlignmentFile(bam_file, "r") as samfile:
        for alignment in samfile:
            
            if alignment.is_unmapped:
                discarded_reads[alignment.query_name] = "unmapped read"
                continue

            signature = signatures.get(alignment.query_name)
            if not signature:
                discarded_reads[alignment.query_name] = "signature not found"
                continue
            
            total_reads += 1
            read_id = alignment.query_name
            if (read_id in blacklist):
                continue
            read_ids.add(read_id)
            read_intervals = sam_to_interval(alignment)

            # print(f"{read_id}: {read_intervals}")

            if not read_intervals:
                discarded_reads[read_id] = "no intervals found"
                continue

            read_splice_pattern = intervals_to_splice_pattern(read_intervals)

            # best_jacard = {'val':0.0, 'isoform':None, 'intervals':[]}
            for isoform in isoforms:

                if isoform["isoform"] == "":
                    print(isoform)
                    assert False


                
                # if (jacard):
                #     jacard_value = get_jacard_similarity(isoform["intervals"], read_intervals)
                #     if (jacard_value > best_jacard['val']):
                #         best_jacard['val'] = jacard_value
                #         best_jacard['isoform'] = isoform["isoform"]
                #         best_jacard['intervals'] = isoform["intervals"]

                match = compare_intervals(isoform["splice_pattern"], read_splice_pattern)
                if len(isoform['intervals']) == 1:
                   jacard_value = get_jacard_similarity(isoform['intervals'], read_intervals)
                   if jacard_value > 0.9:
                        if reads_to_isoform[read_id] != "":
                            if reads_to_isoform[read_id] != isoform["isoform"]:
                                print(f"Warning: Multiple isoforms match read {read_id}: {reads_to_isoform[read_id]} vs {isoform['isoform']}\n{read_splice_pattern} v {isoform['splice_pattern']}")
                                reads_to_isoform[read_id] = ""
                                matched_reads -= 1
                                read_ids.remove(read_id)
                                blacklist.append(read_id)
                                # raise Exception(f"Multiple isoforms match read {read_id}: {reads_to_isoform[read_id]} vs {isoform['isoform']}\n{read_splice_pattern} v {isoform["splice_pattern"]}")
                                continue   
   

                        reads_to_isoform[read_id] = isoform["isoform"]
                        matched_reads += 1
                        continue 
                elif match == True:
                    if reads_to_isoform[read_id] != "":
                        
                        if reads_to_isoform[read_id] != isoform["isoform"]:
                            print(f"Warning: Multiple isoforms match read {read_id}: {reads_to_isoform[read_id]} vs {isoform['isoform']}\n{read_splice_pattern} v {isoform['splice_pattern']}")
                            reads_to_isoform[read_id] = ""
                            matched_reads -= 1
                            read_ids.remove(read_id)
                            blacklist.append(read_id)
                            # raise Exception(f"Multiple isoforms match read {read_id}: {reads_to_isoform[read_id]} vs {isoform['isoform']}\n{read_splice_pattern} v {isoform["splice_pattern"]}")
                        continue    

                    reads_to_isoform[read_id] = isoform["isoform"]
                    matched_reads += 1
                    continue 
                else:
                    discarded_reads[read_id] = "no mathing transcirpt pattern"

            # assert False
               


    print(f"Matched {matched_reads} reads out of {len(read_ids)} total reads")
    print(f"Discarded {len(discarded_reads)} reads")

    return(reads_to_isoform)


def process_multiple_read_files(bam_files, signatures, isoforms, log_file, jacard):

    all_reads_to_isoforms = []

    for i, bam in enumerate(bam_files):
        print(f"\n=== Processing read file {i+1}/{len(bam_files)}: {bam} ===")
        reads_to_isoform = process_read_alignments( bam, signatures, isoforms, log_file, jacard)
        all_reads_to_isoforms.append(reads_to_isoform)

    return all_reads_to_isoforms

def process_reference_alignments(reference_file, isoforms, signatures):
    """Process reference alignments and match to transcript intervals."""
    print(f"Processing reference alignments from {reference_file}")

    transcript_to_refs = defaultdict(list)
    matched_refs = 0
    total_refs = 0
    discarded_refs = {}

    with pysam.AlignmentFile(reference_file, "r") as samfile:
        for alignment in samfile:

            if alignment.is_unmapped:
                discarded_refs[alignment.query_name] = "unmapped reference"
                continue
            
            signature = signatures.get(alignment.query_name)
            if not signature:
                discarded_refs[alignment.query_name] = "signature not found"
                continue

            total_refs += 1
            ref_id = alignment.query_name
            ref_intervals = sam_to_interval(alignment)

            if not ref_intervals:
                discarded_refs[ref_id] = "no intervals found"
                continue

            ref_splice_pattern = intervals_to_splice_pattern(ref_intervals)

            matched = False

            for isoform  in isoforms:
                if compare_intervals(ref_splice_pattern, isoform['splice_pattern']):
                    transcript_to_refs[isoform["isoform"]].append(ref_id)
                    matched_refs += 1
                    matched = True
                    break

            if not matched:
                discarded_refs[ref_id] = "no matching transcript splice pattern"

    print(f"Matched {matched_refs} references out of {total_refs} total references")
    print(f"Discarded {len(discarded_refs)} references")
    print(f"Found matches for {len(transcript_to_refs)} transcripts")

    return transcript_to_refs



def create_signature_transcript_table(all_reads_to_isoforms, transcript_to_refs, signatures, bam_files, output):
    """Create a table with signatures as rows and transcripts as columns, with counts per file.
    Additionally log read - signature - isoform information."""
    

    uniq_signatures = set(signatures.values())
    uniq_signatures_counter = defaultdict(int)
    

    all_isoforms = set()
    for transcripts_to_reads in all_reads_to_isoforms:

        all_isoforms.update(transcripts_to_reads.values())

    uniq_isoforms = sorted(all_isoforms)

    # Create signature-transcript count matrix per file
    # Structure: signature_counts[signature][transcript_id] = [count_file1, count_file2, ...]
    signature_counts = defaultdict(lambda: defaultdict(lambda: [0]* len(bam_files)))
    signature_id_list = defaultdict(lambda: defaultdict(lambda: [[]]* len(bam_files)))
    
    with open(output, "w") as table_file:

        for file_idx, isoforms in enumerate(all_reads_to_isoforms):
            for read_id, iso_id in isoforms.items():
                if read_id in signatures:
                    signature = signatures[read_id]
                    signature_counts[signature][iso_id][file_idx] += 1
                    signature_id_list[signature][iso_id][file_idx] = signature_id_list[signature][iso_id][file_idx] + [read_id]
                    uniq_signatures_counter[signature] += 1
                    print(f"{signature}\t{read_id}\t{iso_id}", file=table_file)
                else:
                    print(f"Warning: Read {read_id} not found in signatures file")

    uniq_signatures_counter_sorted = sorted(uniq_signatures, key = lambda sig: uniq_signatures_counter[sig], reverse = True)

    # Create signture-transcript reference count matrix
    # Structure: signature_ref_counts[signature][transcript_id] = count
    signature_ref_counts = defaultdict(lambda: defaultdict(int))

    for transcript_id, ref_ids in transcript_to_refs.items():
        for ref_id in ref_ids:
            if ref_id in signatures:
                signature = signatures[ref_id]
                signature_ref_counts[signature][transcript_id] += 1
            else:
                print(f"Warning: Reference {ref_id} not found in signatures file")

    return uniq_signatures_counter_sorted, uniq_isoforms, signature_counts, signature_ref_counts, signature_id_list
        # all_isoforms.update()

def main():
    parser = argparse.ArgumentParser(description="Match read to isoform.")
    parser.add_argument("bam_files", nargs="+", help="SAM/BAM files with read alignments")
    parser.add_argument("-i", "--isoforms", help="StringTie identified isoforms", required = True)
    parser.add_argument("-s", "--signatures", help="Tab-delimited file with sequence signatures", required = True)
    parser.add_argument("-r", "--reference", help="SAM/BAM file containing reference sequence alignments", required = True)
    parser.add_argument("-o", "--output", help="Output file for resutls (default: stdout)", default = None)
    parser.add_argument("-u", "--used", help="Collect (used) read ID's with specific isoform and signature combination", default = None)
    parser.add_argument("-d", "--discarded", help="log discarded reads")
    parser.add_argument("-l", "--log", help="Base filename for log files.", default=None)
    # parser.add_argument("--jacard", action="store_true", help="use jacard distance to allow imprecise match of isoforms (set to 0.99)")

    args = parser.parse_args()

    print(f"Extract isoforms from: {args.isoforms}")
    isoforms = read_isoforms(args.isoforms)

    print(f"Extract signatures from: {args.signatures}")
    signatures = read_signature_file(args.signatures)

    print(f"Processing {len(args.bam_files)} read files:")
    for i, bam_file in enumerate(args.bam_files, 1):
        print(f" {i}, {bam_file}")

    try: 
        all_reads_to_isoforms = process_multiple_read_files(
            args.bam_files,signatures, isoforms, args.log, False
        )

        transcript_to_refs = defaultdict(list)
        if args.reference:
            transcript_to_refs = process_reference_alignments(
                args.reference, isoforms, signatures
            )

        # print(all_reads_to_isoforms)
        uniq_signatures, uniq_isoforms, signature_counts, signature_ref_counts, signature_id_list = create_signature_transcript_table(
            all_reads_to_isoforms, transcript_to_refs, signatures, args.bam_files, args.used
        )
    
         # Calculate column totals first to sort transcript_ids by descending totals
        column_totals_dict = {}
        for transcript_id in uniq_isoforms:
            column_total = 0
            for signature in uniq_signatures:
                counts = signature_counts[signature][transcript_id]
                x_value = sum(counts)
                column_total += x_value
            column_totals_dict[transcript_id] = column_total

        # Sort transcript_ids by column totals in descending order
        sorted_isoforms = sorted(uniq_isoforms, key=lambda tid: column_totals_dict[tid], reverse=True)
    

        # Create sorted column_totals array
        sorted_column_totals = [column_totals_dict[tid] for tid in sorted_isoforms]

        output_file = open(args.output, 'w') if args.output else sys.stdout
        try:

            output_file.write("Signature\t" + "\t".join(sorted_isoforms) + "\tRow_Total\n")

            # Write column totals row
            column_total_row = ["Column_Total"] + [str(total) for total in sorted_column_totals]
            column_total_row.append(str(sum(sorted_column_totals)))  # Grand total
            output_file.write("\t".join(column_total_row) + "\n")

            for signature in uniq_signatures:
                    row = [signature]
                    row_total = 0
                    for iso in sorted_isoforms:
                        counts = signature_counts[signature][iso]
                        ref_count = signature_ref_counts[signature][iso]
                        # formatted_count = format_count_cell(counts, ref_count)
                        ref_tag = ""
                        if ref_count > 0:
                            ref_tag = f"+{ref_count}"
                        row.append(f"{sum(counts)}({','.join([str(x) for x in counts])}){ref_tag}")
                        # Add the X value (sum of counts) to row total
                        x_value = sum(counts)
                        row_total += x_value
                    row.append(str(row_total))
                    output_file.write("\t".join(row) + "\n")

        finally:
            output_file.close()

        if args.log:
            read_id_matrix = f"{args.log}_read_id_matrix.json"
            with open(read_id_matrix, 'w') as out_json:
                json.dump(signature_id_list, out_json, indent = 4)



    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()