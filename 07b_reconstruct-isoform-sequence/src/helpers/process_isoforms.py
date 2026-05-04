#!/usr/bin/env python3
"""
Script to process transcript alignments and match reads to transcript intervals.
Modified to accept multiple read files and output in X(y1,..,yn) format.
"""

import argparse
import sys
from collections import defaultdict
import pysam
import traceback
import json


def sam_to_intervals(alignment):
    """Convert a SAM alignment to a list of [start, end] intervals."""
    intervals = []
    if alignment.is_unmapped:
        return intervals

    # Parse CIGAR to get spliced intervals
    ref_pos = alignment.reference_start
    exon_length = 0
    for op, length in alignment.cigartuples or []:
        if op == 0:  # Match/mismatch
            exon_length += length
        elif op == 2:  # Deletion
            exon_length += length
        elif op == 3:  # Skipped region (splice junction)
            intervals.append([ref_pos, ref_pos + exon_length])
            ref_pos += length + exon_length
            exon_length = 0

    intervals.append([ref_pos, ref_pos + exon_length])
    return intervals

def intervals_to_tuple(intervals):
    """Convert intervals list to tuple for hashing."""
    return tuple(tuple(interval) for interval in intervals)

def intervals_to_splice_pattern(intervals):
    """Convert intervals to splice pattern, ignoring first start and last end positions."""
    if not intervals:
        return tuple()

    if len(intervals) == 1:
        # Single interval - no internal structure to compare
        modified_intervals = [(intervals[0][0], intervals[0][1])]
        return tuple(modified_intervals)

    # Create modified intervals: keep first interval's end, keep last interval's start
    # Keep all middle intervals unchanged
    modified_intervals = []

    for i, interval in enumerate(intervals):
        if i == 0:
            # First interval: ignore start, keep end
            modified_intervals.append((None, interval[1]))
        elif i == len(intervals) - 1:
            # Last interval: keep start, ignore end
            modified_intervals.append((interval[0], None))
        else:
            # Middle intervals: keep both start and end
            modified_intervals.append(tuple(interval))

    return tuple(modified_intervals)

def process_transcript_alignments(transcript_file, signature_json, log_file=None):
    """Process transcript alignments and return unique interval sets."""
    print(f"Processing transcript alignments from {transcript_file}")

    transcript_intervals = {}
    interval_counts = defaultdict(int)
    discarded_transcripts = {}
    tuple_to_id = {}
    first_coord = None
    last_coord = None

    # Load signature JSON
    with open(signature_json, 'r') as f:
        signature_positions = json.load(f)
        first_coord = signature_positions[0]['pos']
        last_coord = signature_positions[-1]['pos']

    with pysam.AlignmentFile(transcript_file, "r") as samfile:
        for alignment in samfile:
            if alignment.is_unmapped:
                continue

            transcript_id = alignment.query_name
            intervals = sam_to_intervals(alignment)

            if not intervals:
                reason = "no intervals found"
                discarded_transcripts[transcript_id] = reason
                continue

            # Store both full intervals and splice pattern
            full_intervals = intervals_to_tuple(intervals)
            splice_pattern = intervals_to_splice_pattern(intervals)
            interval_counts[splice_pattern] += 1

            # Only keep first occurrence of each unique splice pattern
            if splice_pattern not in tuple_to_id:
                if full_intervals[0][0] > first_coord:
                    reason = "transcript does not cover full signature - late start"
                    discarded_transcripts[transcript_id] = reason
                elif full_intervals[-1][1] < last_coord:
                    reason = "transcript does not cover full signature - early end"
                    discarded_transcripts[transcript_id] = reason
                else:
                    transcript_intervals[transcript_id] = {
                        'full_intervals': full_intervals,
                        'splice_pattern': splice_pattern
                    }
                    tuple_to_id[splice_pattern] = transcript_id
            else:
                reason = f"same as {tuple_to_id[splice_pattern]}"
                discarded_transcripts[transcript_id] = reason

    print(f"Found {len(transcript_intervals)} unique transcript splice patterns")
    print(f"Discarded {len(discarded_transcripts)} non-unique transcripts")

    # Write discarded isoforms to log file if specified
    if log_file and discarded_transcripts:
        try:
            discarded_log = f"{log_file}_discarded.log"
            with open(discarded_log, 'w') as f:
                f.write("Transcript_ID\tReason\n")
                for transcript_id, reason in discarded_transcripts.items():
                    f.write(f"{transcript_id}\t{reason}\n")
            print(f"Collapsed isoforms written to {discarded_log}")
        except Exception as e:
            print(f"Warning: Could not write collapsed isoforms log file: {e}")

    return transcript_intervals

def read_signature_file(signature_file):
    """Read tab-delimited signature file and return dictionary mapping first column to third column."""
    signatures = {}

    try:
        with open(signature_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    print(f"Warning: Line {line_num} has fewer than 3 columns, skipping: {line}")
                    continue

                key = parts[0]
                value = parts[2]
                signatures[key] = value

        print(f"Loaded {len(signatures)} signatures from {signature_file}")
        return signatures

    except FileNotFoundError:
        print(f"Error: Signature file {signature_file} not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading signature file: {e}")
        sys.exit(1)

def process_read_alignments(read_file, transcript_intervals, signatures, file_index, discard_file=None, log_file = None):
    """Process read alignments and match to transcript intervals."""
    print(f"Processing read alignments from {read_file}")

    transcript_to_reads = defaultdict(list)
    matched_reads = 0
    total_reads = 0
    discarded_reads = {}

    with pysam.AlignmentFile(read_file, "r") as samfile:
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
            read_intervals = sam_to_intervals(alignment)

            if not read_intervals:
                discarded_reads[read_id] = "no intervals found"
                continue

            # Use splice pattern for comparison
            read_splice_pattern = intervals_to_splice_pattern(read_intervals)

            # if alignment.query_name == "SRR22838397.16207":
            #     print(f"Read {read_id} has splice pattern {read_splice_pattern}")

            # Check if read splice pattern matches any transcript splice pattern
            matched = False
            for transcript_id, transcript_data in transcript_intervals.items():
                # if alignment.query_name == "SRR22838397.16207":
                #     print(transcript_data['full_intervals'])
                #     print(transcript_data['splice_pattern'])

                if read_splice_pattern == transcript_data['splice_pattern']:
                    transcript_to_reads[transcript_id].append(read_id)
                    matched_reads += 1
                    matched = True
                    break

            if not matched:
                discarded_reads[read_id] = "no matching transcript splice pattern"

    print(f"Matched {matched_reads} reads out of {total_reads} total reads")
    print(f"Discarded {len(discarded_reads)} reads")
    print(f"Found matches for {len(transcript_to_reads)} transcripts")

    # Write discarded reads to file if specified
    if discard_file:
        try:
            discard_filename = f"{discard_file}_file{file_index+1}"
            with open(discard_filename, 'w') as f:
                f.write("Read_ID\tReason\n")
                for read_id, reason in discarded_reads.items():
                    f.write(f"{read_id}\t{reason}\n")
            print(f"Discarded reads written to {discard_filename}")
        except Exception as e:
            print(f"Warning: Could not write discard file: {e}")

    if log_file:
        try:
            used_log = f"{log_file}_used_read_ids_file{file_index+1}.log"
            with open(used_log, 'w') as f:
                f.write("read_id\tassigned\n")
                for transcript in transcript_to_reads:
                    for read in transcript_to_reads[transcript]:
                        f.write(f"{read}\t{transcript}\n")
        except Exception as e:
            print(f"Warning: Error while logging used reads: {e}")

    return transcript_to_reads

def process_reference_alignments(reference_file, transcript_intervals, signatures):
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
            ref_intervals = sam_to_intervals(alignment)

            if not ref_intervals:
                discarded_refs[ref_id] = "no intervals found"
                continue

            # Use splice pattern for comparison
            ref_splice_pattern = intervals_to_splice_pattern(ref_intervals)

            # Check if reference splice pattern matches any transcript splice pattern
            matched = False
            for transcript_id, transcript_data in transcript_intervals.items():
                if ref_splice_pattern == transcript_data['splice_pattern']:
                    transcript_to_refs[transcript_id].append(ref_id)
                    matched_refs += 1
                    matched = True
                    break

            if not matched:
                discarded_refs[ref_id] = "no matching transcript splice pattern"

    print(f"Matched {matched_refs} references out of {total_refs} total references")
    print(f"Discarded {len(discarded_refs)} references")
    print(f"Found matches for {len(transcript_to_refs)} transcripts")

    return transcript_to_refs

def process_multiple_read_files(read_files, transcript_intervals, signatures, discard_file=None, log_file=None):
    """Process multiple read files and return counts per file."""
    all_transcript_to_reads = []

    for i, read_file in enumerate(read_files):
        print(f"\n=== Processing read file {i+1}/{len(read_files)}: {read_file} ===")
        transcript_to_reads = process_read_alignments(
            read_file, transcript_intervals, signatures, i, discard_file, log_file
        )
        all_transcript_to_reads.append(transcript_to_reads)

    return all_transcript_to_reads

def create_signature_transcript_table(all_transcript_to_reads, transcript_to_refs, signatures, read_files):
    """Create a table with signatures as rows and transcripts as columns, with counts per file and reference counts."""

    # Get unique signatures from the signatures dictionary
    unique_signatures = sorted(set(signatures.values()))
    print(f"Found {len(unique_signatures)} unique signatures")

    # Get all transcript IDs from all files and reference
    all_transcript_ids = set()
    for transcript_to_reads in all_transcript_to_reads:
        all_transcript_ids.update(transcript_to_reads.keys())
    all_transcript_ids.update(transcript_to_refs.keys())
    transcript_ids = sorted(all_transcript_ids)
    print(f"Found {len(transcript_ids)} transcripts with matched reads/references across all files")

    # Create signature-transcript count matrix per file
    # Structure: signature_counts[signature][transcript_id] = [count_file1, count_file2, ...]
    signature_counts = defaultdict(lambda: defaultdict(lambda: [0] * len(read_files)))

    for file_idx, transcript_to_reads in enumerate(all_transcript_to_reads):
        for transcript_id, read_ids in transcript_to_reads.items():
            for read_id in read_ids:
                if read_id in signatures:
                    signature = signatures[read_id]
                    signature_counts[signature][transcript_id][file_idx] += 1
                else:
                    print(f"Warning: Read {read_id} not found in signatures file")

    # Create signature-transcript reference count matrix
    # Structure: signature_ref_counts[signature][transcript_id] = count
    signature_ref_counts = defaultdict(lambda: defaultdict(int))

    for transcript_id, ref_ids in transcript_to_refs.items():
        for ref_id in ref_ids:
            if ref_id in signatures:
                signature = signatures[ref_id]
                signature_ref_counts[signature][transcript_id] += 1
            else:
                print(f"Warning: Reference {ref_id} not found in signatures file")

    return unique_signatures, transcript_ids, signature_counts, signature_ref_counts

def format_count_cell(counts, ref_count=0):
    """Format count cell as X(y1,y2,...,yn)+n where X is sum, y1,y2,...,yn are individual counts, and n is reference count."""
    total = sum(counts)
    if total == 0 and ref_count == 0:
        return "0"

    if total == 0:
        return f"0+{ref_count}" if ref_count > 0 else "0"

    individual_counts = ",".join(str(count) for count in counts)
    if ref_count > 0:
        return f"{total}({individual_counts})+{ref_count}"
    else:
        return f"{total}({individual_counts})"

def join_tab_intervals(interval):
    return "\t".join([str(i) for i in interval])

def main():
    parser = argparse.ArgumentParser(description="Match read alignments to transcript interval patterns")
    parser.add_argument("transcript_file", help="SAM/BAM file containing transcript alignments")
    parser.add_argument("read_files", nargs='+', help="SAM/BAM files containing read alignments")
    parser.add_argument("-s", "--signatures", help="Tab-delimited file with sequence signatures")
    parser.add_argument("-j","--signature_json", help="JSON file with signature structure", required=True)
    parser.add_argument("-r", "--reference", help="SAM/BAM file containing reference sequence alignments")
    parser.add_argument("-o", "--output", help="Output file for results (default: stdout)", default=None)
    parser.add_argument("-d", "--discard", help="Base filename for discarded read IDs and reasons (will append _file1, _file2, etc.)", default=None)
    parser.add_argument("-l", "--log", help="Base filename for isoform log files (will create _collapsed.log and _unused.log)", default=None)

    args = parser.parse_args()

    print(f"Processing {len(args.read_files)} read files:")
    for i, read_file in enumerate(args.read_files, 1):
        print(f"  {i}. {read_file}")

    # Read signature file
    signatures = read_signature_file(args.signatures)

    try:
        # Process transcript alignments to get unique interval sets
        transcript_intervals = process_transcript_alignments(args.transcript_file, args.signature_json, args.log)

        if not transcript_intervals:
            print("No unique transcript intervals found!")
            sys.exit(1)

        # Process multiple read files
        all_transcript_to_reads = process_multiple_read_files(
            args.read_files, transcript_intervals, signatures, args.discard, args.log
        )

        # Process reference alignments if provided
        transcript_to_refs = defaultdict(list)
        if args.reference:
            transcript_to_refs = process_reference_alignments(
                args.reference, transcript_intervals, signatures
            )

        # Create signature-transcript count table
        unique_signatures, transcript_ids, signature_counts, signature_ref_counts = create_signature_transcript_table(
            all_transcript_to_reads, transcript_to_refs, signatures, args.read_files
        )

        # Calculate column totals first to sort transcript_ids by descending totals
        column_totals_dict = {}
        for transcript_id in transcript_ids:
            column_total = 0
            for signature in unique_signatures:
                counts = signature_counts[signature][transcript_id]
                x_value = sum(counts)
                column_total += x_value
            column_totals_dict[transcript_id] = column_total

        # Sort transcript_ids by column totals in descending order
        sorted_transcript_ids = sorted(transcript_ids, key=lambda tid: column_totals_dict[tid], reverse=True)

        # Create sorted column_totals array
        sorted_column_totals = [column_totals_dict[tid] for tid in sorted_transcript_ids]

        # Output results
        output_file = open(args.output, 'w') if args.output else sys.stdout

        try:
            # Write header with sorted transcript IDs
            output_file.write("Signature\t" + "\t".join(sorted_transcript_ids) + "\tRow_Total\n")

            # Write data rows
            for signature in unique_signatures:
                row = [signature]
                row_total = 0
                for transcript_id in sorted_transcript_ids:
                    counts = signature_counts[signature][transcript_id]
                    ref_count = signature_ref_counts[signature][transcript_id]
                    formatted_count = format_count_cell(counts, ref_count)
                    row.append(formatted_count)
                    # Add the X value (sum of counts) to row total
                    x_value = sum(counts)
                    row_total += x_value
                row.append(str(row_total))
                output_file.write("\t".join(row) + "\n")

            # Write column totals row
            column_total_row = ["Column_Total"] + [str(total) for total in sorted_column_totals]
            column_total_row.append(str(sum(sorted_column_totals)))  # Grand total
            output_file.write("\t".join(column_total_row) + "\n")

        finally:
            if args.output:
                output_file.close()

        # Log unused isoforms (transcripts with no matching reads)
        if args.log:
            unused_transcripts = []
            all_used_transcript_ids = set()

            # Collect all transcript IDs that had matching reads
            for transcript_to_reads in all_transcript_to_reads:
                all_used_transcript_ids.update(transcript_to_reads.keys())

            # Find transcripts that were processed but had no matching reads
            for transcript_id in transcript_intervals.keys():
                if transcript_id not in all_used_transcript_ids:
                    unused_transcripts.append(transcript_id)

            if all_used_transcript_ids != set():
                try:
                    used_log = f"{args.log}_used.log"
                    with open(used_log, 'w') as f:
                        f.write("Transcript_ID\tReason\tFull_Intervals\n")
                        for transcript_id in all_used_transcript_ids:
                            interval_data = transcript_intervals[transcript_id]
                            full_intervals_str = "\t".join([join_tab_intervals(interval) for interval in interval_data['full_intervals']])
                            # splice_pattern_str = str(interval_data['splice_pattern'])
                            f.write(f"{transcript_id}\tmatching reads found\t{full_intervals_str}\n")
                    print(f"Used isoforms written to {used_log}")
                    print(f"Found {len(all_used_transcript_ids)} used isoforms")
                except Exception as e:
                    print(f"Warning: Could not write used isoforms log file: {e}")

            if unused_transcripts:
                try:
                    unused_log = f"{args.log}_unused.log"
                    with open(unused_log, 'w') as f:
                        f.write("Transcript_ID\tReason\tFull_Intervals\n")
                        for transcript_id in unused_transcripts:
                            interval_data = transcript_intervals[transcript_id]
                            full_intervals_str = "\t".join([join_tab_intervals(interval) for interval in interval_data['full_intervals']])
                            # full_intervals_str = str(interval_data['full_intervals'])
                            # splice_pattern_str = str(interval_data['splice_pattern'])
                            f.write(f"{transcript_id}\tno matching reads found\t{full_intervals_str}\n")
                    print(f"Unused isoforms written to {unused_log}")
                    print(f"Found {len(unused_transcripts)} unused isoforms")
                except Exception as e:
                    print(f"Warning: Could not write unused isoforms log file: {e}")
            else:
                print("All processed isoforms had matching reads")

        print(f"\nSignature-transcript count table written to {args.output or 'stdout'}")
        print(f"Table dimensions: {len(unique_signatures)} signatures × {len(transcript_ids)} transcripts")
        print(f"Count format: X(y1,y2,...,yn) where X is total and y1,y2,...,yn are counts per file")

    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
