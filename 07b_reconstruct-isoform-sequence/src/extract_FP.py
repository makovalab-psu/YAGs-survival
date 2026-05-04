#!/usr/bin/env python3

import pysam
import vcf
import json
import argparse
import sys
import traceback

def count_unique_strings(list_of_sequences):
  unique_strings = set()
  for sequence in list_of_sequences:
    unique_strings.add(sequence)
  return len(unique_strings)

def find_char_in_string(string,char):
    return [index for index,letter in enumerate(string) if letter == char]

def remove_positions_from_str_list(list_of_sequences, remove_positions):
    """
    Positions specified in 'remove positions' will be removed from all strings
    in the list_of_strings. Individual positions are removed from greatest to
    smallest.
    """

    return_list = []
    for sequence in list_of_sequences:
        new_string = "".join([x for i,x in enumerate(sequence) if i not in remove_positions])
        return_list.append(new_string)

    return return_list

def extract_differences(bam_file, bed_file, optimize = False, species = "", gene = ""):
    """
    Extract positions where reads differ from the reference in a BAM file.

    Args:
        bam_file (str): Path to the input BAM file.
        output_json (str): Path to the output JSON file.
    """
    intervals = []
    with open(bed_file, "r") as bed:
        for line in bed.readlines():
            line=line.strip()
            fields = line.split("\t")
            intervals.append((int(fields[1]), int(fields[2])))

    differences = []
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:  # Skip unmapped reads
                continue

            # print(read.get_aligned_pairs(matches_only=True, with_seq=True))
            for query_pos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                if query_pos is None or ref_base is None:  # Skip gaps
                    continue

                read_base = read.query_sequence[query_pos]
                if read_base != ref_base:  # Check if the read base differs from the reference base
                    pos = ref_pos + 1  # 1-based position
                    if not any(d['pos'] == pos for d in differences):
                        for start, stop in intervals:
                            if start <= pos <= stop:

                                differences.append({
                                    "pos": pos,
                                    "source": "ref_assembly",
                                    "type": "SNV",
                                    "found": ".",
                                    "alt1": ref_base.upper(),
                                    "alt2": read_base

                                })


    if species == "HomSap" and gene == "TSPY":
        differences.append({
            "pos": 816,
            "source": "manual",
            "type": "INS",
            "found": "-",
            "alt1": "T",
            "alt2": "TGGCGGAGGTGGAGGTGG"
        })

    if species == "HomSap" and gene == "HSFY":
        differences.append({
            "pos": 2000,
            "source": "manual",
            "type": "SNV",
            "found": ".",
            "alt1": "T",
            "alt2": "A"
        })

    if species == "HomSap" and gene == "BPY2":
        differences.append({
            "pos": 3700,
            "source": "manual",
            "type": "SNV",
            "found": ".",
            "alt1": "A",
            "alt2": "T"
        })

    if species == "PanPan" and gene == "VCY":
        differences.append({
            "pos": 1303,
            "source": "manual",
            "type": "SNV",
            "found": ".",
            "alt1": "C",
            "alt2": "G"
        })

    differences = sorted(differences, key=lambda x: x['pos'])

    print(differences)
    print(f"{species} {gene}")

    remove_positions_coord = []
    if optimize:
        positions_list = [x['pos'] for x in differences]

        #recunstruct fingerprints for all reads in bam file
        all_sigs = []
        with pysam.AlignmentFile(bam_file, "rb") as bam:

            for read in bam:
                sig = '.'*len(differences)
                if read.is_unmapped:  # Skip unmapped reads
                    continue

                for query_pos, ref_pos, ref_base, cigar in read.get_aligned_pairs(with_seq=True, with_cigar=True):
                    if query_pos is None or ref_base is None:
                        continue

                    read_base = read.query_sequence[query_pos]
                    if ref_pos + 1 in positions_list:
                        ind = positions_list.index(ref_pos+1)
                        sig = sig[:ind] + read_base.upper() + sig[ind+1:]
                all_sigs.append(sig)

        all_sigs = list(set(all_sigs))
        # try to remove dotted positions first
    
        dot_position_lists = [find_char_in_string(x,'.') for x in all_sigs] # find the positions of dots    
        dot_positions = [inner for outer in dot_position_lists for inner in outer] #flatten list of positions per sequence
        dot_positions = list(set(dot_positions)) # make unique
        dot_positions.sort()
        
        remove_positions = []
  
        # this has to be iterative, otherwise I'm losing resolutions
        for i in dot_positions:
            sub_sigs = remove_positions_from_str_list(all_sigs, remove_positions + [i])
            if (len(all_sigs) == count_unique_strings(sub_sigs)):
                remove_positions.append(i)
 

        for i in range(len(all_sigs[0])):
            if i in dot_positions:
                continue
            # I don't want to remove positions that were called from the samples because I'm not testing signatures from all the samples
            if (differences[i]['source'] != "ref_assembly"):
                continue
            sub_sigs = remove_positions_from_str_list(all_sigs, remove_positions + [i])
            if (len(all_sigs) == count_unique_strings(sub_sigs)):
                remove_positions.append(i)


        sub_sigs = remove_positions_from_str_list(all_sigs, remove_positions)
        remove_positions_coord = [differences[i]['pos'] for i in remove_positions ]
        differences = [differences[i] for i in range(len(differences)) if i not in remove_positions]

    # print(remove_positions)

    return [intervals, differences, remove_positions_coord]

def filter_vcf_file(vcf_file, intervals, differences, output_json, remove_positions_coord):
    """
    Load VCF file and filter variants one by one
    """

    vcf_reader = vcf.Reader(open(vcf_file), 'r')

    for variant in vcf_reader:
        for start, stop in intervals:
            if start <= variant.POS <= stop:
                if variant.is_snp:
                    # print(variant.samples)
                    sample = variant.samples[0]
                    if hasattr(sample.data, 'DP') and sample.data.DP is not None:
                        if sample.data.DP > 30:
                            # print(f"{variant.POS}, {variant.REF, variant.ALT[0]}")

                            if not any(d['pos'] == variant.POS for d in differences):
                                # if variant.POS in remove_positions_coord:
                                #     print(f"skipping variant {variant.POS}")
                                #     continue
                                # print(f"{sample.data.AD}")
                                # if sum(sample.data.AD)/sample.data.AD[1] <10: 
                                #     pass
                                differences.append({
                                    "pos": variant.POS,
                                    "source": "var_call",
                                    "type": "SNV",
                                    "found": ".",
                                    "alt1": variant.REF,
                                    "alt2": str(variant.ALT[0])
                                })
    # print(differences)
    # sort differences
    differences = sorted(differences, key=lambda x: x['pos'])
    # Write the differences to a JSON file
    with open(output_json, "w") as json_file:
        json.dump(differences, json_file, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Extract positions where reads differ from the reference in a BAM file.")
    parser.add_argument("--optimize", action='store_true')
    parser.add_argument("--species", help="Specify which spicies we're analyzing.", default = "")
    parser.add_argument("--gene", help="Which gene.", default = "")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("bed_file", help="Path to the input BED file.")
    parser.add_argument("vcf_file", help="Path to the input VCF file")
    parser.add_argument("output_json", help="Path to the output JSON file.")
    args = parser.parse_args()

    try:
        [intervals, differences, remove_positions_coord] = extract_differences(args.bam_file, args.bed_file, args.optimize, args.species, args.gene)
        filter_vcf_file(args.vcf_file, intervals, differences, args.output_json, remove_positions_coord)
        print(f"Differences extracted and saved to {args.output_json}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
