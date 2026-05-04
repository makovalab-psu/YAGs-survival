#!/usr/bin/env python3

import argparse
import json
import sys
import subprocess
import matplotlib.pyplot as plt
from helpers.consensus_sequences_spoa import find_orfs_biopython

from Bio import SeqIO
from Bio.Align import PairwiseAligner

from Bio.Seq import Seq
from collections import defaultdict, Counter

known = [
    {"name":"red"   , "seq": "AYSAYPHSPGQVITGCQLLVYNYQ", "nt":""},
    {"name":"orange", "seq": "EYPTYPDSAFQVTTGYQLPVYNYQ", "nt":""},
    {"name":"blue"  , "seq": "PFPAYPRSPFQVTAGYQLPVYNYQ", "nt":""},
    {"name":"yellow", "seq": "AFPAYPNSPFQVATGYQFPVYNYQ", "nt":""},
    {"name":"pink"  , "seq": "PFPAYPSSPFQVTAGYQLPVYNYQ", "nt":""},
    {"name":"green" , "seq": "AFPAYPNSPVQVTTGYQLPVYNYQ", "nt":""},
    {"name":"gray"  , "seq": "AFPAYPSSPFQVTTGYQLPVYNYQ", "nt":""},
    {"name":"brown" , "seq": "AFPAYPNSAVQVTTGYQFHVYNYQ", "nt":""}
]

# known = [
#     {"name":"yellow", "seq": "AFPAYPNSPFQVATGYQFPVYNYQ", "nt":""},
#     {"name":"gray"  , "seq": "AFPAYPSSPFQVTTGYQLPVYNYQ", "nt":""},
#     {"name":"green" , "seq": "AFPAYPNSPVQVTTGYQLPVYNYQ", "nt":""},
#     {"name":"blue"  , "seq": "PFPAYPRSPFQVTAGYQLPVYNYQ", "nt":""},
#     {"name":"pink"  , "seq": "PFPAYPSSPFQVTAGYQLPVYNYQ", "nt":""},
#     {"name":"brown" , "seq": "AFPAYPNSAVQVTTGYQFHVYNYQ", "nt":""},
#     {"name":"orange", "seq": "EYPTYPDSAFQVTTGYQLPVYNYQ", "nt":""},
#     {"name":"red"   , "seq": "AYSAYPHSPGQVITGCQLLVYNYQ", "nt":""},

# ]

# the cutoff for labeling of DAZ repeat sequence based on frequency of occurence
species_defined_cutoffs = {
    "HomSap":10,
    "PanPan":20,
    "PonAbe":3,
    "PonPyg":10,
    "PanTro":20,
    "GorGor":10
}

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input table produced by hmmsearch", required=True)
    parser.add_argument("-f", "--fasta", help="preprocessed fasta file with DAZ tails", required=True)
    parser.add_argument("-s", "--species", help="Which species are we analyzing", required=True)
    parser.add_argument("-o", "--output", help="output file prefix", required=True)
    parser.add_argument("sample_files", help="Input FASTQ file", nargs="*")

    args = parser.parse_args()

    fasta_file = args.fasta

    sequences_by_id = SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta"))

    with open(args.input, 'r') as table_file:

        read_name_index = 0
        start_index = 6
        end_index = 7

        hits = defaultdict(list)

        # read in the table of hmmsearch hits
        for line in table_file:
            if line.startswith("#"):
                continue
            line = line.strip()
            fields = line.split()
            align_start = fields[start_index]
            align_end = fields[end_index]
            read_name = fields[read_name_index]

            hits[read_name].append({"from": int(align_start), "to": int(align_end), "length" : int(align_end) - int(align_start) + 1 })


        sizes = []
        list_of_all_seqs = []

        # post process hits, the motif recognition struggles slightly in the first few nucleotides
        # in the human all DAZ repeats are 72BP long, adding nucleotides if first exon is shorter
        # in consecutive exomes if there is a gap between motifs adding the difference in the the start
        # of the following exon
        # Edit: Orangutans can have 90bp DAZ repeats.
        for key,value in hits.items():
            # print(value)
            # print(len(value))
            # assert False
            if len(value) >= 1:
                # print(key)
                sorted_positions = sorted(value, key=lambda d: d['from'])

                for index in range(len(sorted_positions)):
                    if index == 0:
                        current_pos = sorted_positions[index]
                        if current_pos['length'] < 72:
                            current_pos['from'] = current_pos['to'] -72 +1
                        elif current_pos['length'] > 80:
                            if current_pos['length'] < 90:
                                current_pos['from'] =current_pos['to'] -90 +1

                        if current_pos['from'] < 1:
                            current_pos['from'] = 1;
                            current_pos['length'] = current_pos['to'] - current_pos['from']


                    if index > 0:
                        current_pos = sorted_positions[index]
                        previous_pos = sorted_positions[index-1]

                        if(current_pos['from'] -  previous_pos['to'] > 1):
                            current_pos['from']  = previous_pos['to'] + 1
                            current_pos['length'] = current_pos['to'] - current_pos['from'] + 1

                    if key in sequences_by_id:
                        selected_record = sequences_by_id[key]
                        print_sequence = selected_record.seq[current_pos['from']-1:current_pos['to']]
                        current_pos['sequence']=str(print_sequence)
                        list_of_all_seqs.append(str(print_sequence))
                        if str(print_sequence) == "":
                            print(f"Found empty sequence: {key}: {value}")
                            assert False

            sizes.append(len(value))

        counts_seqs = Counter(list_of_all_seqs)


        counts = Counter(sizes)
        values = list(counts.keys())
        frequencies = list(counts.values())

        plt.bar(values, frequencies, edgecolor = "black")
        plt.xticks(list(range(1,21)))
        plt.savefig(f"DAZ-{args.species}-repeat_counts-hist.png")


        labeled = [] #list of labeled sequences
        label = 'A'
        counts_seqs_sorted = dict(sorted(counts_seqs.items(), key= lambda item: item[1] , reverse=True ))

        for key in counts_seqs_sorted:
            # print(counts_seqs_sorted[key])
            if counts_seqs_sorted[key] > species_defined_cutoffs[args.species]:
                dna_sequence = Seq(key)
                aa_sequence = str(dna_sequence.translate())
                if(aa_sequence in [x['seq'] for x in labeled]):
                    continue
                labeled.append({"seq":aa_sequence, "name":label, "count":counts_seqs_sorted[key], 'nt':""})
                label = chr(ord(label) + 1)



        labeled_sequences = [x['seq'] for x in labeled]
        known_sequences = [x['seq'] for x in known]

        sequence_dict = {}

        aligner = PairwiseAligner()
        aligner.match_score = 2
        aligner.mismatch_score  = -1
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        aligner.mode = "global"

        with open(f"{args.output}.seq_frequencies.tsv", "w") as out_f:
            # print found sequences and their labels - or closest possible label
            for item in counts_seqs_sorted.items():

                daz_sequence = item[0]
                if daz_sequence == "":
                    print("No DAZ sequence found")
                    print(item)
                    assert False
                sequence_dict[daz_sequence] = ""
                count = item[1]
                s = Seq(daz_sequence)
                tr = s.translate()

                color = "not_assigned"

                if tr in labeled_sequences:
                    ref = [x for x in labeled if x['seq'] == tr][0]
                    if (ref["nt"] == ""):
                        ref["nt"] = daz_sequence

                    color  = [x["name"] for x in labeled if x['seq'] == tr][0]

                best_distance = 100000000
                best_color = ""
                ambiguous = False

                if count <=  species_defined_cutoffs[args.species]:
                    distances = []
                    for ref in labeled:
                        seq1 = ref["nt"]
                        try:
                            alignments = aligner.align(seq1, daz_sequence)
                            best_alignment = alignments[0]
                            max_possible_score = len(max(seq1, daz_sequence, key=len)) * 2
                            distance = max_possible_score  - best_alignment.score
                            distances.append({"ref":seq1, "color": ref["name"], "distance": distance})
                        except Exception as e:
                            print(f"Error when processin alignment of {seq1} vs {daz_sequence}")
                            print(f"Unexpected error: {e}")

                    distances = sorted(distances, key=lambda d: d["distance"])
                    best_color = distances[0]["color"]
                    if distances[0]["distance"] ==  distances[1]["distance"]:
                        ambiguous = True

                    if((not ambiguous) and (best_color !="")):
                        sequence_dict[daz_sequence] = best_color

                else:
                    sequence_dict[daz_sequence] = color
                print(f"{daz_sequence}\t{count}\t{tr}\t{color}\t{best_color}\t{best_distance}\t{ambiguous}", file = out_f)

        color_patterns = defaultdict(list)
        sorted_hits = dict(sorted(hits.items(), key=lambda x : len(x[1]), reverse=True))

        uniq_combo_dict = defaultdict(int)
        truly_uniq_combos = []
        for item in sorted_hits:
            # print(item)
            # print(hits[item])
            # assert False
            exons = sorted(hits[item], key=lambda d: d['from'])
            pattern = []
            if (len(exons) >=1):
                skip = False
                for exon in exons:
                    print(exon)
                    # print(f"{exon['sequence']}\t{sequence_dict[exon['sequence']]}")
                    sequence = exon['sequence']
                    if (not sequence in sequence_dict):
                        skip = True
                        continue
                    color = sequence_dict[sequence]
                    if color == '':
                        skip = True
                    pattern.append(color)

                if(skip):
                    continue

                pattern_tag = "-".join(pattern)

                uniq_combo_dict[pattern_tag] += 1
                color_patterns[pattern_tag].append(item)

                is_uniq = True
                for combo in truly_uniq_combos:
                    if( len(pattern)> len(combo)):
                        continue
                    i = j = 0
                    while i < len(pattern) and j < len(combo):
                        if pattern[i] == combo[j]:
                            i+=1
                        j+=1

                    if i == len(pattern):
                        is_uniq = False

                if is_uniq:
                    truly_uniq_combos.append(pattern)
                    # print(f"uniq combo len: {len(truly_uniq_combos)}")


        with open(f"{args.output}.tail_to_sequenceID_dict.json","w") as out_json:
            json.dump(color_patterns, out_json, indent = 4)

        # print(f"Found {len(truly_uniq_combos)} uniq combos.")
        # for combo in truly_uniq_combos:
        #     print(",".join(combo))

        truly_string = ["-".join(x) for x in truly_uniq_combos]

        # print(uniq_combo_dict)
        with open(f"{args.output}.uniq_tails.tsv", 'w') as out_tsv:
            for key in uniq_combo_dict:
                u = "F"
                if key in truly_string:
                    u = "T"

                print(f"{key}\t{uniq_combo_dict[key]}\t{u}", file=out_tsv)

        sorted_uniq_combo_dict = dict(sorted(uniq_combo_dict.items(), key=lambda item: item[1], reverse=True))

        counts = [value for key,value in uniq_combo_dict.items()]
        sum_counts = sum(counts)
        limit = sum_counts /100 *2
        plt.figure(figsize=(10,6))
        plt.plot(list(sorted_uniq_combo_dict.keys()), list(sorted_uniq_combo_dict.values()), linewidth=2)
        plt.axhline(y=5, color='red', linestyle='--', linewidth=1)
        plt.axhline(y=limit, color='orange', linestyle='-.', linewidth=1)
        plt.xlabel('tail tags')
        plt.ylabel('occurence')
        plt.title(f"{args.species} DAZ tail occurence")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(f"{args.output}.hist.png")


        with open(f"{args.output}.hits.json", 'w') as out_json:
            json.dump(hits, out_json, indent=4)

        #finally create a fastq file for each group passing the limit threshold
        #first read in all the DAZ fastq files
        all_records_dict = {}
        for file in args.sample_files:
            records_from_file =SeqIO.to_dict(SeqIO.parse(file, "fastq"))
            all_records_dict.update(records_from_file)


        with open(f"{args.output}.aa.fa", "w") as outfile:
            for key, value in uniq_combo_dict.items():
                if value > limit and value > 5:
                    fq_file = f"{args.output}.{key}.fq"
                    with open(fq_file,"w") as out_fq:

                        for seqID in color_patterns[key]:
                            record = all_records_dict[seqID]
                            # print(f">{record.id}", file=out_fq)
                            # print(f"{record.seq}", file=out_fq)
                            # print("+", file=out_fq)
                            # print(f"{record.letter_annotations['phred_quality']}", file=out_fq)
                            SeqIO.write(record, out_fq, "fastq")


                    #call SPOA
                    result = subprocess.run(['spoa', fq_file],
                                        capture_output=True, text=True, timeout=300)
                    consensus = result.stdout.split("\n")[1] # remove fasta header >Consensus...
                    print(f"{key} consensus:")
                    print(consensus)

                    #translate consensus
                    orfs = find_orfs_biopython(consensus, 80)

                    print(f"found {len(orfs)} orfs")

                    if orfs:
                        # Sort by length (longest first)
                        orfs.sort(key=lambda x: x['length'], reverse=True)
                        longest_orf = orfs[0]

                        print(longest_orf['protein'])

                        # print(f"  Found {len(orfs)} ORFs (≥{args.min_orf_length} aa)", file = outfile)
                        # print(f"  Longest ORF:", file = outfile)
                        # print(f"    Frame: {longest_orf['frame']}", file = outfile)
                        # print(f"    Position: {longest_orf['start']}-{longest_orf['end']}", file = outfile)
                        # print(f"    Length: {longest_orf['length']} amino acids", file = outfile)
                        # print(f"    Protein: {longest_orf['protein'][:50]}{'...' if len(longest_orf['protein']) > 50 else ''}", file = outfile)

                        print(f">{args.species}::{key}::{value}_reads", file = outfile)
                        print(f"{longest_orf['protein']}", file = outfile)

                    else:
                        print(f"  No ORFs found (≥80 amino acids)")



if __name__ == "__main__":
    main()
