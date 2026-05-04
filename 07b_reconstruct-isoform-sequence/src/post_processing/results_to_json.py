#!/usr/bin/env python3
import glob
import sys
from pathlib import Path
from collections import defaultdict
import json
import argparse

# Add parent directory to Python path so I can import from helpers
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

from helpers.consensus_sequences_spoa import parse_cell_format

result = []


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", help="name of output file", required=True)

    args = parser.parse_args()

    files = glob.glob("data/results/*/*/isoforms/*isoforms_processed.txt")
    global_file = []
    for file_name in files:
        species, gene, *rest = file_name.split("/")[-1].split("_")
        entry = {}
        entry["species"] = species
        entry["gene"] = gene
        isoforms_list = []
        individual_counts = defaultdict(list)
        with open(file_name, "r") as input_file:
            isoforms = []
            isoform_counts = []
            counter = 1
            total = 0
            for line in input_file:
                used = False
                line = line.strip()
                if line.startswith("Signature"):
                    _, *isoforms, total = line.split("\t")
                    # print(signatures)
                    continue

                elif line.startswith("Column_Total"):
                    total = line.split("\t")[-1]
                    total = int(total)
                    if total == 0:
                        break
                    limit = total / 100 * 2
                    if limit < 5:
                        limit = 5
                    # print(f"Species: {species}, gene: {gene}, total: {total}, limit: {limit}")
                    # print(line)
                else:
                    _, *counts, row_total = line.split("\t")
                    # print(line)

                    isoform_counts = [parse_cell_format(x) for x in counts]
                    # print(isoform_countsa
                    # print(len(isoforms))
                    for i in range(len(isoforms)):
                        # print(isoform_counts[i])
                        if (
                            isoform_counts[i]["n"] > limit
                            and isoform_counts[i]["at_least_two_positive"]
                        ):
                            # print(f"Isoform: {isoforms[i]}, position: {counter}, count: {isoform_counts[i]['n']}")
                            individual_counts[isoforms[i]].append(
                                {"ord": counter, "count": isoform_counts[i]["n"]}
                            )
                            used = True
                    if used:
                        counter += 1

            for isoform in isoforms:
                entry = {}
                entry["isoform"] = isoform

                if individual_counts[isoform] != []:
                    entry["counts"] = individual_counts[isoform]
                    isoforms_list.append(entry)

            global_entry = {}
            global_entry["species"] = species
            global_entry["gene"] = gene
            global_entry["total"] = total
            global_entry["isoforms"] = isoforms_list
            global_file.append(global_entry)

    with open(args.output, "w") as outfile:
        json.dump(global_file, outfile, indent=4)


if __name__ == "__main__":
    main()
