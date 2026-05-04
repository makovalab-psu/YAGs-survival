import json
import sys
import csv
from collections import defaultdict


def main(filepath: str) -> None:
    with open(filepath) as f:
        data = json.load(f)

    # {(species, gene): {isoform_name: num_ord_blocks}}
    counts: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)

    for entry in data:
        species = entry["species"]
        gene = entry["gene"]
        for iso in entry["isoforms"]:
            counts[(species, gene)][iso["isoform"]] = len(iso["counts"])

    species_list = sorted({s for s, _ in counts})
    gene_list = sorted({g for _, g in counts})

    writer = csv.writer(sys.stdout)
    writer.writerow(["gene"] + species_list)

    for gene in gene_list:
        row = [gene]
        for species in species_list:
            n = 0
            for key, value in counts.get((species, gene), {}).items():
                n += value
            row.append(n if n > 0 else "")
        writer.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <results.json>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
