#!/usr/env/bin bash
targeted=("SRR22838391" "SRR22838392" "SRR22838393" "SRR22838394" "SRR22838395" "SRR22838396" "SRR22838397" "SRR22838398" "SRR22838399" "SRR22838400" "SRR22838401" "SRR22838402" "SRR22838403" "SRR22838404" "SRR22838405" "SRR22838406")
untargeted=("SRR22452368" "SRR22452369" "SRR22452370" "SRR22306523" "SRR22306524" "SRR23254336")

output_dir="./reads/raw"
mkdir -p "$output_dir"

process_sra(){
    local sra_id=$1
    echo "Processing $sra_id..."
    ~/tools/sratoolkit.3.1.1-ubuntu64/bin/prefetch "$sra_id" --output-directory "$output_dir"
    ~/tools/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump "$output_dir/$sra_id" -O "$output_dir" --split-files

    rm -rf "$output_dir/$sra_id"

    echo "Completed $sra_id"
}

# for sra_id in "${targeted[@]}"; do
#   process_sra "$sra_id"
# done

for sra_id in "${untargeted[@]}"; do
  process_sra "$sra_id"
done

echo "All downloads and conversions are complete."