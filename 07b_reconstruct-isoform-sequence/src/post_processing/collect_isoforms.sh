genes=("BPY2" "CDY" "HSFY" "RBMY" "RBMYB" "TSPY" "VCY")

mkdir -p ../data/collect

# for gene in ${genes[@]};do
#  cat ../data/results/*/${gene}/isoforms/*_ORFs.fa > ./../data/collect/collect_all_iso_${gene}.fa.orig
# done

## For DAZ

cat ../data/results/*/DAZ/*aa.fa > ./../data/collect/collect_all_iso_${gene}.fa.orig
# mkdir -p ../data/collect/dist_plot

# for gene in ${genes[@]};do
#   cp ../data/results/*/${gene}/isoforms/*iso.png ./../data/collect/dist_plot/
# done

# mkdir -p ../data/collect/dist_plot/DAZ/
# cp ../data/results/*/DAZ/*.hist.png ./../data/collect/dist_plot/DAZ/
