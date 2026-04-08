Y_PATH="/Users/kxp5629/proj/Y"

mkdir -p $Y_PATH/data/blast

cd $Y_PATH/data/blast

# Download all human protein files
echo "Downloading human protein sequences..."
for i in {1..12}; do  # Currently there are 12 files, but this might change over time
    wget "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.${i}.protein.faa.gz"
done


# Decompress all files
echo "Decompressing files..."
gunzip human.*.protein.faa.gz

# Concatenate all files into one
echo "Concatenating files..."
cat human.*.protein.faa > human_refseq_protein.faa

rm human.[0-9].protein.faa
rm human.1[0-9].protein.faa


# Make BLAST database from concatenated file
echo "Creating BLAST database..."
makeblastdb -in human_refseq_protein.faa \
    -dbtype prot \
    -parse_seqids \
    -out human_protein

# Optional: remove the concatenated FASTA file to save space
# rm human_refseq_protein.faa

echo "BLAST database created. Use with: -db $(pwd)/human_protein"