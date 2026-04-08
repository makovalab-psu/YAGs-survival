import os
import subprocess
import gzip

# reference sequences are softmasked!
REFS = {
    "HomSap": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/HomSap/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    "GorGor": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/GorGor/data/GCF_029281585.2/GCF_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna",
    "PanPan": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PanPan/data/GCF_029289425.2/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna",
    "PanTro": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PanTro/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna",
    "PonAbe": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PonAbe/data/GCF_028885655.2/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna",
    "PonPyg": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PonPyg/data/GCF_028885625.2/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna",
    "SymSyn": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/SymSyn/data/GCF_028878055.3/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna"
}

REPEATS = {
   "HomSap": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/CHM13.GCF_009914755.1_T2T-CHM13v2.0_rm.bed.gz",
   "GorGor": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/mGorGor1.pri.cur.20231122.combo3.bed.gz",
   "PanPan": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/mPanPan1.pri.cur.20231122.combo3.bed.gz",
   "PanTro": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/mPanTro3.pri.cur.20231122.combo.bed.gz",
   "PonAbe": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/mPonAbe1.pri.cur.20231205.combo.bed.gz",
   "PonPyg": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/mPonPyg2.pri.cur.20231122.combo.bed.gz",
   "SymSyn": "/storage/home/kxp5629/group_storage/shared/T2Tv2.repeats/mSymSyn1.pri.cur.20231205.combo.bed.gz"
}

x_chr = {
    "HomSap": "NC_060947.1",
    "GorGor": "NC_073247.2",
    "PanPan": "NC_073272.2",
    "PanTro": "NC_072421.2",
    "PonAbe": "NC_072008.2",
    "PonPyg": "NC_072396.2",
    "SymSyn": "NC_072447.2"
}

y_chr = {
    "HomSap": "NC_060948.1",
    "GorGor": "NC_073248.2",
    "PanPan": "NC_073273.2",
    "PanTro": "NC_072422.2",
    "PonAbe": "NC_072009.2",
    "PonPyg": "NC_072397.2",
    "SymSyn": "NC_072448.2"
}

SPECIES = ["HomSap", "GorGor", "PanPan", "PanTro", "PonAbe", "PonPyg", "SymSyn"]
#SPECIES = ["GorGor"]

def run_pipeline(chromosome,chr_dict):

    for species in SPECIES:
        
        # 1. use samtools to subset only chrY
        print(REFS[species])

        if not os.path.isdir("./tmp"):
            os.mkdir("./tmp")

        chr_file = f"./tmp/{species}_{chromosome}.fa"
        
        # samtools faidx my_genome.fa chromosome_1 > chromosome_1.fa
        command = ["samtools", "faidx", f"{REFS[species]}", f"{chr_dict[species]}"]

        if species == "HomSap" and chromosome == "chr_y": 
            # in human the heterochromatin region is taking forever to selfalign subsetting the fasta to 27449931 
            # where the heterochromatin region should start according to Supplemntary table 21 in Rhie et a. 2023
            command = ["samtools", "faidx", f"{REFS[species]}", f"{chr_dict[species]}:1-27449931"] 
        with open(chr_file, 'w') as outfile:
            subprocess.run(command, stdout = outfile, text=True, check=True)

        # 2. subset repeats bed file and rename first column
        repeats_file = f"./tmp/{species}_{chromosome}_repeats.bed"
        with open(repeats_file, "w") as outfile:
            with gzip.open(f"{REPEATS[species]}", 'rt') as infile:
                for line in infile:
                    line = line.strip()
                    fields = line.split("\t")
                    if (chromosome in fields[0]) or (chr_dict[species] in fields[0]):
                        if fields[3] == "Low_complexity" or fields[3] == "Simple_repeat" or "Satellite" in fields[3]:
                            print(f"{chr_dict[species]}\t{fields[1]}\t{fields[2]}\t{fields[3]}", file=outfile)


        # collecting the full repeat file in case I need it later
        if not os.path.isdir("./repeats"):
            os.mkdir("./repeats")

        repeats_file_full = f"./repeats/{species}_{chromosome}_repeats_full.bed"
        with open(repeats_file_full, "w") as outfile:
            with gzip.open(f"{REPEATS[species]}", 'rt') as infile:
                for line in infile:
                    line = line.strip()
                    fields = line.split("\t")
                    if "chrY" in fields[0] or chr_dict[species] in fields[0]:
                        print(f"{chr_dict[species]}\t{fields[1]}\t{fields[2]}\t{fields[3]}", file=outfile)


        # 3. create sequence file for lastz alignment 

        sequence_file = f"./tmp/{species}_{chromosome}_squence_file"
        with open(sequence_file, 'w') as outfile:
            print(chr_dict[species], file = outfile)

        # 4. create lastz command
        lastz_file = f"./tmp/{species}_{chromosome}.lastz"
        command = "lastz --scores=lastz_scoring \\\n"
        command += "    --allocate:traceback=2024M \\\n"
        command += "    --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \\\n"
        command +=f"    {chr_file}[subset={sequence_file}] \\\n"
        command +=f"    {chr_file}[subset={sequence_file}] \\\n"
        command +=f"    --strand=minus > {lastz_file}"

        lastz_bash = f"./tmp/{species}_{chromosome}_run_lastz.sh"
        with open(lastz_bash, 'w') as outfile:
            print(command, file = outfile)
        
        # these alignments run up to an hour each, could paralelized if want to reuse
        # subprocess.run(["bash", f"{lastz_bash}"])

        # 5. run palindrover
        bash_script = f"./tmp/{species}_{chromosome}_run_palindrover.sh"

        palindrome_file = f"./tmp/{species}_{chromosome}.pal"

        with open(bash_script, "w") as outfile:
           text = "python3 /storage/home/kxp5629/proj/10_HPRC_R2_Y/src/palindrover/palindrover.py \\\n"
           text = text + "    --minlength=8K \\\n"
           text = text + "    --minidentity=98% \\\n"
           text = text + "    --maxspacer=500K \\\n"
           text = text +f"    --blacklist:80%={repeats_file} \\\n"
           text = text + "    --column:palname \\\n"
           text = text + "    --column:blacklisted \\\n"
           text = text + "    --group:overlaps \\\n"
           text = text + "    --debug=blacklisted% \\\n"
           text = text +f"    --blacklisted=./tmp/{species}_{chromosome}.bl \\\n"
           text = text +f"    < {lastz_file} > {palindrome_file} "

           print(text, file = outfile)

        subprocess.run(["bash", f"{bash_script}"])

        # 6. convert palindrome output to bed file

        if not os.path.isdir("./palindromes/"):
            os.mkdir("./palindromes/")
        with open(f"./palindromes/{species}_{chr_dict[species]}.bed", 'w') as outfile:
            with open(palindrome_file, 'r') as infile:
                for line in infile:
                    line = line.strip()
                    if line.startswith("#"):
                        continue
                    line = line.split("\t")
                    print(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[9]}.A\t100\t-\tsize:{int(line[2]) - int(line[1])};repeat:{line[8]}", file = outfile)
                    print(f"{line[0]}\t{line[5]}\t{line[6]}\t{line[9]}.B\t100\t+\tsize:{int(line[6]) - int(line[5])};repeat:{line[8]}", file = outfile)


def main():

    #run_pipeline("chrY",y_chr)               
    run_pipeline("chrX",x_chr)

if __name__ == "__main__":
    main()
