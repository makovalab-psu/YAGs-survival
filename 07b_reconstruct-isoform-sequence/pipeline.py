#!/usr/bin/env python3
import argparse
import os
import subprocess
import yaml
import logging
import math
import sys
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

ALLOWED_SPECIES = ["HomSap", "PanPan", "PanTro", "PonPyg", "PonAbe", "GorGor"]
ALLOWED_GENES = ["TSPY", "RBMY", "CDY", "BPY2", "HSFY", "VCY", "RBMYB", "DAZ"]

SRC_PATH="/storage/home/kxp5629/proj/09_sequence-diversity/src"

def parse_arguments():
    """Parse command line arguments with restricted choices."""
    parser = argparse.ArgumentParser(description='Genomic alignment pipeline')
    parser.add_argument('--species', required=True, choices=ALLOWED_SPECIES,
                        help=f'Species name (one of: {"," .join(ALLOWED_SPECIES)})')
    parser.add_argument('--gene', required=True, choices=ALLOWED_GENES,
                        help=f'Gene name (one of: {",".join(ALLOWED_GENES)})')
    parser.add_argument('--config', required=True, help='Path to config file (YAML format)')
    parser.add_argument('--isoform', action='store_true', default=False, help='Also run isoform analysis')
    return parser.parse_args()

def read_config(config_path):
    """Read and parse the YAML configuration file."""
    try:
        with open(config_path, 'r') as config_file:
            config = yaml.safe_load(config_file)
        return config

    except FileNotFoundError:
        logger.error(f"Config file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        logger.error(f"Invalid YAML in config file: {config_path}\n{e}")
        sys.exit(1)

def check_output_exists(output_path):
    """Check if the output file already exists."""
    if os.path.exists(output_path):
        logger.info(f"Output file already exists: {output_path}")
        return True
    return False

#output has to be first
def run_minimap(output_path, reference_path, reads_path, threads=24):
    """Run minimap alignment."""
    command_map = [
        "minimap2",
        "-ax",
        "splice:hq",
        "--MD",
        "-t", str(math.floor(threads/2) or 1), # Number of threads
        reference_path,
        reads_path,
    ]

    command_convert = [
        "samtools",
        "view",
        "-@", str(math.floor(threads/4) or 1),
        "-b",
    ]

    command_sort = [
        "samtools",
        "sort",
        "-@", str(math.floor(threads/4) or 1),
        "-o", output_path
    ]

    logger.info(f"Running commands: {' '.join(command_map)}\n\t{' '.join(command_convert)}\n\t{' '.join(command_sort)}")

    try:
        with subprocess.Popen(command_map, stdout=subprocess.PIPE) as map_proc:
            with subprocess.Popen(command_convert, stdin=map_proc.stdout, stdout=subprocess.PIPE ) as convert_proc:
                with subprocess.Popen(command_sort, stdin=convert_proc.stdout, stdout=subprocess.PIPE) as sort_proc:
                    map_proc.stdout.close()
                    convert_proc.stdout.close()
                    sort_proc.communicate()


        logger.info(f"Alignment completed successfully: {output_path}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Alignment failed with error: {e}")
        return False
    except FileNotFoundError:
        logger.error("minimap2 command not found. Please ensure it's installed in your PATH")
        return False

def sam_to_fastq(output_path, bam_path, threads=4):
    """Run samtools fastq."""
    command = [
        "samtools",
        "fastq",
        "-@", str(threads), # Number of threads
        bam_path,
        "-0", output_path
    ]

    logger.info(f"Running: {' '.join(command)}")

    try:
        subprocess.run(command, check=True)
        logger.info(f"Conversion completed successfully: {output_path}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Conversion failed with error: {e}")
        return False
    except FileNotFoundError:
        logger.error("samtools command not found. Please ensure it's installed in your PATH")
        return False

def subset_single_gene_reads(out_file, bam_path, gene_bed):
    """Subset aligned reads to single genome"""
    command_subset = [
        "bedtools",
        "intersect",
        "-abam", f"{bam_path}",
        "-b", f"{gene_bed}",
        "-wa"
    ]
    command_convert_sam = [
        "samtools",
        "fastq",
        "-0", out_file
    ]

    try:
        with subprocess.Popen(command_subset, stdout=subprocess.PIPE) as subset_proc:
            with subprocess.Popen(command_convert_sam, stdin=subset_proc.stdout, stdout=subprocess.PIPE) as convert_proc:
                subset_proc.stdout.close()
                convert_proc.communicate()

        logger.info("Subsetting complete")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Subsetting failed with error: {e}")
        return False
    except FileNotFoundError:
        logger.error("samtools command not found. Please ensure it's installed in your PATH")
        return False

def extract_fingerprint(output_path, bam_genomic_path, bed_file, vcf_file, species, gene, threads=4):
    """Run extraction step."""
    command = [
        f"{SRC_PATH}/extract_FP.py",
        f"{bam_genomic_path}",
        f"{bed_file}",
        f"{vcf_file}",
        f"{output_path}",
        "--optimize",
        "--species",f"{species}",
        "--gene", f"{gene}"
    ]

    logger.info(f"Running: {' '.join(command)}")

    try:
        subprocess.run(command, check=True)
        logger.info(f"Extraction completed successfully: {output_path}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Extraction failed with error: {e}")
        return False
    except FileNotFoundError:
        logger.error("Problem occured while accessing the extract_FP.py script. Please ensure the script is accessible and executable.")
        return False

def merge_bam_files(output_bam, list_of_bams):
    """Combine all bam files from a list into one single bam file"""
    command = [
        "samtools", "merge",
        f"{output_bam}"
    ] + list_of_bams

    logger.info(f"Merging bam files: {' '.join(command)}")

    try:
        subprocess.run(command, check=True)
        logger.info("Bam files merged successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Bam files merge failed, {e}")
        return False
    except FileNotFoundError:
        logger.error("Samtools commandnot found. Please ensure it's installed in your PATH")
        return False

def call_variants(variants_file, reference, bam_file):
    command_call= [
        "freebayes",
        "-f", f"{reference}",
        f"{bam_file}"
    ]
    command_norm =  [
        "bcftools",
        "norm", "-m-any",
        "-f", f"{reference}",
        "-Ov",
        "-o", f"{variants_file}"
    ]

    try:

        logger.info(f"caling '{" ".join(command_call)}' which is piped into '{" ".join(command_norm)}'")
        with subprocess.Popen(command_call, stdout=subprocess.PIPE) as call_proc:
            with subprocess.Popen(command_norm, stdin=call_proc.stdout, stdout=subprocess.PIPE) as norm_proc:
                call_proc.stdout.close()
                norm_proc.communicate()

        logger.info("Subsetting complete")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Conversion failed with error: {e}")
        return False
    except FileNotFoundError:
        logger.error("freebayes or bcftools command not found. Please ensure it's installed in your PATH")
        return False

def run_fingerprint(output_path, bam_file, variant_positions, tab_out_file = None, report_all = False):
    command = [
        f"{SRC_PATH}/signature.py",
        "-i", f"{bam_file}",
        "-v", f"{variant_positions}",
        "-o", f"{output_path}"
    ]
    if report_all:
        command.append("-r")
    if tab_out_file != None:
        command += ["-t", f"{tab_out_file}"]
    try:
        subprocess.run(command, check=True)
        logger.info("fingerprint analysis finished.")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run fingerprint analysis, {e}")
        return False
    except FileNotFoundError:
        logger.error("Failed to find fingerprint command")
        return False

def count_FP(output, output_files):
    command = [
        f"{SRC_PATH}/count_FPs.py",
        f"{output}"
    ] + output_files

    try:
        subprocess.run(command, check=True)
        logger.info("final results count finished.")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Falied to run final count, {e}")
        return False
    except FileNotFoundError:
        logger.error("Failed to find counts command")
        return False

def get_gene_isoforms(output_file, gff_file, gene_id):

    subset_gene = gene_id
    if gene_id == "RBMYB":
        subset_gene = "RBMY"
    command = [
        "grep",
        f"{subset_gene}",
        gff_file
    ]
    try:
        #stdout to output_file
        subprocess.run(command, check=True, stdout=open(output_file,'w'))
        logger.info(f"Subset GFF {gff_file} to {gene_id}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to subset GFF , {e}")
        return False
    except FileNotFoundError:
        logger.error("Failed to subset GFF")
        return False

def concatenate_tab_files(output_file, *tab_files):
    with open(output_file, 'w') as outfile:
        for tab_file in tab_files:
            with open(tab_file, 'r') as infile:
                outfile.write(infile.read())
    logger.info(f"Concatenated tab files to {output_file}")
    return True

def gff_to_fasta(fasta_file, gff_file, reference_path):
    command = [
        "gffread",
        "-w",
        fasta_file,
        "-g",
        reference_path,
        gff_file
    ]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Converted GFF to FASTA: {fasta_file}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to convert GFF to FASTA, {e}")
        return False
    except FileNotFoundError:
        logger.error("Failed to find gffread command")
        return False

# python3
# 	../../src/helpers/process_isoforms.py homsap_TSPY_isoforms.sam
# 	../mapped/HomSap/HomSap_human2_targeted_t1_TSPY.bam
# 	../mapped/HomSap/HomSap_human2_targeted_t2_TSPY.bam
# 	../mapped/HomSap/HomSap_human2_untargeted_t1_TSPY.bam
# 	-s ../results/HomSap_TSPY_human2_fingerprints.tab
# 	-r ../../data/mapped/HomSap/HomSap_TSPY-genomic_alignment.bam \
# 	-j ../variants/signatures/HomSap_TSPY.json \
# 	-o HomSap_human2_output.txt -d HomSap_human2_discarded -l HomSap_human2
def run_isoforms(output, isoforms_file, bam_files, signature_tab, genomic_alignment_file, signature_json, discarded_prefix, logfile_prefix, used_reads):
    command = [
        f"{SRC_PATH}/helpers/process_isoforms_refGenome.py",
        *bam_files,
        "-s", signature_tab,
        "-i", isoforms_file,
        "-r", genomic_alignment_file,
        "-o", output,
        "-u", used_reads,
        "-d", discarded_prefix,
        "-l", logfile_prefix
    ]
    try:
        logger.info(f"Running {" ".join(command)}")
        subprocess.run(command, check=True)
        logger.info(f"Processed isoforms: {output}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to process isoforms, {" ".join(command)} \n {e}")
        return False

def extract_tail(output, outpus_full_sequences, input_files):
    command = [
        f"{SRC_PATH}/extract_DAZ_tail.py",
        *input_files,
        "-o", output,
        "-s", outpus_full_sequences

    ]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Extract DAZ tail: {output}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to extract DAZ tails, {" ".join(command)} \n {e}")
        return False

def hmmbuild(output, msa):
    command = [
        "hmmbuild",
        output,
        msa
    ]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Build hmm model: {output}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build hmm model, {" ".join(command)} \n {e}")
        return False

def nhmmer(output, hmm_file, extracted_sequences):
    command = [
        "nhmmer",
        "--tblout", output,
        hmm_file,
        extracted_sequences

    ]
    try:
        subprocess.run(command, check=True)
        logger.info(f"Build hmm model: {output}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build hmm model, {" ".join(command)} \n {e}")
        return False

def process_daz_table(output, table, species, fasta, sample_list):
    command = [
        f"{SRC_PATH}/process_DAZ_table.py",
        "-o", output,
        "-i", table,
        "-f", fasta,
        "-s", species,
        *sample_list
    ]

    try:
        subprocess.run(command, check=True)
        logger.info(f"Processing table of identified DAZ repeats {output}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to proces DAZ repeats table, {" ".join(command)} \n {e}")
        return False

def run_consensus(output, reads_file, mapping_file, counts_table, plot_name, species):
    command = [
        f"{SRC_PATH}/helpers/consensus_sequences_spoa.py",
        "-o", output,
        "-c", counts_table,
        "-p", plot_name,
        "-s", species,
        reads_file,
        mapping_file
    ]

    try:
        logger.info(f"Processing consensus sequences...")
        subprocess.run(command, check=True)
        logger.info(f"Processing consensus sequences Done.")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Falied to call consensus sequences {e}")
        return False


def check_output_and_launch(output, command, args):
    if check_output_exists(output):
        logger.info(f"Skipping command {command}, output file exists: {output}")
    else:
        logger.info(f"Running command {command}")

        success = command(output, *args)
        if not success:
            logger.error("Error executing command")
            sys.exit(1)

def main():
    # Parse command line arguments
    args = parse_arguments()

    #Read configuration
    logger.info(f"Read YAML config from {args.config}")
    config = read_config(args.config)

    # Extract paths from config
    try:
        reference_path =config['species_genomes'][args.species]
        samples_all = config['sample_species']
        samples = [x for x in samples_all if samples_all[x] == args.species]
        output_dir_mapped = f"data/mapped/{args.species}"
        output_dir_variants = "data/variants"
        # create dir if does not exist
        output_dir_finger = "data/variants/signatures"
        Path(output_dir_finger).mkdir(parents=True, exist_ok=True)
        output_dir_results = f"data/results/{args.species}/{args.gene}"
        Path(output_dir_results).mkdir(parents=True, exist_ok=True)
        gene_family_bed = config["genomic_bed"][args.gene][args.species]
        single_gene_bed = config["single_gene_ref_bed"][args.gene][args.species]
    except KeyError as e:
        logger.error(f"Missing configuration: {e}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    Path(output_dir_mapped).mkdir(parents=True, exist_ok=True)
    read_dir="data/reads"
    Path(read_dir).mkdir(parents=True, exist_ok=True)

    i = 1
    # list of bam files which will be merged into one bam file for variant calling (to get
    # fingerprint info) and then analyzed separately by the fingerprint analysis.
    sample_bam_files_with_sample_name_dict = {}
    sample_bam_file_with_sample_name_dict_algned_to_ref = {}
    list_of_single_gene_bam_files = [] # is redundant use list dict from above and remove
    single_gene_ref = config.get("single_gene_reference","").get(args.gene,"").get(args.species,"")
    if single_gene_ref =="":
        logger.error(f"Didn't find single gene ref for {args.species}/{args.gene}")
        sys.exit(1)
    for sample in samples:

        # # after migrating to roar I downloaded the files from SRA and the input is now fq, conversion not needed
        # reads_path = f"data/reads/{args.species}_{sample}.fq"
        # sample_file = config["sample_files"][sample]
        # check_output_and_launch(reads_path, sam_to_fastq, [sample_file, 24])

        reads_path = config["sample_files"][sample]

        alignment_output = os.path.join(output_dir_mapped, f"{args.species}_{sample}_alignment.bam")
        check_output_and_launch(alignment_output, run_minimap, [reference_path, reads_path])



        single_gene_sequencing_data = f"{read_dir}/{args.species}_{args.gene}_{sample}.fq"
        check_output_and_launch(single_gene_sequencing_data, subset_single_gene_reads, [alignment_output, gene_family_bed ])


        single_gene_alignment_output = alignment_output.replace("alignment.bam", f"{args.gene}.bam")
        list_of_single_gene_bam_files.append(single_gene_alignment_output)
        sample_bam_files_with_sample_name_dict[f"{sample}"] = single_gene_alignment_output
        check_output_and_launch(single_gene_alignment_output, run_minimap, [single_gene_ref, single_gene_sequencing_data])

        whole_ref_alignment_subset_reads = alignment_output.replace("alignment.bam", f"{args.gene}_to_ref.bam")
        check_output_and_launch(whole_ref_alignment_subset_reads, run_minimap, [reference_path, single_gene_sequencing_data])
        sample_bam_file_with_sample_name_dict_algned_to_ref[f"{sample}"] = whole_ref_alignment_subset_reads

        logger.info(f"Alignment completed for {i}/{len(samples)} samples")
        i+=1

    genomic_seq_aligned = os.path.join(output_dir_mapped, f"{args.species}_{args.gene}-genomic_alignment.bam")
    genomic_sequences = config.get("genomic_sequences", {}).get(args.gene, {}).get(args.species, "")
    if genomic_sequences == "":
        logger.error(f"Didn't find genomic sequences for {args.species}/{args.gene}")
    check_output_and_launch(genomic_seq_aligned, run_minimap, [single_gene_ref, genomic_sequences])


    if args.gene == "DAZ":
        print("processing DAZ")
        #extract tail
        logger.info(f"Extract DAZ tail in {args.species} samples")
        sample_list = []
        for sample in samples:
            single_gene_sequencing_data = f"{read_dir}/{args.species}_{args.gene}_{sample}.fq"
            sample_list.append(single_gene_sequencing_data)

        extracted_sequences = f"{read_dir}/{args.species}_{args.gene}_extracted-tail.fq"
        # this might be unneccessary duplication, maybe I could use the full sequences for extracting DAZ repeats as well?
        extracted_sequences_full_reads = f"{read_dir}/{args.species}_{args.gene}_extracted-tail_full-sequence.fq"
        check_output_and_launch(extracted_sequences, extract_tail, [extracted_sequences_full_reads, sample_list])

        daz_msa_file = config["DAZ_msa"][args.species]
        daz_hmm_file = config["DAZ_msa"][args.species].replace(".msa",".hmm")

        check_output_and_launch(daz_hmm_file, hmmbuild, [daz_msa_file])

        daz_hmm_table = extracted_sequences.replace(".fq", ".tbl")
        check_output_and_launch(daz_hmm_table, nhmmer, [daz_hmm_file, extracted_sequences])

        daz_analysis_out = os.path.join(output_dir_results, f"{args.species}_{args.gene}-tail")
        check_output_and_launch(daz_analysis_out, process_daz_table, [daz_hmm_table ,args.species, extracted_sequences, sample_list])


    if args.gene != "DAZ":
        merged_bam_file = os.path.join(output_dir_mapped, f"{args.species}_{args.gene}-merged.bam")
        check_output_and_launch(merged_bam_file, merge_bam_files,[list_of_single_gene_bam_files])

        variants_file = os.path.join(output_dir_variants, f"{args.species}_{args.gene}-variants.vcf")
        check_output_and_launch(variants_file, call_variants, [single_gene_ref, merged_bam_file])

        fingerprint_file = os.path.join(output_dir_finger, f"{args.species}_{args.gene}.json")
        check_output_and_launch(fingerprint_file, extract_fingerprint, [genomic_seq_aligned, single_gene_bed, variants_file, args.species, args.gene])
        # logger.info("Pipeline completed")

        json_results = []
        fingerprint_tabs = []
        for sample in sample_bam_files_with_sample_name_dict:
            bam_file = sample_bam_files_with_sample_name_dict[sample]
            result_path = os.path.join(output_dir_results, f"{args.species}_{args.gene}_{sample}-fingerprint_results.json")
            json_results.append(result_path)
            tab_out_file = os.path.join(output_dir_results, f"{args.species}_{args.gene}_{sample}_fingerprints.tab")
            fingerprint_tabs.append(tab_out_file)
            check_output_and_launch(result_path, run_fingerprint, [bam_file, fingerprint_file, tab_out_file],)

        result_path = os.path.join(output_dir_results, f"{args.species}_{args.gene}_ref-sequences.json")
        json_results.append(result_path)
        tab_out_file = os.path.join(output_dir_results, f"{args.species}_{args.gene}_ref-sequences.tab")
        fingerprint_tabs.append(tab_out_file)
        check_output_and_launch(result_path, run_fingerprint, [genomic_seq_aligned, fingerprint_file, tab_out_file, True])
        counts_table = os.path.join(output_dir_results, f"{args.species}_{args.gene}_counts_table.tab")
        check_output_and_launch(counts_table, count_FP, [json_results])

        if args.isoform:
            isoforms_gff = config.get("isoform_file",{}).get(args.species,"")
            parent_dir = os.path.dirname(isoforms_gff)
            # isoforms_analysis_dir =f"{parent_dir}/{args.species}_{args.gene}/"
            isoforms_analysis_dir = f"{output_dir_results}/isoforms/"
            Path(isoforms_analysis_dir).mkdir(parents=True, exist_ok=True)
            isoforms_gene_gff = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms.gff")
            check_output_and_launch(isoforms_gene_gff, get_gene_isoforms, [isoforms_gff, args.gene])
            isoforms_gene_fasta = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms.fasta")
            check_output_and_launch(isoforms_gene_fasta, gff_to_fasta, [isoforms_gene_gff, reference_path.replace(".gz_splice-hq.mmi", "")])

            # align isoforms to reference
            isoforms_gene_bam = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms.bam")
            check_output_and_launch(isoforms_gene_bam, run_minimap, [single_gene_ref, isoforms_gene_fasta])

            # concatenate fingerprint tab files
            fingerprint_single_tab = os.path.join(output_dir_results, f"{args.species}_{args.gene}_fingerprint.tab")
            check_output_and_launch(fingerprint_single_tab, concatenate_tab_files, fingerprint_tabs)

            ## for human subset for human2 samples also chimp2 for chimp
            sample_list = list(sample_bam_files_with_sample_name_dict.values())
            sample_list_ref = list(sample_bam_file_with_sample_name_dict_algned_to_ref.values())
            print(sample_bam_files_with_sample_name_dict)
            if args.species == "HomSap":
                samples_h2 = {sample: sample_bam_files_with_sample_name_dict[sample] for sample in sample_bam_files_with_sample_name_dict if sample.startswith("human2")}
                sample_list = list(samples_h2.values())

                samples_h2_ref = {sample:sample_bam_file_with_sample_name_dict_algned_to_ref[sample] for sample in sample_bam_file_with_sample_name_dict_algned_to_ref if sample.startswith("human2")}
                sample_list_ref =list(samples_h2_ref.values())

            if args.species == "PanTro":
                samples_c2 = {sample: sample_bam_files_with_sample_name_dict[sample] for sample in sample_bam_files_with_sample_name_dict if sample.startswith("chimp2")}
                sample_list = list(samples_c2.values())

                samples_c2_ref = {sample:sample_bam_file_with_sample_name_dict_algned_to_ref[sample] for sample in sample_bam_file_with_sample_name_dict_algned_to_ref if sample.startswith("chimp2")}
                sample_list_ref =list(samples_c2_ref.values())

            # switch to whole genome alignment instead of single gene reference sequence
            genomic_seq_aligned_to_ref = os.path.join(output_dir_mapped, f"{args.species}_{args.gene}-genomic_alignment_to_ref.bam")
            check_output_and_launch(genomic_seq_aligned_to_ref, run_minimap, [reference_path, genomic_sequences])

            isoforms_ref_bam = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms_ref.bam")
            check_output_and_launch(isoforms_ref_bam, run_minimap, [reference_path, isoforms_gene_fasta])

            #also have to align the reads again, this time to whole genome

            # # Process isoforms
            isoforms_output = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms_processed.txt")
            isoforms_discarded = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms_discarded")
            isoforms_logfile = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_isoforms_log")
            isoforms_used_reads = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_iso_sig_read.tsv")
            # # run_isoforms(                output,                    isoforms_file,  bam_files,        signature_tab, genomic_alignment_file, signature_json, discarded_prefix, logfile_prefix):
            check_output_and_launch(isoforms_output, run_isoforms, [isoforms_gene_gff, sample_list_ref, fingerprint_single_tab, genomic_seq_aligned_to_ref, fingerprint_file, isoforms_discarded, isoforms_logfile, isoforms_used_reads])

            protein_sequences_output = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_ORFs.fa")
            plot_name = os.path.join(isoforms_analysis_dir, f"{args.species}_{args.gene}_iso.png")
            check_output_and_launch(protein_sequences_output, run_consensus, [fingerprint_single_tab, isoforms_used_reads, isoforms_output, plot_name, args.species])

if __name__ == "__main__":
    main()
