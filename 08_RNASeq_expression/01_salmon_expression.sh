#!/bin/bash
#SBATCH --job-name=salmon_yag
#SBATCH --output=logs/salmon_%A_%a.out
#SBATCH --error=logs/salmon_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00

# ==============================================================================
# YAG RNA-seq Expression Analysis with Salmon
# ==============================================================================
# Quantifies gene expression using Salmon pseudoalignment
# Compares standard CDS reference vs enriched reference (CDS + novel YAG transcripts)
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

# Base directories
PROJECT_DIR="/storage/group/kdm16/default/kxp5629/proj/15_RNASeq"
REF_BASE="${HOME}/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq"

# Output directories
INDEX_DIR="${PROJECT_DIR}/salmon_indices"
QUANT_DIR="${PROJECT_DIR}/salmon_quant"

# Novel YAG transcripts directory (species-specific files: {species}_YAG_cds.fasta)
NOVEL_YAG_DIR="/storage/home/kxp5629/proj/15_RNASeq/data/novel_transcripts"

# RNA-seq data directory
FASTQ_DIR="/storage/home/kxp5629/proj/15_RNASeq/data/fastq"

# ------------------------------------------------------------------------------
# Species Configuration
# ------------------------------------------------------------------------------

declare -A SPECIES_CDS=(
    ["HomSap"]="${REF_BASE}/HomSap/data/GCF_009914755.1/cds_from_genomic.fna"
    ["PanTro"]="${REF_BASE}/PanTro/data/GCF_028858775.2/cds_from_genomic.fna"
    ["PanPan"]="${REF_BASE}/PanPan/data/GCF_029289425.2/cds_from_genomic.fna"
    ["GorGor"]="${REF_BASE}/GorGor/data/GCF_029281585.2/cds_from_genomic.fna"
    ["PonAbe"]="${REF_BASE}/PonAbe/data/GCF_028885655.2/cds_from_genomic.fna"
    ["PonPyg"]="${REF_BASE}/PonPyg/data/GCF_028885625.2/cds_from_genomic.fna"
)

declare -A SPECIES_NAMES=(
    ["HomSap"]="Human"
    ["PanTro"]="Chimpanzee"
    ["PanPan"]="Bonobo"
    ["GorGor"]="Gorilla"
    ["PonAbe"]="Sumatran_orangutan"
    ["PonPyg"]="Bornean_orangutan"
)

# Library layout: PAIRED or SINGLE
declare -A LIBRARY_LAYOUT=(
    ["HomSap"]="PAIRED"
    ["PanTro"]="PAIRED"
    ["PanPan"]="SINGLE"
    ["GorGor"]="PAIRED"
    ["PonAbe"]="PAIRED"
    ["PonPyg"]="PAIRED"
)

# FASTQ file prefixes (handles different naming conventions)
# For PAIRED: expects {prefix}_R1.fastq.gz and {prefix}_R2.fastq.gz
# For SINGLE: expects {prefix}.fastq.gz
declare -A FASTQ_PREFIX=(
    ["HomSap"]="Hsapien_16_07_A119b_testis_S1"
    ["PanTro"]="PanTro"
    ["PanPan"]="PanPan_SE"
    ["GorGor"]="GorGor"
    ["PonAbe"]="PonAbe"
    ["PonPyg"]="PonPyg"
)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_dependencies() {
    log "Checking dependencies..."
    if ! command -v salmon &> /dev/null; then
        echo "Error: salmon not found. Please load the salmon module."
        echo "  module load salmon"
        exit 1
    fi
    log "Salmon version: $(salmon --version 2>&1 | head -1)"
}

create_directories() {
    log "Creating output directories..."
    mkdir -p "${INDEX_DIR}/standard"
    mkdir -p "${INDEX_DIR}/enriched"
    mkdir -p "${QUANT_DIR}/standard"
    mkdir -p "${QUANT_DIR}/enriched"
    mkdir -p "${PROJECT_DIR}/logs"
}

# ------------------------------------------------------------------------------
# Phase 1: Build Salmon Indices
# ------------------------------------------------------------------------------

build_index_standard() {
    local species=$1
    local cds_fasta="${SPECIES_CDS[$species]}"
    local index_out="${INDEX_DIR}/standard/${species}"

    if [[ -d "${index_out}" ]]; then
        log "Index already exists for ${species} (standard), skipping..."
        return
    fi

    log "Building standard index for ${species}..."
    salmon index \
        -t "${cds_fasta}" \
        -i "${index_out}" \
        -k 31 \
        -p ${SLURM_CPUS_PER_TASK:-8}
}

build_index_enriched() {
    local species=$1
    local cds_fasta="${SPECIES_CDS[$species]}"
    local novel_fasta="${NOVEL_YAG_DIR}/${species}_YAG_cds.fasta"
    local index_out="${INDEX_DIR}/enriched/${species}"
    local combined_fasta="${INDEX_DIR}/enriched/${species}_combined.fna"

    if [[ ! -f "${novel_fasta}" ]]; then
        log "Warning: Novel YAG file not found for ${species}: ${novel_fasta}"
        return 1
    fi

    if [[ -d "${index_out}" ]]; then
        log "Index already exists for ${species} (enriched), skipping..."
        return
    fi

    log "Building enriched index for ${species}..."

    # Combine CDS with species-specific novel YAG transcripts
    cat "${cds_fasta}" "${novel_fasta}" > "${combined_fasta}"

    salmon index \
        -t "${combined_fasta}" \
        -i "${index_out}" \
        -k 31 \
        -p ${SLURM_CPUS_PER_TASK:-8}
}

build_all_indices() {
    log "=== Building Salmon Indices ==="

    for species in "${!SPECIES_CDS[@]}"; do
        build_index_standard "${species}"
    done

    log "=== Building Enriched Indices ==="
    for species in "${!SPECIES_CDS[@]}"; do
        build_index_enriched "${species}"
    done
}

# ------------------------------------------------------------------------------
# Phase 2: Quantification
# ------------------------------------------------------------------------------

run_quantification_paired() {
    local species=$1
    local ref_type=$2  # "standard" or "enriched"
    local r1=$3
    local r2=$4

    local index="${INDEX_DIR}/${ref_type}/${species}"
    local output="${QUANT_DIR}/${ref_type}/${species}"

    if [[ -f "${output}/quant.sf" ]]; then
        log "Quantification already exists for ${species} (${ref_type}), skipping..."
        return
    fi

    log "Running Salmon quantification for ${species} (${ref_type}) [PAIRED]..."
    salmon quant \
        -i "${index}" \
        -l A \
        -1 "${r1}" \
        -2 "${r2}" \
        -o "${output}" \
        -p ${SLURM_CPUS_PER_TASK:-8} \
        --validateMappings \
        --gcBias \
        --seqBias
}

run_quantification_single() {
    local species=$1
    local ref_type=$2  # "standard" or "enriched"
    local reads=$3

    local index="${INDEX_DIR}/${ref_type}/${species}"
    local output="${QUANT_DIR}/${ref_type}/${species}"

    if [[ -f "${output}/quant.sf" ]]; then
        log "Quantification already exists for ${species} (${ref_type}), skipping..."
        return
    fi

    log "Running Salmon quantification for ${species} (${ref_type}) [SINGLE]..."
    salmon quant \
        -i "${index}" \
        -l A \
        -r "${reads}" \
        -o "${output}" \
        -p ${SLURM_CPUS_PER_TASK:-8} \
        --validateMappings \
        --gcBias \
        --seqBias
}

quantify_all_samples() {
    log "=== Running Salmon Quantification ==="

    if [[ -z "${FASTQ_DIR}" ]]; then
        log "Error: FASTQ_DIR not set. Please set the path to FASTQ files."
        return 1
    fi

    for species in "${!SPECIES_CDS[@]}"; do
        local layout="${LIBRARY_LAYOUT[$species]}"
        local prefix="${FASTQ_PREFIX[$species]}"

        if [[ "${layout}" == "PAIRED" ]]; then
            local r1="${FASTQ_DIR}/${prefix}_R1_001.fastq.gz"
            local r2="${FASTQ_DIR}/${prefix}_R2_001.fastq.gz"

            # Try alternative naming pattern if _001 suffix not found
            if [[ ! -f "${r1}" ]]; then
                r1="${FASTQ_DIR}/${prefix}_R1.fastq.gz"
                r2="${FASTQ_DIR}/${prefix}_R2.fastq.gz"
            fi

            if [[ ! -f "${r1}" ]] || [[ ! -f "${r2}" ]]; then
                log "Warning: FASTQ files not found for ${species} (prefix: ${prefix}), skipping..."
                continue
            fi

            run_quantification_paired "${species}" "standard" "${r1}" "${r2}"

            if [[ -d "${INDEX_DIR}/enriched/${species}" ]]; then
                run_quantification_paired "${species}" "enriched" "${r1}" "${r2}"
            fi
        else
            # Single-end
            local se="${FASTQ_DIR}/${prefix}.fastq.gz"

            if [[ ! -f "${se}" ]]; then
                log "Warning: FASTQ file not found for ${species} (prefix: ${prefix}), skipping..."
                continue
            fi

            run_quantification_single "${species}" "standard" "${se}"

            if [[ -d "${INDEX_DIR}/enriched/${species}" ]]; then
                run_quantification_single "${species}" "enriched" "${se}"
            fi
        fi
    done
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

usage() {
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  index     Build Salmon indices (standard and enriched)"
    echo "  quant     Run Salmon quantification on all samples"
    echo "  all       Run complete pipeline (index + quant)"
}

main() {
    local cmd="${1:-}"

    check_dependencies
    create_directories

    case "${cmd}" in
        index)
            build_all_indices
            ;;
        quant)
            quantify_all_samples
            ;;
        all)
            build_all_indices
            quantify_all_samples
            ;;
        *)
            usage
            exit 1
            ;;
    esac

    log "Done!"
}

main "$@"
