#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_download_%A_%a.out
#SBATCH --error=logs/sra_download_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=8:00:00

# ==============================================================================
# Download RNA-seq data from SRA
# ==============================================================================
# Downloads paired-end FASTQ files for YAG expression analysis
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

PROJECT_DIR="/storage/group/kdm16/default/kxp5629/proj/15_RNASeq"
FASTQ_DIR="${PROJECT_DIR}/fastq"

# SRA toolkit path
SRA_TOOLKIT="${HOME}/tools/sratoolkit.3.1.1-ubuntu64/bin"
export PATH="${SRA_TOOLKIT}:${PATH}"

# SRA accessions mapped to species codes
declare -A SRA_ACCESSIONS=(
    ["GorGor"]="SRR3053573"
    ["PanTro"]="SRR2040590"
    ["PanPan"]="SRR306837"
    ["PonAbe"]="SRR10393299"
    ["PonPyg"]="SRR2176207"
    # ["HomSap"]=""  # TODO: Add human accession
)

# Library layout: PAIRED or SINGLE
declare -A LIBRARY_LAYOUT=(
    ["GorGor"]="PAIRED"
    ["PanTro"]="PAIRED"
    ["PanPan"]="SINGLE"
    ["PonAbe"]="PAIRED"
    ["PonPyg"]="PAIRED"
)

declare -A SPECIES_NAMES=(
    ["GorGor"]="Gorilla"
    ["PanTro"]="Chimpanzee"
    ["PanPan"]="Bonobo"
    ["PonAbe"]="Sumatran_orangutan"
    ["PonPyg"]="Bornean_orangutan"
    # ["HomSap"]="Human"
)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_dependencies() {
    log "Checking dependencies..."

    if ! command -v fasterq-dump &> /dev/null; then
        echo "Error: fasterq-dump not found at ${SRA_TOOLKIT}"
        exit 1
    fi

    log "SRA toolkit: $(fasterq-dump --version 2>&1 | head -1)"
}

create_directories() {
    log "Creating output directories..."
    mkdir -p "${FASTQ_DIR}"
    mkdir -p "${PROJECT_DIR}/logs"
}

download_sample() {
    local species=$1
    local accession="${SRA_ACCESSIONS[$species]}"
    local species_name="${SPECIES_NAMES[$species]}"
    local layout="${LIBRARY_LAYOUT[$species]}"

    log "Downloading ${accession} for ${species} (${species_name}) [${layout}]..."

    if [[ "${layout}" == "PAIRED" ]]; then
        local r1="${FASTQ_DIR}/${species}_R1.fastq.gz"
        local r2="${FASTQ_DIR}/${species}_R2.fastq.gz"

        # Skip if already downloaded
        if [[ -f "${r1}" ]] && [[ -f "${r2}" ]]; then
            log "FASTQ files already exist for ${species}, skipping..."
            return
        fi

        # Download paired-end
        fasterq-dump "${accession}" \
            --outdir "${FASTQ_DIR}" \
            --threads ${SLURM_CPUS_PER_TASK:-8} \
            --split-files \
            --progress

        # Rename to match expected naming convention
        log "Renaming files for ${species}..."
        mv "${FASTQ_DIR}/${accession}_1.fastq" "${FASTQ_DIR}/${species}_R1.fastq"
        mv "${FASTQ_DIR}/${accession}_2.fastq" "${FASTQ_DIR}/${species}_R2.fastq"

        # Compress
        log "Compressing FASTQ files for ${species}..."
        gzip "${FASTQ_DIR}/${species}_R1.fastq"
        gzip "${FASTQ_DIR}/${species}_R2.fastq"

    else
        # Single-end
        local se="${FASTQ_DIR}/${species}_SE.fastq.gz"

        # Skip if already downloaded
        if [[ -f "${se}" ]]; then
            log "FASTQ file already exists for ${species}, skipping..."
            return
        fi

        # Download single-end
        fasterq-dump "${accession}" \
            --outdir "${FASTQ_DIR}" \
            --threads ${SLURM_CPUS_PER_TASK:-8} \
            --progress

        # Rename
        log "Renaming file for ${species}..."
        mv "${FASTQ_DIR}/${accession}.fastq" "${FASTQ_DIR}/${species}_SE.fastq"

        # Compress
        log "Compressing FASTQ file for ${species}..."
        gzip "${FASTQ_DIR}/${species}_SE.fastq"
    fi

    log "Completed download for ${species}"
}

download_all() {
    log "=== Downloading SRA data ==="

    for species in "${!SRA_ACCESSIONS[@]}"; do
        download_sample "${species}"
    done
}

verify_downloads() {
    log "=== Verifying downloads ==="

    for species in "${!SRA_ACCESSIONS[@]}"; do
        local layout="${LIBRARY_LAYOUT[$species]}"

        if [[ "${layout}" == "PAIRED" ]]; then
            local r1="${FASTQ_DIR}/${species}_R1.fastq.gz"
            local r2="${FASTQ_DIR}/${species}_R2.fastq.gz"

            if [[ -f "${r1}" ]] && [[ -f "${r2}" ]]; then
                local r1_size=$(du -h "${r1}" | cut -f1)
                local r2_size=$(du -h "${r2}" | cut -f1)
                log "${species} [PAIRED]: R1=${r1_size}, R2=${r2_size}"
            else
                log "${species} [PAIRED]: MISSING"
            fi
        else
            local se="${FASTQ_DIR}/${species}_SE.fastq.gz"

            if [[ -f "${se}" ]]; then
                local se_size=$(du -h "${se}" | cut -f1)
                log "${species} [SINGLE]: SE=${se_size}"
            else
                log "${species} [SINGLE]: MISSING"
            fi
        fi
    done
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

usage() {
    echo "Usage: $0 [command]"
    echo ""
    echo "Commands:"
    echo "  download    Download all SRA samples (default)"
    echo "  verify      Verify downloaded files"
    echo "  single      Download single sample: $0 single <species_code>"
    echo ""
    echo "Species codes: GorGor, PanTro, PanPan, PonAbe, PonPyg"
}

main() {
    local cmd="${1:-download}"

    check_dependencies
    create_directories

    case "${cmd}" in
        download)
            download_all
            verify_downloads
            ;;
        verify)
            verify_downloads
            ;;
        single)
            local species="${2:-}"
            if [[ -z "${species}" ]] || [[ -z "${SRA_ACCESSIONS[$species]:-}" ]]; then
                echo "Error: Invalid species code"
                usage
                exit 1
            fi
            download_sample "${species}"
            ;;
        *)
            usage
            exit 1
            ;;
    esac

    log "Done!"
}

main "$@"
