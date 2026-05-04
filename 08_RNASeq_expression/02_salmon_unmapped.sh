#!/bin/bash
#SBATCH --job-name=salmon_unmapped
#SBATCH --output=logs/salmon_unmapped_%A_%a.out
#SBATCH --error=logs/salmon_unmapped_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=6:00:00
#SBATCH --array=0-5

# ==============================================================================
# Re-run Salmon with --writeUnmappedNames to capture unmapped read IDs
# Enables comparison of which reads gain mapping in the enriched reference
# ==============================================================================

set -euo pipefail

PROJECT_DIR="/storage/group/kdm16/default/kxp5629/proj/15_RNASeq"
INDEX_DIR="${PROJECT_DIR}/salmon_indices"
QUANT_DIR="${PROJECT_DIR}/salmon_quant_unmapped"
FASTQ_DIR="/storage/home/kxp5629/proj/15_RNASeq/data/fastq"

SPECIES_LIST=(HomSap PanTro PanPan GorGor PonAbe PonPyg)

declare -A LIBRARY_LAYOUT=(
    ["HomSap"]="PAIRED"
    ["PanTro"]="PAIRED"
    ["PanPan"]="SINGLE"
    ["GorGor"]="PAIRED"
    ["PonAbe"]="PAIRED"
    ["PonPyg"]="PAIRED"
)

declare -A FASTQ_PREFIX=(
    ["HomSap"]="Hsapien_16_07_A119b_testis_S1"
    ["PanTro"]="PanTro"
    ["PanPan"]="PanPan_SE"
    ["GorGor"]="GorGor"
    ["PonAbe"]="PonAbe"
    ["PonPyg"]="PonPyg"
)

SPECIES="${SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
LAYOUT="${LIBRARY_LAYOUT[$SPECIES]}"
PREFIX="${FASTQ_PREFIX[$SPECIES]}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

mkdir -p "${QUANT_DIR}/standard/${SPECIES}"
mkdir -p "${QUANT_DIR}/enriched/${SPECIES}"
mkdir -p "${PROJECT_DIR}/logs"

run_salmon() {
    local ref_type=$1
    local index="${INDEX_DIR}/${ref_type}/${SPECIES}"
    local output="${QUANT_DIR}/${ref_type}/${SPECIES}"
    shift 1
    local reads=("$@")

    log "Running Salmon (${ref_type}) for ${SPECIES}..."
    salmon quant \
        -i "${index}" \
        -l A \
        "${reads[@]}" \
        -o "${output}" \
        -p "${SLURM_CPUS_PER_TASK:-8}" \
        --validateMappings \
        --gcBias \
        --seqBias \
        --writeUnmappedNames
}

if [[ "${LAYOUT}" == "PAIRED" ]]; then
    R1="${FASTQ_DIR}/${PREFIX}_R1_001.fastq.gz"
    R2="${FASTQ_DIR}/${PREFIX}_R2_001.fastq.gz"
    if [[ ! -f "${R1}" ]]; then
        R1="${FASTQ_DIR}/${PREFIX}_R1.fastq.gz"
        R2="${FASTQ_DIR}/${PREFIX}_R2.fastq.gz"
    fi
    READS=(-1 "${R1}" -2 "${R2}")
else
    SE="${FASTQ_DIR}/${PREFIX}.fastq.gz"
    READS=(-r "${SE}")
fi

run_salmon "standard" "${READS[@]}"
run_salmon "enriched" "${READS[@]}"

log "Done: ${SPECIES}"
