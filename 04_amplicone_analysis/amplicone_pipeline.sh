#!/usr/bin/env bash
# =============================================================================
# AmpliCoNE Y-chromosome Copy Number Estimation Pipeline
# =============================================================================
#
# PURPOSE:
#   Estimates copy numbers of multicopy Y-chromosome gene families (TSPY, VCY,
#   CDY, HSFY, RBMY, BPY2, DAZ) in primate species using short-read whole-genome
#   sequencing data and the AmpliCoNE-count tool.
#
# SPECIES ANALYZED (final run Nov 2024):
#   - Gorilla gorilla (GorGor) — 18 individuals
#   - Pan troglodytes (PanTro) — 21 individuals
#
# PIPELINE OVERVIEW:
#   Step 1  — Decompress reference genomes if only .gz is present
#   Step 1b — Build GEM mappability index for each reference (run once)
#   Step 2  — Extract chrY FASTA and run Tandem Repeat Finder
#   Step 2.1— Convert TRF .dat output to BED format
#   Step 3  — Extract chrY rows from RepeatMasker .out files
#   Step 4  — Build AmpliCoNE annotation .tab files (AmpliCoNE-build.py)
#   Step 5  — Batch-run AmpliCoNE-count on all individuals (03_run_amplicone.py)
#   Step 6  — Aggregate per-individual results into summary TSV (04_collect_results.py)
#
# CONFIGURATION:
#   Edit config.yaml (in this directory) to set all input paths and tool
#   locations before running this script.  No paths need to be changed here.
#
# PREREQUISITE TOOLS (must be on PATH or configured in config.yaml):
#   - samtools    https://www.htslib.org/
#   - bgzip       (part of htslib)
#   - TRF         https://tandem.bu.edu/trf/trf.html
#   - GEM suite   https://sourceforge.net/projects/gemlibrary/
#   - wig2bed     (part of BEDOPS) https://bedops.readthedocs.io/
#   - Python 3 + pyyaml  (pip install pyyaml)
#   - AmpliCoNE-tool in AmpliCoNE-tool/ subdirectory
#
# PREREQUISITE PYTHON PACKAGES:
#   pip install pyyaml bcbio-gff
#
# OUTPUT:
#   results/ampl_results.YYYYMMDD.tsv  — copy numbers per individual per gene
#
# NOTES:
#   - GEM mappability computed with 101 bp reads, -m 2 -e 2, 32 threads.
#   - TRF parameters: 2 5 7 80 10 50 2000 -l 10 -h -d
#   - GorGor uses RepeatMasker RM1; RM2 was too aggressive (masked BPY2).
#   - Gorilla VCY is absent — column is empty in results.
#   - Gene definitions required manual pseudogene additions (see gene_definitions/).
#   - AmpliCoNE-count.py requires Python 2.7; all wrappers use Python 3.
#
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve pipeline root (directory containing this script)
# ---------------------------------------------------------------------------
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"
CONFIG="${PIPELINE_DIR}/config.yaml"

# Relative working directories (created automatically as needed)
GENE_DEFS_DIR="${PIPELINE_DIR}/gene_definitions/targets"
AMPLICONE_TOOL="${PIPELINE_DIR}/AmpliCoNE-tool"
Y_CHR_DIR="${PIPELINE_DIR}/y_chromosomes"
TRF_DIR="${PIPELINE_DIR}/trf_output"
RM_DIR="${PIPELINE_DIR}/repeat_masker"
PANELS_DIR="${PIPELINE_DIR}/panels"
RESULTS_DIR="${PIPELINE_DIR}/results"
TMP_DIR="${PIPELINE_DIR}/tmp"

# ---------------------------------------------------------------------------
# Load all external paths and species parameters from config.yaml
# ---------------------------------------------------------------------------
# This eval sources variables of the form:
#   TRF_BIN, GEM_BIN_DIR, BEDOPS_BIN_DIR, ALIGNMENTS_ROOT
#   SPECIES  (space-separated list)
#   {SPECIES}_REFERENCE, {SPECIES}_GFF, {SPECIES}_REPEATMASKER_OUT,
#   {SPECIES}_CHRY_RM_TAG, {SPECIES}_CHRY_ACCESSION, {SPECIES}_CHRY_LENGTH,
#   {SPECIES}_GENE_DEFINITIONS, {SPECIES}_ANNOTATION_TAB
eval "$(python3 "${SCRIPTS_DIR}/config_to_env.py" "${CONFIG}")"

# Add GEM and BEDOPS to PATH (values come from config.yaml)
export PATH="${GEM_BIN_DIR}:${BEDOPS_BIN_DIR}:${PATH}"

# ---------------------------------------------------------------------------
# Helper: get a per-species config variable via indirect expansion
#   e.g.  sp_var GorGor reference  →  value of GORGOR_REFERENCE
# ---------------------------------------------------------------------------
sp_var() {
    local sp="$1" field="$2"
    local varname="${sp^^}_${field^^}"   # bash 4+: ^^ = uppercase
    echo "${!varname}"
}

# ---------------------------------------------------------------------------
# STEP 1 — DECOMPRESS REFERENCES (if only .gz exists)
# ---------------------------------------------------------------------------
# GEM indexer and some samtools operations require plain (non-bgzipped) FASTA.
# Skip if the uncompressed file already exists.

module load htslib   # provides bgzip; remove if bgzip is already on PATH

decompress_ref() {
    local ref="$1"
    if [ ! -f "${ref}" ] && [ -f "${ref}.gz" ]; then
        echo "[Step 1] Decompressing ${ref}.gz ..."
        bgzip -d "${ref}.gz" -c > "${ref}"
    fi
}

for sp in $SPECIES; do
    decompress_ref "$(sp_var "$sp" reference)"
done

# ---------------------------------------------------------------------------
# STEP 1b — GEM MAPPABILITY INDEX (run once per reference)
# ---------------------------------------------------------------------------
# Produces a per-base mappability BED file ({ref}.bed) used by AmpliCoNE-build.
# Parameters: read length 101 bp, up to 2 mismatches and 2 indels, 32 threads.
# Runtime: ~50–75 min per genome even with 32 threads.

build_mappability() {
    local ref="$1"
    if [ ! -f "${ref}.bed" ]; then
        echo "[Step 1b] Building GEM mappability index for $(basename "${ref}") ..."
        gem-indexer      -i "${ref}" -o "${ref}" -T 32
        gem-mappability  -I "${ref}.gem" -l 101 -o "${ref}" -m 2 -e 2 -T 32
        gem-2-wig        -I "${ref}.gem" -i "${ref}.mappability" -o "${ref}"
        wig2bed < "${ref}.wig" > "${ref}.bed"
    else
        echo "[Step 1b] Mappability BED exists for $(basename "${ref}"), skipping."
    fi
}

for sp in $SPECIES; do
    build_mappability "$(sp_var "$sp" reference)"
done

# ---------------------------------------------------------------------------
# STEP 2 — EXTRACT chrY FASTA & RUN TANDEM REPEAT FINDER
# ---------------------------------------------------------------------------
# TRF identifies simple tandem repeats that are later masked in the annotation.
# TRF parameters: 2 5 7 80 10 50 2000 -l 10 -h -d
#   2    = match weight          80   = match probability (%)
#   5    = mismatch penalty      10   = indel probability (%)
#   7    = indel penalty         50   = minimum alignment score
#   2000 = maximum period size   -l 10 = max TR array length (Mbp)
#   -h   = suppress HTML         -d   = write .dat data file

module load samtools   # remove if samtools is already on PATH

mkdir -p "${Y_CHR_DIR}" "${TRF_DIR}"

run_trf() {
    local sp="$1"
    local ref accession fa dat
    ref="$(sp_var "$sp" reference)"
    accession="$(sp_var "$sp" chry_accession)"
    fa="${Y_CHR_DIR}/${sp}_${accession}.fa"
    # TRF appends its parameters to the input filename to form the .dat name
    dat="${fa}.2.5.7.80.10.50.2000.dat"

    if [ ! -f "${fa}" ]; then
        echo "[Step 2] Extracting ${accession} from ${sp} reference ..."
        samtools faidx "${ref}" "${accession}" > "${fa}"
    fi

    if [ ! -f "${TRF_DIR}/$(basename "${dat}")" ]; then
        echo "[Step 2] Running TRF on ${fa} ..."
        # TRF writes output next to the input file; move to trf_output/ afterwards
        "${TRF_BIN}" "${fa}" 2 5 7 80 10 50 2000 -l 10 -h -d
        mv "${dat}" "${TRF_DIR}/"
    else
        echo "[Step 2] TRF output exists for ${sp}, skipping."
    fi
}

for sp in $SPECIES; do
    run_trf "$sp"
done

# ---------------------------------------------------------------------------
# STEP 2.1 — CONVERT TRF .dat TO BED FORMAT
# ---------------------------------------------------------------------------

convert_trf() {
    local sp="$1"
    local accession dat bed
    accession="$(sp_var "$sp" chry_accession)"
    dat="${TRF_DIR}/${sp}_${accession}.fa.2.5.7.80.10.50.2000.dat"
    bed="${dat%.dat}.bed"

    if [ ! -f "${bed}" ]; then
        echo "[Step 2.1] Converting TRF .dat to BED for ${sp} ..."
        python3 "${SCRIPTS_DIR}/01_convert_trf_to_bed.py" "${dat}" "${bed}"
    else
        echo "[Step 2.1] TRF BED exists for ${sp}, skipping."
    fi
}

for sp in $SPECIES; do
    convert_trf "$sp"
done

# ---------------------------------------------------------------------------
# STEP 3 — EXTRACT chrY ROWS FROM REPEATMASKER OUTPUT
# ---------------------------------------------------------------------------
# Filters the genome-wide RepeatMasker .out to chrY only, renaming the
# chromosome label to the NCBI accession used in this pipeline.
#
# Note: the chrY label in the RM file differs between species (see config.yaml).
# GorGor RM1: chrY_pat_hsaY   PanTro RM1: chrY_hap2_hsaY
# WARNING: RM2 was found to mask coding regions too aggressively (e.g. BPY2);
# RM1 outputs were used for the final GorGor and PanTro runs.

mkdir -p "${RM_DIR}"

extract_repeatmasker() {
    local sp="$1"
    local rm_src rm_tag accession out
    rm_src="$(sp_var "$sp" repeatmasker_out)"
    rm_tag="$(sp_var "$sp" chry_rm_tag)"
    accession="$(sp_var "$sp" chry_accession)"
    out="${RM_DIR}/${sp}_chrY_rm.out"

    if [ ! -f "${out}" ]; then
        echo "[Step 3] Extracting ${sp} chrY RepeatMasker rows ..."
        grep "${rm_tag}" "${rm_src}" \
            | sed "s/${rm_tag}/${accession}/g" \
            > "${out}"
    else
        echo "[Step 3] RepeatMasker output exists for ${sp}, skipping."
    fi
}

for sp in $SPECIES; do
    extract_repeatmasker "$sp"
done

# ---------------------------------------------------------------------------
# STEP 4 — BUILD AMPLICONE ANNOTATION FILES
# ---------------------------------------------------------------------------
# AmpliCoNE-build.py integrates mappability, TRF, RepeatMasker, and gene
# definitions to label each chrY window as ampliconic, single-copy control,
# or excluded (repeat / low-mappability).
#
# Gene definition TSVs (gene_definitions/targets/) were manually curated to
# include pseudogene loci — CDY/RBMY in GorGor, CDY/RBMY/TSPY in PanTro —
# to prevent AmpliCoNE-build from returning NaN for those families.

mkdir -p "${PANELS_DIR}" "${TMP_DIR}"

build_annotation() {
    local sp="$1"
    local ref accession gene_defs annotation_tab out_tab
    ref="$(sp_var "$sp" reference)"
    accession="$(sp_var "$sp" chry_accession)"
    gene_defs="${GENE_DEFS_DIR}/$(sp_var "$sp" gene_definitions)"
    annotation_tab="$(sp_var "$sp" annotation_tab)"
    out_tab="${PANELS_DIR}/${annotation_tab}"

    if [ ! -f "${out_tab}" ]; then
        echo "[Step 4] Building AmpliCoNE annotation for ${sp} ..."
        python3 "${AMPLICONE_TOOL}/AmpliCoNE-build.py" \
            -c "${accession}" \
            -i "${ref}" \
            -m "${ref}.bed" \
            -r "${RM_DIR}/${sp}_chrY_rm.out" \
            -t "${TRF_DIR}/${sp}_${accession}.fa.2.5.7.80.10.50.2000.bed" \
            -g "${gene_defs}" \
            -o "${out_tab}" \
            -d "${TMP_DIR}" \
            --debug
    else
        echo "[Step 4] Annotation exists for ${sp}, skipping."
    fi
}

for sp in $SPECIES; do
    build_annotation "$sp"
done

# ---------------------------------------------------------------------------
# STEP 5 — BATCH RUN AMPLICONE-COUNT (per individual)
# ---------------------------------------------------------------------------
# For each individual BAM in {alignments_root}/{species}/{individual}/merged.100.bam:
#   - auto-indexes the BAM if .bai is missing
#   - runs AmpliCoNE-count.py, writing output to results/{species}/{individual}/
# AmpliCoNE-count.py uses read depth in annotated windows and normalises
# against single-copy control regions to estimate copy numbers.

echo "[Step 5] Running AmpliCoNE batch estimation ..."
python3 "${SCRIPTS_DIR}/03_run_amplicone.py" \
    --config  "${CONFIG}" \
    --panels_dir  "${PANELS_DIR}" \
    --gene_defs_dir "${GENE_DEFS_DIR}" \
    --results_root  "${RESULTS_DIR}" \
    --amplicone_tool "${AMPLICONE_TOOL}"

# ---------------------------------------------------------------------------
# STEP 6 — COLLECT RESULTS INTO SUMMARY TABLE
# ---------------------------------------------------------------------------
# Reads OutputAmpliconic_Summary.txt from each individual's results directory
# and writes one TSV row per individual.
# Genes: TSPY  VCY  CDY  HSFY  RBMY  BPY2  DAZ

echo "[Step 6] Collecting results ..."
python3 "${SCRIPTS_DIR}/04_collect_results.py" \
    --config     "${CONFIG}" \
    --results_root "${RESULTS_DIR}" \
    --output_dir   "${RESULTS_DIR}"

echo "Pipeline complete. Results in ${RESULTS_DIR}/ampl_results.*.tsv"
