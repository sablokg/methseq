#!/usr/bin/env bash
###############################################################################
# methylation_pipeline.sh Gaurav Sablok gsablok@proton.me
#
# Complete DNA methylation analysis pipeline (WGBS / RRBS)
#
# Steps:
#   1. Raw read QC                (FastQC)
#   2. Adapter/quality trimming   (Trim Galore, bisulfite-aware)
#   3. Post-trim QC               (FastQC)
#   4. Genome preparation         (Bismark, run once, cached)
#   5. Bisulfite alignment        (Bismark / Bowtie2)
#   6. Deduplication              (deduplicate_bismark)
#   7. Methylation extraction     (bismark_methylation_extractor)
#   8. Genome-wide report         (bismark2report / bismark2summary)
#   9. Coverage -> CpG report     (coverage2cytosine)
#  10. Differential methylation   (methylKit, via Rscript)
#
# Requirements (install via conda/mamba is easiest):
#   fastqc, trim_galore, cutadapt, bismark, bowtie2 (or hisat2), samtools,
#   multiqc, R (methylKit, DSS, GenomicRanges)
#
# Usage:
#   ./methylation_pipeline.sh -c config.sh
#
# Or edit the CONFIG block below directly and run with no arguments.
###############################################################################

set -euo pipefail
trap 'echo "[ERROR] Pipeline failed at line $LINENO. See log: $LOGFILE" >&2' ERR

###############################################################################
# 0. CONFIG — edit these, or pass "-c your_config.sh" to override
###############################################################################

# --- Input / output ---
RAW_DIR="${RAW_DIR:-./raw_fastq}"          # folder with *_R1.fastq.gz / *_R2.fastq.gz
GENOME_DIR="${GENOME_DIR:-./genome}"       # folder with reference FASTA (one file)
OUT_DIR="${OUT_DIR:-./methylation_results}"
THREADS="${THREADS:-8}"

# --- Library type ---
# "paired"  -> expects <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
# "single"  -> expects <sample>.fastq.gz
LIBRARY_TYPE="${LIBRARY_TYPE:-paired}"

# --- Assay type ---
# "wgbs" -> whole genome bisulfite seq
# "rrbs" -> reduced representation (adds --rrbs trimming mode)
ASSAY_TYPE="${ASSAY_TYPE:-wgbs}"

# --- Differential methylation (optional, group-based) ---
RUN_DMR="${RUN_DMR:-false}"                # true/false
GROUP_A_SAMPLES="${GROUP_A_SAMPLES:-}"     # space-separated sample names
GROUP_B_SAMPLES="${GROUP_B_SAMPLES:-}"
DMR_CONTEXT="${DMR_CONTEXT:-CpG}"          # CpG, CHG, CHH

###############################################################################
# Parse optional external config file
###############################################################################
while getopts "c:" opt; do
  case $opt in
    c) source "$OPTARG" ;;
    *) echo "Usage: $0 [-c config.sh]"; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR"
LOGDIR="$OUT_DIR/logs"
mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "=============================================================="
echo " Methylation Pipeline started: $(date)"
echo " RAW_DIR=$RAW_DIR  GENOME_DIR=$GENOME_DIR  OUT_DIR=$OUT_DIR"
echo " LIBRARY_TYPE=$LIBRARY_TYPE  ASSAY_TYPE=$ASSAY_TYPE  THREADS=$THREADS"
echo "=============================================================="

###############################################################################
# Helper: check required tools exist
###############################################################################
require_tool () {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Required tool not found in PATH: $1" >&2; exit 1; }
}
for t in fastqc trim_galore bismark bismark_genome_preparation \
         deduplicate_bismark bismark_methylation_extractor \
         coverage2cytosine bismark2report bismark2summary bowtie2 samtools multiqc; do
  require_tool "$t"
done

###############################################################################
# 1. Discover samples
###############################################################################
mkdir -p "$OUT_DIR"/{01_fastqc_raw,02_trimmed,03_fastqc_trimmed,04_bam,05_dedup,06_methyl_extract,07_reports,08_cytosine_reports}

if [[ "$LIBRARY_TYPE" == "paired" ]]; then
  mapfile -t SAMPLES < <(find "$RAW_DIR" -maxdepth 1 -name "*_R1.fastq.gz" -printf "%f\n" \
                          | sed 's/_R1\.fastq\.gz$//' | sort)
else
  mapfile -t SAMPLES < <(find "$RAW_DIR" -maxdepth 1 -name "*.fastq.gz" -printf "%f\n" \
                          | sed 's/\.fastq\.gz$//' | sort)
fi

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "[ERROR] No samples found in $RAW_DIR" >&2
  exit 1
fi
echo "Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

###############################################################################
# 2. Raw QC
###############################################################################
echo ">>> Step 1/9: FastQC on raw reads"
fastqc -t "$THREADS" -o "$OUT_DIR/01_fastqc_raw" "$RAW_DIR"/*.fastq.gz

###############################################################################
# 3. Trimming (bisulfite-aware)
###############################################################################
echo ">>> Step 2/9: Trim Galore"
TRIM_EXTRA=""
[[ "$ASSAY_TYPE" == "rrbs" ]] && TRIM_EXTRA="--rrbs"

for s in "${SAMPLES[@]}"; do
  echo "  Trimming sample: $s"
  if [[ "$LIBRARY_TYPE" == "paired" ]]; then
    trim_galore --paired --cores "$THREADS" $TRIM_EXTRA \
      --output_dir "$OUT_DIR/02_trimmed" \
      "$RAW_DIR/${s}_R1.fastq.gz" "$RAW_DIR/${s}_R2.fastq.gz"
  else
    trim_galore --cores "$THREADS" $TRIM_EXTRA \
      --output_dir "$OUT_DIR/02_trimmed" \
      "$RAW_DIR/${s}.fastq.gz"
  fi
done

###############################################################################
# 4. Post-trim QC
###############################################################################
echo ">>> Step 3/9: FastQC on trimmed reads"
fastqc -t "$THREADS" -o "$OUT_DIR/03_fastqc_trimmed" "$OUT_DIR/02_trimmed"/*.fq.gz

###############################################################################
# 5. Genome preparation (cached — only runs if not already converted)
###############################################################################
echo ">>> Step 4/9: Bismark genome preparation"
if [[ ! -d "$GENOME_DIR/Bisulfite_Genome" ]]; then
  bismark_genome_preparation --bowtie2 "$GENOME_DIR"
else
  echo "  Bisulfite genome index already exists, skipping."
fi

###############################################################################
# 6. Bisulfite alignment
###############################################################################
echo ">>> Step 5/9: Bismark alignment"
for s in "${SAMPLES[@]}"; do
  echo "  Aligning sample: $s"
  if [[ "$LIBRARY_TYPE" == "paired" ]]; then
    R1="$OUT_DIR/02_trimmed/${s}_R1_val_1.fq.gz"
    R2="$OUT_DIR/02_trimmed/${s}_R2_val_2.fq.gz"
    bismark --bowtie2 --multicore $(( THREADS / 2 > 0 ? THREADS / 2 : 1 )) \
      --genome "$GENOME_DIR" -1 "$R1" -2 "$R2" \
      -o "$OUT_DIR/04_bam" --temp_dir "$OUT_DIR/04_bam/tmp_${s}"
  else
    R1="$OUT_DIR/02_trimmed/${s}_trimmed.fq.gz"
    bismark --bowtie2 --multicore $(( THREADS / 2 > 0 ? THREADS / 2 : 1 )) \
      --genome "$GENOME_DIR" "$R1" \
      -o "$OUT_DIR/04_bam" --temp_dir "$OUT_DIR/04_bam/tmp_${s}"
  fi
done

###############################################################################
# 7. Deduplication (skip for RRBS — amplicon-like, dedup not recommended)
###############################################################################
echo ">>> Step 6/9: Deduplication"
for s in "${SAMPLES[@]}"; do
  BAM=$(find "$OUT_DIR/04_bam" -name "${s}*bismark_bt2*.bam" ! -name "*dedup*" | head -n1)
  if [[ "$ASSAY_TYPE" == "rrbs" ]]; then
    echo "  RRBS detected — skipping deduplication for $s (copying BAM through)."
    cp "$BAM" "$OUT_DIR/05_dedup/$(basename "$BAM")"
  else
    echo "  Deduplicating: $BAM"
    if [[ "$LIBRARY_TYPE" == "paired" ]]; then
      deduplicate_bismark --paired --bam --output_dir "$OUT_DIR/05_dedup" "$BAM"
    else
      deduplicate_bismark --single --bam --output_dir "$OUT_DIR/05_dedup" "$BAM"
    fi
  fi
done

###############################################################################
# 8. Methylation extraction
###############################################################################
echo ">>> Step 7/9: Methylation extraction"
PAIR_FLAG="--single-end"
[[ "$LIBRARY_TYPE" == "paired" ]] && PAIR_FLAG="--paired-end"

for f in "$OUT_DIR"/05_dedup/*.bam; do
  echo "  Extracting: $f"
  bismark_methylation_extractor $PAIR_FLAG --comprehensive --multicore "$THREADS" \
    --cytosine_report --CX --genome_folder "$GENOME_DIR" \
    -o "$OUT_DIR/06_methyl_extract" "$f"
done

###############################################################################
# 9. Reports (per-sample + combined summary)
###############################################################################
echo ">>> Step 8/9: Bismark reports"
cd "$OUT_DIR/04_bam" && bismark2report --dir "$OUT_DIR/07_reports" || true
cd - >/dev/null
bismark2summary "$OUT_DIR"/04_bam/*.bam -o "$OUT_DIR/07_reports/bismark_summary" || true

echo ">>> MultiQC aggregate report"
multiqc "$OUT_DIR" -o "$OUT_DIR/07_reports" -n multiqc_report

###############################################################################
# 10. Optional: differential methylation with methylKit (CpG-context DMRs)
###############################################################################
if [[ "$RUN_DMR" == "true" ]]; then
  echo ">>> Step 9/9: Differential methylation (methylKit)"
  if [[ -z "$GROUP_A_SAMPLES" || -z "$GROUP_B_SAMPLES" ]]; then
    echo "[WARN] RUN_DMR=true but GROUP_A_SAMPLES/GROUP_B_SAMPLES not set. Skipping DMR step."
  else
    Rscript "$(dirname "$0")/run_methylkit.R" \
      --cov_dir "$OUT_DIR/06_methyl_extract" \
      --group_a "$GROUP_A_SAMPLES" \
      --group_b "$GROUP_B_SAMPLES" \
      --context "$DMR_CONTEXT" \
      --out_dir "$OUT_DIR/09_dmr"
  fi
else
  echo ">>> Step 9/9: Skipped (RUN_DMR=false)"
fi

echo "=============================================================="
echo " Pipeline finished: $(date)"
echo " Results in: $OUT_DIR"
echo " Full log:   $LOGFILE"
echo "=============================================================="
