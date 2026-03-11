#!/usr/bin/env bash
# ==============================================================================
# run_hic_pipeline.sh
# Automated HiC Pipeline: Trimming → Alignment → Tag Directory → PCA →
#                          TAD/Loop Calling → Comparison → Annotation
#
# INPUT: data.txt  (tab-separated, no header)
#   Col 1: sample name   (e.g.  HiC-P0)
#   Col 2: comparison    (e.g.  1vs2   means row-1 sample vs row-2 sample)
#   Col 3: R1 fastq path
#   Col 4: R2 fastq path
#
# EXAMPLE data.txt:
#   HiC-P0        1vs2   /data/HiC_P0_R1.fastq.gz   /data/HiC_P0_R2.fastq.gz
#   HiC-7dpi      1vs2   /data/HiC_7dpi_R1.fastq.gz /data/HiC_7dpi_R2.fastq.gz
#   HiC-Uninj     2vs3   /data/HiC_Uninj_R1.fastq.gz /data/HiC_Uninj_R2.fastq.gz
#
# USAGE:
#   bash run_hic_pipeline.sh [OPTIONS]
#
# OPTIONS:
#   -d  data.txt          Input sample sheet         [required]
#   -g  mm10              Genome build                [default: mm10]
#   -b  /path/bowtie2idx  Bowtie2 index prefix        [required]
#   -r  gene_split.csv    Gene annotation for R script[default: gene_split.csv]
#   -o  results/          Output root directory       [default: HiC_Results]
#   -t  15                CPU threads                 [default: 15]
#   -j  /path/juicer.jar  Juicer tools jar            [optional]
#   -s  skip_align        Skip trimming+alignment if tag dirs exist [flag]
#   -h                    Show this help
#
# ==============================================================================

set -euo pipefail

# ---------- defaults ----------------------------------------------------------
GENOME="mm10"
BOWTIE2_IDX=""
DATA_FILE=""
GENE_FILE="gene_split.csv"
OUT_ROOT="HiC_Results"
THREADS=15
JUICER_JAR=""
SKIP_ALIGN=false
TAD_RES=25000
LOOP_RES=5000
PCA_RES=50000
PCA_WIN=100000
RESTRICTION_SITE="GATC"
TAD_SCRIPT="Tad_ann_CLI.R"

# ---------- colours -----------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[WARN] $*${NC}"; }
err()  { echo -e "${RED}[ERROR] $*${NC}" >&2; exit 1; }

# ---------- usage -------------------------------------------------------------
usage() {
cat <<EOF
Usage: bash run_hic_pipeline.sh -d data.txt -b /path/to/bowtie2_index [OPTIONS]

  -d  FILE    Sample sheet (tab-separated): name | comparison | R1.fastq.gz | R2.fastq.gz
  -g  STR     Genome build                  [default: mm10]
  -b  PREFIX  Bowtie2 index prefix          [REQUIRED]
  -r  FILE    Gene CSV for R annotation     [default: gene_split.csv]
  -o  DIR     Output root directory         [default: HiC_Results]
  -t  INT     CPU threads                   [default: 15]
  -j  JAR     Juicer tools jar path         [optional, enables .hic export]
  -s          Skip trimming + alignment (use existing tag directories)
  -h          Show this help message
EOF
exit 0
}

# ---------- argument parsing --------------------------------------------------
while getopts "d:g:b:r:o:t:j:sh" opt; do
  case $opt in
    d) DATA_FILE="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    b) BOWTIE2_IDX="$OPTARG" ;;
    r) GENE_FILE="$OPTARG" ;;
    o) OUT_ROOT="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    j) JUICER_JAR="$OPTARG" ;;
    s) SKIP_ALIGN=true ;;
    h) usage ;;
    *) usage ;;
  esac
done

# ---------- validate required -------------------------------------------------
[[ -z "$DATA_FILE" ]]   && err "-d data.txt is required"
[[ ! -f "$DATA_FILE" ]] && err "data.txt not found: $DATA_FILE"
[[ -z "$BOWTIE2_IDX" ]] && err "-b bowtie2 index prefix is required"

mkdir -p "$OUT_ROOT"
SAM_DIR="${OUT_ROOT}/SAM"
TAGDIR_ROOT="${OUT_ROOT}/TagDirs"
COMP_ROOT="${OUT_ROOT}/Comparisons"
mkdir -p "$SAM_DIR" "$TAGDIR_ROOT" "$COMP_ROOT"

LOG_FILE="${OUT_ROOT}/pipeline_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1

log "==========================================================="
log "  HiC Automated Pipeline"
log "  Genome   : $GENOME"
log "  Threads  : $THREADS"
log "  Data     : $DATA_FILE"
log "  Output   : $OUT_ROOT"
log "==========================================================="

# ---------- read sample sheet into arrays ------------------------------------
declare -a NAMES COMPARISONS R1S R2S
while IFS=$'\t' read -r name comp r1 r2 || [[ -n "$name" ]]; do
  [[ "$name" =~ ^#.*$ || -z "$name" ]] && continue   # skip comments/blanks
  NAMES+=("$name")
  COMPARISONS+=("$comp")
  R1S+=("$r1")
  R2S+=("$r2")
done < "$DATA_FILE"

NSAMP=${#NAMES[@]}
log "Found $NSAMP samples in $DATA_FILE"

# ==============================================================================
# STEP 1 – TRIM + ALIGN + MAKE TAG DIRECTORY (per sample)
# ==============================================================================
for (( i=0; i<NSAMP; i++ )); do
  NAME="${NAMES[$i]}"
  R1="${R1S[$i]}"
  R2="${R2S[$i]}"
  TAGDIR="${TAGDIR_ROOT}/${NAME}"

  if $SKIP_ALIGN && [[ -d "$TAGDIR" ]]; then
    warn "Skipping trim/align for $NAME (tag directory exists)"
    continue
  fi

  log "---------- Processing sample: $NAME ----------"

  # --- trim ---
  log "[$NAME] Trimming R1..."
  homerTools trim -3 "$RESTRICTION_SITE" -mis 0 -matchStart 20 -min 20 "$R1"

  log "[$NAME] Trimming R2..."
  homerTools trim -3 "$RESTRICTION_SITE" -mis 0 -matchStart 20 -min 20 "$R2"

  R1_TRIM="${R1}.trimmed"
  R2_TRIM="${R2}.trimmed"

  # --- align ---
  R1_SAM="${SAM_DIR}/${NAME}_R1.sam"
  R2_SAM="${SAM_DIR}/${NAME}_R2.sam"

  log "[$NAME] Aligning R1..."
  bowtie2 -p "$THREADS" -x "$BOWTIE2_IDX" -U "$R1_TRIM" > "$R1_SAM"

  log "[$NAME] Aligning R2..."
  bowtie2 -p "$THREADS" -x "$BOWTIE2_IDX" -U "$R2_TRIM" > "$R2_SAM"

  # --- tag directory ---
  log "[$NAME] Making tag directory..."
  makeTagDirectory "$TAGDIR" \
    "${R1_SAM},${R2_SAM}" \
    -tbp 1 \
    -genome "$GENOME" \
    -checkGC \
    -restrictionSite "$RESTRICTION_SITE" \
    -removePEbg \
    -removeSelfLigation \
    -removeRestrictionEnds \
    -removeSpikes 10000 5

  log "[$NAME] Tag directory created: $TAGDIR"
done

# ==============================================================================
# STEP 2 – HiC PCA (compartments) + compaction stats (per sample)
# ==============================================================================
for (( i=0; i<NSAMP; i++ )); do
  NAME="${NAMES[$i]}"
  TAGDIR="${TAGDIR_ROOT}/${NAME}"

  log "[$NAME] Running HiC PCA (res=${PCA_RES}, win=${PCA_WIN})..."
  runHiCpca.pl auto "$TAGDIR" \
    -res "$PCA_RES" \
    -window "$PCA_WIN" \
    -genome "$GENOME" \
    -cpu "$THREADS"

  log "[$NAME] Running analyzeHiC (compaction stats)..."
  analyzeHiC "$TAGDIR" \
    -res "$PCA_RES" \
    -window "$PCA_WIN" \
    -nomatrix \
    -compactionStats auto \
    -cpu "$THREADS" \
    -ifc auto

  # --- optional: export .hic file for Juicebox ---
  if [[ -n "$JUICER_JAR" && -f "$JUICER_JAR" ]]; then
    log "[$NAME] Exporting .hic file for Juicebox..."
    tagDir2hicFile.pl "$TAGDIR" \
      -juicer auto \
      -genome "$GENOME" \
      -p "$THREADS" \
      -juicerExe "java -Xms4G -Xmx200G -jar ${JUICER_JAR}"
  fi
done

# ==============================================================================
# STEP 3 – TAD + LOOP CALLING (per sample, individual)
# ==============================================================================
for (( i=0; i<NSAMP; i++ )); do
  NAME="${NAMES[$i]}"
  TAGDIR="${TAGDIR_ROOT}/${NAME}"

  log "[$NAME] Calling TADs (res=${TAD_RES})..."
  findTADsAndLoops.pl find "$TAGDIR" \
    -genome "$GENOME" \
    -res "$TAD_RES" \
    -cpu "$THREADS"

  log "[$NAME] Calling Loops (res=${LOOP_RES})..."
  findTADsAndLoops.pl find "$TAGDIR" \
    -genome "$GENOME" \
    -res "$LOOP_RES" \
    -cpu "$THREADS"
done

# ==============================================================================
# STEP 4 – PAIRWISE COMPARISONS (driven by comparison column in data.txt)
# ==============================================================================

# Build a lookup: comparison_id -> list of sample indices
declare -A COMP_MAP
for (( i=0; i<NSAMP; i++ )); do
  COMP="${COMPARISONS[$i]}"
  if [[ -z "${COMP_MAP[$COMP]+_}" ]]; then
    COMP_MAP["$COMP"]="$i"
  else
    COMP_MAP["$COMP"]="${COMP_MAP[$COMP]} $i"
  fi
done

for COMP_ID in "${!COMP_MAP[@]}"; do
  IDX_LIST=(${COMP_MAP[$COMP_ID]})
  NIDX=${#IDX_LIST[@]}

  if [[ $NIDX -lt 2 ]]; then
    warn "Comparison '$COMP_ID' has only 1 sample – skipping"
    continue
  fi

  # Do all pairwise combinations within this comparison group
  for (( a=0; a<NIDX-1; a++ )); do
    for (( b=a+1; b<NIDX; b++ )); do
      IDX_A=${IDX_LIST[$a]}
      IDX_B=${IDX_LIST[$b]}
      NAME_A="${NAMES[$IDX_A]}"
      NAME_B="${NAMES[$IDX_B]}"
      TAGDIR_A="${TAGDIR_ROOT}/${NAME_A}"
      TAGDIR_B="${TAGDIR_ROOT}/${NAME_B}"
      PAIR_DIR="${COMP_ROOT}/${NAME_A}_vs_${NAME_B}"

      log "========== Comparison: $NAME_A  vs  $NAME_B =========="
      mkdir -p "$PAIR_DIR"

      # ---- merge TAD beds ----
      TAD_A="${TAGDIR_A}/${NAME_A}.tad.2D.bed"
      TAD_B="${TAGDIR_B}/${NAME_B}.tad.2D.bed"
      MERGED_TAD="${PAIR_DIR}/merged_TADs.bed"

      log "Merging TAD beds..."
      merge2Dbed.pl "$TAD_A" "$TAD_B" -tad > "$MERGED_TAD"

      # ---- merge Loop beds ----
      LOOP_A="${TAGDIR_A}/${NAME_A}.loop.2D.bed"
      LOOP_B="${TAGDIR_B}/${NAME_B}.loop.2D.bed"
      MERGED_LOOP="${PAIR_DIR}/merged_Loops.bed"

      log "Merging Loop beds..."
      merge2Dbed.pl "$LOOP_A" "$LOOP_B" -loop > "$MERGED_LOOP"

      # ---- score TADs/Loops across both tag dirs ----
      SCORE_PREFIX="${PAIR_DIR}/TADLoop_score"
      log "Scoring TADs and Loops..."
      findTADsAndLoops.pl score \
        -tad "$MERGED_TAD" \
        -loop "$MERGED_LOOP" \
        -o "$SCORE_PREFIX" \
        -d "$TAGDIR_A" "$TAGDIR_B" \
        -cpu "$THREADS" \
        -res "$TAD_RES" \
        -window "$((TAD_RES * 5))"

      # ---- compartment comparison (PCA bedGraphs) ----
      PC1_A="${TAGDIR_A}/${NAME_A}.${PCA_RES}x${PCA_WIN}.PC1.txt"
      PC1_BG_A="${TAGDIR_A}/${NAME_A}.${PCA_RES}x${PCA_WIN}.PC1.bedGraph"
      PC1_BG_B="${TAGDIR_B}/${NAME_B}.${PCA_RES}x${PCA_WIN}.PC1.bedGraph"
      COMP_OUT="${PAIR_DIR}/Chromatin_compartment_${NAME_A}_vs_${NAME_B}.txt"

      if [[ -f "$PC1_A" && -f "$PC1_BG_A" && -f "$PC1_BG_B" ]]; then
        log "Comparing chromatin compartments..."
        annotatePeaks.pl "$PC1_A" "$GENOME" \
          -noblanks \
          -bedGraph "$PC1_BG_A" "$PC1_BG_B" \
          > "$COMP_OUT"
      else
        warn "PC1 files not found for $NAME_A – skipping compartment comparison"
      fi

      # ---- R annotation (TAD_ann_CLI.R) ----
      if [[ -f "$TAD_SCRIPT" ]]; then
        log "Running TAD/Loop R annotation..."
        Rscript "$TAD_SCRIPT" \
          --mode both \
          --n1  "$NAME_A" \
          --n2  "$NAME_B" \
          --t1  "$TAD_A" \
          --t2  "$TAD_B" \
          --l1  "$LOOP_A" \
          --l2  "$LOOP_B" \
          --gene "$GENE_FILE" \
          --out  "${PAIR_DIR}/R_annotation"
      else
        warn "Tad_ann_CLI.R not found at '$TAD_SCRIPT' – skipping R annotation"
      fi

      log "Comparison $NAME_A vs $NAME_B complete → $PAIR_DIR"
    done
  done
done

# ==============================================================================
log "==========================================================="
log "  ALL DONE.  Results in: $OUT_ROOT"
log "  Log saved to: $LOG_FILE"
log "==========================================================="
