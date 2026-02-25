#!/usr/bin/env bash
# ============================================================
#  Purge_dups → Chromeister → RagTag pipeline (get_seqs-enabled)
#
#  Requirements implemented:
#   - Chromeister outputs written in $OUTDIR
#   - Initial step: create a length-sorted copy of input assembly (long→short)
#   - If REF has >20 sequences: order PURGED ragtag scaffolds by length (desc)
#     else: reorder PURGED scaffolds by REF order, KEEPING unmatched at end
#   - FINAL_$SAMPLE/: rename sequences to ${SAMPLE}_scaff1..N (+ mapping TSV)
#
#  Core behavior:
#   - split_fa + minimap2 -xasm5 -DP (like your Snakemake rule)
#   - purge_dups -2 -T cutoffs (no PB.base.cov)
#   - get_seqs -e -c with BED cleanup + oneline FASTA
# ============================================================
# Usage:
#   bash run_purge_ragtag.sh <assembly.fasta[.gz]> <reference.fasta[.gz]> <outdir> <sample_prefix> [threads]
# ============================================================

set -euo pipefail

# -------------------- USER-EDITABLE CUT-OFFS --------------------
LOW=5
MID=35
HIGH=90

# -------------------- Parameters --------------------
ASM=${1:?Provide assembly FASTA (.fa/.fasta) optionally gzipped}
REF=${2:?Provide reference FASTA (.fna/.fa) optionally gzipped}
OUTDIR=${3:?Provide output directory path}
SAMPLE=${4:?Provide sample/output prefix}
THREADS=${5:-32}

# -------------------- Resolve input paths to absolute (before cd) --------------------
abs_path() {
  local p="$1"
  if command -v readlink >/dev/null 2>&1; then
    readlink -f "$p"
  else
    python - <<'PY' "$p"
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
  fi
}

ASM="$(abs_path "${ASM}")"
REF="$(abs_path "${REF}")"
OUTDIR="$(abs_path "${OUTDIR}")"

# -------------------- Setup output dir --------------------
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# -------------------- Logging setup --------------------
LOG="${SAMPLE}_pipeline.log"
exec > >(tee -a "${LOG}") 2>&1

echo "============================================================"
echo ">>> Purge_dups–Chromeister–RagTag pipeline (get_seqs-enabled)"
echo "Sample/prefix: ${SAMPLE}"
echo "Assembly: ${ASM}"
echo "Reference: ${REF}"
echo "Output dir: ${OUTDIR}"
echo "Threads: ${THREADS}"
echo "Cutoffs: low=${LOW} mid=${MID} high=${HIGH}"
echo "Start time: $(date)"
echo "============================================================"
echo

# -------------------- Check required tools --------------------
req_tools=(split_fa minimap2 purge_dups get_seqs samtools ragtag.py Rscript CHROMEISTER gunzip gzip tee awk grep head wc zcat sort cut sed dirname zgrep comm mktemp)
for t in "${req_tools[@]}"; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: Required tool not found in PATH: $t"; exit 1; }
done

# -------------------- Normalize inputs (handle gz) --------------------
ASM_USED="${ASM}"
REF_USED="${REF}"
ASM_TMP=""
REF_TMP=""

if [[ "${ASM}" == *.gz ]]; then
  echo ">>> [$(date)] Assembly is gzipped; creating temporary copy..."
  ASM_TMP="${SAMPLE}_assembly_tmp.fasta"
  gunzip -c "${ASM}" > "${ASM_TMP}"
  ASM_USED="${ASM_TMP}"
fi

if [[ "${REF}" == *.gz ]]; then
  echo ">>> [$(date)] Reference is gzipped; creating temporary copy..."
  REF_TMP="${SAMPLE}_ref_tmp.fasta"
  gunzip -c "${REF}" > "${REF_TMP}"
  REF_USED="${REF_TMP}"
fi

cleanup() {
  if [[ -n "${ASM_TMP}" && -f "${ASM_TMP}" ]]; then rm -f "${ASM_TMP}"; fi
  if [[ -n "${REF_TMP}" && -f "${REF_TMP}" ]]; then rm -f "${REF_TMP}"; fi
}
trap cleanup EXIT

# -------------------- Count REF sequences (for conditional ordering) --------------------
if [[ "${REF}" == *.gz ]]; then
  REF_NSEQ=$(zgrep -c '^>' "${REF}" || true)
else
  REF_NSEQ=$(grep -c '^>' "${REF}" || true)
fi
echo ">>> [$(date)] Reference sequences: ${REF_NSEQ}"

# -------------------- Sanity check input FASTA --------------------
if [[ ! -f "${ASM_USED}" ]]; then
  echo "ERROR: Assembly file not found: ${ASM_USED}"
  exit 1
fi
if ! grep -q '^>' "${ASM_USED}"; then
  echo "ERROR: Assembly file does not look like FASTA (no '>' headers): ${ASM_USED}"
  head -n 5 "${ASM_USED}" || true
  exit 1
fi

# ============================================================
# Helpers
# ============================================================
fasta_to_oneline() {
  local in="$1" out="$2"
  awk '
    BEGIN{seq="";hdr=""}
    /^>/{
      if(hdr!=""){print hdr; print seq}
      hdr=$0; seq=""; next
    }
    {
      gsub(/\r$/, "", $0);
      if($0!="") seq=seq $0
    }
    END{ if(hdr!=""){print hdr; print seq} }
  ' "$in" > "$out"
}

sort_fasta_by_len_desc() {
  local in="$1" out="$2"
  local tmp_oneline
  tmp_oneline="$(mktemp -p . "${SAMPLE}.tmp.oneline.XXXXXX.fa")"
  fasta_to_oneline "$in" "$tmp_oneline"
  awk '
    NR%2==1 {hdr=$0; next}
    NR%2==0 {print hdr "\t" length($0) "\t" $0}
  ' "$tmp_oneline" \
  | sort -k2,2nr \
  | awk '{print $1; print $3}' > "$out"
  rm -f "$tmp_oneline"
}

reorder_scaffolds_by_ref_keep_unmatched() {
  local ref="$1" agp="$2" fasta="$3" out_sorted="$4"
  local dir
  dir="$(dirname "${out_sorted}")"
  mkdir -p "${dir}"

  local reforder="${dir}/ref.order.txt"
  local matched_order="${dir}/scaffold.order.matched.txt"
  local fasta_all="${dir}/scaffold.order.all.txt"
  local unmatched_set="${dir}/scaffold.order.unmatched.set.txt"
  local final="${dir}/scaffold.order.final.txt"

  echo ">>> [$(date)] Reordering scaffolds to match reference order (keeping unmatched at end)..."

  if [[ "${ref}" == *.gz ]]; then
    zgrep '^>' "${ref}" | cut -d' ' -f1 | sed 's/^>//' > "${reforder}"
  else
    grep '^>' "${ref}" | cut -d' ' -f1 | sed 's/^>//' > "${reforder}"
  fi

  awk '
    NR==FNR { ord[$1]=NR; next }
    /^#/ { next }
    $1!=p {
      x=$1
      sub(/_RagTag$/, "", x)
      if (x in ord) print ord[x] "\t" $1
      p=$1
    }
  ' "${reforder}" "${agp}" \
  | sort -k1,1n | cut -f2 > "${matched_order}"

  grep '^>' "${fasta}" | sed 's/^>//' | cut -d' ' -f1 > "${fasta_all}"

  comm -23 <(sort "${fasta_all}") <(sort "${matched_order}") > "${unmatched_set}"

  cp -f "${matched_order}" "${final}"
  awk 'NR==FNR{u[$1]=1; next} ($1 in u){print $1}' "${unmatched_set}" "${fasta_all}" >> "${final}"

  if [[ ! -s "${final}" ]]; then
    echo "WARNING: Could not build final scaffold order. Copying FASTA unchanged."
    cp -f "${fasta}" "${out_sorted}"
    samtools faidx "${out_sorted}" || true
    return 0
  fi

  samtools faidx "${fasta}"
  # shellcheck disable=SC2046
  samtools faidx "${fasta}" $(cat "${final}") > "${out_sorted}"
  samtools faidx "${out_sorted}"

  echo ">>> [$(date)] Reordering summary:"
  echo ">>> [$(date)]   matched:   $(wc -l < "${matched_order}" | tr -d " ")"
  echo ">>> [$(date)]   unmatched: $(wc -l < "${unmatched_set}" | tr -d " ")"
  echo ">>> [$(date)]   total:     $(wc -l < "${final}" | tr -d " ")"
}

rename_fasta_to_sample_scaffs() {
  local in="$1" out="$2" map="$3" prefix="$4"
  local tmp_oneline
  tmp_oneline="$(mktemp -p . "${SAMPLE}.tmp.rename.oneline.XXXXXX.fa")"
  fasta_to_oneline "$in" "$tmp_oneline"

  : > "$out"
  : > "$map"

  awk -v P="${prefix}" -v OUT="${out}" -v MAP="${map}" '
    function clean(h,   x){
      x=h; sub(/^>/,"",x); sub(/[ \t].*$/,"",x); return x
    }
    NR%2==1 { old=clean($0); next }
    NR%2==0 {
      n++
      new=P "_scaff" n
      print old "\t" new >> MAP
      print ">" new >> OUT
      print $0 >> OUT
    }
  ' "$tmp_oneline"

  rm -f "$tmp_oneline"
}

# ============================================================
# 0) Length-sorted copy of input assembly (longest → shortest)
# ============================================================
echo ">>> [$(date)] Creating length-sorted copy of input assembly (longest → shortest)..."
sort_fasta_by_len_desc "${ASM_USED}" "${SAMPLE}.assembly.len_sorted.fasta"
echo ">>> [$(date)] Length-sorted assembly: ${SAMPLE}.assembly.len_sorted.fasta"

# ============================================================
# 1) cutoffs + split_fa
# ============================================================
echo ">>> [$(date)] Writing cutoffs..."
cat > cutoffs <<EOF
low=${LOW}
mid=${MID}
high=${HIGH}
EOF
cat cutoffs

echo ">>> [$(date)] Running split_fa..."
split_fa "${ASM_USED}" > "${SAMPLE}.split.fasta"
nseq=$(grep -c '^>' "${SAMPLE}.split.fasta" || true)
if [[ "${nseq}" -eq 0 ]]; then
  echo "ERROR: split_fa produced empty FASTA: ${SAMPLE}.split.fasta"
  exit 1
fi
echo ">>> [$(date)] Split FASTA contains ${nseq} sequences."

# ============================================================
# 2) minimap2 self-alignment
# ============================================================
echo ">>> [$(date)] Running minimap2 self-alignment..."
minimap2 -xasm5 -DP -t "${THREADS}" "${SAMPLE}.split.fasta" "${SAMPLE}.split.fasta" \
  | gzip -c > "${SAMPLE}.split.self.paf.gz"
paf_n=$(zcat "${SAMPLE}.split.self.paf.gz" | wc -l || true)
echo ">>> [$(date)] PAF lines: ${paf_n}"

# ============================================================
# 3) purge_dups (no PB.base.cov)
# ============================================================
echo ">>> [$(date)] Running purge_dups..."
purge_dups -2 -T cutoffs "${SAMPLE}.split.self.paf.gz" > "dups_${SAMPLE}.bed"
bed_n=$(wc -l < "dups_${SAMPLE}.bed" | tr -d ' ')
echo ">>> [$(date)] dups BED lines: ${bed_n}"
if [[ "${bed_n}" -gt 0 ]]; then
  echo ">>> [$(date)] dups BED (first 5 lines):"
  head -n 5 "dups_${SAMPLE}.bed"
fi

# ============================================================
# 4) get_seqs with safe preprocessing
# ============================================================
echo ">>> [$(date)] Preparing safe inputs for get_seqs..."

awk 'BEGIN{OFS="\t"} $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $3>$2{print}' \
  "dups_${SAMPLE}.bed" > "dups_${SAMPLE}.clean.bed"
clean_bed_n=$(wc -l < "dups_${SAMPLE}.clean.bed" | tr -d ' ')
echo ">>> [$(date)] Clean BED lines: ${clean_bed_n}"

fasta_to_oneline "${ASM_USED}" "${SAMPLE}.asm.clean.oneline.fa"
samtools faidx "${SAMPLE}.asm.clean.oneline.fa" >/dev/null 2>&1 || true

if [[ "${clean_bed_n}" -eq 0 ]]; then
  echo "[INFO] ${SAMPLE}: BED empty after cleaning; copying original assembly to purged.fa"
  cp -f "${ASM_USED}" purged.fa
  echo -e ">hap_empty\nN" > hap.fa
else
  echo ">>> [$(date)] Running get_seqs (conservative: -e -c)..."
  set +e
  get_seqs -e -c "dups_${SAMPLE}.clean.bed" "${SAMPLE}.asm.clean.oneline.fa"
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "[WARN] ${SAMPLE}: get_seqs failed (exit ${rc}). Falling back to original assembly."
    cp -f "${ASM_USED}" purged.fa
    echo -e ">hap_empty\nN" > hap.fa
  elif [[ ! -f purged.fa || ! -s purged.fa ]]; then
    echo "[WARN] ${SAMPLE}: purged.fa missing/empty after get_seqs. Falling back."
    cp -f "${ASM_USED}" purged.fa
    echo -e ">hap_empty\nN" > hap.fa
  else
    if [[ ! -f hap.fa || ! -s hap.fa ]]; then
      echo -e ">hap_empty\nN" > hap.fa
    fi
  fi
fi

samtools faidx purged.fa >/dev/null 2>&1 || true
samtools faidx hap.fa   >/dev/null 2>&1 || true

purged_n=$(grep -c '^>' purged.fa 2>/dev/null || true)
hap_n=$(grep -c '^>' hap.fa 2>/dev/null || true)
echo ">>> [$(date)] purged.fa contigs: ${purged_n}"
echo ">>> [$(date)] hap.fa contigs:    ${hap_n}"

# ============================================================
# 5) Chromeister (pre-scaffold)  ---> outputs in $OUTDIR
# ============================================================
echo ">>> [$(date)] Running Chromeister & score calculation (pre-scaffold)..."
CHROMEISTER -query "${REF_USED}" -db "${ASM_USED}" -out "${SAMPLE}.assembly.mat" -dimension 2000
Rscript ~/anaconda3/envs/chromeister/bin/compute_score.R "${SAMPLE}.assembly.mat" 2000

# ============================================================
# 6) RagTag scaffolding
# ============================================================
echo ">>> [$(date)] Running RagTag scaffolding (purged)..."
ragtag.py scaffold -t "${THREADS}" -o "purged_${SAMPLE}_ragtag" "${REF_USED}" purged.fa

HAP_PLACEHOLDER=0
if grep -q '^>hap_empty' hap.fa; then
  echo ">>> [$(date)] hap.fa is placeholder; skipping hap ragtag."
  HAP_PLACEHOLDER=1
else
  echo ">>> [$(date)] Running RagTag scaffolding (hap)..."
  ragtag.py scaffold -t "${THREADS}" -o "hap_${SAMPLE}_ragtag" "${REF_USED}" hap.fa
fi

# ============================================================
# 7) Ordering of PURGED scaffolds:
#    - if REF_NSEQ > 20 : order by scaffold length (desc)
#    - else             : reorder by reference order (keeping unmatched)
# ============================================================
PURGED_AGP="purged_${SAMPLE}_ragtag/ragtag.scaffold.agp"
PURGED_FASTA="purged_${SAMPLE}_ragtag/ragtag.scaffold.fasta"

if [[ ! -s "${PURGED_FASTA}" ]]; then
  echo "ERROR: Purged RagTag scaffold FASTA not found: ${PURGED_FASTA}"
  exit 1
fi

if [[ "${REF_NSEQ}" -gt 20 ]]; then
  PURGED_ORDERED="purged_${SAMPLE}_ragtag/ragtag.scaffold.len_sorted.fasta"
  echo ">>> [$(date)] REF has >80 sequences; ordering PURGED scaffolds by size (desc)."
  sort_fasta_by_len_desc "${PURGED_FASTA}" "${PURGED_ORDERED}"
  samtools faidx "${PURGED_ORDERED}" >/dev/null 2>&1 || true
else
  PURGED_ORDERED="purged_${SAMPLE}_ragtag/ragtag.scaffold.ref_sorted.fasta"
  if [[ -s "${PURGED_AGP}" ]]; then
    reorder_scaffolds_by_ref_keep_unmatched "${REF}" "${PURGED_AGP}" "${PURGED_FASTA}" "${PURGED_ORDERED}"
  else
    echo "WARNING: ${PURGED_AGP} missing; leaving PURGED scaffolds as-is."
    PURGED_ORDERED="${PURGED_FASTA}"
  fi
fi

# HAP scaffolds: reorder by ref (keep unmatched) if produced
HAP_FASTA="hap_${SAMPLE}_ragtag/ragtag.scaffold.fasta"
HAP_AGP="hap_${SAMPLE}_ragtag/ragtag.scaffold.agp"
HAP_ORDERED=""
if [[ "${HAP_PLACEHOLDER}" -eq 0 && -s "${HAP_FASTA}" ]]; then
  HAP_ORDERED="hap_${SAMPLE}_ragtag/ragtag.scaffold.ref_sorted.fasta"
  if [[ -s "${HAP_AGP}" ]]; then
    reorder_scaffolds_by_ref_keep_unmatched "${REF}" "${HAP_AGP}" "${HAP_FASTA}" "${HAP_ORDERED}"
  else
    echo "WARNING: ${HAP_AGP} missing; leaving HAP scaffolds as-is."
    HAP_ORDERED="${HAP_FASTA}"
  fi
fi

# ============================================================
# 8) FINAL folder: rename to SAMPLE_scaff1..N
# ============================================================
FINAL_DIR="final_${SAMPLE}"
mkdir -p "${FINAL_DIR}"

PURGED_FINAL="${FINAL_DIR}/${SAMPLE}.purged.final.fasta"
PURGED_MAP="${FINAL_DIR}/${SAMPLE}.purged.name_map.tsv"
echo ">>> [$(date)] Writing FINAL renamed PURGED scaffolds: ${PURGED_FINAL}"
rename_fasta_to_sample_scaffs "${PURGED_ORDERED}" "${PURGED_FINAL}" "${PURGED_MAP}" "${SAMPLE}"
samtools faidx "${PURGED_FINAL}" >/dev/null 2>&1 || true

if [[ -n "${HAP_ORDERED}" && -s "${HAP_ORDERED}" ]]; then
  HAP_FINAL="${FINAL_DIR}/${SAMPLE}.hap.final.fasta"
  HAP_MAP="${FINAL_DIR}/${SAMPLE}.hap.name_map.tsv"
  echo ">>> [$(date)] Writing FINAL renamed HAP scaffolds: ${HAP_FINAL}"
  rename_fasta_to_sample_scaffs "${HAP_ORDERED}" "${HAP_FINAL}" "${HAP_MAP}" "${SAMPLE}"
  samtools faidx "${HAP_FINAL}" >/dev/null 2>&1 || true
else
  echo ">>> [$(date)] No hap scaffolds to finalize (placeholder or missing)."
fi

# ============================================================
# 9) Post-scaffolding Chromeister scores  ---> outputs in $OUTDIR
# ============================================================
echo ">>> [$(date)] Scoring FINAL scaffolds (Chromeister)..."
CHROMEISTER -query "${REF_USED}" -db "${PURGED_FINAL}" -out "${SAMPLE}.purged.final.mat" -dimension 2000
Rscript ~/anaconda3/envs/chromeister/bin/compute_score.R "${SAMPLE}.purged.final.mat" 2000

if [[ -f "${FINAL_DIR}/${SAMPLE}.hap.final.fasta" && -s "${FINAL_DIR}/${SAMPLE}.hap.final.fasta" ]]; then
  CHROMEISTER -query "${REF_USED}" -db "${FINAL_DIR}/${SAMPLE}.hap.final.fasta" -out "${SAMPLE}.hap.final.mat" -dimension 2000
  Rscript ~/anaconda3/envs/chromeister/bin/compute_score.R "${SAMPLE}.hap.final.mat" 2000
fi

echo
echo "============================================================"
echo "Pipeline completed successfully for ${SAMPLE}"
echo "End time: $(date)"
echo "Log saved to: ${OUTDIR}/${LOG}"
echo "FINAL outputs in: ${OUTDIR}/${FINAL_DIR}/"
echo "Chromeister outputs in: ${OUTDIR}/"
echo "============================================================"
