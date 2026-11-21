#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "$0")/common.sh"

require_cmd roary
mapfile -t GFFS < <(find "${OUTPUT_DIR}/prokka" -mindepth 2 -maxdepth 2 -type f -name '*.gff' -print | sort)
((${#GFFS[@]} >= 2)) || die "Need Prokka GFFs for at least two genomes. Run 02_annotation_quality.sh first."

ROARY_OUT="${OUTPUT_DIR}/roary"
mkdir -p "$ROARY_OUT"
log "Roary pan-genome analysis on ${#GFFS[@]} Prokka GFF files"
roary -e -n -p "$THREADS" -i "$ROARY_MIN_IDENTITY" -cd "$ROARY_CORE_CUTOFF" \
  -f "$ROARY_OUT" "${GFFS[@]}"

ALIGNMENT="${ROARY_OUT}/core_gene_alignment.aln"
[[ -s "$ALIGNMENT" ]] || die "Roary core_gene_alignment.aln not found."

RAXML_BIN="$(command -v raxmlHPC-PTHREADS-AVX2 || command -v raxmlHPC-PTHREADS || true)"
[[ -n "$RAXML_BIN" ]] || die "RAxML v8 PTHREADS executable not found."
RAXML_OUT="${OUTPUT_DIR}/phylogeny/roary_raxml"
mkdir -p "$RAXML_OUT"
(
  cd "$RAXML_OUT"
  log "RAxML maximum-likelihood tree from Roary core-gene alignment"
  "$RAXML_BIN" -T "$THREADS" -f a -x 12345 -p 12345 -# "$RAXML_BOOTSTRAPS" \
    -m GTRGAMMA -s "$ALIGNMENT" -n Kp_Roary_Core
)

log "Roary/RAxML analysis completed. Tree: ${RAXML_OUT}/RAxML_bipartitions.Kp_Roary_Core"
