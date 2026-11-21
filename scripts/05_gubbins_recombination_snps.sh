#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "$0")/common.sh"

require_cmd generate_ska_alignment.py
require_cmd run_gubbins.py
require_cmd mask_gubbins_aln.py
require_cmd snp-dists

[[ -s "$REFERENCE_FASTA" ]] || die "Set REFERENCE_FASTA in config/config.sh to one selected complete reference genome."
mapfile -t GENOMES < <(assembly_files)
((${#GENOMES[@]} >= 3)) || die "At least three assemblies are recommended for this comparative analysis."

GUBBINS_DIR="${OUTPUT_DIR}/gubbins"
mkdir -p "$GUBBINS_DIR"
INPUT_LIST="${GUBBINS_DIR}/ska_input.tsv"
: > "$INPUT_LIST"
REF_BASE="$(basename "$REFERENCE_FASTA")"
REF_ID="${REF_BASE%.*}"

# SKA2/Gubbins input: isolate name, then FASTA path. The designated reference
# is supplied separately to generate_ska_alignment.py and is excluded here.
for GENOME in "${GENOMES[@]}"; do
  GENOME_ID="$(basename "$GENOME" .fa)"
  if [[ "$GENOME_ID" == "$REF_ID" || "$(readlink -f "$GENOME")" == "$(readlink -f "$REFERENCE_FASTA")" ]]; then
    continue
  fi
  SAMPLE="$(basename "$GENOME" .fa)"
  printf '%s\t%s\n' "$SAMPLE" "$(readlink -f "$GENOME")" >> "$INPUT_LIST"
done

ALIGNMENT="${GUBBINS_DIR}/KP_whole_genome_alignment.aln"
log "Generating SKA2 whole-genome alignment using reference: $REFERENCE_FASTA"
generate_ska_alignment.py --reference "$REFERENCE_FASTA" --input "$INPUT_LIST" --out "$ALIGNMENT"
[[ -s "$ALIGNMENT" ]] || die "SKA2/Gubbins alignment was not generated."

PREFIX="${GUBBINS_DIR}/KP_Gubbins"
log "Running Gubbins recombination detection"
run_gubbins.py --prefix "$PREFIX" --threads "$THREADS" \
  --first-tree-builder rapidnj --first-model JC \
  --tree-builder "$GUBBINS_TREE_BUILDER" --model "$GUBBINS_MODEL" \
  --bootstrap "$GUBBINS_BOOTSTRAPS" "$ALIGNMENT"

GFF="${PREFIX}.recombination_predictions.gff"
[[ -s "$GFF" ]] || die "Gubbins recombination_predictions.gff not found."
MASKED="${GUBBINS_DIR}/KP_recombination_masked_alignment.aln"
mask_gubbins_aln.py --aln "$ALIGNMENT" --gff "$GFF" --out "$MASKED"

log "Calculating pairwise SNP distances from recombination-masked whole-genome alignment"
snp-dists -b -j "$THREADS" "$MASKED" > "${GUBBINS_DIR}/KP_pairwise_SNP_distances.tsv"

log "Gubbins, recombination masking and SNP-distance matrix completed."
