#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "$0")/common.sh"

require_cmd prokka
require_cmd prodigal
require_cmd busco
require_cmd mlst

mapfile -t GENOMES < <(assembly_files)
((${#GENOMES[@]} > 0)) || die "No assemblies found. Run 01_qc_assembly.sh first."
mkdir -p "${OUTPUT_DIR}/prokka" "${OUTPUT_DIR}/prodigal" "${OUTPUT_DIR}/busco" "${OUTPUT_DIR}/mlst"

INDEX=0
for GENOME in "${GENOMES[@]}"; do
  SAMPLE="$(basename "$GENOME" .fa)"
  INDEX=$((INDEX + 1))
  LOCUSTAG="$(printf 'KP%03d' "$INDEX")"
  PROKKA_OUT="${OUTPUT_DIR}/prokka/${SAMPLE}"

  log "Prokka annotation: $SAMPLE"
  PROKKA_ARGS=(--outdir "$PROKKA_OUT" --prefix "$SAMPLE" --locustag "$LOCUSTAG" \
    --kingdom Bacteria --genus Klebsiella --species pneumoniae --cpus "$THREADS" --force)
  if [[ -n "$PROKKA_PROTEINS" ]]; then
    [[ -s "$PROKKA_PROTEINS" ]] || die "PROKKA_PROTEINS not found: $PROKKA_PROTEINS"
    PROKKA_ARGS+=(--proteins "$PROKKA_PROTEINS")
  fi
  prokka "${PROKKA_ARGS[@]}" "$GENOME"

  log "Prodigal gene prediction: $SAMPLE"
  PRODIGAL_OUT="${OUTPUT_DIR}/prodigal/${SAMPLE}"
  mkdir -p "$PRODIGAL_OUT"
  prodigal -i "$GENOME" -p single -f gff \
    -o "${PRODIGAL_OUT}/${SAMPLE}.gff" \
    -a "${PRODIGAL_OUT}/${SAMPLE}.faa" \
    -d "${PRODIGAL_OUT}/${SAMPLE}.ffn" \
    -s "${PRODIGAL_OUT}/${SAMPLE}.genes.txt"

  if [[ -n "$HMM_PROFILE" ]]; then
    require_cmd hmmsearch
    [[ -s "$HMM_PROFILE" ]] || die "HMM_PROFILE not found: $HMM_PROFILE"
    hmmsearch --cpu "$THREADS" --tblout "${PRODIGAL_OUT}/${SAMPLE}.hmmsearch.tbl" \
      "$HMM_PROFILE" "${PRODIGAL_OUT}/${SAMPLE}.faa" \
      > "${PRODIGAL_OUT}/${SAMPLE}.hmmsearch.txt"
  fi

  log "BUSCO genome completeness: $SAMPLE"
  busco -i "$GENOME" -m genome -l "$BUSCO_LINEAGE" -c "$THREADS" \
    -o "${SAMPLE}_BUSCO" --out_path "${OUTPUT_DIR}/busco"
done

log "Pasteur K. pneumoniae MLST typing"
mlst --scheme kpneumoniae "${GENOMES[@]}" > "${OUTPUT_DIR}/mlst/Kp_MLST_results.tsv"

if [[ "$RUN_CHECKM" == "1" ]]; then
  require_cmd checkm
  log "CheckM quality assessment for all assemblies (CheckM 1.x syntax)"
  checkm lineage_wf -t "$THREADS" -x fa "$OUTPUT_DIR/assemblies" "${OUTPUT_DIR}/checkm"
fi

log "Annotation, gene prediction, BUSCO, MLST and optional CheckM completed."
