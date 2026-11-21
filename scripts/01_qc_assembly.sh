#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "$0")/common.sh"

require_cmd fastqc
require_cmd trimmomatic
require_cmd multiqc
require_cmd shovill
require_cmd quast.py

[[ -f "${ADAPTERS}" ]] || die "Set ADAPTERS in config/config.sh to the correct Trimmomatic adapter FASTA."

RAW_QC="${OUTPUT_DIR}/qc/raw_fastqc"
TRIM_DIR="${OUTPUT_DIR}/trimmed"
TRIM_QC="${OUTPUT_DIR}/qc/trimmed_fastqc"
ASSEMBLY_DIR="${OUTPUT_DIR}/assemblies"
mkdir -p "$RAW_QC" "$TRIM_DIR" "$TRIM_QC" "$ASSEMBLY_DIR" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/plasmidspades"

while IFS=$'\t' read -r SAMPLE R1 R2; do
  [[ -s "$R1" ]] || die "Missing R1 FASTQ: $R1"
  [[ -s "$R2" ]] || die "Missing R2 FASTQ: $R2"
  log "Raw-read FastQC: $SAMPLE"
  fastqc --threads "$THREADS" --outdir "$RAW_QC" "$R1" "$R2"
done < <(read_samples)

log "MultiQC summary of raw-read FastQC"
multiqc "$RAW_QC" --outdir "${OUTPUT_DIR}/qc/multiqc_raw" --force

while IFS=$'\t' read -r SAMPLE R1 R2; do
  P1="${TRIM_DIR}/${SAMPLE}_R1.paired.fastq.gz"
  U1="${TRIM_DIR}/${SAMPLE}_R1.unpaired.fastq.gz"
  P2="${TRIM_DIR}/${SAMPLE}_R2.paired.fastq.gz"
  U2="${TRIM_DIR}/${SAMPLE}_R2.unpaired.fastq.gz"

  log "Trimmomatic PE: $SAMPLE"
  trimmomatic PE -threads "$THREADS" -phred33 \
    "$R1" "$R2" "$P1" "$U1" "$P2" "$U2" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10" \
    "${TRIMMOMATIC_EXTRA_ARGS[@]}" \
    2> "${OUTPUT_DIR}/logs/${SAMPLE}.trimmomatic.log"

  log "Trimmed-read FastQC: $SAMPLE"
  fastqc --threads "$THREADS" --outdir "$TRIM_QC" "$P1" "$P2"

done < <(read_samples)

multiqc "$TRIM_QC" --outdir "${OUTPUT_DIR}/qc/multiqc_trimmed" --force

while IFS=$'\t' read -r SAMPLE R1 R2; do
  P1="${TRIM_DIR}/${SAMPLE}_R1.paired.fastq.gz"
  P2="${TRIM_DIR}/${SAMPLE}_R2.paired.fastq.gz"
  SHOVILL_OUT="${OUTPUT_DIR}/shovill/${SAMPLE}"
  mkdir -p "${OUTPUT_DIR}/shovill"

  log "Shovill/SPAdes assembly: $SAMPLE"
  shovill --R1 "$P1" --R2 "$P2" --outdir "$SHOVILL_OUT" \
    --assembler spades --cpus "$THREADS" --ram "$RAM_GB" \
    --depth "$SHOVILL_DEPTH" --force

  [[ -s "${SHOVILL_OUT}/contigs.fa" ]] || die "Shovill did not produce contigs.fa for $SAMPLE"
  cp "${SHOVILL_OUT}/contigs.fa" "${ASSEMBLY_DIR}/${SAMPLE}.fa"

  if [[ "$RUN_PILON" == "1" ]]; then
    require_cmd bwa
    require_cmd samtools
    [[ -n "$PILON_JAR" && -f "$PILON_JAR" ]] || require_cmd pilon
    PILON_DIR="${OUTPUT_DIR}/pilon/${SAMPLE}"
    mkdir -p "$PILON_DIR"
    DRAFT="${ASSEMBLY_DIR}/${SAMPLE}.fa"
    BAM="${PILON_DIR}/${SAMPLE}.sorted.bam"

    log "Pilon read mapping and polishing: $SAMPLE"
    bwa index "$DRAFT"
    bwa mem -t "$THREADS" "$DRAFT" "$P1" "$P2" \
      | samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index "$BAM"

    if [[ -n "$PILON_JAR" ]]; then
      java -Xmx"${RAM_GB}"G -jar "$PILON_JAR" \
        --genome "$DRAFT" --frags "$BAM" --output "${SAMPLE}.pilon" \
        --outdir "$PILON_DIR" --threads "$THREADS" --changes --vcf
    else
      pilon --genome "$DRAFT" --frags "$BAM" --output "${SAMPLE}.pilon" \
        --outdir "$PILON_DIR" --threads "$THREADS" --changes --vcf
    fi

    [[ -s "${PILON_DIR}/${SAMPLE}.pilon.fasta" ]] || die "Pilon output FASTA missing for $SAMPLE"
    cp "${PILON_DIR}/${SAMPLE}.pilon.fasta" "${ASSEMBLY_DIR}/${SAMPLE}.fa"
  fi

  if [[ "$RUN_PLASMIDSPADES" == "1" ]]; then
    require_cmd "$SPADES4_BIN"
    PLASMID_OUT="${OUTPUT_DIR}/plasmidspades/${SAMPLE}"
    log "Separate plasmidSPAdes assembly (SPAdes v4.x): $SAMPLE"
    "$SPADES4_BIN" --plasmid -1 "$P1" -2 "$P2" \
      -o "$PLASMID_OUT" -t "$THREADS" -m "$RAM_GB" \
      > "${OUTPUT_DIR}/logs/${SAMPLE}.plasmidspades.log" 2>&1
  fi
done < <(read_samples)

# Stage the public reference assemblies; keep their names unique and informative.
shopt -s nullglob
REFS=("${REFERENCES_DIR}"/*.fa "${REFERENCES_DIR}"/*.fna "${REFERENCES_DIR}"/*.fasta)
for REF in "${REFS[@]}"; do
  [[ -f "$REF" ]] || continue
  REF_BASE="$(basename "$REF")"
  REF_NAME="${REF_BASE%.*}"
  cp "$REF" "${ASSEMBLY_DIR}/${REF_NAME}.fa"
done
shopt -u nullglob

mapfile -t GENOMES < <(find "$ASSEMBLY_DIR" -maxdepth 1 -type f -name '*.fa' -print | sort)
((${#GENOMES[@]} > 0)) || die "No genome FASTA files available for QUAST."
quast.py --threads "$THREADS" --output-dir "${OUTPUT_DIR}/quast" "${GENOMES[@]}"

log "QC, trimming, assembly, optional Pilon polishing, plasmidSPAdes and QUAST completed."
