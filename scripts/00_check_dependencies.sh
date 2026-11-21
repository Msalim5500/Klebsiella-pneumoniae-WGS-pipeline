#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "$0")/common.sh"

mkdir -p "${OUTPUT_DIR}/logs"
REPORT="${OUTPUT_DIR}/logs/software_versions.txt"
: > "$REPORT"

TOOLS=(fastqc trimmomatic multiqc shovill quast.py bwa samtools pilon "$SPADES4_BIN" busco prodigal hmmsearch prokka mlst abricate amrfinder staramr plasmidfinder.py mob_recon pathogenfinder2 roary raxmlHPC-PTHREADS-AVX2 raxmlHPC-PTHREADS run_gubbins.py generate_ska_alignment.py mask_gubbins_aln.py snp-dists)

for tool in "${TOOLS[@]}"; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '\n### %s\n' "$tool" | tee -a "$REPORT"
    case "$(basename "$tool")" in
      trimmomatic) trimmomatic -version 2>&1 | tee -a "$REPORT" || true ;;
      quast.py) quast.py --version 2>&1 | tee -a "$REPORT" || true ;;
      abricate) abricate --version 2>&1 | tee -a "$REPORT" || true; abricate --list 2>&1 | tee -a "$REPORT" || true ;;
      mlst) mlst --version 2>&1 | tee -a "$REPORT" || true ;;
      *) "$tool" --version 2>&1 | head -n 3 | tee -a "$REPORT" || true ;;
    esac
  else
    printf 'MISSING\t%s\n' "$tool" | tee -a "$REPORT"
  fi
done

log "Version report written to ${REPORT}"
log "Review missing tools and save database versions/access dates before running the full analysis."
