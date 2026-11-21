#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname "$0")/common.sh"

require_cmd abricate
mapfile -t GENOMES < <(assembly_files)
((${#GENOMES[@]} > 0)) || die "No assemblies found. Run 01_qc_assembly.sh first."
mkdir -p "${OUTPUT_DIR}/abricate" "${OUTPUT_DIR}/amrfinder" "${OUTPUT_DIR}/staramr" \
  "${OUTPUT_DIR}/plasmidfinder" "${OUTPUT_DIR}/mob_suite" "${OUTPUT_DIR}/virulencefinder" \
  "${OUTPUT_DIR}/pathogenfinder2"

# Run only databases actually installed in this ABRicate installation.
INSTALLED_DBS="$(abricate --list | awk 'NR>1 {print $1}')"
for DB in "${ABRICATE_DATABASES[@]}"; do
  if ! grep -Fxq "$DB" <<< "$INSTALLED_DBS"; then
    log "ABRicate database '$DB' is not installed; skipping it. Use abricate --list to inspect available DBs."
    continue
  fi
  log "ABRicate database: $DB"
  abricate --db "$DB" --minid "$ABRICATE_MINID" --mincov "$ABRICATE_MINCOV" \
    "${GENOMES[@]}" > "${OUTPUT_DIR}/abricate/${DB}_all_genomes.tsv"
done

# Full NCBI AMRFinderPlus is run independently from ABRicate's optional 'ncbi' DB.
if command -v amrfinder >/dev/null 2>&1; then
  for GENOME in "${GENOMES[@]}"; do
    SAMPLE="$(basename "$GENOME" .fa)"
    log "NCBI AMRFinderPlus: $SAMPLE"
    amrfinder --nucleotide "$GENOME" --organism "$AMRFINDER_ORGANISM" --plus \
      --output "${OUTPUT_DIR}/amrfinder/${SAMPLE}_AMRFinderPlus.tsv"
  done
else
  log "AMRFinderPlus not found; skipping. Install it separately and update its database before analysis."
fi

# StarAMR search is run without PointFinder organism flag. Its documented
# PointFinder species list does not currently include K. pneumoniae; do not
# claim K. pneumoniae PointFinder mutations from this command.
if [[ "$RUN_STARAMR" == "1" ]]; then
  if command -v staramr >/dev/null 2>&1; then
    for GENOME in "${GENOMES[@]}"; do
      SAMPLE="$(basename "$GENOME" .fa)"
      log "StarAMR gene/plasmid screen (no PointFinder organism specified): $SAMPLE"
      staramr search -o "${OUTPUT_DIR}/staramr/${SAMPLE}" "$GENOME"
    done
  else
    log "StarAMR not found; skipping."
  fi
fi

# VirulenceFinder is disabled by default because the current CGE documentation
# lists E. coli, Enterococcus, S. aureus and Listeria, not K. pneumoniae.
if [[ "$RUN_VIRULENCEFINDER" == "1" ]]; then
  command -v virulencefinder.py >/dev/null 2>&1 || die "virulencefinder.py not found."
  [[ -d "$VIRULENCEFINDER_DB" ]] || die "Set VIRULENCEFINDER_DB to the installed database directory."
  for GENOME in "${GENOMES[@]}"; do
    SAMPLE="$(basename "$GENOME" .fa)"
    virulencefinder.py -i "$GENOME" -o "${OUTPUT_DIR}/virulencefinder/${SAMPLE}" \
      -p "$VIRULENCEFINDER_DB"
  done
fi

# Dedicated PlasmidFinder; thresholds correspond to 95% identity and 60% coverage.
if command -v plasmidfinder.py >/dev/null 2>&1 && [[ -d "$PLASMIDFINDER_DB" ]]; then
  for GENOME in "${GENOMES[@]}"; do
    SAMPLE="$(basename "$GENOME" .fa)"
    log "PlasmidFinder: $SAMPLE"
    plasmidfinder.py -i "$GENOME" -o "${OUTPUT_DIR}/plasmidfinder/${SAMPLE}" \
      -p "$PLASMIDFINDER_DB" -l 60 -t 95 -x -q
  done
else
  log "Standalone PlasmidFinder or its database not found; skipping dedicated PlasmidFinder."
fi

if command -v mob_recon >/dev/null 2>&1; then
  for GENOME in "${GENOMES[@]}"; do
    SAMPLE="$(basename "$GENOME" .fa)"
    log "MOB-suite: $SAMPLE"
    mob_recon -i "$GENOME" -o "${OUTPUT_DIR}/mob_suite/${SAMPLE}" \
      -n "$THREADS" -f
  done
else
  log "MOB-suite not found; skipping."
fi

if [[ "$RUN_PATHOGENFINDER2" == "1" ]]; then
  require_cmd pathogenfinder2
  LIST="${OUTPUT_DIR}/pathogenfinder2/genomes.txt"
  printf '%s\n' "${GENOMES[@]}" > "$LIST"
  pathogenfinder2 predict --multipleFiles "$LIST" --formatSeq genome \
    --outputFolder "${OUTPUT_DIR}/pathogenfinder2/output"
fi

log "Resistome, optional virulence, plasmid, mobility and pathogenicity steps completed."
