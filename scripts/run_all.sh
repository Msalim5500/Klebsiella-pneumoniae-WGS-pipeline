#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

"${SCRIPT_DIR}/00_check_dependencies.sh"
"${SCRIPT_DIR}/01_qc_assembly.sh"
"${SCRIPT_DIR}/02_annotation_quality.sh"
"${SCRIPT_DIR}/03_resistome_virulence_plasmids.sh"
"${SCRIPT_DIR}/04_roary_core_tree.sh"
"${SCRIPT_DIR}/05_gubbins_recombination_snps.sh"

echo "Main pipeline completed. Review all logs, software versions, and QC outputs before interpreting results."
