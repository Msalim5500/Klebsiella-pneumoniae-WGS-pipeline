#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=../config/config.sh
source "${PROJECT_ROOT}/config/config.sh"

mkdir -p "${OUTPUT_DIR}"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required executable not found on PATH: $1"; }
resolve_path() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *)  printf '%s/%s\n' "$PROJECT_ROOT" "$1" ;;
  esac
}

read_samples() {
  [[ -s "${SAMPLES_TSV}" ]] || die "Sample sheet not found: ${SAMPLES_TSV}"
  while IFS=$'\t' read -r sample_id r1 r2; do
    [[ -z "${sample_id:-}" || "$sample_id" == "sample_id" ]] && continue
    [[ -n "${r1:-}" && -n "${r2:-}" ]] || die "Missing R1/R2 in sample sheet row for ${sample_id}"
    printf '%s\t%s\t%s\n' "$sample_id" "$(resolve_path "$r1")" "$(resolve_path "$r2")"
  done < "$SAMPLES_TSV"
}

assembly_files() {
  find "${OUTPUT_DIR}/assemblies" -maxdepth 1 -type f -name '*.fa' -print 2>/dev/null | sort
}
