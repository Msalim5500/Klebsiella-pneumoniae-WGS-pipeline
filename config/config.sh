#!/usr/bin/env bash
# User-editable configuration for the K. pneumoniae WGS workflow.
# Paths in metadata/samples.tsv are interpreted relative to PROJECT_ROOT.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
READS_DIR="${DATA_DIR}/reads"
REFERENCES_DIR="${DATA_DIR}/references"
SAMPLES_TSV="${PROJECT_ROOT}/metadata/samples.tsv"
OUTPUT_DIR="${PROJECT_ROOT}/results"

# Compute resources
THREADS=16
RAM_GB=32

# Read processing
ADAPTERS="/EDIT/THIS/PATH/Trimmomatic/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_EXTRA_ARGS=(LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50)
SHOVILL_DEPTH=0                  # 0 = do not downsample reads
RUN_PILON=1                      # Set to 0 to skip Pilon polishing
PILON_JAR=""                    # Optional: /path/to/pilon.jar; otherwise `pilon` must be on PATH
RUN_PLASMIDSPADES=1              # Separate plasmidSPAdes run; not merged into the chromosome assembly
SPADES4_BIN="/EDIT/THIS/PATH/SPAdes-4.0.0/bin/spades.py"  # Separate SPAdes v4.0.0 executable for plasmidSPAdes

# Annotation and quality assessment
PROKKA_PROTEINS=""              # Optional K. pneumoniae reference proteome FASTA
BUSCO_LINEAGE="bacteria_odb10"  # Confirm with: busco --list-datasets
HMM_PROFILE=""                  # Optional profile HMM file for standalone hmmsearch
RUN_CHECKM=1                    # CheckM 1.x command; adjust if using CheckM2

# Database directories (edit these to match local installations)
PLASMIDFINDER_DB="/EDIT/THIS/PATH/plasmidfinder_db"
VIRULENCEFINDER_DB="/EDIT/THIS/PATH/virulencefinder_db"

# ABRicate database names. Run `abricate --list` and remove names not installed.
# 'ncbi' in ABRicate is NOT equivalent to a full AMRFinderPlus run.
ABRICATE_DATABASES=(card resfinder ncbi megares vfdb victors plasmidfinder bacmet2 argannot)
ABRICATE_MINID=90
ABRICATE_MINCOV=80

# Separate AMRFinderPlus run
AMRFINDER_ORGANISM="Klebsiella_pneumoniae"

# Optional tools. VirulenceFinder and StarAMR limitations for K. pneumoniae
# are documented in README.md. Leave VirulenceFinder disabled unless the exact
# version/database you used explicitly supports K. pneumoniae.
RUN_STARAMR=1
RUN_VIRULENCEFINDER=0
RUN_PATHOGENFINDER2=1

# Comparative genomics
ROARY_MIN_IDENTITY=95
ROARY_CORE_CUTOFF=99
RAXML_BOOTSTRAPS=1000
GUBBINS_BOOTSTRAPS=1000
GUBBINS_TREE_BUILDER="raxml"
GUBBINS_MODEL="GTRGAMMA"
REFERENCE_FASTA="${REFERENCES_DIR}/EDIT_THIS_TO_THE_SELECTED_COMPLETE_REFERENCE.fa"

# Optional plotting; edit metadata/sankey_metadata.tsv and run scripts/06_sankey.R manually.
