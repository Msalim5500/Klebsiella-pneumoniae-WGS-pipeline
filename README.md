# Klebsiella pneumoniae WGS analysis pipeline

A reproducible, command-line workflow template for paired-end Illumina/DNBSEQ PE150 *Klebsiella pneumoniae* whole-genome sequencing. The repository is configured for seven local NOHU isolates and fourteen externally selected reference genomes, but the sample sheet and reference directory can be edited for other datasets.

**This is a configurable research workflow, not a pre-installed software environment.** Install the named tools, configure the database paths and versions, test each stage on one isolate, and retain the software/database version report with the analysis.

## Workflow

1. Raw-read QC: FastQC v0.12.1; MultiQC v1.33 summary.
2. Adapter/quality trimming: Trimmomatic v0.40.
3. Draft assemblies: Shovill v1.4.2 using SPAdes v3.15.5.
4. Optional short-read polishing: BWA + SAMtools + Pilon.
5. Separate plasmid-focused assembly: plasmidSPAdes mode in SPAdes v4.0.0.
6. Assembly statistics: QUAST v5.2.0.
7. Genome quality: BUSCO v6.0.0; optional CheckM for genome quality screening.
8. Gene prediction and annotation: Prodigal v2.6.3, optional HMMER `hmmsearch` v3.4, and Prokka v1.14.6.
9. Sequence typing: `mlst` v2.23.0 using the Pasteur *K. pneumoniae* scheme (verify the scheme name installed locally).
10. Resistance/virulence screening: ABRicate v1.2.0, NCBI AMRFinderPlus (separate executable), StarAMR v0.12.1 (optional), VFDB and VICTORS.
11. Plasmid characterization: PlasmidFinder v2.1.6 and MOB-suite v2.0.1.
12. Pathogenicity prediction: PathogenFinder2.
13. Pan-genome and core-gene tree: Roary + RAxML, 1,000 bootstrap replicates.
14. Recombination and SNP analysis: SKA2 whole-genome alignment → Gubbins v3.4.3 → masked alignment → `snp-dists` pairwise matrix.
15. Optional Sankey diagram: R + ggplot2 + ggalluvial.

## Repository layout

```text
Kp_WGS_GitHub_Pipeline/
├── config/
│   └── config.sh
├── metadata/
│   ├── samples.tsv
│   ├── reference_genomes.tsv
│   └── sankey_metadata_template.tsv
├── scripts/
│   ├── common.sh
│   ├── 00_check_dependencies.sh
│   ├── 01_qc_assembly.sh
│   ├── 02_annotation_quality.sh
│   ├── 03_resistome_virulence_plasmids.sh
│   ├── 04_roary_core_tree.sh
│   ├── 05_gubbins_recombination_snps.sh
│   ├── 06_sankey.R
│   ├── filter_checkm_refs.py
│   └── run_all.sh
├── docs/
└── README.md
```

Raw reads, reference FASTA files and generated results are deliberately excluded by `.gitignore`. Do not commit patient-derived data or unpublished sequence files to a public repository unless you have permission and have checked the relevant data-sharing requirements.

## 1. Prepare the data

Create these folders:

```bash
mkdir -p data/reads data/references
```

Place each paired FASTQ file in `data/reads/` using the naming convention shown in `metadata/samples.tsv`, for example:

```text
NOHU9_2KP_R1.fastq.gz
NOHU9_2KP_R2.fastq.gz
NOHU38_7KP_R1.fastq.gz
NOHU38_7KP_R2.fastq.gz
...
```

Place the fourteen downloaded NCBI reference genome FASTA files in `data/references/`. Use clear, unique file names, preferably including the accession, e.g. `GCF_XXXXXXXXX.1.fa`. Fill out `metadata/reference_genomes.tsv` with accession, source, host, country, year, ST, CheckM completeness and contamination.

Edit `metadata/samples.tsv` if your filenames or sample IDs differ. It must be a tab-separated file with exactly these first three columns:

```text
sample_id    R1    R2
```

Paths can be absolute or relative to the repository root.

## 1A. Select the fourteen public reference genomes

The selection of the fourteen comparator genomes is a curated step and is not fully automated by the main scripts. Use BLASTn to identify close *K. pneumoniae* assemblies, then retain only human clinical genomes that meet the study criteria. For a local NCBI `nt` database, an example command is:

```bash
blastn -query results/assemblies/NOHU9_2KP.fa -db /path/to/nt \
  -num_threads 16 -max_target_seqs 100 \
  -outfmt "6 qseqid saccver pident length evalue bitscore stitle" \
  > results/NOHU9_2KP_BLASTn_hits.tsv
```

Repeat for the local isolates and use the hits to shortlist candidate accessions. An NCBI remote BLAST search is also possible, but it sends the query sequences to NCBI; use it only when that is acceptable for your unpublished data. Download the selected public assemblies into `data/references/`, then record their metadata in `metadata/reference_genomes.tsv`.

After running CheckM, create a table showing which genomes pass the stated quality criteria:

```bash
python scripts/filter_checkm_refs.py \
  --report results/checkm/quality_report.tsv \
  --output results/checkm/reference_genome_filter.tsv \
  --min-completeness 90 --max-contamination 5
```

Review the PASS/FAIL table and retain the fourteen eligible human clinical genomes for the final comparative analysis. The script reports the threshold decision; it does not automatically remove or replace any FASTA files.

## 2. Configure the software paths and analysis parameters

Open:

```bash
nano config/config.sh
```

At minimum, edit:

- `ADAPTERS`: full path to the correct Trimmomatic adapter FASTA for your library preparation.
- `SPADES4_BIN`: full path to the SPAdes v4.0.0 executable used for the separate plasmidSPAdes run. Keep the SPAdes v3.15.5 installation required by Shovill separate from this SPAdes v4 installation, e.g. use separate Conda/micromamba environments.
- `PLASMIDFINDER_DB`: local PlasmidFinder database directory.
- `REFERENCE_FASTA`: one selected complete, high-quality *K. pneumoniae* reference genome from `data/references/`.
- `PROKKA_PROTEINS`: optional *K. pneumoniae* proteome FASTA used by Prokka, if you have the exact file.
- `BUSCO_LINEAGE`: confirm that the requested lineage exists in your BUSCO installation. The manuscript draft says `bacteria_odb10`; check with `busco --list-datasets` and record the exact dataset/version.
- `HMM_PROFILE`: optional HMM database/profile file used by the stand-alone `hmmsearch` step. Do not run `hmmsearch` without a biologically relevant profile HMM.
- `THREADS` and `RAM_GB`: set according to your Linux workstation/server.

Use the database snapshot and tool versions actually used in the analysis. Avoid silently updating databases between isolates.

## 3. Check installed software

Activate the relevant micromamba/Conda environments, then run:

```bash
chmod +x scripts/*.sh
bash scripts/00_check_dependencies.sh
```

The report is written to `results/logs/software_versions.txt`. The script reports missing tools but does not install them. Some tools have separate environments or database dependencies; PathogenFinder2 should be installed in its own Python environment as recommended by its developers.

Before running all genomes, test the relevant commands independently:

```bash
fastqc --version
trimmomatic -version
shovill --version
quast.py --version
prokka --version
mlst --list | grep -i klebsiella
abricate --list
busco --list-datasets | grep -i bacteria
run_gubbins.py --help
pathogenfinder2 predict -h
```

If the MLST scheme is named differently in your installed database, update `scripts/02_annotation_quality.sh` from `kpneumoniae` to the exact scheme name shown by `mlst --list`. Confirm that the alleles correspond to the Pasteur scheme (`gapA`, `infB`, `mdh`, `pgi`, `phoE`, `rpoB`, `tonB`).

## 4. Run the workflow

Run one stage at a time first:

```bash
bash scripts/01_qc_assembly.sh
bash scripts/02_annotation_quality.sh
bash scripts/03_resistome_virulence_plasmids.sh
bash scripts/04_roary_core_tree.sh
bash scripts/05_gubbins_recombination_snps.sh
```

Or run all main stages sequentially:

```bash
bash scripts/run_all.sh
```

The full workflow can take substantial CPU time, memory and disk space. Check the log/output for each stage before continuing. For a first test, reduce `metadata/samples.tsv` to one or two isolates and use a small test reference set.

## 5. Key output files

| Analysis | Output location / file |
|---|---|
| Raw FastQC | `results/qc/raw_fastqc/` |
| Raw MultiQC | `results/qc/multiqc_raw/` |
| Trimmed FASTQ | `results/trimmed/` |
| Trimmed FastQC and MultiQC | `results/qc/trimmed_fastqc/`, `results/qc/multiqc_trimmed/` |
| Shovill assemblies | `results/shovill/<sample>/` |
| Polished assemblies | `results/assemblies/<sample>.fa` (when Pilon is enabled) |
| plasmidSPAdes | `results/plasmidspades/<sample>/` |
| QUAST | `results/quast/` |
| Prokka | `results/prokka/<genome>/` |
| Prodigal and optional hmmsearch | `results/prodigal/<genome>/` |
| BUSCO | `results/busco/` |
| MLST | `results/mlst/Kp_MLST_results.tsv` |
| ABRicate | `results/abricate/` |
| AMRFinderPlus | `results/amrfinder/` |
| StarAMR | `results/staramr/` |
| PlasmidFinder | `results/plasmidfinder/` |
| MOB-suite | `results/mob_suite/` |
| PathogenFinder2 | `results/pathogenfinder2/output/` |
| Roary | `results/roary/` |
| RAxML tree | `results/phylogeny/roary_raxml/RAxML_bipartitions.Kp_Roary_Core` |
| Gubbins | `results/gubbins/KP_Gubbins.*` |
| Pairwise SNP matrix | `results/gubbins/KP_pairwise_SNP_distances.tsv` |

## 6. Important interpretation and reproducibility notes

### Roary and Gubbins are separate alignments

Roary's `core_gene_alignment.aln` is used for the Roary-based core-gene tree. Gubbins is run on a separate whole-genome alignment created with SKA2. Do **not** feed the concatenated Roary core-gene alignment into Gubbins: the gene concatenation does not retain the original genomic coordinates needed for recombination mapping. Gubbins requires a whole-genome alignment and its `mask_gubbins_aln.py` helper can mask inferred recombinant regions. Pairwise SNP distances are then calculated from that masked alignment.

### Pilon and plasmidSPAdes outputs are distinct

Pilon polishing uses paired reads mapped back to the Shovill draft assembly (BWA + SAMtools). The plasmidSPAdes run is a separate plasmid-focused assembly from the same trimmed paired reads. Do not merge plasmidSPAdes contigs into the chromosome assembly without a justified validation workflow, and do not describe plasmidSPAdes contigs as closed plasmids unless circularity/assembly evidence supports that conclusion.

### AMR database terminology

ABRicate's `ncbi` database is not the same analysis as a full NCBI AMRFinderPlus run. This repository runs `amrfinder` separately. The manuscript should distinguish them. ABRicate supports databases such as CARD, ResFinder, NCBI, MEGARes, VFDB, VICTORS, PlasmidFinder and BacMet; the exact local database names and dates are shown by `abricate --list`.

The Methods draft says “eleven databases” but enumerates eight named databases. Verify the exact list actually run, list each database once, and do not report eleven unless eleven distinct databases were used. ABRicate has other possible databases, but availability does not automatically make all of them appropriate for *K. pneumoniae*. Avoid including organism-specific databases (for example, *E. coli*-specific virulence databases) as evidence for *K. pneumoniae* without a reasoned justification.

### StarAMR / PointFinder for *K. pneumoniae*

The current StarAMR documentation lists PointFinder organism options for a limited set of species and does not list *K. pneumoniae*. Therefore this template runs StarAMR without `--pointfinder-organism`. Do not report StarAMR as having screened *K. pneumoniae* PointFinder mutations in `gyrA`, `parC`, `ompK35` or `ompK36` unless the precise StarAMR/PointFinder database version you used explicitly supports those targets. Porin changes such as OmpK35/OmpK36 require a separate, validated sequence comparison and interpretation workflow.

### VirulenceFinder for *K. pneumoniae*

The current official VirulenceFinder documentation lists *E. coli*, *Enterococcus*, *S. aureus* and *Listeria* as available species, not *K. pneumoniae*. Accordingly, `RUN_VIRULENCEFINDER=0` by default. For *K. pneumoniae*, use appropriately curated VFDB/VICTORS screening and/or a validated species-specific database. Only enable VirulenceFinder if the exact tool/database version used for the study supports *K. pneumoniae* and you can document that support.

### CheckM, BUSCO and reference filtering

BUSCO and CheckM provide different completeness/contamination assessments. For the fourteen NCBI comparison genomes, apply the study's stated filter (completeness >=90% and contamination <=5%) using the relevant CheckM output before finalizing the comparison set. Keep the accepted/rejected accession list and quality report in your analysis records. BUSCO lineage names change across releases; confirm `bacteria_odb10` exists in the installed BUSCO data, otherwise choose and report the most appropriate current bacterial lineage. Do not silently substitute lineages after analysis.

### Prodigal and hmmsearch

Prokka performs gene prediction/annotation and uses Prodigal internally. The additional Prodigal run in this repository records explicit nucleotide/protein gene predictions. The optional `hmmsearch` command only runs when `HMM_PROFILE` is set. The Methods must name the HMM profile/database, version and cutoffs used; `hmmsearch` alone does not identify orthologues without a defined profile set.

### RAxML vs Gubbins trees

The RAxML tree from Roary and the Gubbins tree from a recombination-aware whole-genome alignment are different analyses. Report them separately, and state which tree is shown in each figure. The RAxML command in this repository targets the RAxML v8 PTHREADS executable. If using RAxML-NG, use its distinct command-line syntax.

### NCBI GenBank submissions

The workflow does not automatically submit genomes to GenBank. Submit assemblies and metadata through the NCBI submission process, then add the returned accession numbers to Table 1 and `metadata/reference_genomes.tsv` as appropriate. Do not upload controlled or patient-identifiable metadata to GitHub.

## 7. Optional Sankey diagram

Fill `metadata/sankey_metadata_template.tsv` with the three categorical variables to be displayed (columns `stage1`, `stage2`, `stage3`) and a numeric `value`. Then run:

```bash
Rscript scripts/06_sankey.R metadata/sankey_metadata_template.tsv results/figures/Kp_Sankey_diagram.pdf
```

The stage labels in the script are generic placeholders. Edit them to match the actual variables shown in Fig. 4 and make sure the data structure is appropriate for an alluvial diagram.

## 8. Upload this repository to GitHub

First, create a new empty repository on GitHub, for example `Klebsiella-pneumoniae-WGS-pipeline`. Do not initialize it with another README because this folder already contains one.

On your Linux machine, open a terminal in this downloaded repository folder and run:

```bash
cd Kp_WGS_GitHub_Pipeline

git init
git add .
git status
# Confirm that data/reads, data/references and results are NOT staged.
git commit -m "Add reproducible K. pneumoniae WGS analysis pipeline"
git branch -M main
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/Klebsiella-pneumoniae-WGS-pipeline.git
git push -u origin main
```

Replace `YOUR_GITHUB_USERNAME` with your GitHub username and use the exact repository URL shown by GitHub. GitHub may ask you to authenticate. For later updates:

```bash
git add README.md config/ metadata/ scripts/ docs/
git commit -m "Update WGS pipeline"
git push
```

Before publishing, inspect `git status` and `.gitignore` carefully. Do not commit FASTQ files, patient metadata, raw sequencing reads, reference sequence data with redistribution restrictions, or large analysis outputs unless you have checked permission and repository limits.

## 9. Main software documentation

- FastQC: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
- Trimmomatic: <https://github.com/usadellab/Trimmomatic>
- MultiQC: <https://multiqc.info/>
- Shovill: <https://github.com/tseemann/shovill>
- SPAdes / plasmidSPAdes: <https://ablab.github.io/spades/running.html>
- Pilon: <https://github.com/broadinstitute/pilon>
- QUAST: <https://quast.sourceforge.net/>
- BUSCO: <https://busco.ezlab.org/busco_userguide>
- Prokka: <https://github.com/tseemann/prokka>
- MLST: <https://github.com/tseemann/mlst>
- ABRicate: <https://github.com/tseemann/abricate>
- NCBI AMRFinderPlus: <https://github.com/ncbi/amr>
- StarAMR: <https://github.com/phac-nml/staramr>
- PlasmidFinder: <https://github.com/genomicepidemiology/plasmidfinder>
- MOB-suite: <https://github.com/phac-nml/mob-suite>
- PathogenFinder2: <https://github.com/genomicepidemiology/PathogenFinder2>
- Roary: <https://github.com/sanger-pathogens/Roary>
- Gubbins: <https://github.com/nickjcroucher/gubbins>
- snp-dists: <https://github.com/tseemann/snp-dists>
- iTOL: <https://itol.embl.de/>

### Visualize the phylogeny in iTOL

After the tree has been generated, open <https://itol.embl.de/> and upload the Newick tree you intend to present (for example, `RAxML_bipartitions.Kp_Roary_Core` for the Roary tree or `KP_Gubbins.final_tree.tre` for the Gubbins tree). Upload or create iTOL annotation datasets separately using the sample metadata (ST, country, year, source, resistance phenotype, etc.). Do not combine annotations from different tree tip labels; the isolate names in the annotation file must match the Newick tip labels exactly.
