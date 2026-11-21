#!/usr/bin/env python3
"""Classify genome rows in a CheckM quality_report.tsv by completeness/contamination."""
import argparse
import csv
import re
import sys
from pathlib import Path


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, help="CheckM quality_report.tsv")
    parser.add_argument("--output", required=True, help="Output TSV with PASS/FAIL labels")
    parser.add_argument("--min-completeness", type=float, default=90.0)
    parser.add_argument("--max-contamination", type=float, default=5.0)
    args = parser.parse_args()

    report = Path(args.report)
    if not report.is_file():
        parser.error(f"CheckM report not found: {report}")

    with report.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            parser.error("Could not read a tab-separated header from the CheckM report.")
        columns = {norm(name): name for name in reader.fieldnames}
        bin_col = next((columns[k] for k in columns if k in {"binid", "bin", "genome", "genomeid"}), None)
        comp_col = next((columns[k] for k in columns if "completeness" in k), None)
        cont_col = next((columns[k] for k in columns if "contamination" in k), None)
        if not all((bin_col, comp_col, cont_col)):
            parser.error(
                "Could not identify Bin Id, Completeness and Contamination columns. "
                f"Found headers: {reader.fieldnames}"
            )

        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t")
            writer.writerow(["genome_id", "completeness_percent", "contamination_percent", "filter", "reason"])
            for row in reader:
                try:
                    completeness = float(row[comp_col])
                    contamination = float(row[cont_col])
                except (TypeError, ValueError):
                    continue
                passed = completeness >= args.min_completeness and contamination <= args.max_contamination
                reason = "meets_thresholds" if passed else "below_completeness_or_above_contamination"
                writer.writerow([row[bin_col], completeness, contamination, "PASS" if passed else "FAIL", reason])

    print(f"Wrote reference-genome QC filter: {output}")
    print(f"Criteria: completeness >= {args.min_completeness}% and contamination <= {args.max_contamination}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
