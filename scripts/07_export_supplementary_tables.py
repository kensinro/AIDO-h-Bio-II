#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AIDO-h-Biology II supplementary table exporter.

Purpose
-------
This utility organizes CSV result tables produced by the main analysis scripts into a
GitHub- and submission-friendly supplementary-table directory.

It does not rerun the analysis. It copies or converts existing result tables into a
numbered supplementary-table pack, and writes a manifest.

Typical inputs
--------------
1) Main AIDO-h multidatabase output directory, usually containing:
   - Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv
   - Table_1_cancer_database_observation_readiness_summary.csv
   - Table_2_database_overall_observation_readiness_summary.csv
   - Table_3_before_after_AIDOh_filter_summary.csv
   - Table_4_low_resolution_and_near_unobservable_terms.csv
   - Table_5_top_BP_cancer_by_D_survival.csv

2) Size-bin random baseline output directory, usually containing:
   - Table_SizeBinRandom_T300_bin_level_baselines.csv
   - Table_SizeBinRandom_T300_term_level_mapped_results.csv
   - Table_SizeBinRandom_T300_database_summary.csv
   - Table_SizeBinRandom_T300_cancer_database_summary.csv
   - Table_SizeBinRandom_T300_top_structured_favored_terms.csv
   - Table_SizeBinRandom_T300_highD_but_not_above_random_q95.csv

3) Reviewer-analysis output directory, if available.

Output
------
A directory with copied CSV files named as Supplementary_Table_S*.csv plus
SUPPLEMENTARY_TABLE_MANIFEST.csv.

Large files can optionally be split row-wise using --max-mb.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import pandas as pd


DEFAULT_MAPPING = [
    ("Supplementary_Table_S1_cancer_labels_input_file_status.csv", [
        "Run_status_all_cancers_MultiDB.csv",
        "cancer_labels_input_file_status.csv",
    ]),
    ("Supplementary_Table_S2_cancer_database_observation_readiness_summary.csv", [
        "Table_1_cancer_database_observation_readiness_summary.csv",
        "Table_1_cancer_database_observation_readiness_summary.xlsx",
    ]),
    ("Supplementary_Table_S3_database_observation_readiness_summary.csv", [
        "Table_2_database_overall_observation_readiness_summary.csv",
        "Table_2_database_overall_observation_readiness_summary.xlsx",
    ]),
    ("Supplementary_Table_S4_before_after_filter_survival_summary.csv", [
        "Table_3_before_after_AIDOh_filter_summary.csv",
        "Table_3_before_after_AIDOh_filter_summary.xlsx",
    ]),
    ("Supplementary_Table_S5_low_resolution_near_unobservable_terms.csv", [
        "Table_4_low_resolution_and_near_unobservable_terms.csv",
        "Table_4_low_resolution_and_near_unobservable_terms.xlsx",
    ]),
    ("Supplementary_Table_S6_all_cancer_database_term_full_results.csv", [
        "Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv",
    ]),
    ("Supplementary_Table_S7_top_BP_cancer_by_D_survival.csv", [
        "Table_5_top_BP_cancer_by_D_survival.csv",
    ]),
    ("Supplementary_Table_S8_random_baseline_term_level_mapped_results.csv", [
        "Table_SizeBinRandom_T300_term_level_mapped_results.csv",
    ]),
    ("Supplementary_Table_S9_random_baseline_bin_level_baselines.csv", [
        "Table_SizeBinRandom_T300_bin_level_baselines.csv",
    ]),
    ("Supplementary_Table_S10_random_baseline_database_summary.csv", [
        "Table_SizeBinRandom_T300_database_summary.csv",
    ]),
    ("Supplementary_Table_S11_random_baseline_cancer_database_summary.csv", [
        "Table_SizeBinRandom_T300_cancer_database_summary.csv",
    ]),
    ("Supplementary_Table_S12_top_structured_favored_terms.csv", [
        "Table_SizeBinRandom_T300_top_structured_favored_terms.csv",
    ]),
    ("Supplementary_Table_S13_highD_but_not_above_random_q95.csv", [
        "Table_SizeBinRandom_T300_highD_but_not_above_random_q95.csv",
    ]),
    ("Supplementary_Table_S14_threshold_sensitivity.csv", [
        "Table_threshold_sensitivity.csv",
        "threshold_sensitivity_summary.csv",
    ]),
    ("Supplementary_Table_S15_event_count_sensitivity.csv", [
        "Table_event_count_sensitivity.csv",
        "event_count_sensitivity.csv",
    ]),
    ("Supplementary_Table_S16_redundancy_or_overlap_analysis.csv", [
        "Table_GO_Reactome_redundancy_top_highD.csv",
        "Table_GO_Reactome_to_Hallmark_overlap.csv",
        "go_reactome_hallmark_overlap.csv",
    ]),
    ("Supplementary_Table_S17_cancer_gene_overlap_results.csv", [
        "Table_cancer_gene_overlap_results.csv",
        "cancer_gene_overlap_results.csv",
    ]),
    ("Supplementary_Table_S18_combined_cohort_sensitivity.csv", [
        "Table_combined_cohort_sensitivity.csv",
        "combined_cohort_sensitivity.csv",
    ]),
]


def find_candidate(filename_options: List[str], search_roots: List[Path]) -> Optional[Path]:
    for root in search_roots:
        if not root.exists():
            continue
        for name in filename_options:
            direct = root / name
            if direct.exists():
                return direct
        # recursive fallback, useful when tables are nested under /tables
        for name in filename_options:
            matches = list(root.rglob(name))
            if matches:
                return matches[0]
    return None


def split_csv_by_size(input_csv: Path, output_dir: Path, output_stem: str, max_mb: float) -> List[Path]:
    """Split a CSV into row-wise parts, preserving header in every part."""
    max_bytes = int(max_mb * 1024 * 1024)
    if input_csv.stat().st_size <= max_bytes:
        dst = output_dir / f"{output_stem}.csv"
        shutil.copy2(input_csv, dst)
        return [dst]

    parts = []
    with open(input_csv, "r", newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader)
        part_idx = 1
        current_rows = []
        current_size = len(",".join(header).encode("utf-8")) + 1

        def flush(rows, idx):
            out = output_dir / f"{output_stem}_part{idx:02d}.csv"
            with open(out, "w", newline="", encoding="utf-8") as g:
                writer = csv.writer(g)
                writer.writerow(header)
                writer.writerows(rows)
            parts.append(out)

        for row in reader:
            row_size = len(",".join(row).encode("utf-8")) + 1
            if current_rows and current_size + row_size > max_bytes:
                flush(current_rows, part_idx)
                part_idx += 1
                current_rows = []
                current_size = len(",".join(header).encode("utf-8")) + 1
            current_rows.append(row)
            current_size += row_size

        if current_rows:
            flush(current_rows, part_idx)

    return parts


def convert_to_csv(src: Path, temp_dir: Path) -> Path:
    if src.suffix.lower() == ".csv":
        return src
    if src.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(src)
        out = temp_dir / (src.stem + ".csv")
        df.to_csv(out, index=False)
        return out
    raise ValueError(f"Unsupported table format: {src}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main-root", type=Path, required=True, help="Main multidatabase AIDO-h output root or /tables folder.")
    parser.add_argument("--random-root", type=Path, default=None, help="Random-baseline output root or /tables folder.")
    parser.add_argument("--reviewer-root", type=Path, default=None, help="Reviewer-analysis output root or /tables folder.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output supplementary-table pack directory.")
    parser.add_argument("--max-mb", type=float, default=24.0, help="Split output CSVs larger than this size.")
    args = parser.parse_args()

    search_roots = [args.main_root]
    if args.random_root:
        search_roots.append(args.random_root)
    if args.reviewer_root:
        search_roots.append(args.reviewer_root)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.out_dir / "_tmp_converted"
    temp_dir.mkdir(exist_ok=True)

    manifest_rows = []
    for supp_name, candidates in DEFAULT_MAPPING:
        src = find_candidate(candidates, search_roots)
        if src is None:
            manifest_rows.append({
                "supplementary_table": supp_name,
                "source_file": "",
                "output_file": "",
                "status": "missing_source",
                "note": "No matching source table found in provided roots.",
            })
            continue

        csv_src = convert_to_csv(src, temp_dir)
        stem = Path(supp_name).stem
        outputs = split_csv_by_size(csv_src, args.out_dir, stem, args.max_mb)
        for out in outputs:
            manifest_rows.append({
                "supplementary_table": supp_name,
                "source_file": str(src),
                "output_file": out.name,
                "status": "exported",
                "note": "Split row-wise if multiple parts exist. Header is preserved in each part.",
            })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(args.out_dir / "SUPPLEMENTARY_TABLE_MANIFEST.csv", index=False)

    readme = args.out_dir / "README.md"
    readme.write_text(
        "# Supplementary Tables\n\n"
        "This directory was generated by `scripts/07_export_supplementary_tables.py`.\n\n"
        "Large supplementary tables may be split into row-wise parts to comply with file-size limits. "
        "Each split file preserves the same header as the original table.\n\n"
        "The manifest file `SUPPLEMENTARY_TABLE_MANIFEST.csv` records source files, exported files, "
        "and missing optional tables.\n",
        encoding="utf-8",
    )

    shutil.rmtree(temp_dir, ignore_errors=True)
    print(f"Done. Supplementary table pack written to: {args.out_dir}")


if __name__ == "__main__":
    main()
