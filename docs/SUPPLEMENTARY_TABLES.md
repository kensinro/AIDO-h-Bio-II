# Supplementary tables

Large supplementary tables are not included directly in this GitHub repository.

Reason:
- GitHub web upload has file-size limits.
- Generated all-cancer term-level tables can be large.
- The repository should primarily serve as the reproducible code package.

Use:

```bash
python scripts/07_export_supplementary_tables.py --main-root <main_output_root> --random-root <random_output_root> --reviewer-root <reviewer_output_root> --out-dir <supplementary_table_output_root>
```

This creates:
- `Supplementary_Table_S*.csv`
- row-wise split parts for large files
- `SUPPLEMENTARY_TABLE_MANIFEST.csv`
- `README.md`

Important wording for manuscript/revision:
"The GitHub repository contains the analysis code used to generate the source result
tables underlying the supplementary tables. The full formatted supplementary tables
were provided with the manuscript submission."
