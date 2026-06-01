# ============================================================
# AIDO-h-Biology II
# Reviewer Analysis Pack
#
# Purpose:
#   Complete supplementary / reviewer-defense analyses:
#
#   1. AIDO-h threshold sensitivity
#   2. Top-D before vs after AIDO-h filtering
#   3. Database-level and cancer-level regime summaries
#   4. Event-count sensitivity
#   5. Low-resolution / near-unobservable characterization
#   6. GO/Reactome redundancy among top high-D terms
#   7. Combined-cohort sensitivity
#   8. Optional exact Nmatched random validation for top candidates
#
# Input:
#   Main MultiDB result:
#     D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/
#       Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv
#
#   Optional size-bin random results:
#     D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300/tables/
#       Table_SizeBinRandom_T300_term_level_mapped_results.csv
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-ReviewerAnalyses/
# ============================================================

import os
import re
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# ============================================================
# 0. CONFIG
# ============================================================

BASE_DATA_DIR = Path(r"D:/AIDO-Data/UCSC_XENA")

MAIN_RESULT_FILE = Path(
    r"D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/"
    r"Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv"
)

SIZEBIN_RANDOM_FILE = Path(
    r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300/tables/"
    r"Table_SizeBinRandom_T300_term_level_mapped_results.csv"
)

OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-ReviewerAnalyses")

TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
LOG_DIR = OUT_DIR / "logs"
EXACT_RANDOM_DIR = OUT_DIR / "exact_random_validation"

for d in [OUT_DIR, TABLE_DIR, FIG_DIR, LOG_DIR, EXACT_RANDOM_DIR]:
    d.mkdir(parents=True, exist_ok=True)

DATABASES = ["Hallmark", "GO_BP", "Reactome"]

D_THRESHOLD = 1.301

# AIDO-h sensitivity thresholds
AIDOH_THRESHOLDS = [5, 10, 15, 20]

# Near-unobservable threshold is kept fixed.
NEAR_UNOBSERVABLE_LT = 3

# Combined / overlapping cohort labels.
COMBINED_COHORTS = ["COADREAD", "LUNG"]

# For redundancy analysis
REDUNDANCY_TOP_N_PER_CANCER_DATABASE = 100
REDUNDANCY_ONLY_HIGH_D = True
REDUNDANCY_DATABASES = ["GO_BP", "Reactome"]

# Exact random validation
RUN_EXACT_RANDOM_VALIDATION = True

EXACT_RANDOM_T = 1000

# Candidate selection for exact validation:
# safer default: not too huge
EXACT_TOP_N_PER_DATABASE_GLOBAL = 20
EXACT_TOP_N_PER_CANCER_DATABASE = 5

# If too slow, reduce this:
MAX_EXACT_RANDOM_TERMS = 250
# MAX_EXACT_RANDOM_TERMS = 100

MIN_SURVIVAL_PATIENTS = 30
MIN_SURVIVAL_EVENTS = 5

RANDOM_SEED = 20260525
np.random.seed(RANDOM_SEED)


# ============================================================
# 1. UTILITIES
# ============================================================

def log_message(msg):
    print(msg)
    with open(LOG_DIR / "run_log.txt", "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")


def sanitize_filename(x):
    x = str(x)
    x = re.sub(r"[^\w\-_\.]+", "_", x)
    return x[:180]


def clean_gene_symbol(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def normalize_tcga_barcode_to_patient(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip().replace(".", "-")
    return x[:12]


def parse_gene_list(x):
    if pd.isna(x):
        return []
    genes = str(x).split(";")
    genes = [clean_gene_symbol(g) for g in genes if clean_gene_symbol(g) != ""]
    return sorted(list(set(genes)))


def safe_neglog10_p(p):
    if p is None or pd.isna(p):
        return np.nan
    if p <= 0:
        return 300.0
    return -np.log10(p)


def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


# ============================================================
# 2. LOAD INPUTS
# ============================================================

def load_main_results():
    if not MAIN_RESULT_FILE.exists():
        raise FileNotFoundError(f"Main result file not found: {MAIN_RESULT_FILE}")

    df = pd.read_csv(MAIN_RESULT_FILE)

    required = [
        "cancer", "database", "bp_name", "N_defined", "N_matched",
        "matched_fraction", "observation_class", "D_survival",
        "matched_genes"
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in main result table: {missing}")

    df = df[df["database"].isin(DATABASES)].copy()
    df["N_matched"] = pd.to_numeric(df["N_matched"], errors="coerce")
    df["D_survival"] = pd.to_numeric(df["D_survival"], errors="coerce")
    df["high_D"] = df["D_survival"] >= D_THRESHOLD

    return df


def load_sizebin_random_optional():
    if SIZEBIN_RANDOM_FILE.exists():
        df = pd.read_csv(SIZEBIN_RANDOM_FILE)
        log_message(f"Loaded size-bin random mapped results: {SIZEBIN_RANDOM_FILE}")
        return df
    else:
        log_message("Size-bin random mapped result not found. Random-linked summaries will be skipped.")
        return None


# ============================================================
# 3. AIDO-h THRESHOLD SENSITIVITY
# ============================================================

def classify_by_cutoff(n_matched, cutoff):
    if pd.isna(n_matched):
        return np.nan
    if n_matched >= cutoff:
        return "observation_ready"
    if n_matched >= NEAR_UNOBSERVABLE_LT:
        return "low_resolution"
    return "near_unobservable"


def threshold_sensitivity(df, random_df=None):
    rows = []

    for cutoff in AIDOH_THRESHOLDS:
        tmp = df.copy()
        tmp["class_at_cutoff"] = tmp["N_matched"].apply(lambda x: classify_by_cutoff(x, cutoff))
        tmp["ready_at_cutoff"] = tmp["class_at_cutoff"] == "observation_ready"

        # If size-bin random exists, merge structured_favored from random_df
        tmp_random = None
        if random_df is not None and "structured_favored" in random_df.columns:
            cols = ["cancer", "database", "bp_name", "structured_favored"]
            tmp_random = random_df[cols].drop_duplicates(["cancer", "database", "bp_name"])
            tmp = tmp.merge(tmp_random, on=["cancer", "database", "bp_name"], how="left")

        for (db, cancer), sub in tmp.groupby(["database", "cancer"]):
            ready = sub[sub["ready_at_cutoff"]].copy()

            row = {
                "cutoff_N_matched": cutoff,
                "database": db,
                "cancer": cancer,
                "n_total": sub.shape[0],
                "n_observation_ready": ready.shape[0],
                "fraction_observation_ready": ready.shape[0] / sub.shape[0] if sub.shape[0] > 0 else np.nan,
                "n_low_resolution": int((sub["class_at_cutoff"] == "low_resolution").sum()),
                "n_near_unobservable": int((sub["class_at_cutoff"] == "near_unobservable").sum()),
                "median_D_all": sub["D_survival"].median(),
                "median_D_ready": ready["D_survival"].median(),
                "n_high_D_all": int((sub["D_survival"] >= D_THRESHOLD).sum()),
                "n_high_D_ready": int((ready["D_survival"] >= D_THRESHOLD).sum()),
                "fraction_high_D_ready": float((ready["D_survival"] >= D_THRESHOLD).mean()) if ready.shape[0] > 0 else np.nan
            }

            if "structured_favored" in sub.columns:
                ready_sf = ready.dropna(subset=["structured_favored"])
                row["structured_favored_fraction_ready"] = ready_sf["structured_favored"].mean() if ready_sf.shape[0] > 0 else np.nan
                row["n_structured_favored_ready"] = int(ready_sf["structured_favored"].sum()) if ready_sf.shape[0] > 0 else np.nan

            rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "Threshold_sensitivity_by_cancer_database.csv", index=False)

    overall = (
        out.groupby(["cutoff_N_matched", "database"])
        .agg(
            n_total=("n_total", "sum"),
            n_observation_ready=("n_observation_ready", "sum"),
            n_low_resolution=("n_low_resolution", "sum"),
            n_near_unobservable=("n_near_unobservable", "sum"),
            median_fraction_observation_ready=("fraction_observation_ready", "median"),
            median_D_ready=("median_D_ready", "median"),
            total_high_D_ready=("n_high_D_ready", "sum"),
            median_fraction_high_D_ready=("fraction_high_D_ready", "median"),
            median_structured_favored_fraction_ready=("structured_favored_fraction_ready", "median") if "structured_favored_fraction_ready" in out.columns else ("fraction_high_D_ready", "median")
        )
        .reset_index()
    )

    overall["fraction_observation_ready_global"] = (
        overall["n_observation_ready"] / overall["n_total"]
    )

    overall.to_csv(TABLE_DIR / "Threshold_sensitivity_overall_by_database.csv", index=False)

    # Figures
    for metric, filename, ylabel in [
        ("fraction_observation_ready_global", "Figure_threshold_fraction_ready_by_database.png", "Global fraction observation-ready"),
        ("median_D_ready", "Figure_threshold_median_D_ready_by_database.png", "Median D among ready terms"),
        ("median_fraction_high_D_ready", "Figure_threshold_highD_fraction_ready_by_database.png", "Median high-D fraction among ready terms")
    ]:
        plt.figure(figsize=(7, 4))
        for db in sorted(overall["database"].unique()):
            sub = overall[overall["database"] == db].sort_values("cutoff_N_matched")
            plt.plot(sub["cutoff_N_matched"], sub[metric], marker="o", label=db)
        plt.xlabel("AIDO-h matched-gene cutoff")
        plt.ylabel(ylabel)
        plt.title(ylabel + " vs AIDO-h cutoff")
        plt.legend()
        save_fig(FIG_DIR / filename)

    return out, overall


# ============================================================
# 4. TOP-D BEFORE VS AFTER AIDO-h FILTERING
# ============================================================

def topD_before_after_filtering(df):
    records = []

    for db in DATABASES:
        sub = df[df["database"] == db].copy()
        sub = sub.dropna(subset=["D_survival"])

        before = sub.sort_values("D_survival", ascending=False).head(50).copy()
        after = sub[sub["observation_class"] == "observation_ready"].sort_values("D_survival", ascending=False).head(50).copy()

        before["ranking_context"] = "before_AIDOh_filtering"
        after["ranking_context"] = "after_AIDOh_filtering"

        records.append(before)
        records.append(after)

    out = pd.concat(records, axis=0, ignore_index=True)
    out.to_csv(TABLE_DIR / "TopD_before_vs_after_AIDOh_filtering_by_database.csv", index=False)

    # Compact top 20 table for manuscript
    compact = out.groupby(["database", "ranking_context"]).head(20).copy()
    compact_cols = [
        "database", "ranking_context", "cancer", "bp_name",
        "N_defined", "N_matched", "matched_fraction",
        "observation_class", "D_survival", "p_value_survival"
    ]
    compact_cols = [c for c in compact_cols if c in compact.columns]
    compact[compact_cols].to_csv(TABLE_DIR / "Top20_before_vs_after_AIDOh_filtering_compact.csv", index=False)

    # Figure: top D distributions before vs after
    for db in DATABASES:
        plot_df = compact[compact["database"] == db].copy()
        if plot_df.shape[0] == 0:
            continue

        labels = []
        values = []
        colors_group = []

        for context in ["before_AIDOh_filtering", "after_AIDOh_filtering"]:
            sub = plot_df[plot_df["ranking_context"] == context].head(20)
            labels.extend([f"{context.replace('_AIDOh_filtering','')}_{i+1}" for i in range(sub.shape[0])])
            values.extend(sub["D_survival"].values)
            colors_group.extend([context] * sub.shape[0])

        plt.figure(figsize=(10, 4))
        plt.bar(range(len(values)), values)
        plt.axhline(D_THRESHOLD, linestyle="--", linewidth=1.5)
        plt.xticks(range(len(values)), labels, rotation=90, fontsize=7)
        plt.ylabel("D survival")
        plt.title(f"Top-D before vs after AIDO-h filtering: {db}")
        save_fig(FIG_DIR / f"Figure_topD_before_after_{sanitize_filename(db)}.png")

    return out


# ============================================================
# 5. DATABASE AND CANCER REGIME SUMMARIES
# ============================================================

def database_cancer_regime_summaries(df, random_df=None):
    tmp = df.copy()

    if random_df is not None and "structured_favored" in random_df.columns:
        cols = [
            "cancer", "database", "bp_name",
            "structured_favored", "D_excess_vs_random_median",
            "D_excess_vs_random_q95", "random_D_median", "random_D_q95"
        ]
        cols = [c for c in cols if c in random_df.columns]
        rand = random_df[cols].drop_duplicates(["cancer", "database", "bp_name"])
        tmp = tmp.merge(rand, on=["cancer", "database", "bp_name"], how="left")

    rows = []

    for (cancer, db), sub in tmp.groupby(["cancer", "database"]):
        ready = sub[sub["observation_class"] == "observation_ready"].copy()
        low = sub[sub["observation_class"] == "low_resolution"].copy()
        near = sub[sub["observation_class"] == "near_unobservable"].copy()

        row = {
            "cancer": cancer,
            "database": db,
            "n_total": sub.shape[0],
            "n_ready": ready.shape[0],
            "n_low_resolution": low.shape[0],
            "n_near_unobservable": near.shape[0],
            "fraction_ready": ready.shape[0] / sub.shape[0] if sub.shape[0] else np.nan,
            "fraction_low_resolution": low.shape[0] / sub.shape[0] if sub.shape[0] else np.nan,
            "fraction_near_unobservable": near.shape[0] / sub.shape[0] if sub.shape[0] else np.nan,
            "median_N_matched": sub["N_matched"].median(),
            "median_N_matched_ready": ready["N_matched"].median(),
            "median_D_ready": ready["D_survival"].median(),
            "max_D_ready": ready["D_survival"].max(),
            "n_high_D_ready": int((ready["D_survival"] >= D_THRESHOLD).sum()),
            "fraction_high_D_ready": float((ready["D_survival"] >= D_THRESHOLD).mean()) if ready.shape[0] else np.nan
        }

        if "structured_favored" in ready.columns:
            ready_rand = ready.dropna(subset=["structured_favored"])
            row["structured_favored_fraction_ready"] = ready_rand["structured_favored"].mean() if ready_rand.shape[0] else np.nan
            row["median_D_excess_vs_random_median_ready"] = ready_rand["D_excess_vs_random_median"].median() if "D_excess_vs_random_median" in ready_rand.columns else np.nan
            row["n_highD_and_structured_favored_ready"] = int(((ready_rand["D_survival"] >= D_THRESHOLD) & (ready_rand["structured_favored"])).sum())

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "Cancer_database_regime_summary.csv", index=False)

    db_summary = (
        out.groupby("database")
        .agg(
            n_cancers=("cancer", "nunique"),
            total_terms=("n_total", "sum"),
            total_ready=("n_ready", "sum"),
            total_low_resolution=("n_low_resolution", "sum"),
            total_near_unobservable=("n_near_unobservable", "sum"),
            median_fraction_ready=("fraction_ready", "median"),
            median_fraction_low_resolution=("fraction_low_resolution", "median"),
            median_D_ready=("median_D_ready", "median"),
            median_fraction_high_D_ready=("fraction_high_D_ready", "median"),
            median_structured_favored_fraction_ready=("structured_favored_fraction_ready", "median") if "structured_favored_fraction_ready" in out.columns else ("fraction_high_D_ready", "median")
        )
        .reset_index()
    )

    db_summary["global_fraction_ready"] = db_summary["total_ready"] / db_summary["total_terms"]
    db_summary["global_fraction_low_resolution"] = db_summary["total_low_resolution"] / db_summary["total_terms"]
    db_summary["global_fraction_near_unobservable"] = db_summary["total_near_unobservable"] / db_summary["total_terms"]

    db_summary.to_csv(TABLE_DIR / "Database_regime_summary.csv", index=False)

    # Heatmaps
    for metric in ["fraction_ready", "fraction_high_D_ready"]:
        pivot = out.pivot(index="cancer", columns="database", values=metric).sort_index()
        plt.figure(figsize=(8, max(5, 0.35 * pivot.shape[0])))
        plt.imshow(pivot.values, aspect="auto")
        plt.colorbar(label=metric)
        plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=30, ha="right")
        plt.yticks(range(len(pivot.index)), pivot.index)
        plt.title(metric.replace("_", " "))
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                val = pivot.values[i, j]
                if pd.notna(val):
                    plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
        save_fig(FIG_DIR / f"Figure_heatmap_{metric}.png")

    if "structured_favored_fraction_ready" in out.columns:
        metric = "structured_favored_fraction_ready"
        pivot = out.pivot(index="cancer", columns="database", values=metric).sort_index()
        plt.figure(figsize=(8, max(5, 0.35 * pivot.shape[0])))
        plt.imshow(pivot.values, aspect="auto")
        plt.colorbar(label=metric)
        plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=30, ha="right")
        plt.yticks(range(len(pivot.index)), pivot.index)
        plt.title("structured-favored fraction among ready terms")
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                val = pivot.values[i, j]
                if pd.notna(val):
                    plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
        save_fig(FIG_DIR / "Figure_heatmap_structured_favored_fraction_ready.png")

    return out, db_summary


# ============================================================
# 6. EVENT-COUNT SENSITIVITY
# ============================================================

def event_count_sensitivity(df, random_df=None):
    tmp = df[df["observation_class"] == "observation_ready"].copy()

    if random_df is not None and "structured_favored" in random_df.columns:
        cols = ["cancer", "database", "bp_name", "structured_favored"]
        rand = random_df[cols].drop_duplicates(["cancer", "database", "bp_name"])
        tmp = tmp.merge(rand, on=["cancer", "database", "bp_name"], how="left")

    rows = []

    for (cancer, db), sub in tmp.groupby(["cancer", "database"]):
        n_events = sub["n_events"].dropna().median() if "n_events" in sub.columns else np.nan
        n_patients = sub["n_patients_survival"].dropna().median() if "n_patients_survival" in sub.columns else np.nan

        row = {
            "cancer": cancer,
            "database": db,
            "n_events": n_events,
            "n_patients_survival": n_patients,
            "median_D": sub["D_survival"].median(),
            "max_D": sub["D_survival"].max(),
            "n_high_D": int((sub["D_survival"] >= D_THRESHOLD).sum()),
            "fraction_high_D": float((sub["D_survival"] >= D_THRESHOLD).mean()),
        }

        if "structured_favored" in sub.columns:
            rsub = sub.dropna(subset=["structured_favored"])
            row["structured_favored_fraction"] = rsub["structured_favored"].mean() if rsub.shape[0] else np.nan
            row["n_structured_favored"] = int(rsub["structured_favored"].sum()) if rsub.shape[0] else np.nan

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "Event_count_sensitivity_by_cancer_database.csv", index=False)

    corr_rows = []

    for db, sub in out.groupby("database"):
        for y in ["median_D", "max_D", "fraction_high_D", "structured_favored_fraction"]:
            if y not in sub.columns:
                continue
            use = sub.dropna(subset=["n_events", y])
            if use.shape[0] < 3:
                continue

            spearman = use[["n_events", y]].corr(method="spearman").iloc[0, 1]
            pearson = use[["n_events", y]].corr(method="pearson").iloc[0, 1]

            corr_rows.append({
                "database": db,
                "x": "n_events",
                "y": y,
                "n_cancers": use.shape[0],
                "spearman_r": spearman,
                "pearson_r": pearson
            })

    corr = pd.DataFrame(corr_rows)
    corr.to_csv(TABLE_DIR / "Event_count_correlation_summary.csv", index=False)

    # Figures
    for db in sorted(out["database"].unique()):
        sub = out[out["database"] == db].copy()
        for y in ["median_D", "fraction_high_D"]:
            plt.figure(figsize=(6, 4))
            plt.scatter(sub["n_events"], sub[y], alpha=0.8)
            for _, r in sub.iterrows():
                if pd.notna(r["n_events"]) and pd.notna(r[y]):
                    plt.text(r["n_events"], r[y], str(r["cancer"]), fontsize=7)
            plt.xlabel("Number of OS events")
            plt.ylabel(y)
            plt.title(f"Event-count sensitivity: {db}")
            save_fig(FIG_DIR / f"Figure_event_count_{y}_{sanitize_filename(db)}.png")

    return out, corr


# ============================================================
# 7. LOW-RESOLUTION CHARACTERIZATION
# ============================================================

def low_resolution_characterization(df):
    low = df[df["observation_class"].isin(["low_resolution", "near_unobservable"])].copy()
    low.to_csv(TABLE_DIR / "Low_resolution_and_near_unobservable_full_table.csv", index=False)

    rows = []

    for (db, cls), sub in low.groupby(["database", "observation_class"]):
        rows.append({
            "database": db,
            "observation_class": cls,
            "n_terms_cancer_observations": sub.shape[0],
            "n_unique_bp": sub["bp_name"].nunique(),
            "n_cancers": sub["cancer"].nunique(),
            "median_N_defined": sub["N_defined"].median(),
            "median_N_matched": sub["N_matched"].median(),
            "median_matched_fraction": sub["matched_fraction"].median(),
            "median_D_survival": sub["D_survival"].median(),
            "max_D_survival": sub["D_survival"].max(),
            "n_high_D": int((sub["D_survival"] >= D_THRESHOLD).sum())
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(TABLE_DIR / "Low_resolution_characterization_summary.csv", index=False)

    # Top high-D low-resolution examples
    top_low = low.sort_values("D_survival", ascending=False).head(200)
    top_low.to_csv(TABLE_DIR / "Top_highD_low_resolution_examples.csv", index=False)

    # Figures
    for db in sorted(low["database"].unique()):
        sub = low[low["database"] == db]
        if sub.shape[0] == 0:
            continue

        plt.figure(figsize=(7, 4))
        plt.hist(sub["N_matched"].dropna(), bins=range(0, 11))
        plt.xlabel("N matched")
        plt.ylabel("Number of BP-cancer observations")
        plt.title(f"Low-resolution / near-unobservable Nmatched: {db}")
        save_fig(FIG_DIR / f"Figure_lowres_Nmatched_distribution_{sanitize_filename(db)}.png")

        plt.figure(figsize=(7, 4))
        plt.scatter(sub["N_matched"], sub["D_survival"], alpha=0.3)
        plt.axhline(D_THRESHOLD, linestyle="--", linewidth=1.5)
        plt.xlabel("N matched")
        plt.ylabel("D survival")
        plt.title(f"Low-resolution D distribution: {db}")
        save_fig(FIG_DIR / f"Figure_lowres_D_vs_Nmatched_{sanitize_filename(db)}.png")

    return low, summary


# ============================================================
# 8. REDUNDANCY / JACCARD OVERLAP
# ============================================================

def jaccard(a, b):
    a = set(a)
    b = set(b)
    if len(a) == 0 and len(b) == 0:
        return np.nan
    inter = len(a.intersection(b))
    union = len(a.union(b))
    return inter / union if union > 0 else np.nan


def redundancy_analysis(df):
    rows = []
    pair_rows = []

    tmp = df[
        (df["database"].isin(REDUNDANCY_DATABASES)) &
        (df["observation_class"] == "observation_ready")
    ].copy()

    if REDUNDANCY_ONLY_HIGH_D:
        tmp = tmp[tmp["D_survival"] >= D_THRESHOLD].copy()

    for (cancer, db), sub in tmp.groupby(["cancer", "database"]):
        sub = sub.sort_values("D_survival", ascending=False).head(REDUNDANCY_TOP_N_PER_CANCER_DATABASE).copy()
        sub["gene_list"] = sub["matched_genes"].apply(parse_gene_list)

        n = sub.shape[0]
        if n < 2:
            continue

        vals = []
        names = sub["bp_name"].tolist()
        genes = sub["gene_list"].tolist()
        Ds = sub["D_survival"].tolist()

        for i in range(n):
            for j in range(i + 1, n):
                jac = jaccard(genes[i], genes[j])
                vals.append(jac)

                pair_rows.append({
                    "cancer": cancer,
                    "database": db,
                    "bp1": names[i],
                    "bp2": names[j],
                    "D1": Ds[i],
                    "D2": Ds[j],
                    "jaccard": jac
                })

        vals = np.array(vals, dtype=float)
        vals = vals[~np.isnan(vals)]

        rows.append({
            "cancer": cancer,
            "database": db,
            "n_terms": n,
            "n_pairs": len(vals),
            "median_jaccard": float(np.median(vals)) if len(vals) else np.nan,
            "mean_jaccard": float(np.mean(vals)) if len(vals) else np.nan,
            "q90_jaccard": float(np.quantile(vals, 0.90)) if len(vals) else np.nan,
            "fraction_pairs_jaccard_ge_0_25": float(np.mean(vals >= 0.25)) if len(vals) else np.nan,
            "fraction_pairs_jaccard_ge_0_50": float(np.mean(vals >= 0.50)) if len(vals) else np.nan
        })

    summary = pd.DataFrame(rows)
    pairs = pd.DataFrame(pair_rows)

    summary.to_csv(TABLE_DIR / "Redundancy_Jaccard_summary_top_terms.csv", index=False)
    pairs.to_csv(TABLE_DIR / "Redundancy_Jaccard_pairwise_top_terms.csv", index=False)

    # Figures
    for db in sorted(summary["database"].unique()) if summary.shape[0] else []:
        sub = summary[summary["database"] == db]
        plt.figure(figsize=(8, max(4, 0.35 * sub.shape[0])))
        plt.barh(sub["cancer"], sub["median_jaccard"])
        plt.xlabel("Median pairwise Jaccard")
        plt.ylabel("Cancer")
        plt.title(f"Redundancy among top high-D terms: {db}")
        save_fig(FIG_DIR / f"Figure_redundancy_median_jaccard_{sanitize_filename(db)}.png")

    return summary, pairs


# ============================================================
# 9. COMBINED COHORT SENSITIVITY
# ============================================================

def combined_cohort_sensitivity(df, random_df=None):
    tmp = df.copy()
    tmp["analysis_set"] = np.where(tmp["cancer"].isin(COMBINED_COHORTS), "combined_or_overlapping", "primary_noncombined")

    if random_df is not None and "structured_favored" in random_df.columns:
        cols = ["cancer", "database", "bp_name", "structured_favored", "D_excess_vs_random_median"]
        cols = [c for c in cols if c in random_df.columns]
        rand = random_df[cols].drop_duplicates(["cancer", "database", "bp_name"])
        tmp = tmp.merge(rand, on=["cancer", "database", "bp_name"], how="left")

    rows = []

    for (analysis_set, db), sub in tmp.groupby(["analysis_set", "database"]):
        ready = sub[sub["observation_class"] == "observation_ready"].copy()

        row = {
            "analysis_set": analysis_set,
            "database": db,
            "n_cancers": sub["cancer"].nunique(),
            "n_total": sub.shape[0],
            "n_ready": ready.shape[0],
            "fraction_ready": ready.shape[0] / sub.shape[0] if sub.shape[0] else np.nan,
            "median_D_ready": ready["D_survival"].median(),
            "max_D_ready": ready["D_survival"].max(),
            "n_high_D_ready": int((ready["D_survival"] >= D_THRESHOLD).sum()),
            "fraction_high_D_ready": float((ready["D_survival"] >= D_THRESHOLD).mean()) if ready.shape[0] else np.nan
        }

        if "structured_favored" in ready.columns:
            rr = ready.dropna(subset=["structured_favored"])
            row["structured_favored_fraction_ready"] = rr["structured_favored"].mean() if rr.shape[0] else np.nan
            row["median_D_excess_vs_random_median_ready"] = rr["D_excess_vs_random_median"].median() if "D_excess_vs_random_median" in rr.columns else np.nan

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "Combined_cohort_sensitivity_summary.csv", index=False)

    primary = tmp[tmp["analysis_set"] == "primary_noncombined"].copy()
    primary.to_csv(TABLE_DIR / "Main_results_excluding_combined_cohorts.csv", index=False)

    return out, primary


# ============================================================
# 10. EXACT RANDOM VALIDATION FOR TOP CANDIDATES
# ============================================================

def find_cancer_dir(cancer_code):
    for p in BASE_DATA_DIR.iterdir():
        if not p.is_dir():
            continue
        if f"({cancer_code})" in p.name:
            return p
        if p.name.upper() == cancer_code.upper():
            return p
    return None


def find_ge_file(cancer_dir):
    candidate = cancer_dir / "GE.tsv"
    if candidate.exists():
        return candidate

    patterns = [
        r"^GE\.tsv$",
        r"HiSeqV2.*\.tsv$",
        r".*expression.*\.tsv$",
        r".*gene.*expression.*\.tsv$"
    ]

    for f in cancer_dir.iterdir():
        if not f.is_file():
            continue
        for pat in patterns:
            if re.search(pat, f.name, flags=re.IGNORECASE):
                return f

    return None


def find_survival_or_clinical_file(cancer_dir, cancer_code):
    preferred_names = [
        f"TCGA-{cancer_code}.survival.tsv",
        f"TCGA-{cancer_code}.survival.txt",
        f"{cancer_code}.survival.tsv",
        f"{cancer_code}.survival.txt",
        "survival.tsv",
        "survival.txt",
        f"TCGA.{cancer_code}.sampleMap_{cancer_code}_clinicalMatrix",
        f"TCGA.{cancer_code}.sampleMap_{cancer_code}_clinicalMatrix.tsv",
        f"TCGA.{cancer_code}.sampleMap_{cancer_code}_clinicalMatrix.txt",
        "Phenotype.tsv",
        "Phenotype.txt",
        "phenotype.tsv",
        "phenotype.txt"
    ]

    for name in preferred_names:
        p = cancer_dir / name
        if p.exists():
            return p

    files = [f for f in cancer_dir.iterdir() if f.is_file()]

    for f in files:
        name = f.name.lower()
        if "survival" in name and f.suffix.lower() in ["", ".tsv", ".txt", ".csv"]:
            return f

    for f in files:
        name = f.name.lower()
        if "clinicalmatrix" in name or "clinical_matrix" in name:
            return f

    for f in files:
        name = f.name.lower()
        if "samplemap" in name and "clinical" in name:
            return f

    for f in files:
        name = f.name.lower()
        if "phenotype" in name and f.suffix.lower() in ["", ".tsv", ".txt", ".csv"]:
            return f

    return None


def read_table_flexible(path):
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path, sep=None, engine="python")


def find_column_by_candidates(df, candidates, exact_first=True):
    lower_map = {str(c).lower(): c for c in df.columns}

    if exact_first:
        for cand in candidates:
            if cand.lower() in lower_map:
                return lower_map[cand.lower()]

    for c in df.columns:
        cl = str(c).lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c

    return None


def parse_event_value(x):
    if pd.isna(x):
        return np.nan
    xs = str(x).strip().lower()

    dead_values = ["1", "1.0", "true", "yes", "dead", "deceased", "death", "event", "died"]
    alive_values = ["0", "0.0", "false", "no", "alive", "living", "censored", "not dead"]

    if xs in dead_values:
        return 1
    if xs in alive_values:
        return 0

    try:
        v = float(xs)
        if v == 1:
            return 1
        if v == 0:
            return 0
    except Exception:
        pass

    return np.nan


def load_ge_matrix(ge_path):
    df = pd.read_csv(ge_path, sep="\t", index_col=0)
    df.index = [clean_gene_symbol(x) for x in df.index]
    df = df.loc[df.index != ""]
    df = df.apply(pd.to_numeric, errors="coerce")

    if df.index.duplicated().any():
        df = df.groupby(df.index).mean(numeric_only=True)

    df.columns = [normalize_tcga_barcode_to_patient(c) for c in df.columns]

    if pd.Index(df.columns).duplicated().any():
        df = df.T.groupby(level=0).mean(numeric_only=True).T

    return df


def zscore_by_gene(ge):
    mu = ge.mean(axis=1, skipna=True)
    sd = ge.std(axis=1, skipna=True).replace(0, np.nan)
    return ge.sub(mu, axis=0).div(sd, axis=0)


def load_survival_from_any_clinical_file(path):
    df = read_table_flexible(path)

    id_candidates = [
        "sample", "Sample", "sample_id", "SampleID", "sampleID",
        "_PATIENT", "patient", "Patient", "patient_id", "PatientID",
        "bcr_patient_barcode", "submitter_id", "case_submitter_id"
    ]

    id_col = find_column_by_candidates(df, id_candidates)
    if id_col is None:
        id_col = df.columns[0]

    df["patient_id"] = df[id_col].apply(normalize_tcga_barcode_to_patient)

    time_col = find_column_by_candidates(
        df,
        [
            "OS_Time_nature2012", "OS.time", "OS_Time", "OS.time.days",
            "OS_days", "OS.time_months", "overall_survival",
            "overall_survival_time", "overall_survival_days"
        ]
    )

    event_col = find_column_by_candidates(
        df,
        ["OS_event_nature2012", "OS_event", "OS", "event", "death_event", "vital_status"]
    )

    death_time_col = find_column_by_candidates(
        df,
        ["days_to_death", "days.death", "days_to_death.diagnoses", "days_to_death.demographic"]
    )

    followup_col = find_column_by_candidates(
        df,
        [
            "days_to_last_followup", "days_to_last_follow_up",
            "days_to_last_known_alive", "days_to_last_followup.diagnoses",
            "days_to_last_follow_up.diagnoses",
            "days_to_last_known_alive.diagnoses"
        ]
    )

    vital_col = find_column_by_candidates(
        df,
        ["vital_status", "vital.status", "patient.vital_status", "demographic.vital_status"]
    )

    if time_col is not None:
        df["OS_time"] = pd.to_numeric(df[time_col], errors="coerce")

    elif death_time_col is not None or followup_col is not None:
        death_time = pd.to_numeric(df[death_time_col], errors="coerce") if death_time_col is not None else pd.Series(np.nan, index=df.index)
        follow_time = pd.to_numeric(df[followup_col], errors="coerce") if followup_col is not None else pd.Series(np.nan, index=df.index)
        df["OS_time"] = death_time.copy()
        missing = df["OS_time"].isna()
        df.loc[missing, "OS_time"] = follow_time[missing]

    else:
        raise ValueError(f"No usable OS time column found in {path}")

    if event_col is not None:
        df["OS_event"] = df[event_col].apply(parse_event_value)
    elif vital_col is not None:
        df["OS_event"] = df[vital_col].apply(parse_event_value)
    elif death_time_col is not None:
        death_time = pd.to_numeric(df[death_time_col], errors="coerce")
        df["OS_event"] = death_time.notna().astype(int)
    else:
        raise ValueError(f"No usable OS event column found in {path}")

    out = df[["patient_id", "OS_time", "OS_event"]].copy()
    out = out.dropna(subset=["patient_id", "OS_time", "OS_event"])
    out = out[out["OS_time"] > 0]
    out = out.sort_values("OS_time", ascending=False)
    out = out.drop_duplicates("patient_id", keep="first")
    out["OS_event"] = out["OS_event"].astype(int)

    return out


def compute_logrank_D(score, survival_df):
    try:
        from lifelines.statistics import logrank_test
    except Exception:
        raise ImportError("Please install lifelines: pip install lifelines")

    temp = survival_df.copy()
    temp["score"] = score.reindex(temp["patient_id"]).values
    temp = temp.dropna(subset=["score", "OS_time", "OS_event"])

    if temp.shape[0] < MIN_SURVIVAL_PATIENTS or temp["OS_event"].sum() < MIN_SURVIVAL_EVENTS:
        return np.nan

    median_val = temp["score"].median()
    low = temp[temp["score"] <= median_val]
    high = temp[temp["score"] > median_val]

    if low.shape[0] < 10 or high.shape[0] < 10:
        return np.nan

    try:
        result = logrank_test(
            low["OS_time"],
            high["OS_time"],
            event_observed_A=low["OS_event"],
            event_observed_B=high["OS_event"]
        )
        return safe_neglog10_p(result.p_value)
    except Exception:
        return np.nan


def select_exact_random_candidates(df):
    ready = df[
        (df["observation_class"] == "observation_ready") &
        (df["D_survival"].notna())
    ].copy()

    candidates = []

    # Global top per database
    for db, sub in ready.groupby("database"):
        candidates.append(sub.sort_values("D_survival", ascending=False).head(EXACT_TOP_N_PER_DATABASE_GLOBAL))

    # Top per cancer/database
    for (cancer, db), sub in ready.groupby(["cancer", "database"]):
        candidates.append(sub.sort_values("D_survival", ascending=False).head(EXACT_TOP_N_PER_CANCER_DATABASE))

    cand = pd.concat(candidates, axis=0, ignore_index=True)
    cand = cand.drop_duplicates(["cancer", "database", "bp_name"])
    cand = cand.sort_values("D_survival", ascending=False)

    if MAX_EXACT_RANDOM_TERMS is not None and cand.shape[0] > MAX_EXACT_RANDOM_TERMS:
        cand = cand.head(MAX_EXACT_RANDOM_TERMS)

    cand.to_csv(TABLE_DIR / "Exact_random_validation_selected_candidates.csv", index=False)
    return cand


def summarize_random(vals, d_real):
    vals = np.array(vals, dtype=float)
    vals = vals[~np.isnan(vals)]

    if vals.size == 0:
        return {
            "exact_random_D_median": np.nan,
            "exact_random_D_q95": np.nan,
            "exact_random_D_q99": np.nan,
            "exact_D_excess_vs_random_median": np.nan,
            "exact_D_excess_vs_random_q95": np.nan,
            "exact_frac_random_ge_real": np.nan,
            "exact_empirical_p_random": np.nan,
            "exact_structured_favored": np.nan,
            "exact_n_valid_random": 0
        }

    q95 = np.quantile(vals, 0.95)
    med = np.median(vals)

    return {
        "exact_random_D_median": float(med),
        "exact_random_D_q95": float(q95),
        "exact_random_D_q99": float(np.quantile(vals, 0.99)),
        "exact_D_excess_vs_random_median": float(d_real - med),
        "exact_D_excess_vs_random_q95": float(d_real - q95),
        "exact_frac_random_ge_real": float(np.mean(vals >= d_real)),
        "exact_empirical_p_random": float((np.sum(vals >= d_real) + 1) / (vals.size + 1)),
        "exact_structured_favored": bool(d_real > q95),
        "exact_n_valid_random": int(vals.size)
    }


def run_exact_random_validation(df):
    candidates = select_exact_random_candidates(df)

    log_message(f"Exact random validation candidates: {candidates.shape[0]}")
    log_message(f"EXACT_RANDOM_T = {EXACT_RANDOM_T}")

    all_records = []
    status_rows = []

    rng_global = np.random.default_rng(RANDOM_SEED)

    for cancer, cancer_terms in candidates.groupby("cancer"):
        start_cancer = time.time()

        cancer_dir = find_cancer_dir(cancer)
        if cancer_dir is None:
            log_message(f"Exact random skip {cancer}: cancer dir not found")
            continue

        ge_file = find_ge_file(cancer_dir)
        surv_file = find_survival_or_clinical_file(cancer_dir, cancer)

        if ge_file is None or surv_file is None:
            log_message(f"Exact random skip {cancer}: missing GE or survival")
            continue

        ge = load_ge_matrix(ge_file)
        survival_df = load_survival_from_any_clinical_file(surv_file)

        common = sorted(list(set(ge.columns).intersection(set(survival_df["patient_id"]))))
        ge = ge[common]
        survival_df = survival_df[survival_df["patient_id"].isin(common)].copy()
        ge_z = zscore_by_gene(ge)

        all_genes = np.array(ge_z.index.tolist())

        log_message("")
        log_message(f"Exact random cancer: {cancer}, terms={cancer_terms.shape[0]}, patients={ge_z.shape[1]}, events={int(survival_df['OS_event'].sum())}")

        for i, row in cancer_terms.reset_index(drop=True).iterrows():
            term_start = time.time()

            bp_name = row["bp_name"]
            db = row["database"]
            d_real = float(row["D_survival"])
            matched_genes = parse_gene_list(row["matched_genes"])
            matched_genes = [g for g in matched_genes if g in ge_z.index]
            n_matched = len(matched_genes)

            if n_matched < 10:
                continue

            matched_set = set(matched_genes)
            pool = np.array([g for g in all_genes if g not in matched_set])

            if len(pool) < n_matched:
                pool = all_genes

            random_D = []

            rng_seed = abs(hash((cancer, db, bp_name, RANDOM_SEED))) % (2**32)
            rng = np.random.default_rng(rng_seed)

            for t in range(EXACT_RANDOM_T):
                sampled = rng.choice(pool, size=n_matched, replace=False)
                score = ge_z.loc[sampled].mean(axis=0, skipna=True)
                random_D.append(compute_logrank_D(score, survival_df))

            rec = {
                "cancer": cancer,
                "database": db,
                "bp_name": bp_name,
                "N_matched": n_matched,
                "D_real": d_real,
                "T_exact_random": EXACT_RANDOM_T
            }
            rec.update(summarize_random(random_D, d_real))
            all_records.append(rec)

            if (i + 1) % 10 == 0 or (i + 1) == cancer_terms.shape[0]:
                elapsed = time.time() - term_start
                log_message(f"  {cancer}: {i+1}/{cancer_terms.shape[0]} terms, last={elapsed:.2f}s")

            # checkpoint
            pd.DataFrame(all_records).to_csv(
                EXACT_RANDOM_DIR / "Exact_random_validation_checkpoint.csv",
                index=False
            )

        elapsed_cancer = time.time() - start_cancer
        status_rows.append({
            "cancer": cancer,
            "n_terms": cancer_terms.shape[0],
            "seconds": elapsed_cancer,
            "minutes": elapsed_cancer / 60
        })

        pd.DataFrame(status_rows).to_csv(
            EXACT_RANDOM_DIR / "Exact_random_validation_status.csv",
            index=False
        )

    out = pd.DataFrame(all_records)
    out.to_csv(TABLE_DIR / "Exact_random_validation_T1000_results.csv", index=False)

    return out


# ============================================================
# 11. MAIN
# ============================================================

def main():
    start = time.time()

    log_file = LOG_DIR / "run_log.txt"
    if log_file.exists():
        log_file.unlink()

    log_message("============================================================")
    log_message("AIDO-h-Biology II | Reviewer Analysis Pack")
    log_message("============================================================")
    log_message(f"Main result file: {MAIN_RESULT_FILE}")
    log_message(f"Size-bin random file: {SIZEBIN_RANDOM_FILE}")
    log_message(f"Output: {OUT_DIR}")
    log_message("============================================================")

    df = load_main_results()
    random_df = load_sizebin_random_optional()

    log_message(f"Loaded main results: {df.shape}")
    log_message(f"Databases: {sorted(df['database'].unique())}")
    log_message(f"Cancers: {df['cancer'].nunique()}")

    # 1. Threshold sensitivity
    log_message("Running threshold sensitivity...")
    threshold_sensitivity(df, random_df=random_df)

    # 2. Top-D before vs after AIDO-h filtering
    log_message("Running top-D before/after filtering...")
    topD_before_after_filtering(df)

    # 3. Database and cancer regime summaries
    log_message("Running database/cancer regime summaries...")
    database_cancer_regime_summaries(df, random_df=random_df)

    # 4. Event-count sensitivity
    log_message("Running event-count sensitivity...")
    event_count_sensitivity(df, random_df=random_df)

    # 5. Low-resolution characterization
    log_message("Running low-resolution characterization...")
    low_resolution_characterization(df)

    # 6. Redundancy analysis
    log_message("Running redundancy / Jaccard analysis...")
    redundancy_analysis(df)

    # 7. Combined cohort sensitivity
    log_message("Running combined cohort sensitivity...")
    combined_cohort_sensitivity(df, random_df=random_df)

    # 8. Exact random validation
    if RUN_EXACT_RANDOM_VALIDATION:
        log_message("Running exact Nmatched random validation for top candidates...")
        run_exact_random_validation(df)
    else:
        log_message("Skipping exact random validation because RUN_EXACT_RANDOM_VALIDATION=False")

    elapsed = time.time() - start

    config = {
        "MAIN_RESULT_FILE": str(MAIN_RESULT_FILE),
        "SIZEBIN_RANDOM_FILE": str(SIZEBIN_RANDOM_FILE),
        "OUT_DIR": str(OUT_DIR),
        "DATABASES": DATABASES,
        "D_THRESHOLD": D_THRESHOLD,
        "AIDOH_THRESHOLDS": AIDOH_THRESHOLDS,
        "NEAR_UNOBSERVABLE_LT": NEAR_UNOBSERVABLE_LT,
        "COMBINED_COHORTS": COMBINED_COHORTS,
        "REDUNDANCY_TOP_N_PER_CANCER_DATABASE": REDUNDANCY_TOP_N_PER_CANCER_DATABASE,
        "REDUNDANCY_ONLY_HIGH_D": REDUNDANCY_ONLY_HIGH_D,
        "REDUNDANCY_DATABASES": REDUNDANCY_DATABASES,
        "RUN_EXACT_RANDOM_VALIDATION": RUN_EXACT_RANDOM_VALIDATION,
        "EXACT_RANDOM_T": EXACT_RANDOM_T,
        "EXACT_TOP_N_PER_DATABASE_GLOBAL": EXACT_TOP_N_PER_DATABASE_GLOBAL,
        "EXACT_TOP_N_PER_CANCER_DATABASE": EXACT_TOP_N_PER_CANCER_DATABASE,
        "MAX_EXACT_RANDOM_TERMS": MAX_EXACT_RANDOM_TERMS,
        "elapsed_seconds": elapsed,
        "elapsed_hours": elapsed / 3600,
        "note": (
            "This script generates reviewer-defense analyses for AIDO-h-Biology II, "
            "including threshold sensitivity, top-D filtering effect, event-count sensitivity, "
            "low-resolution characterization, redundancy checks, combined cohort sensitivity, "
            "and optional exact random validation for top candidates."
        )
    }

    with open(LOG_DIR / "run_config_reviewer_analysis_pack.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    log_message("============================================================")
    log_message("DONE.")
    log_message(f"Elapsed: {elapsed/3600:.2f} hours")
    log_message("Key outputs:")
    log_message(str(TABLE_DIR / "Threshold_sensitivity_overall_by_database.csv"))
    log_message(str(TABLE_DIR / "Top20_before_vs_after_AIDOh_filtering_compact.csv"))
    log_message(str(TABLE_DIR / "Cancer_database_regime_summary.csv"))
    log_message(str(TABLE_DIR / "Event_count_correlation_summary.csv"))
    log_message(str(TABLE_DIR / "Low_resolution_characterization_summary.csv"))
    log_message(str(TABLE_DIR / "Redundancy_Jaccard_summary_top_terms.csv"))
    log_message(str(TABLE_DIR / "Combined_cohort_sensitivity_summary.csv"))
    log_message(str(TABLE_DIR / "Exact_random_validation_T1000_results.csv"))
    log_message("============================================================")


if __name__ == "__main__":
    main()
