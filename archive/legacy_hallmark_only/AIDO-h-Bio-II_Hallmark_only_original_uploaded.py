# ============================================================
# AIDO-h-Biology II
# TRUE All-Cancer Hallmark AIDO-h + Optional Survival-D Pipeline
#
# Core principle:
#   GE.tsv is required.
#   Survival / clinicalMatrix is optional.
#   Every cancer with GE.tsv will be processed for AIDO-h.
#
# Supports clinical files named like:
#   TCGA.UCEC.sampleMap_UCEC_clinicalMatrix
#   TCGA.BRCA.sampleMap_BRCA_clinicalMatrix
#   TCGA.LUAD.sampleMap_LUAD_clinicalMatrix
#
# Input root:
#   D:/AIDO-Data/UCSC_XENA/
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-TRUE/
# ============================================================

import os
import re
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# ============================================================
# 0. CONFIG
# ============================================================

BASE_DIR = Path(r"D:/AIDO-Data/UCSC_XENA")
OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-TRUE")

TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
LOG_DIR = OUT_DIR / "logs"
SCORE_DIR = OUT_DIR / "bp_scores"
PER_CANCER_DIR = OUT_DIR / "per_cancer_results"

for d in [OUT_DIR, TABLE_DIR, FIG_DIR, LOG_DIR, SCORE_DIR, PER_CANCER_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Hallmark GMT fallback.
# If each cancer folder does not contain a GMT file, this one will be used.
FALLBACK_GMT = BASE_DIR / "Breast Cancer (BRCA)" / "h.all.v2026.1.Hs.symbols.gmt"

MIN_OBSERVATION_READY = 10
LOW_RES_MIN = 3

D_THRESHOLD = 1.301  # -log10(0.05)

MIN_SURVIVAL_PATIENTS = 30
MIN_SURVIVAL_EVENTS = 5

RANDOM_SEED = 20260523
np.random.seed(RANDOM_SEED)


# ============================================================
# 1. GENERAL UTILITIES
# ============================================================

def log_message(msg):
    print(msg)
    with open(LOG_DIR / "run_log.txt", "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")


def normalize_tcga_barcode_to_patient(x):
    """
    Convert TCGA sample barcode to patient barcode.
    Example:
        TCGA-AB-1234-01A -> TCGA-AB-1234
    """
    if pd.isna(x):
        return np.nan
    x = str(x).strip().replace(".", "-")
    return x[:12]


def clean_gene_symbol(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def infer_cancer_code(folder_name):
    """
    Extract cancer code from folder name.
    Example:
        Breast Cancer (BRCA) -> BRCA
        Uterine Corpus Endometrial Carcinoma (UCEC) -> UCEC
    """
    m = re.search(r"\(([^()]+)\)", folder_name)
    if m:
        return m.group(1).strip()
    return folder_name.replace(" ", "_").replace("-", "_")


def safe_neglog10_p(p):
    if p is None or pd.isna(p):
        return np.nan
    if p <= 0:
        return 300.0
    return -np.log10(p)


def classify_observation_readiness(n_matched):
    if n_matched >= MIN_OBSERVATION_READY:
        return "observation_ready"
    elif n_matched >= LOW_RES_MIN:
        return "low_resolution"
    else:
        return "near_unobservable"


# ============================================================
# 2. FILE DISCOVERY
# ============================================================

def find_ge_file(cancer_dir):
    """
    Find GE matrix.
    Preferred file name:
        GE.tsv
    Fallback:
        HiSeqV2*.tsv
        *expression*.tsv
        *gene*expression*.tsv
    """
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


def find_gmt_file(cancer_dir):
    """
    Find Hallmark GMT file.
    First search local cancer folder.
    If absent, use FALLBACK_GMT.
    """
    for f in cancer_dir.iterdir():
        if f.is_file() and f.suffix.lower() == ".gmt":
            name = f.name.lower()
            if "h.all" in name or "hallmark" in name:
                return f

    if FALLBACK_GMT.exists():
        return FALLBACK_GMT

    return None


def find_survival_or_clinical_file(cancer_dir, cancer_code):
    """
    Survival / clinical file is optional.

    Search priority:
    1. TCGA-CODE.survival.tsv / txt
    2. any file containing survival
    3. TCGA.CODE.sampleMap_CODE_clinicalMatrix
    4. any file containing clinicalMatrix
    5. Phenotype.tsv / Phenotype.txt
    """

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

    # survival-like files first
    for f in files:
        name = f.name.lower()
        if "survival" in name and f.suffix.lower() in ["", ".tsv", ".txt", ".csv"]:
            return f

    # Xena clinicalMatrix files
    for f in files:
        name = f.name.lower()
        if "clinicalmatrix" in name or "clinical_matrix" in name:
            return f

    # Xena sampleMap clinical files
    for f in files:
        name = f.name.lower()
        if "samplemap" in name and "clinical" in name:
            return f

    # phenotype fallback
    for f in files:
        name = f.name.lower()
        if "phenotype" in name and f.suffix.lower() in ["", ".tsv", ".txt", ".csv"]:
            return f

    return None


# ============================================================
# 3. TABLE READERS
# ============================================================

def read_table_flexible(path):
    """
    Read txt / tsv / csv / extension-less UCSC Xena clinicalMatrix files.
    Most Xena clinicalMatrix files are tab-separated even without extension.
    """
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path, sep=None, engine="python")


def find_column_by_candidates(df, candidates, exact_first=True):
    """
    Flexible column finder.
    First exact lower-case match, then substring match.
    """
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


# ============================================================
# 4. LOAD GMT
# ============================================================

def load_gmt(gmt_path):
    """
    Load GMT file:
        gene_set_name, description, gene1, gene2, ...
    """
    gene_sets = {}

    with open(gmt_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")

            if len(parts) < 3:
                continue

            name = parts[0].strip()
            genes = [clean_gene_symbol(g) for g in parts[2:] if clean_gene_symbol(g) != ""]
            genes = sorted(list(set(genes)))

            gene_sets[name] = genes

    return gene_sets


# ============================================================
# 5. LOAD GE MATRIX
# ============================================================

def load_ge_matrix(ge_path):
    """
    Expected GE.tsv format:
        rows = genes
        columns = samples

    First column should be gene symbol or gene ID.
    """
    df = pd.read_csv(ge_path, sep="\t", index_col=0)

    # Clean gene symbols
    df.index = [clean_gene_symbol(x) for x in df.index]
    df = df.loc[df.index != ""]

    # Convert to numeric
    df = df.apply(pd.to_numeric, errors="coerce")

    # Collapse duplicated genes
    if df.index.duplicated().any():
        df = df.groupby(df.index).mean(numeric_only=True)

    # Normalize TCGA sample barcodes to patient IDs
    df.columns = [normalize_tcga_barcode_to_patient(c) for c in df.columns]

    # Collapse duplicated patient columns
    if pd.Index(df.columns).duplicated().any():
        df = df.T.groupby(level=0).mean(numeric_only=True).T

    return df


# ============================================================
# 6. FLEXIBLE SURVIVAL / CLINICALMATRIX PARSER
# ============================================================

def parse_event_value(x):
    """
    Convert event / vital status to OS_event:
        dead / event = 1
        alive / censored = 0
    """
    if pd.isna(x):
        return np.nan

    xs = str(x).strip().lower()

    dead_values = [
        "1", "1.0", "true", "yes",
        "dead", "deceased", "death", "event", "died"
    ]

    alive_values = [
        "0", "0.0", "false", "no",
        "alive", "living", "censored", "not dead"
    ]

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


def load_survival_from_any_clinical_file(path):
    """
    Parse OS_time and OS_event from:
        - standard TCGA survival tables
        - UCSC Xena clinicalMatrix
        - phenotype tables

    Returns:
        survival_df, used_time_col, used_event_col
    """

    df = read_table_flexible(path)

    # --------------------------------------------------------
    # 1. Patient ID
    # --------------------------------------------------------
    id_candidates = [
        "sample",
        "Sample",
        "sample_id",
        "SampleID",
        "sampleID",
        "_PATIENT",
        "patient",
        "Patient",
        "patient_id",
        "PatientID",
        "bcr_patient_barcode",
        "submitter_id",
        "case_submitter_id"
    ]

    id_col = find_column_by_candidates(df, id_candidates)

    if id_col is None:
        id_col = df.columns[0]

    df["patient_id"] = df[id_col].apply(normalize_tcga_barcode_to_patient)

    # --------------------------------------------------------
    # 2. Direct OS time columns
    # --------------------------------------------------------
    time_col = find_column_by_candidates(
        df,
        [
            "OS_Time_nature2012",
            "OS.time",
            "OS_Time",
            "OS.time.days",
            "OS_days",
            "OS.time_months",
            "overall_survival",
            "overall_survival_time",
            "overall_survival_days"
        ]
    )

    # --------------------------------------------------------
    # 3. Direct OS event columns
    # --------------------------------------------------------
    event_col = find_column_by_candidates(
        df,
        [
            "OS_event_nature2012",
            "OS_event",
            "OS",
            "event",
            "death_event",
            "vital_status"
        ]
    )

    # --------------------------------------------------------
    # 4. UCSC Xena clinicalMatrix fallback columns
    # --------------------------------------------------------
    death_time_col = find_column_by_candidates(
        df,
        [
            "days_to_death",
            "days.death",
            "days_to_death.diagnoses",
            "days_to_death.demographic"
        ]
    )

    followup_col = find_column_by_candidates(
        df,
        [
            "days_to_last_followup",
            "days_to_last_follow_up",
            "days_to_last_known_alive",
            "days_to_last_followup.diagnoses",
            "days_to_last_follow_up.diagnoses",
            "days_to_last_known_alive.diagnoses"
        ]
    )

    vital_col = find_column_by_candidates(
        df,
        [
            "vital_status",
            "vital.status",
            "patient.vital_status",
            "demographic.vital_status"
        ]
    )

    # --------------------------------------------------------
    # 5. Build OS_time
    # --------------------------------------------------------
    if time_col is not None:
        df["OS_time"] = pd.to_numeric(df[time_col], errors="coerce")
        used_time_col = time_col

    elif death_time_col is not None or followup_col is not None:
        if death_time_col is not None:
            death_time = pd.to_numeric(df[death_time_col], errors="coerce")
        else:
            death_time = pd.Series(np.nan, index=df.index)

        if followup_col is not None:
            follow_time = pd.to_numeric(df[followup_col], errors="coerce")
        else:
            follow_time = pd.Series(np.nan, index=df.index)

        df["OS_time"] = death_time.copy()

        missing = df["OS_time"].isna()
        df.loc[missing, "OS_time"] = follow_time[missing]

        used_time_col = f"{death_time_col} + {followup_col}"

    else:
        raise ValueError(
            f"No usable OS time column found in {path}. "
            f"Columns={list(df.columns)}"
        )

    # --------------------------------------------------------
    # 6. Build OS_event
    # --------------------------------------------------------
    if event_col is not None:
        df["OS_event"] = df[event_col].apply(parse_event_value)
        used_event_col = event_col

    elif vital_col is not None:
        df["OS_event"] = df[vital_col].apply(parse_event_value)
        used_event_col = vital_col

    elif death_time_col is not None:
        death_time = pd.to_numeric(df[death_time_col], errors="coerce")
        df["OS_event"] = death_time.notna().astype(int)
        used_event_col = f"inferred_from_{death_time_col}"

    else:
        raise ValueError(
            f"No usable OS event column found in {path}. "
            f"Columns={list(df.columns)}"
        )

    # --------------------------------------------------------
    # 7. Clean survival table
    # --------------------------------------------------------
    out = df[["patient_id", "OS_time", "OS_event"]].copy()

    out = out.dropna(subset=["patient_id", "OS_time", "OS_event"])
    out = out[out["OS_time"] > 0]

    # One row per patient.
    # If duplicated, keep longest follow-up.
    out = out.sort_values("OS_time", ascending=False)
    out = out.drop_duplicates("patient_id", keep="first")

    out["OS_event"] = out["OS_event"].astype(int)

    return out, used_time_col, used_event_col


# ============================================================
# 7. AIDO-h LITE METRICS
# ============================================================

def compute_aido_h_lite_metrics(gene_sets, ge_genes, cancer_code, database_name="Hallmark"):
    records = []
    ge_gene_set = set(ge_genes)

    for bp_name, genes in gene_sets.items():
        defined_genes = sorted(list(set([clean_gene_symbol(g) for g in genes])))
        matched_genes = sorted(list(set(defined_genes).intersection(ge_gene_set)))

        n_defined = len(defined_genes)
        n_matched = len(matched_genes)

        records.append({
            "cancer": cancer_code,
            "database": database_name,
            "bp_name": bp_name,
            "N_defined": n_defined,
            "N_matched": n_matched,
            "matched_fraction": n_matched / n_defined if n_defined > 0 else np.nan,
            "observation_class": classify_observation_readiness(n_matched),
            "matched_genes": ";".join(matched_genes)
        })

    return pd.DataFrame(records)


# ============================================================
# 8. BP SCORE CONSTRUCTION
# ============================================================

def zscore_by_gene(ge):
    """
    Z-score each gene across patients.
    """
    mu = ge.mean(axis=1, skipna=True)
    sd = ge.std(axis=1, skipna=True).replace(0, np.nan)
    z = ge.sub(mu, axis=0).div(sd, axis=0)
    return z


def construct_bp_scores(ge_z, h_metrics):
    """
    BP score = mean z-scored expression across matched genes.

    Construct score if N_matched >= LOW_RES_MIN.
    Observation-ready classification still uses MIN_OBSERVATION_READY.
    """
    score_records = []
    score_matrix = {}

    for _, row in h_metrics.iterrows():
        bp = row["bp_name"]

        matched_genes = str(row["matched_genes"]).split(";")
        matched_genes = [g for g in matched_genes if g in ge_z.index]

        if len(matched_genes) < LOW_RES_MIN:
            continue

        bp_score = ge_z.loc[matched_genes].mean(axis=0, skipna=True)
        score_matrix[bp] = bp_score

        score_records.append({
            "cancer": row["cancer"],
            "database": row["database"],
            "bp_name": bp,
            "N_defined": row["N_defined"],
            "N_matched": row["N_matched"],
            "matched_fraction": row["matched_fraction"],
            "observation_class": row["observation_class"],
            "score_mean": float(bp_score.mean(skipna=True)),
            "score_sd": float(bp_score.std(skipna=True)),
            "score_iqr": float(bp_score.quantile(0.75) - bp_score.quantile(0.25)),
            "score_missing_fraction": float(bp_score.isna().mean())
        })

    score_df = pd.DataFrame(score_records)

    if len(score_matrix) > 0:
        score_mat = pd.DataFrame(score_matrix)
        score_mat.index.name = "patient_id"
    else:
        score_mat = pd.DataFrame()

    return score_df, score_mat


# ============================================================
# 9. SURVIVAL-D
# ============================================================

def compute_survival_discriminability(bp_scores, survival_df, cancer_code):
    """
    Survival D:
        median split BP score
        log-rank test
        D_survival = -log10(p)
    """
    try:
        from lifelines.statistics import logrank_test
    except Exception:
        log_message("WARNING: lifelines not installed. Survival D skipped.")
        log_message("Install with: pip install lifelines")
        return pd.DataFrame()

    merged = survival_df.merge(
        bp_scores.reset_index(),
        on="patient_id",
        how="inner"
    )

    records = []

    for bp in bp_scores.columns:
        temp = merged[["patient_id", "OS_time", "OS_event", bp]].dropna()

        if temp.shape[0] < MIN_SURVIVAL_PATIENTS or temp["OS_event"].sum() < MIN_SURVIVAL_EVENTS:
            records.append({
                "cancer": cancer_code,
                "database": "Hallmark",
                "bp_name": bp,
                "n_patients_survival": temp.shape[0],
                "n_events": int(temp["OS_event"].sum()) if temp.shape[0] > 0 else 0,
                "grouping_method": "median_split",
                "median_score": np.nan,
                "p_value_survival": np.nan,
                "D_survival": np.nan,
                "event_rate_low": np.nan,
                "event_rate_high": np.nan,
                "event_rate_ratio_high_vs_low": np.nan
            })
            continue

        median_val = temp[bp].median()

        low = temp[temp[bp] <= median_val]
        high = temp[temp[bp] > median_val]

        if low.shape[0] < 10 or high.shape[0] < 10:
            p = np.nan
            d = np.nan
            event_rate_low = np.nan
            event_rate_high = np.nan
            event_rate_ratio = np.nan

        else:
            try:
                result = logrank_test(
                    low["OS_time"],
                    high["OS_time"],
                    event_observed_A=low["OS_event"],
                    event_observed_B=high["OS_event"]
                )

                p = result.p_value
                d = safe_neglog10_p(p)

                event_rate_low = low["OS_event"].sum() / low.shape[0]
                event_rate_high = high["OS_event"].sum() / high.shape[0]
                event_rate_ratio = (
                    event_rate_high / event_rate_low
                    if event_rate_low > 0 else np.nan
                )

            except Exception:
                p = np.nan
                d = np.nan
                event_rate_low = np.nan
                event_rate_high = np.nan
                event_rate_ratio = np.nan

        records.append({
            "cancer": cancer_code,
            "database": "Hallmark",
            "bp_name": bp,
            "n_patients_survival": temp.shape[0],
            "n_events": int(temp["OS_event"].sum()),
            "grouping_method": "median_split",
            "median_score": median_val,
            "p_value_survival": p,
            "D_survival": d,
            "event_rate_low": event_rate_low,
            "event_rate_high": event_rate_high,
            "event_rate_ratio_high_vs_low": event_rate_ratio
        })

    return pd.DataFrame(records)


# ============================================================
# 10. MERGE TABLES
# ============================================================

def merge_result_tables(h_metrics, score_metrics, survival_d):
    out = h_metrics.copy()

    if score_metrics is not None and score_metrics.shape[0] > 0:
        drop_cols = [
            "cancer", "database", "N_defined", "N_matched",
            "matched_fraction", "observation_class"
        ]

        out = out.merge(
            score_metrics.drop(columns=drop_cols, errors="ignore"),
            on="bp_name",
            how="left"
        )

    if survival_d is not None and survival_d.shape[0] > 0:
        out = out.merge(
            survival_d.drop(columns=["database"], errors="ignore"),
            on=["cancer", "bp_name"],
            how="left"
        )
    else:
        out["n_patients_survival"] = np.nan
        out["n_events"] = np.nan
        out["grouping_method"] = np.nan
        out["median_score"] = np.nan
        out["p_value_survival"] = np.nan
        out["D_survival"] = np.nan
        out["event_rate_low"] = np.nan
        out["event_rate_high"] = np.nan
        out["event_rate_ratio_high_vs_low"] = np.nan

    return out


# ============================================================
# 11. SUMMARY TABLES
# ============================================================

def make_cancer_summary(all_results):
    rows = []

    for (cancer, db), sub in all_results.groupby(["cancer", "database"]):
        n_total = sub.shape[0]
        n_ready = int((sub["observation_class"] == "observation_ready").sum())
        n_low = int((sub["observation_class"] == "low_resolution").sum())
        n_near = int((sub["observation_class"] == "near_unobservable").sum())

        rows.append({
            "cancer": cancer,
            "database": db,
            "n_gene_sets_defined": n_total,
            "n_gene_sets_with_matched_genes": int((sub["N_matched"] > 0).sum()),
            "n_observation_ready": n_ready,
            "n_low_resolution": n_low,
            "n_near_unobservable": n_near,
            "fraction_observation_ready": n_ready / n_total if n_total > 0 else np.nan,
            "fraction_low_resolution": n_low / n_total if n_total > 0 else np.nan,
            "fraction_near_unobservable": n_near / n_total if n_total > 0 else np.nan,
            "median_N_defined": sub["N_defined"].median(),
            "median_N_matched": sub["N_matched"].median(),
            "median_matched_fraction": sub["matched_fraction"].median(),
            "median_score_iqr": sub["score_iqr"].median() if "score_iqr" in sub.columns else np.nan,
            "max_D_survival": sub["D_survival"].max() if "D_survival" in sub.columns else np.nan,
            "median_D_survival": sub["D_survival"].median() if "D_survival" in sub.columns else np.nan,
            "n_high_D_survival": int((sub["D_survival"] >= D_THRESHOLD).sum()) if "D_survival" in sub.columns else 0
        })

    return pd.DataFrame(rows)


def make_overall_summary(all_results):
    rows = []

    for db, sub in all_results.groupby("database"):
        n_total = sub.shape[0]
        n_ready = int((sub["observation_class"] == "observation_ready").sum())
        n_low = int((sub["observation_class"] == "low_resolution").sum())
        n_near = int((sub["observation_class"] == "near_unobservable").sum())

        rows.append({
            "database": db,
            "n_cancers": sub["cancer"].nunique(),
            "n_bp_cancer_observations": n_total,
            "n_observation_ready": n_ready,
            "n_low_resolution": n_low,
            "n_near_unobservable": n_near,
            "fraction_observation_ready": n_ready / n_total if n_total > 0 else np.nan,
            "fraction_low_resolution": n_low / n_total if n_total > 0 else np.nan,
            "fraction_near_unobservable": n_near / n_total if n_total > 0 else np.nan,
            "median_N_matched": sub["N_matched"].median(),
            "median_matched_fraction": sub["matched_fraction"].median(),
            "median_D_survival": sub["D_survival"].median() if "D_survival" in sub.columns else np.nan,
            "max_D_survival": sub["D_survival"].max() if "D_survival" in sub.columns else np.nan,
            "n_high_D_survival": int((sub["D_survival"] >= D_THRESHOLD).sum()) if "D_survival" in sub.columns else 0
        })

    return pd.DataFrame(rows)


def make_before_after_filter_summary(all_results):
    rows = []

    for (cancer, db), sub in all_results.groupby(["cancer", "database"]):
        before = sub.copy()
        after = sub[sub["observation_class"] == "observation_ready"].copy()

        row = {
            "cancer": cancer,
            "database": db,
            "n_bp_before_filter": before.shape[0],
            "n_bp_after_filter": after.shape[0],
            "fraction_retained": after.shape[0] / before.shape[0] if before.shape[0] > 0 else np.nan,
            "D_survival_median_before": before["D_survival"].median() if "D_survival" in before.columns else np.nan,
            "D_survival_median_after": after["D_survival"].median() if "D_survival" in after.columns else np.nan,
            "D_survival_max_before": before["D_survival"].max() if "D_survival" in before.columns else np.nan,
            "D_survival_max_after": after["D_survival"].max() if "D_survival" in after.columns else np.nan,
            "D_survival_n_high_before": int((before["D_survival"] >= D_THRESHOLD).sum()) if "D_survival" in before.columns else 0,
            "D_survival_n_high_after": int((after["D_survival"] >= D_THRESHOLD).sum()) if "D_survival" in after.columns else 0
        }

        rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# 12. FIGURES
# ============================================================

def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def figure_readiness_by_cancer(summary):
    plot_df = summary.sort_values("fraction_observation_ready", ascending=True)

    plt.figure(figsize=(8, max(5, 0.35 * plot_df.shape[0])))
    plt.barh(plot_df["cancer"], plot_df["fraction_observation_ready"])
    plt.xlim(0, 1.05)
    plt.xlabel("Fraction observation-ready")
    plt.ylabel("Cancer")
    plt.title("AIDO-h observation-readiness across cancers")
    save_fig(FIG_DIR / "Figure_2A_fraction_observation_ready_by_cancer.png")


def figure_median_matched_by_cancer(summary):
    plot_df = summary.sort_values("median_N_matched", ascending=True)

    plt.figure(figsize=(8, max(5, 0.35 * plot_df.shape[0])))
    plt.barh(plot_df["cancer"], plot_df["median_N_matched"])
    plt.axvline(MIN_OBSERVATION_READY, linestyle="--", linewidth=1.5)
    plt.xlabel("Median matched genes per Hallmark BP")
    plt.ylabel("Cancer")
    plt.title("Effective Hallmark representation across cancers")
    save_fig(FIG_DIR / "Figure_2B_median_matched_genes_by_cancer.png")


def figure_Dmax_by_cancer(summary):
    plot_df = summary.dropna(subset=["max_D_survival"]).copy()

    if plot_df.shape[0] == 0:
        return

    plot_df = plot_df.sort_values("max_D_survival", ascending=True)

    plt.figure(figsize=(8, max(5, 0.35 * plot_df.shape[0])))
    plt.barh(plot_df["cancer"], plot_df["max_D_survival"])
    plt.axvline(D_THRESHOLD, linestyle="--", linewidth=1.5)
    plt.xlabel("Maximum survival D")
    plt.ylabel("Cancer")
    plt.title("Maximum Hallmark survival discriminability by cancer")
    save_fig(FIG_DIR / "Figure_3A_max_D_survival_by_cancer.png")


def figure_highD_count_by_cancer(summary):
    plot_df = summary.copy()
    plot_df = plot_df.sort_values("n_high_D_survival", ascending=True)

    plt.figure(figsize=(8, max(5, 0.35 * plot_df.shape[0])))
    plt.barh(plot_df["cancer"], plot_df["n_high_D_survival"])
    plt.xlabel(f"Number of BP with D >= {D_THRESHOLD}")
    plt.ylabel("Cancer")
    plt.title("Number of survival-discriminative Hallmark BP by cancer")
    save_fig(FIG_DIR / "Figure_3B_highD_BP_count_by_cancer.png")


def figure_D_distribution(all_results):
    vals = all_results["D_survival"].dropna()

    if vals.shape[0] == 0:
        return

    plt.figure(figsize=(7, 4))
    plt.hist(vals, bins=40)
    plt.axvline(D_THRESHOLD, linestyle="--", linewidth=1.5)
    plt.xlabel("Survival D")
    plt.ylabel("Number of BP-cancer observations")
    plt.title("Distribution of survival discriminability")
    save_fig(FIG_DIR / "Figure_4_D_survival_distribution.png")


def figure_readiness_vs_D(all_results):
    plot_df = all_results.dropna(subset=["N_matched", "D_survival"]).copy()

    if plot_df.shape[0] == 0:
        return

    plt.figure(figsize=(7, 5))
    plt.scatter(plot_df["N_matched"], plot_df["D_survival"], alpha=0.5)
    plt.axvline(MIN_OBSERVATION_READY, linestyle="--", linewidth=1.5)
    plt.axhline(D_THRESHOLD, linestyle="--", linewidth=1.5)
    plt.xlabel("Matched genes per BP")
    plt.ylabel("Survival D")
    plt.title("Observation readiness vs survival discriminability")
    save_fig(FIG_DIR / "Figure_5_readiness_vs_D_survival.png")


def figure_top_bp_global(all_results, top_n=30):
    plot_df = all_results.dropna(subset=["D_survival"]).copy()

    if plot_df.shape[0] == 0:
        return

    plot_df = plot_df.sort_values("D_survival", ascending=False).head(top_n)

    labels = (
        plot_df["cancer"].astype(str)
        + " | "
        + plot_df["bp_name"].str.replace("HALLMARK_", "", regex=False)
    )

    plt.figure(figsize=(10, max(6, 0.30 * plot_df.shape[0])))
    plt.barh(labels[::-1], plot_df["D_survival"].values[::-1])
    plt.axvline(D_THRESHOLD, linestyle="--", linewidth=1.5)
    plt.xlabel("Survival D")
    plt.ylabel("Cancer | Hallmark BP")
    plt.title(f"Top {top_n} BP-cancer observations by survival D")
    save_fig(FIG_DIR / f"Figure_top_{top_n}_BP_cancer_by_D_survival.png")


# ============================================================
# 13. SINGLE-CANCER PIPELINE
# ============================================================

def run_one_cancer(cancer_dir):
    folder_name = cancer_dir.name
    cancer_code = infer_cancer_code(folder_name)

    log_message("")
    log_message("------------------------------------------------------------")
    log_message(f"Processing {folder_name} [{cancer_code}]")
    log_message("------------------------------------------------------------")

    ge_file = find_ge_file(cancer_dir)

    if ge_file is None:
        log_message(f"SKIP {cancer_code}: no GE file found.")
        return None, {
            "cancer": cancer_code,
            "folder": str(cancer_dir),
            "status": "skipped_no_GE"
        }

    gmt_file = find_gmt_file(cancer_dir)

    if gmt_file is None:
        log_message(f"SKIP {cancer_code}: no GMT file and no fallback GMT.")
        return None, {
            "cancer": cancer_code,
            "folder": str(cancer_dir),
            "status": "skipped_no_GMT"
        }

    survival_file = find_survival_or_clinical_file(cancer_dir, cancer_code)

    try:
        log_message(f"GE file: {ge_file}")
        log_message(f"GMT file: {gmt_file}")

        if survival_file is not None:
            log_message(f"Optional survival/clinical file: {survival_file}")
        else:
            log_message("Optional survival/clinical file: NOT FOUND")

        # Load inputs
        gene_sets = load_gmt(gmt_file)
        ge = load_ge_matrix(ge_file)

        log_message(f"GE loaded: genes={ge.shape[0]}, patients={ge.shape[1]}")
        log_message(f"Gene sets loaded: {len(gene_sets)}")

        # AIDO-h always runs
        h_metrics = compute_aido_h_lite_metrics(
            gene_sets=gene_sets,
            ge_genes=ge.index,
            cancer_code=cancer_code,
            database_name="Hallmark"
        )

        # BP score construction
        ge_z = zscore_by_gene(ge)
        score_metrics, bp_scores = construct_bp_scores(ge_z, h_metrics)

        if bp_scores.shape[1] > 0:
            bp_scores.to_csv(SCORE_DIR / f"{cancer_code}_Hallmark_BP_scores.csv")

        # Optional survival-D
        survival_d = pd.DataFrame()
        survival_status = "not_available"

        n_survival_patients = np.nan
        n_common_patients = np.nan
        n_events = np.nan
        survival_time_col = ""
        survival_event_col = ""

        if survival_file is not None:
            try:
                survival_df, survival_time_col, survival_event_col = load_survival_from_any_clinical_file(survival_file)

                common_patients = sorted(
                    list(set(bp_scores.index).intersection(set(survival_df["patient_id"])))
                )

                n_common_patients = len(common_patients)
                n_survival_patients = survival_df.shape[0]
                n_events = int(survival_df["OS_event"].sum())

                log_message(
                    f"Survival loaded: patients={n_survival_patients}, "
                    f"events={n_events}, common={n_common_patients}, "
                    f"time_col={survival_time_col}, event_col={survival_event_col}"
                )

                if (
                    n_common_patients >= MIN_SURVIVAL_PATIENTS
                    and n_events >= MIN_SURVIVAL_EVENTS
                    and bp_scores.shape[1] > 0
                ):
                    survival_df = survival_df[survival_df["patient_id"].isin(common_patients)].copy()

                    survival_d = compute_survival_discriminability(
                        bp_scores=bp_scores,
                        survival_df=survival_df,
                        cancer_code=cancer_code
                    )

                    survival_status = "computed"
                else:
                    survival_status = "insufficient_common_patients_or_events"

            except Exception as e:
                log_message(f"WARNING {cancer_code}: survival parsing failed: {repr(e)}")
                survival_status = f"failed: {repr(e)}"

        # Merge result
        result = merge_result_tables(h_metrics, score_metrics, survival_d)

        result.to_csv(
            PER_CANCER_DIR / f"{cancer_code}_AIDOh_Biology_II_Hallmark_full_results.csv",
            index=False
        )

        status = {
            "cancer": cancer_code,
            "folder": str(cancer_dir),
            "status": "success_AIDOh",
            "survival_status": survival_status,
            "n_GE_genes": ge.shape[0],
            "n_GE_patients": ge.shape[1],
            "n_survival_patients": n_survival_patients,
            "n_common_patients": n_common_patients,
            "n_events": n_events,
            "n_gene_sets": len(gene_sets),
            "n_BP_scores": bp_scores.shape[1],
            "n_observation_ready": int((result["observation_class"] == "observation_ready").sum()),
            "n_low_resolution": int((result["observation_class"] == "low_resolution").sum()),
            "n_near_unobservable": int((result["observation_class"] == "near_unobservable").sum()),
            "max_D_survival": result["D_survival"].max() if "D_survival" in result.columns else np.nan,
            "median_D_survival": result["D_survival"].median() if "D_survival" in result.columns else np.nan,
            "n_high_D_survival": int((result["D_survival"] >= D_THRESHOLD).sum()) if "D_survival" in result.columns else 0,
            "GE_file": str(ge_file),
            "GMT_file": str(gmt_file),
            "survival_or_clinical_file": str(survival_file) if survival_file is not None else "",
            "survival_time_col": survival_time_col,
            "survival_event_col": survival_event_col
        }

        log_message(
            f"DONE {cancer_code}: "
            f"AIDO-h success, ready={status['n_observation_ready']}, "
            f"low={status['n_low_resolution']}, "
            f"near={status['n_near_unobservable']}, "
            f"survival_status={survival_status}, "
            f"maxD={status['max_D_survival']}"
        )

        return result, status

    except Exception as e:
        log_message(f"ERROR {cancer_code}: {repr(e)}")
        return None, {
            "cancer": cancer_code,
            "folder": str(cancer_dir),
            "status": "error",
            "error": repr(e)
        }


# ============================================================
# 14. MAIN
# ============================================================

def main():
    log_file = LOG_DIR / "run_log.txt"

    if log_file.exists():
        log_file.unlink()

    log_message("============================================================")
    log_message("AIDO-h-Biology II | TRUE All-Cancer GE-first Pipeline")
    log_message("With automatic UCSC Xena clinicalMatrix survival parsing")
    log_message("============================================================")
    log_message(f"Input root: {BASE_DIR}")
    log_message(f"Output root: {OUT_DIR}")
    log_message(f"Fallback GMT: {FALLBACK_GMT}")
    log_message("============================================================")

    if not BASE_DIR.exists():
        raise FileNotFoundError(f"BASE_DIR not found: {BASE_DIR}")

    cancer_dirs = sorted(
        [p for p in BASE_DIR.iterdir() if p.is_dir()],
        key=lambda x: x.name
    )

    all_results = []
    status_records = []

    for cancer_dir in cancer_dirs:
        result, status = run_one_cancer(cancer_dir)
        status_records.append(status)

        if result is not None and result.shape[0] > 0:
            all_results.append(result)

    status_df = pd.DataFrame(status_records)
    status_df.to_csv(TABLE_DIR / "Run_status_all_cancers_TRUE.csv", index=False)

    if len(all_results) == 0:
        log_message("No successful AIDO-h results.")
        return

    all_results = pd.concat(all_results, axis=0, ignore_index=True)

    # Full tables
    all_results.to_csv(
        TABLE_DIR / "Table_Main_AIDOh_Biology_II_AllCancers_TRUE_Hallmark_full_results.csv",
        index=False
    )

    cancer_summary = make_cancer_summary(all_results)
    cancer_summary.to_csv(
        TABLE_DIR / "Table_1_cancer_level_observation_readiness_summary_TRUE.csv",
        index=False
    )

    overall_summary = make_overall_summary(all_results)
    overall_summary.to_csv(
        TABLE_DIR / "Table_2_overall_observation_readiness_summary_TRUE.csv",
        index=False
    )

    before_after = make_before_after_filter_summary(all_results)
    before_after.to_csv(
        TABLE_DIR / "Table_3_before_after_AIDOh_filter_summary_TRUE.csv",
        index=False
    )

    top_d = all_results.sort_values(
        "D_survival",
        ascending=False,
        na_position="last"
    )

    top_d.to_csv(
        TABLE_DIR / "Table_4_top_BP_cancer_by_D_survival_TRUE.csv",
        index=False
    )

    lowres = all_results[all_results["observation_class"] != "observation_ready"].copy()
    lowres.to_csv(
        TABLE_DIR / "Table_S_low_resolution_and_near_unobservable_terms_TRUE.csv",
        index=False
    )

    # Figures
    log_message("Generating figures...")

    figure_readiness_by_cancer(cancer_summary)
    figure_median_matched_by_cancer(cancer_summary)
    figure_Dmax_by_cancer(cancer_summary)
    figure_highD_count_by_cancer(cancer_summary)
    figure_D_distribution(all_results)
    figure_readiness_vs_D(all_results)
    figure_top_bp_global(all_results, top_n=30)

    # Config
    config = {
        "BASE_DIR": str(BASE_DIR),
        "OUT_DIR": str(OUT_DIR),
        "FALLBACK_GMT": str(FALLBACK_GMT),
        "MIN_OBSERVATION_READY": MIN_OBSERVATION_READY,
        "LOW_RES_MIN": LOW_RES_MIN,
        "D_THRESHOLD": D_THRESHOLD,
        "MIN_SURVIVAL_PATIENTS": MIN_SURVIVAL_PATIENTS,
        "MIN_SURVIVAL_EVENTS": MIN_SURVIVAL_EVENTS,
        "RANDOM_SEED": RANDOM_SEED,
        "important_note": (
            "GE.tsv is required. Survival/clinicalMatrix is optional. "
            "All cancers with GE.tsv are processed for AIDO-h. "
            "UCSC Xena clinicalMatrix files are parsed directly without renaming."
        )
    }

    with open(LOG_DIR / "run_config_TRUE.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    log_message("============================================================")
    log_message("DONE.")
    log_message(f"Total cancer folders scanned: {len(cancer_dirs)}")
    log_message(f"Successful AIDO-h cancers: {sum(status_df['status'] == 'success_AIDOh')}")
    log_message(f"Total BP-cancer observations: {all_results.shape[0]}")
    log_message(f"Output directory: {OUT_DIR}")
    log_message("Main files:")
    log_message(str(TABLE_DIR / "Run_status_all_cancers_TRUE.csv"))
    log_message(str(TABLE_DIR / "Table_Main_AIDOh_Biology_II_AllCancers_TRUE_Hallmark_full_results.csv"))
    log_message(str(TABLE_DIR / "Table_1_cancer_level_observation_readiness_summary_TRUE.csv"))
    log_message(str(TABLE_DIR / "Table_4_top_BP_cancer_by_D_survival_TRUE.csv"))
    log_message("============================================================")


if __name__ == "__main__":
    main()
