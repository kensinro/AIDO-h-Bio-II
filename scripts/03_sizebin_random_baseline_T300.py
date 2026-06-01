# ============================================================
# AIDO-h-Biology II
# Global Size-Bin-Matched Random Baseline
#
# Purpose:
#   Fast global random baseline after AIDO-h MultiDB analysis.
#
# Key idea:
#   Instead of exact random baseline for every BP term:
#       each cancer  x  database  x  N_matched size-bin gets one
#       random baseline distribution.
#
# This reduces runtime massively while preserving size-aware control.
#
# Input:
#   D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/
#       Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300/
#
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

MULTIDB_RESULT_FILE = Path(
    r"D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/"
    r"Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv"
)

OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300")

TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
LOG_DIR = OUT_DIR / "logs"
CHECKPOINT_DIR = OUT_DIR / "checkpoints"

for d in [OUT_DIR, TABLE_DIR, FIG_DIR, LOG_DIR, CHECKPOINT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

T_RANDOM = 300

D_THRESHOLD = 1.301
OBSERVATION_CLASS_REQUIRED = "observation_ready"

DATABASES_TO_RUN = ["Hallmark", "GO_BP", "Reactome"]

# Use None for all cancers.
CANCERS_TO_RUN = None
# CANCERS_TO_RUN = ["BRCA", "KIRC", "SKCM"]

MIN_SURVIVAL_PATIENTS = 30
MIN_SURVIVAL_EVENTS = 5

RANDOM_SEED = 20260525

# Size bins based on effective matched genes.
# label, lower inclusive, upper inclusive; upper=None means infinity.
SIZE_BINS = [
    ("10_14", 10, 14),
    ("15_19", 15, 19),
    ("20_29", 20, 29),
    ("30_49", 30, 49),
    ("50_99", 50, 99),
    ("100_199", 100, 199),
    ("200_499", 200, 499),
    ("500_plus", 500, None),
]

# Random set size within each bin.
# "median" uses median N_matched of real terms in that cancer/database/bin.
# "lower" uses bin lower bound.
# "random_within_bin" randomly samples size within observed N_matched values.
RANDOM_SIZE_MODE = "median"

# Runtime test option. Set None for full run.
MAX_BINS_PER_CANCER_DATABASE = None
# MAX_BINS_PER_CANCER_DATABASE = 2


# ============================================================
# 1. UTILITIES
# ============================================================

def log_message(msg):
    print(msg)
    with open(LOG_DIR / "run_log.txt", "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")


def normalize_tcga_barcode_to_patient(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip().replace(".", "-")
    return x[:12]


def clean_gene_symbol(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def infer_cancer_code(folder_name):
    m = re.search(r"\(([^()]+)\)", folder_name)
    if m:
        return m.group(1).strip()
    return folder_name.replace(" ", "_").replace("-", "_")


def sanitize_filename(x):
    x = str(x)
    x = re.sub(r"[^\w\-_\.]+", "_", x)
    return x[:180]


def safe_neglog10_p(p):
    if p is None or pd.isna(p):
        return np.nan
    if p <= 0:
        return 300.0
    return -np.log10(p)


def assign_size_bin(n):
    if pd.isna(n):
        return np.nan

    n = int(n)

    for label, lo, hi in SIZE_BINS:
        if hi is None:
            if n >= lo:
                return label
        else:
            if lo <= n <= hi:
                return label

    return np.nan


def get_bin_bounds(label):
    for lab, lo, hi in SIZE_BINS:
        if lab == label:
            return lo, hi
    return None, None


# ============================================================
# 2. FILE DISCOVERY
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


# ============================================================
# 3. DATA LOADERS
# ============================================================

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
        [
            "OS_event_nature2012", "OS_event", "OS",
            "event", "death_event", "vital_status"
        ]
    )

    death_time_col = find_column_by_candidates(
        df,
        [
            "days_to_death", "days.death", "days_to_death.diagnoses",
            "days_to_death.demographic"
        ]
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
        [
            "vital_status", "vital.status", "patient.vital_status",
            "demographic.vital_status"
        ]
    )

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
        raise ValueError(f"No usable OS time column found in {path}")

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
        raise ValueError(f"No usable OS event column found in {path}")

    out = df[["patient_id", "OS_time", "OS_event"]].copy()
    out = out.dropna(subset=["patient_id", "OS_time", "OS_event"])
    out = out[out["OS_time"] > 0]
    out = out.sort_values("OS_time", ascending=False)
    out = out.drop_duplicates("patient_id", keep="first")
    out["OS_event"] = out["OS_event"].astype(int)

    return out, used_time_col, used_event_col


# ============================================================
# 4. SURVIVAL D
# ============================================================

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


# ============================================================
# 5. LOAD TERMS AND CREATE SIZE BINS
# ============================================================

def load_and_prepare_terms():
    if not MULTIDB_RESULT_FILE.exists():
        raise FileNotFoundError(f"Input table not found: {MULTIDB_RESULT_FILE}")

    df = pd.read_csv(MULTIDB_RESULT_FILE)

    required_cols = [
        "cancer", "database", "bp_name", "N_matched",
        "observation_class", "D_survival"
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if len(missing) > 0:
        raise ValueError(f"Missing required columns: {missing}")

    df = df[df["database"].isin(DATABASES_TO_RUN)].copy()

    if CANCERS_TO_RUN is not None:
        df = df[df["cancer"].isin(CANCERS_TO_RUN)].copy()

    df = df[df["observation_class"] == OBSERVATION_CLASS_REQUIRED].copy()
    df = df.dropna(subset=["N_matched", "D_survival"]).copy()

    df["N_matched"] = df["N_matched"].astype(int)
    df["size_bin"] = df["N_matched"].apply(assign_size_bin)

    df = df.dropna(subset=["size_bin"]).copy()

    df.to_csv(TABLE_DIR / "Selected_observation_ready_terms_with_size_bins.csv", index=False)

    return df


def get_random_size_for_bin(sub_terms, size_bin, rng):
    lo, hi = get_bin_bounds(size_bin)

    if RANDOM_SIZE_MODE == "lower":
        return int(lo)

    if RANDOM_SIZE_MODE == "median":
        return int(np.round(sub_terms["N_matched"].median()))

    if RANDOM_SIZE_MODE == "random_within_bin":
        observed_sizes = sub_terms["N_matched"].dropna().astype(int).values
        if len(observed_sizes) == 0:
            return int(lo)
        return int(rng.choice(observed_sizes))

    return int(np.round(sub_terms["N_matched"].median()))


# ============================================================
# 6. SIZE-BIN RANDOM BASELINE CORE
# ============================================================

def summarize_random_D(random_D):
    vals = np.array(random_D, dtype=float)
    vals = vals[~np.isnan(vals)]

    if vals.size == 0:
        return {
            "random_D_mean": np.nan,
            "random_D_median": np.nan,
            "random_D_q90": np.nan,
            "random_D_q95": np.nan,
            "random_D_q99": np.nan,
            "random_D_max": np.nan,
            "n_valid_random": 0
        }

    return {
        "random_D_mean": float(np.mean(vals)),
        "random_D_median": float(np.median(vals)),
        "random_D_q90": float(np.quantile(vals, 0.90)),
        "random_D_q95": float(np.quantile(vals, 0.95)),
        "random_D_q99": float(np.quantile(vals, 0.99)),
        "random_D_max": float(np.max(vals)),
        "n_valid_random": int(vals.size)
    }


def run_random_for_one_size_bin(
    cancer,
    database,
    size_bin,
    sub_terms,
    ge_z,
    survival_df,
    rng
):
    n_random_size = get_random_size_for_bin(sub_terms, size_bin, rng)

    all_genes = np.array(ge_z.index.tolist())

    if n_random_size < 3:
        return None

    if n_random_size > len(all_genes):
        return None

    random_D = []

    for t in range(T_RANDOM):
        sampled_genes = rng.choice(all_genes, size=n_random_size, replace=False)
        random_score = ge_z.loc[sampled_genes].mean(axis=0, skipna=True)
        d_rand = compute_logrank_D(random_score, survival_df)
        random_D.append(d_rand)

    summary = summarize_random_D(random_D)

    lo, hi = get_bin_bounds(size_bin)

    record = {
        "cancer": cancer,
        "database": database,
        "size_bin": size_bin,
        "bin_lower": lo,
        "bin_upper": hi if hi is not None else "inf",
        "random_size_used": n_random_size,
        "T_random": T_RANDOM,
        "n_real_terms_in_bin": sub_terms.shape[0],
        "median_N_matched_real_terms": float(sub_terms["N_matched"].median()),
        "min_N_matched_real_terms": int(sub_terms["N_matched"].min()),
        "max_N_matched_real_terms": int(sub_terms["N_matched"].max())
    }

    record.update(summary)

    return record


def checkpoint_file_for_cancer_database(cancer, database):
    return CHECKPOINT_DIR / f"{sanitize_filename(cancer)}_{sanitize_filename(database)}_sizebin_random_T{T_RANDOM}.csv"


def run_one_cancer_database(cancer, database, sub_terms):
    checkpoint_file = checkpoint_file_for_cancer_database(cancer, database)

    if checkpoint_file.exists():
        log_message(f"SKIP completed {cancer} | {database}: {checkpoint_file}")
        existing = pd.read_csv(checkpoint_file)
        return existing, {
            "cancer": cancer,
            "database": database,
            "status": "skipped_completed",
            "n_bins": existing.shape[0],
            "checkpoint_file": str(checkpoint_file)
        }

    cancer_dir = find_cancer_dir(cancer)

    if cancer_dir is None:
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "error_no_cancer_dir"
        }

    ge_file = find_ge_file(cancer_dir)
    survival_file = find_survival_or_clinical_file(cancer_dir, cancer)

    if ge_file is None:
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "error_no_GE"
        }

    if survival_file is None:
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "error_no_survival"
        }

    log_message("")
    log_message("------------------------------------------------------------")
    log_message(f"Size-bin random baseline: {cancer} | {database}")
    log_message("------------------------------------------------------------")
    log_message(f"Real terms: {sub_terms.shape[0]}")
    log_message(f"GE file: {ge_file}")
    log_message(f"Survival file: {survival_file}")

    start = time.time()

    ge = load_ge_matrix(ge_file)
    survival_df, used_time_col, used_event_col = load_survival_from_any_clinical_file(survival_file)

    common = sorted(list(set(ge.columns).intersection(set(survival_df["patient_id"]))))
    ge = ge[common]
    survival_df = survival_df[survival_df["patient_id"].isin(common)].copy()

    ge_z = zscore_by_gene(ge)

    n_events = int(survival_df["OS_event"].sum())

    log_message(
        f"{cancer} | {database}: genes={ge_z.shape[0]}, "
        f"patients={ge_z.shape[1]}, survival={survival_df.shape[0]}, events={n_events}"
    )

    if survival_df.shape[0] < MIN_SURVIVAL_PATIENTS or n_events < MIN_SURVIVAL_EVENTS:
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "skipped_insufficient_survival",
            "n_patients": survival_df.shape[0],
            "n_events": n_events
        }

    rng_seed = abs(hash((cancer, database, RANDOM_SEED))) % (2**32)
    rng = np.random.default_rng(rng_seed)

    records = []

    bins = sorted(sub_terms["size_bin"].dropna().unique().tolist())

    if MAX_BINS_PER_CANCER_DATABASE is not None:
        bins = bins[:MAX_BINS_PER_CANCER_DATABASE]

    for i, size_bin in enumerate(bins):
        bin_terms = sub_terms[sub_terms["size_bin"] == size_bin].copy()

        bin_start = time.time()

        rec = run_random_for_one_size_bin(
            cancer=cancer,
            database=database,
            size_bin=size_bin,
            sub_terms=bin_terms,
            ge_z=ge_z,
            survival_df=survival_df,
            rng=rng
        )

        if rec is not None:
            records.append(rec)

        elapsed = time.time() - bin_start

        log_message(
            f"{cancer} | {database}: bin {i+1}/{len(bins)} "
            f"{size_bin}, terms={bin_terms.shape[0]}, "
            f"time={elapsed:.2f}s"
        )

        pd.DataFrame(records).to_csv(checkpoint_file, index=False)

    result = pd.DataFrame(records)
    result.to_csv(checkpoint_file, index=False)

    elapsed_total = time.time() - start

    status = {
        "cancer": cancer,
        "database": database,
        "status": "completed",
        "n_bins": result.shape[0],
        "n_real_terms": sub_terms.shape[0],
        "n_patients": ge_z.shape[1],
        "n_events": n_events,
        "seconds": elapsed_total,
        "minutes": elapsed_total / 60,
        "checkpoint_file": str(checkpoint_file),
        "survival_time_col": used_time_col,
        "survival_event_col": used_event_col
    }

    log_message(
        f"DONE {cancer} | {database}: bins={result.shape[0]}, "
        f"time={elapsed_total/60:.2f} min"
    )

    return result, status


# ============================================================
# 7. MAP SIZE-BIN RANDOM BASELINE BACK TO REAL TERMS
# ============================================================

def map_random_baseline_to_terms(real_terms, bin_random_df):
    merged = real_terms.merge(
        bin_random_df,
        on=["cancer", "database", "size_bin"],
        how="left",
        suffixes=("", "_randombin")
    )

    merged["D_real"] = merged["D_survival"]

    merged["D_excess_vs_random_median"] = (
        merged["D_real"] - merged["random_D_median"]
    )

    merged["D_excess_vs_random_q95"] = (
        merged["D_real"] - merged["random_D_q95"]
    )

    merged["structured_favored"] = (
        merged["D_real"] > merged["random_D_q95"]
    )

    # Since this is size-bin baseline, not exact per-term baseline,
    # we estimate whether real D exceeds random distribution using q95 and excess.
    merged["high_D_real"] = merged["D_real"] >= D_THRESHOLD

    return merged


# ============================================================
# 8. SUMMARY TABLES
# ============================================================

def make_summary_tables(mapped_df, bin_random_df):
    bin_random_df.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_bin_level_baselines.csv", index=False)
    mapped_df.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_term_level_mapped_results.csv", index=False)

    db_rows = []

    for db, sub in mapped_df.groupby("database"):
        db_rows.append({
            "database": db,
            "n_terms": sub.shape[0],
            "n_cancers": sub["cancer"].nunique(),
            "n_size_bins_used": sub[["cancer", "database", "size_bin"]].drop_duplicates().shape[0],
            "median_D_real": sub["D_real"].median(),
            "median_random_D_median": sub["random_D_median"].median(),
            "median_random_D_q95": sub["random_D_q95"].median(),
            "median_D_excess_vs_random_median": sub["D_excess_vs_random_median"].median(),
            "median_D_excess_vs_random_q95": sub["D_excess_vs_random_q95"].median(),
            "structured_favored_fraction": sub["structured_favored"].mean(),
            "n_structured_favored": int(sub["structured_favored"].sum()),
            "n_high_D_real": int(sub["high_D_real"].sum()),
            "n_high_D_and_structured_favored": int((sub["high_D_real"] & sub["structured_favored"]).sum())
        })

    db_summary = pd.DataFrame(db_rows)
    db_summary.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_database_summary.csv", index=False)

    cancer_db_rows = []

    for (cancer, db), sub in mapped_df.groupby(["cancer", "database"]):
        cancer_db_rows.append({
            "cancer": cancer,
            "database": db,
            "n_terms": sub.shape[0],
            "n_size_bins_used": sub["size_bin"].nunique(),
            "median_D_real": sub["D_real"].median(),
            "median_random_D_median": sub["random_D_median"].median(),
            "median_random_D_q95": sub["random_D_q95"].median(),
            "median_D_excess_vs_random_median": sub["D_excess_vs_random_median"].median(),
            "median_D_excess_vs_random_q95": sub["D_excess_vs_random_q95"].median(),
            "structured_favored_fraction": sub["structured_favored"].mean(),
            "n_structured_favored": int(sub["structured_favored"].sum()),
            "n_high_D_real": int(sub["high_D_real"].sum()),
            "n_high_D_and_structured_favored": int((sub["high_D_real"] & sub["structured_favored"]).sum())
        })

    cancer_db_summary = pd.DataFrame(cancer_db_rows)
    cancer_db_summary.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_cancer_database_summary.csv", index=False)

    top_structured = mapped_df.sort_values(
        ["structured_favored", "D_excess_vs_random_q95", "D_real"],
        ascending=[False, False, False]
    )

    top_structured.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_top_structured_favored_terms.csv", index=False)

    highD_not_favored = mapped_df[
        (mapped_df["high_D_real"]) &
        (~mapped_df["structured_favored"])
    ].copy()

    highD_not_favored = highD_not_favored.sort_values("D_real", ascending=False)
    highD_not_favored.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_highD_but_not_above_random_q95.csv", index=False)

    return db_summary, cancer_db_summary


# ============================================================
# 9. FIGURES
# ============================================================

def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def figure_database_structured_fraction(db_summary):
    plot_df = db_summary.sort_values("structured_favored_fraction", ascending=True)

    plt.figure(figsize=(7, 4))
    plt.barh(plot_df["database"], plot_df["structured_favored_fraction"])
    plt.xlim(0, 1.05)
    plt.xlabel("Fraction with D_real > size-bin random D q95")
    plt.ylabel("Database")
    plt.title("Structured BP advantage over size-bin random baseline")
    save_fig(FIG_DIR / "Figure_SizeBinRandom_structured_favored_fraction_by_database.png")


def figure_database_D_excess(db_summary):
    plot_df = db_summary.sort_values("median_D_excess_vs_random_median", ascending=True)

    plt.figure(figsize=(7, 4))
    plt.barh(plot_df["database"], plot_df["median_D_excess_vs_random_median"])
    plt.axvline(0, linestyle="--", linewidth=1.5)
    plt.xlabel("Median D_real - median size-bin random D")
    plt.ylabel("Database")
    plt.title("Excess survival discriminability over size-bin random baseline")
    save_fig(FIG_DIR / "Figure_SizeBinRandom_D_excess_by_database.png")


def figure_D_real_vs_random_q95(mapped_df):
    for db in sorted(mapped_df["database"].dropna().unique()):
        sub = mapped_df[mapped_df["database"] == db].copy()

        if sub.shape[0] == 0:
            continue

        plt.figure(figsize=(6, 6))
        plt.scatter(sub["random_D_q95"], sub["D_real"], alpha=0.25)
        max_val = np.nanmax([sub["random_D_q95"].max(), sub["D_real"].max()])
        plt.plot([0, max_val], [0, max_val], linestyle="--", linewidth=1.5)
        plt.xlabel("Size-bin random D q95")
        plt.ylabel("D_real")
        plt.title(f"D_real vs size-bin random q95: {db}")
        save_fig(FIG_DIR / f"Figure_SizeBinRandom_Dreal_vs_random_q95_{sanitize_filename(db)}.png")


def figure_heatmap(cancer_db_summary, value_col, filename, title):
    pivot = cancer_db_summary.pivot(index="cancer", columns="database", values=value_col)
    pivot = pivot.sort_index()

    plt.figure(figsize=(8, max(5, 0.35 * pivot.shape[0])))
    plt.imshow(pivot.values, aspect="auto")
    plt.colorbar(label=value_col)
    plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=30, ha="right")
    plt.yticks(range(len(pivot.index)), pivot.index)
    plt.title(title)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if pd.notna(val):
                plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)

    save_fig(FIG_DIR / filename)


def make_figures(mapped_df, db_summary, cancer_db_summary):
    figure_database_structured_fraction(db_summary)
    figure_database_D_excess(db_summary)
    figure_D_real_vs_random_q95(mapped_df)

    figure_heatmap(
        cancer_db_summary,
        value_col="structured_favored_fraction",
        filename="Figure_SizeBinRandom_structured_favored_fraction_heatmap.png",
        title="Fraction of terms exceeding size-bin random q95"
    )

    figure_heatmap(
        cancer_db_summary,
        value_col="median_D_excess_vs_random_median",
        filename="Figure_SizeBinRandom_D_excess_heatmap.png",
        title="Median D excess over size-bin random baseline"
    )


# ============================================================
# 10. MAIN
# ============================================================

def main():
    log_file = LOG_DIR / "run_log.txt"

    if log_file.exists():
        log_file.unlink()

    log_message("============================================================")
    log_message("AIDO-h-Biology II | Size-Bin Random Baseline T=300")
    log_message("Fast global screening after AIDO-h filtering")
    log_message("============================================================")
    log_message(f"Input: {MULTIDB_RESULT_FILE}")
    log_message(f"Output: {OUT_DIR}")
    log_message(f"T_RANDOM: {T_RANDOM}")
    log_message(f"RANDOM_SIZE_MODE: {RANDOM_SIZE_MODE}")
    log_message("============================================================")

    start_global = time.time()

    real_terms = load_and_prepare_terms()

    log_message(f"Selected observation-ready terms: {real_terms.shape[0]}")
    log_message("Terms by cancer/database:")
    log_message(real_terms.groupby(["cancer", "database"]).size().to_string())

    all_bin_results = []
    status_records = []

    for (cancer, database), sub in real_terms.groupby(["cancer", "database"]):
        result, status = run_one_cancer_database(cancer, database, sub)
        status_records.append(status)

        if result is not None and result.shape[0] > 0:
            all_bin_results.append(result)

        pd.DataFrame(status_records).to_csv(TABLE_DIR / "Run_status_sizebin_random_T300.csv", index=False)

        if len(all_bin_results) > 0:
            partial = pd.concat(all_bin_results, axis=0, ignore_index=True)
            partial.to_csv(TABLE_DIR / "Table_SizeBinRandom_T300_partial_bin_results.csv", index=False)

    if len(all_bin_results) == 0:
        log_message("No random bin results generated.")
        return

    bin_random_df = pd.concat(all_bin_results, axis=0, ignore_index=True)

    mapped_df = map_random_baseline_to_terms(real_terms, bin_random_df)

    db_summary, cancer_db_summary = make_summary_tables(mapped_df, bin_random_df)

    make_figures(mapped_df, db_summary, cancer_db_summary)

    total_elapsed = time.time() - start_global

    config = {
        "BASE_DATA_DIR": str(BASE_DATA_DIR),
        "MULTIDB_RESULT_FILE": str(MULTIDB_RESULT_FILE),
        "OUT_DIR": str(OUT_DIR),
        "T_RANDOM": T_RANDOM,
        "D_THRESHOLD": D_THRESHOLD,
        "OBSERVATION_CLASS_REQUIRED": OBSERVATION_CLASS_REQUIRED,
        "DATABASES_TO_RUN": DATABASES_TO_RUN,
        "CANCERS_TO_RUN": CANCERS_TO_RUN,
        "SIZE_BINS": SIZE_BINS,
        "RANDOM_SIZE_MODE": RANDOM_SIZE_MODE,
        "MIN_SURVIVAL_PATIENTS": MIN_SURVIVAL_PATIENTS,
        "MIN_SURVIVAL_EVENTS": MIN_SURVIVAL_EVENTS,
        "MAX_BINS_PER_CANCER_DATABASE": MAX_BINS_PER_CANCER_DATABASE,
        "RANDOM_SEED": RANDOM_SEED,
        "total_elapsed_seconds": total_elapsed,
        "total_elapsed_hours": total_elapsed / 3600,
        "important_note": (
            "This is a fast global size-bin-matched random baseline. "
            "AIDO-h observation-ready terms are grouped by effective N_matched bins. "
            "Each cancer  x  database  x  size-bin receives one T=300 random baseline, "
            "which is mapped back to all real terms in the same bin. "
            "Top candidates should still receive exact N_matched validation if needed."
        )
    }

    with open(LOG_DIR / "run_config_sizebin_random_T300.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    log_message("============================================================")
    log_message("DONE.")
    log_message(f"Total elapsed: {total_elapsed/3600:.2f} hours")
    log_message("Main outputs:")
    log_message(str(TABLE_DIR / "Table_SizeBinRandom_T300_bin_level_baselines.csv"))
    log_message(str(TABLE_DIR / "Table_SizeBinRandom_T300_term_level_mapped_results.csv"))
    log_message(str(TABLE_DIR / "Table_SizeBinRandom_T300_database_summary.csv"))
    log_message(str(TABLE_DIR / "Table_SizeBinRandom_T300_cancer_database_summary.csv"))
    log_message(str(TABLE_DIR / "Table_SizeBinRandom_T300_top_structured_favored_terms.csv"))
    log_message(str(TABLE_DIR / "Table_SizeBinRandom_T300_highD_but_not_above_random_q95.csv"))
    log_message("============================================================")


if __name__ == "__main__":
    main()
