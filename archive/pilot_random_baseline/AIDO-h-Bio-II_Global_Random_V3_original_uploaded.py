# ============================================================
# AIDO-h-Biology II
# Global Size-Matched Random Baseline
#
# Purpose:
#   Continue after MultiDB AIDO-h + AIDO-D results.
#   Run formal random baseline for observation-ready terms.
#
# Main design:
#   - Use observation-ready terms only
#   - Size-matched random gene sets
#   - T_RANDOM = 300, based on pilot convergence
#   - All available cancers
#   - Hallmark + GO_BP + Reactome
#   - Checkpoint/resume enabled
#
# Input:
#   D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/
#       Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-T300/
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

OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-T300")

TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
LOG_DIR = OUT_DIR / "logs"
CHECKPOINT_DIR = OUT_DIR / "checkpoints"

for d in [OUT_DIR, TABLE_DIR, FIG_DIR, LOG_DIR, CHECKPOINT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Formal global random baseline repeat number
T_RANDOM = 300

# Optional checkpoints within each term
T_CHECKPOINTS = [50, 100, 200, 300]

# Main filter
OBSERVATION_CLASS_REQUIRED = "observation_ready"

# Use all databases by default
DATABASES_TO_RUN = ["Hallmark", "GO_BP", "Reactome"]

# Use all cancers in the main result table by default.
# If you want to test only some cancers, replace None with a list:
# CANCERS_TO_RUN = ["BRCA", "KIRC", "SKCM"]
CANCERS_TO_RUN = None

# D threshold
D_THRESHOLD = 1.301

# Survival minimum requirements
MIN_SURVIVAL_PATIENTS = 30
MIN_SURVIVAL_EVENTS = 5

# Chunking / checkpoint behavior
SAVE_AFTER_EACH_CANCER_DATABASE = True
SKIP_COMPLETED_CANCER_DATABASE = True

# Runtime safety:
# If True, random baseline is computed for all observation-ready terms.
# If False, it will only run D_real >= D_THRESHOLD terms.
RUN_ALL_OBSERVATION_READY_TERMS = True

# Optional: for testing only. Set None for full run.
MAX_TERMS_PER_CANCER_DATABASE = None
# MAX_TERMS_PER_CANCER_DATABASE = 200

RANDOM_SEED = 20260524
np.random.seed(RANDOM_SEED)


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


def safe_neglog10_p(p):
    if p is None or pd.isna(p):
        return np.nan
    if p <= 0:
        return 300.0
    return -np.log10(p)


def infer_cancer_code(folder_name):
    m = re.search(r"\(([^()]+)\)", folder_name)
    if m:
        return m.group(1).strip()
    return folder_name.replace(" ", "_").replace("-", "_")


def sanitize_filename(x):
    x = str(x)
    x = re.sub(r"[^\w\-_\.]+", "_", x)
    return x[:180]


def parse_matched_genes(x):
    if pd.isna(x):
        return []
    genes = str(x).split(";")
    genes = [clean_gene_symbol(g) for g in genes if clean_gene_symbol(g) != ""]
    return sorted(list(set(genes)))


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
# 5. TERM SELECTION
# ============================================================

def load_and_filter_terms():
    if not MULTIDB_RESULT_FILE.exists():
        raise FileNotFoundError(f"Main result file not found: {MULTIDB_RESULT_FILE}")

    df = pd.read_csv(MULTIDB_RESULT_FILE)

    required_cols = [
        "cancer", "database", "bp_name", "N_matched",
        "observation_class", "D_survival", "matched_genes"
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in main result table: {missing}")

    df = df[df["database"].isin(DATABASES_TO_RUN)].copy()

    if CANCERS_TO_RUN is not None:
        df = df[df["cancer"].isin(CANCERS_TO_RUN)].copy()

    df = df[df["observation_class"] == OBSERVATION_CLASS_REQUIRED].copy()
    df = df.dropna(subset=["N_matched", "D_survival"]).copy()

    df["N_matched"] = df["N_matched"].astype(int)

    if not RUN_ALL_OBSERVATION_READY_TERMS:
        df = df[df["D_survival"] >= D_THRESHOLD].copy()

    if MAX_TERMS_PER_CANCER_DATABASE is not None:
        selected = []
        for (cancer, database), sub in df.groupby(["cancer", "database"]):
            sub = sub.sort_values("D_survival", ascending=False).head(MAX_TERMS_PER_CANCER_DATABASE)
            selected.append(sub)
        df = pd.concat(selected, axis=0, ignore_index=True)

    df = df.sort_values(["cancer", "database", "D_survival"], ascending=[True, True, False])
    df.to_csv(TABLE_DIR / "Selected_terms_for_global_random_T300.csv", index=False)

    return df


# ============================================================
# 6. RANDOM BASELINE CORE
# ============================================================

def summarize_random_D(random_D, d_real):
    random_D = np.array(random_D, dtype=float)
    vals = random_D[~np.isnan(random_D)]

    if vals.size == 0:
        return {
            "random_D_mean": np.nan,
            "random_D_median": np.nan,
            "random_D_q90": np.nan,
            "random_D_q95": np.nan,
            "random_D_q99": np.nan,
            "random_D_max": np.nan,
            "D_excess_vs_random_median": np.nan,
            "D_excess_vs_random_q95": np.nan,
            "frac_random_ge_real": np.nan,
            "empirical_p_random": np.nan,
            "structured_favored": np.nan,
            "n_valid_random": 0
        }

    random_median = float(np.median(vals))
    random_q95 = float(np.quantile(vals, 0.95))
    frac_random_ge_real = float(np.mean(vals >= d_real))

    return {
        "random_D_mean": float(np.mean(vals)),
        "random_D_median": random_median,
        "random_D_q90": float(np.quantile(vals, 0.90)),
        "random_D_q95": random_q95,
        "random_D_q99": float(np.quantile(vals, 0.99)),
        "random_D_max": float(np.max(vals)),
        "D_excess_vs_random_median": float(d_real - random_median),
        "D_excess_vs_random_q95": float(d_real - random_q95),
        "frac_random_ge_real": frac_random_ge_real,
        "empirical_p_random": float((np.sum(vals >= d_real) + 1) / (vals.size + 1)),
        "structured_favored": bool(d_real > random_q95),
        "n_valid_random": int(vals.size)
    }


def run_random_for_one_term(
    cancer,
    database,
    bp_name,
    matched_genes,
    n_matched,
    d_real,
    ge_z,
    survival_df,
    rng
):
    all_genes = np.array(ge_z.index.tolist())

    # Remove true BP genes from random pool if possible.
    # This prevents random sets from reusing the same genes as the real BP.
    matched_genes_set = set(matched_genes)
    pool_genes = np.array([g for g in all_genes if g not in matched_genes_set])

    if len(pool_genes) < n_matched:
        pool_genes = all_genes

    random_D = []

    for t in range(T_RANDOM):
        sampled_genes = rng.choice(pool_genes, size=n_matched, replace=False)
        random_score = ge_z.loc[sampled_genes].mean(axis=0, skipna=True)
        d_rand = compute_logrank_D(random_score, survival_df)
        random_D.append(d_rand)

    summary = summarize_random_D(random_D, d_real)

    record = {
        "cancer": cancer,
        "database": database,
        "bp_name": bp_name,
        "N_matched": n_matched,
        "D_real": d_real,
        "T_random": T_RANDOM
    }

    record.update(summary)

    # Optional checkpoint summaries at lower T values
    for T in T_CHECKPOINTS:
        vals = np.array(random_D[:T], dtype=float)
        vals = vals[~np.isnan(vals)]

        if vals.size == 0:
            record[f"random_D_median_T{T}"] = np.nan
            record[f"random_D_q95_T{T}"] = np.nan
            record[f"structured_favored_T{T}"] = np.nan
        else:
            q95 = float(np.quantile(vals, 0.95))
            med = float(np.median(vals))
            record[f"random_D_median_T{T}"] = med
            record[f"random_D_q95_T{T}"] = q95
            record[f"structured_favored_T{T}"] = bool(d_real > q95)

    return record


# ============================================================
# 7. PROCESS ONE CANCER / DATABASE
# ============================================================

def completed_checkpoint_path(cancer, database):
    return CHECKPOINT_DIR / f"{sanitize_filename(cancer)}_{sanitize_filename(database)}_random_T{T_RANDOM}.csv"


def run_one_cancer_database(cancer, database, sub_terms):
    checkpoint_file = completed_checkpoint_path(cancer, database)

    if SKIP_COMPLETED_CANCER_DATABASE and checkpoint_file.exists():
        log_message(f"SKIP completed {cancer} | {database}: {checkpoint_file}")
        try:
            existing = pd.read_csv(checkpoint_file)
            return existing, {
                "cancer": cancer,
                "database": database,
                "status": "skipped_completed",
                "n_terms": existing.shape[0],
                "checkpoint_file": str(checkpoint_file)
            }
        except Exception:
            log_message(f"Existing checkpoint unreadable; rerunning {cancer} | {database}")

    cancer_dir = find_cancer_dir(cancer)

    if cancer_dir is None:
        log_message(f"ERROR {cancer}: cancer folder not found.")
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "error_no_cancer_dir",
            "n_terms": 0
        }

    ge_file = find_ge_file(cancer_dir)
    survival_file = find_survival_or_clinical_file(cancer_dir, cancer)

    if ge_file is None:
        log_message(f"ERROR {cancer}: GE file not found.")
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "error_no_GE",
            "n_terms": 0
        }

    if survival_file is None:
        log_message(f"ERROR {cancer}: survival/clinical file not found.")
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "error_no_survival",
            "n_terms": 0
        }

    log_message("")
    log_message("------------------------------------------------------------")
    log_message(f"Random baseline: {cancer} | {database}")
    log_message("------------------------------------------------------------")
    log_message(f"Terms: {sub_terms.shape[0]}")
    log_message(f"GE file: {ge_file}")
    log_message(f"Survival file: {survival_file}")

    start_time = time.time()

    ge = load_ge_matrix(ge_file)
    survival_df, used_time_col, used_event_col = load_survival_from_any_clinical_file(survival_file)

    common = sorted(list(set(ge.columns).intersection(set(survival_df["patient_id"]))))
    ge = ge[common]
    survival_df = survival_df[survival_df["patient_id"].isin(common)].copy()

    ge_z = zscore_by_gene(ge)

    log_message(
        f"{cancer} | {database}: genes={ge_z.shape[0]}, patients={ge_z.shape[1]}, "
        f"survival={survival_df.shape[0]}, events={int(survival_df['OS_event'].sum())}"
    )

    if survival_df.shape[0] < MIN_SURVIVAL_PATIENTS or survival_df["OS_event"].sum() < MIN_SURVIVAL_EVENTS:
        log_message(f"WARNING {cancer}: insufficient survival patients/events. Skipping.")
        return pd.DataFrame(), {
            "cancer": cancer,
            "database": database,
            "status": "skipped_insufficient_survival",
            "n_terms": sub_terms.shape[0],
            "n_patients": survival_df.shape[0],
            "n_events": int(survival_df["OS_event"].sum())
        }

    rng_seed = abs(hash((cancer, database, RANDOM_SEED))) % (2**32)
    rng = np.random.default_rng(rng_seed)

    records = []
    timing_records = []

    sub_terms = sub_terms.reset_index(drop=True)

    for i, row in sub_terms.iterrows():
        term_start = time.time()

        bp_name = row["bp_name"]
        n_matched = int(row["N_matched"])
        d_real = float(row["D_survival"])
        matched_genes = parse_matched_genes(row["matched_genes"])

        # Safety: keep only genes present in this GE matrix
        matched_genes = [g for g in matched_genes if g in ge_z.index]

        if len(matched_genes) < n_matched:
            n_matched = len(matched_genes)

        if n_matched < 10:
            # Should not happen for observation_ready, but guard anyway.
            continue

        try:
            rec = run_random_for_one_term(
                cancer=cancer,
                database=database,
                bp_name=bp_name,
                matched_genes=matched_genes,
                n_matched=n_matched,
                d_real=d_real,
                ge_z=ge_z,
                survival_df=survival_df,
                rng=rng
            )
            records.append(rec)

        except Exception as e:
            records.append({
                "cancer": cancer,
                "database": database,
                "bp_name": bp_name,
                "N_matched": n_matched,
                "D_real": d_real,
                "T_random": T_RANDOM,
                "error": repr(e)
            })

        elapsed = time.time() - term_start
        timing_records.append({
            "cancer": cancer,
            "database": database,
            "bp_name": bp_name,
            "N_matched": n_matched,
            "seconds": elapsed,
            "seconds_per_random": elapsed / T_RANDOM
        })

        if (i + 1) % 25 == 0 or (i + 1) == sub_terms.shape[0]:
            log_message(
                f"{cancer} | {database}: {i+1}/{sub_terms.shape[0]} terms done. "
                f"Last={elapsed:.2f}s"
            )

            # Save partial checkpoint every 25 terms
            pd.DataFrame(records).to_csv(checkpoint_file, index=False)

    result_df = pd.DataFrame(records)
    result_df.to_csv(checkpoint_file, index=False)

    timing_df = pd.DataFrame(timing_records)
    timing_file = CHECKPOINT_DIR / f"{sanitize_filename(cancer)}_{sanitize_filename(database)}_timing_T{T_RANDOM}.csv"
    timing_df.to_csv(timing_file, index=False)

    total_elapsed = time.time() - start_time

    status = {
        "cancer": cancer,
        "database": database,
        "status": "completed",
        "n_terms": result_df.shape[0],
        "n_patients": ge_z.shape[1],
        "n_events": int(survival_df["OS_event"].sum()),
        "seconds": total_elapsed,
        "minutes": total_elapsed / 60,
        "checkpoint_file": str(checkpoint_file),
        "timing_file": str(timing_file),
        "survival_time_col": used_time_col,
        "survival_event_col": used_event_col
    }

    log_message(
        f"DONE {cancer} | {database}: {result_df.shape[0]} terms, "
        f"{total_elapsed/60:.2f} min"
    )

    return result_df, status


# ============================================================
# 8. SUMMARY AND FIGURES
# ============================================================

def make_summary_tables(random_df):
    random_df.to_csv(TABLE_DIR / "Table_RandomGlobal_T300_term_level_results.csv", index=False)

    rows = []

    for (database, cancer), sub in random_df.groupby(["database", "cancer"]):
        rows.append({
            "database": database,
            "cancer": cancer,
            "n_terms": sub.shape[0],
            "median_D_real": sub["D_real"].median(),
            "median_random_D_median": sub["random_D_median"].median(),
            "median_random_D_q95": sub["random_D_q95"].median(),
            "median_D_excess_vs_random_median": sub["D_excess_vs_random_median"].median(),
            "median_D_excess_vs_random_q95": sub["D_excess_vs_random_q95"].median(),
            "structured_favored_fraction": sub["structured_favored"].mean(),
            "median_frac_random_ge_real": sub["frac_random_ge_real"].median(),
            "n_structured_favored": int(sub["structured_favored"].sum()),
            "n_high_D_real": int((sub["D_real"] >= D_THRESHOLD).sum())
        })

    cancer_db_summary = pd.DataFrame(rows)
    cancer_db_summary.to_csv(TABLE_DIR / "Table_RandomGlobal_T300_cancer_database_summary.csv", index=False)

    rows = []

    for database, sub in random_df.groupby("database"):
        rows.append({
            "database": database,
            "n_terms": sub.shape[0],
            "n_cancers": sub["cancer"].nunique(),
            "median_D_real": sub["D_real"].median(),
            "median_random_D_median": sub["random_D_median"].median(),
            "median_random_D_q95": sub["random_D_q95"].median(),
            "median_D_excess_vs_random_median": sub["D_excess_vs_random_median"].median(),
            "median_D_excess_vs_random_q95": sub["D_excess_vs_random_q95"].median(),
            "structured_favored_fraction": sub["structured_favored"].mean(),
            "median_frac_random_ge_real": sub["frac_random_ge_real"].median(),
            "n_structured_favored": int(sub["structured_favored"].sum()),
            "n_high_D_real": int((sub["D_real"] >= D_THRESHOLD).sum())
        })

    db_summary = pd.DataFrame(rows)
    db_summary.to_csv(TABLE_DIR / "Table_RandomGlobal_T300_database_summary.csv", index=False)

    top_structured = random_df.sort_values(
        ["structured_favored", "D_excess_vs_random_q95", "D_real"],
        ascending=[False, False, False]
    )
    top_structured.to_csv(TABLE_DIR / "Table_RandomGlobal_T300_top_structured_favored_terms.csv", index=False)

    highD_not_favored = random_df[
        (random_df["D_real"] >= D_THRESHOLD) &
        (random_df["structured_favored"] == False)
    ].copy()

    highD_not_favored = highD_not_favored.sort_values("D_real", ascending=False)
    highD_not_favored.to_csv(TABLE_DIR / "Table_RandomGlobal_T300_highD_but_not_above_random_q95.csv", index=False)

    return cancer_db_summary, db_summary


def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def figure_database_structured_fraction(db_summary):
    plot_df = db_summary.sort_values("structured_favored_fraction", ascending=True)

    plt.figure(figsize=(7, 4))
    plt.barh(plot_df["database"], plot_df["structured_favored_fraction"])
    plt.xlim(0, 1.05)
    plt.xlabel("Fraction with D_real > random D q95")
    plt.ylabel("Database")
    plt.title("Structured BP advantage over size-matched random baseline")
    save_fig(FIG_DIR / "Figure_RandomGlobal_structured_favored_fraction_by_database.png")


def figure_database_D_excess(db_summary):
    plot_df = db_summary.sort_values("median_D_excess_vs_random_median", ascending=True)

    plt.figure(figsize=(7, 4))
    plt.barh(plot_df["database"], plot_df["median_D_excess_vs_random_median"])
    plt.axvline(0, linestyle="--", linewidth=1.5)
    plt.xlabel("Median D_real - median random D")
    plt.ylabel("Database")
    plt.title("Excess survival discriminability over random baseline")
    save_fig(FIG_DIR / "Figure_RandomGlobal_D_excess_by_database.png")


def figure_cancer_database_heatmap(cancer_db_summary, value_col, filename, title):
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


def figure_scatter_D_real_vs_random_q95(random_df):
    for db in sorted(random_df["database"].unique()):
        sub = random_df[random_df["database"] == db].copy()

        plt.figure(figsize=(6, 6))
        plt.scatter(sub["random_D_q95"], sub["D_real"], alpha=0.35)
        max_val = np.nanmax([sub["random_D_q95"].max(), sub["D_real"].max()])
        plt.plot([0, max_val], [0, max_val], linestyle="--", linewidth=1.5)
        plt.xlabel("Random D q95")
        plt.ylabel("D_real")
        plt.title(f"D_real vs size-matched random q95: {db}")
        save_fig(FIG_DIR / f"Figure_RandomGlobal_Dreal_vs_random_q95_{sanitize_filename(db)}.png")


def make_figures(random_df, cancer_db_summary, db_summary):
    figure_database_structured_fraction(db_summary)
    figure_database_D_excess(db_summary)

    figure_cancer_database_heatmap(
        cancer_db_summary,
        value_col="structured_favored_fraction",
        filename="Figure_RandomGlobal_structured_favored_fraction_heatmap.png",
        title="Fraction of BP terms exceeding random q95"
    )

    figure_cancer_database_heatmap(
        cancer_db_summary,
        value_col="median_D_excess_vs_random_median",
        filename="Figure_RandomGlobal_D_excess_heatmap.png",
        title="Median D excess over random baseline"
    )

    figure_scatter_D_real_vs_random_q95(random_df)


# ============================================================
# 9. MAIN
# ============================================================

def main():
    log_file = LOG_DIR / "run_log.txt"
    if log_file.exists():
        log_file.unlink()

    log_message("============================================================")
    log_message("AIDO-h-Biology II | Global Random Baseline T=300")
    log_message("Observation-ready terms only")
    log_message("============================================================")
    log_message(f"Input result table: {MULTIDB_RESULT_FILE}")
    log_message(f"Output directory: {OUT_DIR}")
    log_message(f"T_RANDOM: {T_RANDOM}")
    log_message(f"Databases: {DATABASES_TO_RUN}")
    log_message(f"Cancers: {CANCERS_TO_RUN if CANCERS_TO_RUN is not None else 'ALL'}")
    log_message("============================================================")

    selected_terms = load_and_filter_terms()

    log_message(f"Selected terms for random baseline: {selected_terms.shape[0]}")
    log_message("Selected terms by cancer/database:")
    log_message(selected_terms.groupby(["cancer", "database"]).size().to_string())

    all_random_results = []
    status_records = []

    global_start = time.time()

    grouped = selected_terms.groupby(["cancer", "database"])

    for (cancer, database), sub_terms in grouped:
        result_df, status = run_one_cancer_database(cancer, database, sub_terms)
        status_records.append(status)

        if result_df is not None and result_df.shape[0] > 0:
            all_random_results.append(result_df)

        status_df = pd.DataFrame(status_records)
        status_df.to_csv(TABLE_DIR / "Run_status_random_global_T300.csv", index=False)

        # Merge checkpoints so far
        if len(all_random_results) > 0:
            partial = pd.concat(all_random_results, axis=0, ignore_index=True)
            partial.to_csv(TABLE_DIR / "Table_RandomGlobal_T300_partial_merged_results.csv", index=False)

    if len(all_random_results) == 0:
        log_message("No random baseline results generated.")
        return

    random_df = pd.concat(all_random_results, axis=0, ignore_index=True)

    cancer_db_summary, db_summary = make_summary_tables(random_df)
    make_figures(random_df, cancer_db_summary, db_summary)

    total_elapsed = time.time() - global_start

    config = {
        "BASE_DATA_DIR": str(BASE_DATA_DIR),
        "MULTIDB_RESULT_FILE": str(MULTIDB_RESULT_FILE),
        "OUT_DIR": str(OUT_DIR),
        "T_RANDOM": T_RANDOM,
        "T_CHECKPOINTS": T_CHECKPOINTS,
        "OBSERVATION_CLASS_REQUIRED": OBSERVATION_CLASS_REQUIRED,
        "DATABASES_TO_RUN": DATABASES_TO_RUN,
        "CANCERS_TO_RUN": CANCERS_TO_RUN,
        "D_THRESHOLD": D_THRESHOLD,
        "MIN_SURVIVAL_PATIENTS": MIN_SURVIVAL_PATIENTS,
        "MIN_SURVIVAL_EVENTS": MIN_SURVIVAL_EVENTS,
        "RUN_ALL_OBSERVATION_READY_TERMS": RUN_ALL_OBSERVATION_READY_TERMS,
        "MAX_TERMS_PER_CANCER_DATABASE": MAX_TERMS_PER_CANCER_DATABASE,
        "RANDOM_SEED": RANDOM_SEED,
        "total_elapsed_seconds": total_elapsed,
        "total_elapsed_hours": total_elapsed / 3600,
        "important_note": (
            "Random baseline is applied after AIDO-h observation-readiness filtering. "
            "Each real BP observable is compared against T=300 size-matched random gene sets. "
            "Random genes are sampled from measured genes, excluding the real BP genes when possible."
        )
    }

    with open(LOG_DIR / "run_config_random_global_T300.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    log_message("============================================================")
    log_message("DONE.")
    log_message(f"Total elapsed: {total_elapsed/3600:.2f} hours")
    log_message(f"Total random baseline terms: {random_df.shape[0]}")
    log_message("Main outputs:")
    log_message(str(TABLE_DIR / "Table_RandomGlobal_T300_term_level_results.csv"))
    log_message(str(TABLE_DIR / "Table_RandomGlobal_T300_database_summary.csv"))
    log_message(str(TABLE_DIR / "Table_RandomGlobal_T300_cancer_database_summary.csv"))
    log_message(str(TABLE_DIR / "Table_RandomGlobal_T300_top_structured_favored_terms.csv"))
    log_message(str(TABLE_DIR / "Table_RandomGlobal_T300_highD_but_not_above_random_q95.csv"))
    log_message("============================================================")


if __name__ == "__main__":
    main()
