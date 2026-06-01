# ============================================================
# AIDO-h-Biology II
# Random Baseline Pilot / Convergence Test
#
# Goal:
#   Run size-matched random baselines on selected cancers
#   to determine appropriate random repeat number T.
#
# Recommended pilot cancers:
#   KIRC, SKCM, BRCA, COAD, THCA
#
# Input:
#   D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/
#       Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-RandomPilot/
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

OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomPilot")

TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
LOG_DIR = OUT_DIR / "logs"

for d in [OUT_DIR, TABLE_DIR, FIG_DIR, LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Pilot cancers chosen to represent different D/event regimes
PILOT_CANCERS = ["KIRC", "SKCM", "BRCA", "COAD", "THCA"]

# Databases to test
PILOT_DATABASES = ["Hallmark", "Reactome", "GO_BP"]

# Use only observation-ready terms for primary random baseline
OBSERVATION_CLASS_REQUIRED = "observation_ready"

# To reduce runtime in pilot:
#   "highD_only" = only D_real >= D_THRESHOLD
#   "topN_per_cancer_db" = top N per cancer/database
#   "mixed" = top N high-D + random ready terms
PILOT_MODE = "mixed"

D_THRESHOLD = 1.301

TOP_N_PER_CANCER_DB = 80
RANDOM_READY_N_PER_CANCER_DB = 40

# Random repeat checkpoints
T_CHECKPOINTS = [50, 100, 200, 500, 1000]
T_MAX = max(T_CHECKPOINTS)

# If too slow, set T_MAX=500 and T_CHECKPOINTS=[50,100,200,500]
# T_CHECKPOINTS = [50, 100, 200, 500]
# T_MAX = 500

MIN_SURVIVAL_PATIENTS = 30
MIN_SURVIVAL_EVENTS = 5

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
# 3. LOAD DATA
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
# 5. SELECT PILOT TERMS
# ============================================================

def select_pilot_terms(main_df):
    df = main_df.copy()

    df = df[df["cancer"].isin(PILOT_CANCERS)]
    df = df[df["database"].isin(PILOT_DATABASES)]
    df = df[df["observation_class"] == OBSERVATION_CLASS_REQUIRED]
    df = df.dropna(subset=["D_survival", "N_matched"])

    df["N_matched"] = df["N_matched"].astype(int)

    selected = []

    for (cancer, db), sub in df.groupby(["cancer", "database"]):
        sub = sub.copy()

        if PILOT_MODE == "highD_only":
            pick = sub[sub["D_survival"] >= D_THRESHOLD].copy()

        elif PILOT_MODE == "topN_per_cancer_db":
            pick = sub.sort_values("D_survival", ascending=False).head(TOP_N_PER_CANCER_DB)

        elif PILOT_MODE == "mixed":
            top = sub.sort_values("D_survival", ascending=False).head(TOP_N_PER_CANCER_DB)

            ready_pool = sub.drop(index=top.index, errors="ignore")
            if ready_pool.shape[0] > RANDOM_READY_N_PER_CANCER_DB:
                random_ready = ready_pool.sample(
                    n=RANDOM_READY_N_PER_CANCER_DB,
                    random_state=RANDOM_SEED
                )
            else:
                random_ready = ready_pool

            pick = pd.concat([top, random_ready], axis=0).drop_duplicates(
                subset=["cancer", "database", "bp_name"]
            )

        else:
            raise ValueError(f"Unknown PILOT_MODE: {PILOT_MODE}")

        selected.append(pick)

    if len(selected) == 0:
        return pd.DataFrame()

    selected_df = pd.concat(selected, axis=0, ignore_index=True)

    selected_df.to_csv(TABLE_DIR / "Pilot_selected_terms.csv", index=False)

    return selected_df


# ============================================================
# 6. RANDOM BASELINE
# ============================================================

def run_random_baseline_for_term(
    cancer,
    database,
    bp_name,
    n_matched,
    d_real,
    ge_z,
    survival_df,
    rng
):
    genes = np.array(ge_z.index.tolist())

    if n_matched < 3:
        return []

    if n_matched > len(genes):
        return []

    random_D = []

    for t in range(T_MAX):
        sampled_genes = rng.choice(genes, size=n_matched, replace=False)
        random_score = ge_z.loc[sampled_genes].mean(axis=0, skipna=True)
        d_rand = compute_logrank_D(random_score, survival_df)
        random_D.append(d_rand)

    random_D = np.array(random_D, dtype=float)

    records = []

    for T in T_CHECKPOINTS:
        vals = random_D[:T]
        vals = vals[~np.isnan(vals)]

        if vals.size == 0:
            rec = {
                "cancer": cancer,
                "database": database,
                "bp_name": bp_name,
                "N_matched": n_matched,
                "D_real": d_real,
                "T_random": T,
                "random_D_mean": np.nan,
                "random_D_median": np.nan,
                "random_D_q90": np.nan,
                "random_D_q95": np.nan,
                "random_D_q99": np.nan,
                "random_D_max": np.nan,
                "D_excess_vs_random_median": np.nan,
                "empirical_p_random": np.nan,
                "real_exceeds_random_q95": np.nan,
                "n_valid_random": 0
            }
        else:
            random_q95 = np.quantile(vals, 0.95)
            random_median = np.median(vals)

            rec = {
                "cancer": cancer,
                "database": database,
                "bp_name": bp_name,
                "N_matched": n_matched,
                "D_real": d_real,
                "T_random": T,
                "random_D_mean": float(np.mean(vals)),
                "random_D_median": float(random_median),
                "random_D_q90": float(np.quantile(vals, 0.90)),
                "random_D_q95": float(random_q95),
                "random_D_q99": float(np.quantile(vals, 0.99)),
                "random_D_max": float(np.max(vals)),
                "D_excess_vs_random_median": float(d_real - random_median),
                "empirical_p_random": float((np.sum(vals >= d_real) + 1) / (vals.size + 1)),
                "real_exceeds_random_q95": bool(d_real > random_q95),
                "n_valid_random": int(vals.size)
            }

        records.append(rec)

    return records


def run_random_pilot(selected_terms):
    all_records = []
    timing_records = []

    rng = np.random.default_rng(RANDOM_SEED)

    for cancer in selected_terms["cancer"].unique():
        log_message("")
        log_message("============================================================")
        log_message(f"Random pilot cancer: {cancer}")
        log_message("============================================================")

        cancer_dir = find_cancer_dir(cancer)

        if cancer_dir is None:
            log_message(f"SKIP {cancer}: cancer folder not found.")
            continue

        ge_file = find_ge_file(cancer_dir)
        survival_file = find_survival_or_clinical_file(cancer_dir, cancer)

        if ge_file is None or survival_file is None:
            log_message(f"SKIP {cancer}: missing GE or survival/clinical file.")
            continue

        ge = load_ge_matrix(ge_file)
        ge_z = zscore_by_gene(ge)

        survival_df, time_col, event_col = load_survival_from_any_clinical_file(survival_file)

        common = sorted(list(set(ge_z.columns).intersection(set(survival_df["patient_id"]))))
        survival_df = survival_df[survival_df["patient_id"].isin(common)].copy()
        ge_z = ge_z[common]

        log_message(
            f"{cancer}: GE genes={ge_z.shape[0]}, patients={ge_z.shape[1]}, "
            f"survival patients={survival_df.shape[0]}, events={int(survival_df['OS_event'].sum())}"
        )

        sub_terms = selected_terms[selected_terms["cancer"] == cancer].copy()
        sub_terms = sub_terms.reset_index(drop=True)

        cancer_start = time.time()

        for i, row in sub_terms.iterrows():
            term_start = time.time()

            database = row["database"]
            bp_name = row["bp_name"]
            n_matched = int(row["N_matched"])
            d_real = float(row["D_survival"])

            recs = run_random_baseline_for_term(
                cancer=cancer,
                database=database,
                bp_name=bp_name,
                n_matched=n_matched,
                d_real=d_real,
                ge_z=ge_z,
                survival_df=survival_df,
                rng=rng
            )

            all_records.extend(recs)

            elapsed = time.time() - term_start

            timing_records.append({
                "cancer": cancer,
                "database": database,
                "bp_name": bp_name,
                "N_matched": n_matched,
                "seconds_for_Tmax": elapsed,
                "T_MAX": T_MAX,
                "seconds_per_random": elapsed / T_MAX if T_MAX > 0 else np.nan
            })

            if (i + 1) % 20 == 0 or (i + 1) == sub_terms.shape[0]:
                log_message(
                    f"{cancer}: {i+1}/{sub_terms.shape[0]} terms done. "
                    f"Last term {elapsed:.2f}s"
                )

        cancer_elapsed = time.time() - cancer_start
        log_message(f"{cancer}: completed in {cancer_elapsed/60:.2f} min")

        # Save intermediate after each cancer
        pd.DataFrame(all_records).to_csv(
            TABLE_DIR / "RandomPilot_term_level_checkpoint_results.csv",
            index=False
        )

        pd.DataFrame(timing_records).to_csv(
            TABLE_DIR / "RandomPilot_timing_checkpoint.csv",
            index=False
        )

    result_df = pd.DataFrame(all_records)
    timing_df = pd.DataFrame(timing_records)

    return result_df, timing_df


# ============================================================
# 7. CONVERGENCE SUMMARIES
# ============================================================

def make_convergence_summary(result_df):
    rows = []

    # Compare each T against T_MAX as reference
    ref = result_df[result_df["T_random"] == T_MAX].copy()
    ref_cols = [
        "cancer", "database", "bp_name",
        "random_D_median", "random_D_q95",
        "empirical_p_random", "real_exceeds_random_q95"
    ]

    ref = ref[ref_cols].rename(columns={
        "random_D_median": "ref_random_D_median",
        "random_D_q95": "ref_random_D_q95",
        "empirical_p_random": "ref_empirical_p_random",
        "real_exceeds_random_q95": "ref_real_exceeds_random_q95"
    })

    merged = result_df.merge(
        ref,
        on=["cancer", "database", "bp_name"],
        how="left"
    )

    merged["abs_delta_q95_vs_Tmax"] = (
        merged["random_D_q95"] - merged["ref_random_D_q95"]
    ).abs()

    merged["abs_delta_median_vs_Tmax"] = (
        merged["random_D_median"] - merged["ref_random_D_median"]
    ).abs()

    merged["abs_delta_empirical_p_vs_Tmax"] = (
        merged["empirical_p_random"] - merged["ref_empirical_p_random"]
    ).abs()

    merged["q95_class_agrees_with_Tmax"] = (
        merged["real_exceeds_random_q95"] == merged["ref_real_exceeds_random_q95"]
    )

    for T, sub in merged.groupby("T_random"):
        rows.append({
            "T_random": T,
            "n_terms": sub.shape[0],
            "median_abs_delta_q95_vs_Tmax": sub["abs_delta_q95_vs_Tmax"].median(),
            "mean_abs_delta_q95_vs_Tmax": sub["abs_delta_q95_vs_Tmax"].mean(),
            "p90_abs_delta_q95_vs_Tmax": sub["abs_delta_q95_vs_Tmax"].quantile(0.90),
            "median_abs_delta_median_vs_Tmax": sub["abs_delta_median_vs_Tmax"].median(),
            "median_abs_delta_empirical_p_vs_Tmax": sub["abs_delta_empirical_p_vs_Tmax"].median(),
            "q95_class_agreement_with_Tmax": sub["q95_class_agrees_with_Tmax"].mean(),
            "fraction_real_exceeds_q95": sub["real_exceeds_random_q95"].mean()
        })

    summary = pd.DataFrame(rows).sort_values("T_random")

    merged.to_csv(TABLE_DIR / "RandomPilot_convergence_term_level_vs_Tmax.csv", index=False)
    summary.to_csv(TABLE_DIR / "RandomPilot_convergence_summary.csv", index=False)

    return summary, merged


def make_runtime_summary(timing_df, selected_terms):
    if timing_df.shape[0] == 0:
        return pd.DataFrame()

    seconds_per_random_median = timing_df["seconds_per_random"].median()
    seconds_per_random_mean = timing_df["seconds_per_random"].mean()

    n_pilot_terms = selected_terms.shape[0]

    rows = []

    for T in T_CHECKPOINTS:
        rows.append({
            "T_random": T,
            "n_pilot_terms": n_pilot_terms,
            "median_seconds_per_random_test": seconds_per_random_median,
            "mean_seconds_per_random_test": seconds_per_random_mean,
            "estimated_seconds_for_pilot": seconds_per_random_median * T * n_pilot_terms,
            "estimated_minutes_for_pilot": seconds_per_random_median * T * n_pilot_terms / 60,
            "estimated_hours_for_pilot": seconds_per_random_median * T * n_pilot_terms / 3600
        })

    runtime = pd.DataFrame(rows)
    runtime.to_csv(TABLE_DIR / "RandomPilot_runtime_summary.csv", index=False)

    return runtime


# ============================================================
# 8. FIGURES
# ============================================================

def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def figure_convergence_q95(summary):
    plt.figure(figsize=(7, 4))
    plt.plot(summary["T_random"], summary["median_abs_delta_q95_vs_Tmax"], marker="o")
    plt.xlabel("Number of random repeats T")
    plt.ylabel("Median |Delta random q95| vs Tmax")
    plt.title("Random baseline q95 convergence")
    save_fig(FIG_DIR / "Figure_RandomPilot_q95_convergence.png")


def figure_class_agreement(summary):
    plt.figure(figsize=(7, 4))
    plt.plot(summary["T_random"], summary["q95_class_agreement_with_Tmax"], marker="o")
    plt.ylim(0, 1.05)
    plt.xlabel("Number of random repeats T")
    plt.ylabel("Agreement with Tmax classification")
    plt.title("Stability of real > random q95 classification")
    save_fig(FIG_DIR / "Figure_RandomPilot_q95_class_agreement.png")


def figure_runtime(runtime):
    plt.figure(figsize=(7, 4))
    plt.plot(runtime["T_random"], runtime["estimated_hours_for_pilot"], marker="o")
    plt.xlabel("Number of random repeats T")
    plt.ylabel("Estimated pilot runtime, hours")
    plt.title("Estimated runtime vs random repeats")
    save_fig(FIG_DIR / "Figure_RandomPilot_runtime.png")


# ============================================================
# 9. MAIN
# ============================================================

def main():
    log_file = LOG_DIR / "run_log.txt"
    if log_file.exists():
        log_file.unlink()

    log_message("============================================================")
    log_message("AIDO-h-Biology II | Random Baseline Pilot")
    log_message("============================================================")
    log_message(f"Input result table: {MULTIDB_RESULT_FILE}")
    log_message(f"Output directory: {OUT_DIR}")
    log_message(f"Pilot cancers: {PILOT_CANCERS}")
    log_message(f"Pilot databases: {PILOT_DATABASES}")
    log_message(f"Pilot mode: {PILOT_MODE}")
    log_message(f"T checkpoints: {T_CHECKPOINTS}")
    log_message("============================================================")

    if not MULTIDB_RESULT_FILE.exists():
        raise FileNotFoundError(f"Main result file not found: {MULTIDB_RESULT_FILE}")

    main_df = pd.read_csv(MULTIDB_RESULT_FILE)

    selected_terms = select_pilot_terms(main_df)

    if selected_terms.shape[0] == 0:
        log_message("No selected terms. Stop.")
        return

    log_message(f"Selected pilot terms: {selected_terms.shape[0]}")
    log_message(
        selected_terms.groupby(["cancer", "database"]).size().to_string()
    )

    start_time = time.time()

    result_df, timing_df = run_random_pilot(selected_terms)

    result_df.to_csv(TABLE_DIR / "RandomPilot_term_level_results.csv", index=False)
    timing_df.to_csv(TABLE_DIR / "RandomPilot_timing.csv", index=False)

    convergence_summary, convergence_term = make_convergence_summary(result_df)
    runtime_summary = make_runtime_summary(timing_df, selected_terms)

    figure_convergence_q95(convergence_summary)
    figure_class_agreement(convergence_summary)
    figure_runtime(runtime_summary)

    total_elapsed = time.time() - start_time

    config = {
        "BASE_DATA_DIR": str(BASE_DATA_DIR),
        "MULTIDB_RESULT_FILE": str(MULTIDB_RESULT_FILE),
        "OUT_DIR": str(OUT_DIR),
        "PILOT_CANCERS": PILOT_CANCERS,
        "PILOT_DATABASES": PILOT_DATABASES,
        "OBSERVATION_CLASS_REQUIRED": OBSERVATION_CLASS_REQUIRED,
        "PILOT_MODE": PILOT_MODE,
        "D_THRESHOLD": D_THRESHOLD,
        "TOP_N_PER_CANCER_DB": TOP_N_PER_CANCER_DB,
        "RANDOM_READY_N_PER_CANCER_DB": RANDOM_READY_N_PER_CANCER_DB,
        "T_CHECKPOINTS": T_CHECKPOINTS,
        "T_MAX": T_MAX,
        "MIN_SURVIVAL_PATIENTS": MIN_SURVIVAL_PATIENTS,
        "MIN_SURVIVAL_EVENTS": MIN_SURVIVAL_EVENTS,
        "RANDOM_SEED": RANDOM_SEED,
        "total_elapsed_seconds": total_elapsed,
        "total_elapsed_minutes": total_elapsed / 60
    }

    with open(LOG_DIR / "run_config_random_pilot.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    log_message("============================================================")
    log_message("DONE.")
    log_message(f"Elapsed time: {total_elapsed/60:.2f} min")
    log_message("Main outputs:")
    log_message(str(TABLE_DIR / "Pilot_selected_terms.csv"))
    log_message(str(TABLE_DIR / "RandomPilot_term_level_results.csv"))
    log_message(str(TABLE_DIR / "RandomPilot_convergence_summary.csv"))
    log_message(str(TABLE_DIR / "RandomPilot_runtime_summary.csv"))
    log_message("============================================================")


if __name__ == "__main__":
    main()
