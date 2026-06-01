# ============================================================
# AIDO-h-Biology II
# Post-AIDO-h Random Baseline Pilot / Convergence Test
#
# Purpose:
#   Run random baseline repeat-number pilot:
#       T = 50, 100, 200, 300, 400, 500
#
#   Use the pilot to decide final random repeat number for
#   AIDO-h-Biology II / CSBJ manuscript.
#
# Core principle:
#   Random baseline is performed AFTER AIDO-h filtering.
#   Primary comparison uses observation_ready terms only.
#
# Input:
#   D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/
#       Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv
#
# Required raw data:
#   D:/AIDO-Data/UCSC_XENA/<Cancer folder>/GE.tsv
#   Survival/clinical file in each cancer folder
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-RandomPilot-T50-500/
#
# Recommended run:
#   python AIDO-h-Bio-II-RandomPilot-T50-500.py
#
# Dependencies:
#   pip install pandas numpy matplotlib lifelines
# ============================================================

import os
import re
import time
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

BASE_DATA_DIR = Path(r"D:/AIDO-Data/UCSC_XENA")

MULTIDB_RESULT_FILE = Path(
    r"D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB/tables/"
    r"Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv"
)

OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomPilot-T50-500")

TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
LOG_DIR = OUT_DIR / "logs"
CACHE_DIR = OUT_DIR / "cache"

for d in [OUT_DIR, TABLE_DIR, FIG_DIR, LOG_DIR, CACHE_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Pilot cancers
# ------------------------------------------------------------
# These cancers cover different regimes:
#   BRCA: large cohort / stable
#   SKCM: often strong immune/transcriptomic signal
#   KIRC: often strong survival signal
#   COAD: possible projection / endpoint mismatch
#   THCA: weaker OS/event-limited regime
# ------------------------------------------------------------

PILOT_CANCERS = ["BRCA", "SKCM", "KIRC", "COAD", "THCA"]


# ------------------------------------------------------------
# Databases
# ------------------------------------------------------------

PILOT_DATABASES = ["Hallmark", "GO_BP", "Reactome"]


# ------------------------------------------------------------
# AIDO-h class required for primary random baseline
# ------------------------------------------------------------

OBSERVATION_CLASS_REQUIRED = "observation_ready"


# ------------------------------------------------------------
# Pilot selection mode
# ------------------------------------------------------------
# Options:
#   "highD_only":
#       only terms with D_survival >= D_THRESHOLD
#
#   "topN_per_cancer_db":
#       top N by D_survival for each cancer/database
#
#   "mixed":
#       top N by D_survival + random observation-ready terms
#
# Recommended:
#   "mixed"
# ------------------------------------------------------------

PILOT_MODE = "mixed"

D_THRESHOLD = 1.301  # -log10(0.05)

TOP_N_PER_CANCER_DB = 60
RANDOM_READY_N_PER_CANCER_DB = 30


# ------------------------------------------------------------
# Random repeat checkpoints
# ------------------------------------------------------------

T_CHECKPOINTS = [50, 100, 200, 300, 400, 500]
T_MAX = max(T_CHECKPOINTS)


# ------------------------------------------------------------
# Survival requirements
# ------------------------------------------------------------

MIN_SURVIVAL_PATIENTS = 30
MIN_SURVIVAL_EVENTS = 5
MIN_GROUP_SIZE = 10


# ------------------------------------------------------------
# Random gene universe options
# ------------------------------------------------------------
# If True, random gene sets will exclude the genes of the
# corresponding real BP/pathway.
#
# For primary manuscript baseline, False is usually acceptable
# because random sets are sampled from the measured gene universe
# and matched by effective gene count.
# ------------------------------------------------------------

EXCLUDE_REAL_GENES_FROM_RANDOM = False


# ------------------------------------------------------------
# Save detailed random records?
# ------------------------------------------------------------
# True: saves one row per real term per random repeat.
# Can be large but useful for debugging/reproducibility.
# ------------------------------------------------------------

SAVE_ALL_RANDOM_REPEATS = True


# ------------------------------------------------------------
# Random seed
# ------------------------------------------------------------

RANDOM_SEED = 20260524
rng = np.random.default_rng(RANDOM_SEED)


# ============================================================
# 1. LOGGING AND UTILITIES
# ============================================================

def log_message(msg):
    msg = str(msg)
    print(msg)
    with open(LOG_DIR / "run_log.txt", "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def normalize_tcga_barcode_to_patient(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip().replace(".", "-")
    return x[:12]


def clean_gene_symbol(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    x = re.sub(r"\.\d+$", "", x)
    return x.upper()


def safe_neglog10_p(p):
    if p is None or pd.isna(p):
        return np.nan
    try:
        p = float(p)
    except Exception:
        return np.nan
    if p <= 0:
        return 300.0
    return -np.log10(p)


def infer_cancer_code(folder_name):
    m = re.search(r"\(([^()]+)\)", str(folder_name))
    if m:
        return m.group(1).strip()
    return str(folder_name).replace(" ", "_").replace("-", "_")


def sanitize_filename(x):
    x = str(x)
    x = re.sub(r"[^\w\-_\.]+", "_", x)
    return x[:180]


def parse_matched_genes(x):
    if pd.isna(x):
        return []
    if isinstance(x, list):
        return [clean_gene_symbol(g) for g in x if clean_gene_symbol(g) != ""]
    genes = str(x).replace(",", ";").split(";")
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
# 3. TABLE READERS
# ============================================================

def read_table_flexible(path):
    path = Path(path)

    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)

    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if df.shape[1] > 1:
            return df
    except Exception:
        pass

    try:
        df = pd.read_csv(path, sep=",", low_memory=False)
        if df.shape[1] > 1:
            return df
    except Exception:
        pass

    return pd.read_csv(path, sep=None, engine="python", low_memory=False)


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


# ============================================================
# 4. LOAD GE MATRIX
# ============================================================

def load_ge_matrix(ge_path):
    log_message(f"[LOAD] GE: {ge_path}")

    df = pd.read_csv(ge_path, sep="\t", index_col=0, low_memory=False)

    df.index = [clean_gene_symbol(x) for x in df.index]
    df = df.loc[df.index != ""]

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=0, how="all")
    df = df.dropna(axis=1, how="all")

    if df.index.duplicated().any():
        df = df.groupby(df.index).mean(numeric_only=True)

    df.columns = [normalize_tcga_barcode_to_patient(c) for c in df.columns]
    df = df.loc[:, [c for c in df.columns if isinstance(c, str) and c.startswith("TCGA-")]]

    if pd.Index(df.columns).duplicated().any():
        df = df.T.groupby(level=0).mean(numeric_only=True).T

    log_message(f"[OK] GE loaded: genes={df.shape[0]}, patients={df.shape[1]}")
    return df


def zscore_by_gene(ge):
    mu = ge.mean(axis=1, skipna=True)
    sd = ge.std(axis=1, skipna=True).replace(0, np.nan)
    z = ge.sub(mu, axis=0).div(sd, axis=0)
    z = z.replace([np.inf, -np.inf], np.nan)
    return z


# ============================================================
# 5. LOAD SURVIVAL
# ============================================================

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


def load_survival_from_any_clinical_file(path):
    log_message(f"[LOAD] Survival/clinical: {path}")

    df = read_table_flexible(path)
    df.columns = [str(c).strip() for c in df.columns]

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

    # Build OS_time
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

    # Build OS_event
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

    log_message(
        f"[OK] Survival loaded: patients={out.shape[0]}, "
        f"events={int(out['OS_event'].sum())}, "
        f"time_col={used_time_col}, event_col={used_event_col}"
    )

    return out, used_time_col, used_event_col


# ============================================================
# 6. SURVIVAL D
# ============================================================

def compute_logrank_D(score, survival_df):
    try:
        from lifelines.statistics import logrank_test
    except Exception:
        raise ImportError("Please install lifelines first: pip install lifelines")

    temp = survival_df.copy()
    temp["score"] = score.reindex(temp["patient_id"]).values
    temp = temp.dropna(subset=["score", "OS_time", "OS_event"])

    if temp.shape[0] < MIN_SURVIVAL_PATIENTS:
        return np.nan

    if temp["OS_event"].sum() < MIN_SURVIVAL_EVENTS:
        return np.nan

    median_val = temp["score"].median()

    low = temp[temp["score"] <= median_val]
    high = temp[temp["score"] > median_val]

    if low.shape[0] < MIN_GROUP_SIZE or high.shape[0] < MIN_GROUP_SIZE:
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
# 7. PILOT TERM SELECTION
# ============================================================

def load_main_multidb_results():
    if not MULTIDB_RESULT_FILE.exists():
        raise FileNotFoundError(
            f"Main MultiDB result file not found:\n{MULTIDB_RESULT_FILE}\n\n"
            f"Please run AIDO-h-Bio-II.py first."
        )

    df = pd.read_csv(MULTIDB_RESULT_FILE, low_memory=False)

    required = ["cancer", "database", "bp_name", "N_matched", "observation_class"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Required columns missing from main result file: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    if "D_survival" not in df.columns:
        raise ValueError(
            "D_survival column not found. Random pilot requires survival D results."
        )

    if "matched_genes" not in df.columns:
        raise ValueError(
            "matched_genes column not found. Random pilot requires matched gene lists."
        )

    return df


def select_pilot_terms(main_df):
    df = main_df.copy()

    df = df[df["cancer"].isin(PILOT_CANCERS)]
    df = df[df["database"].isin(PILOT_DATABASES)]
    df = df[df["observation_class"] == OBSERVATION_CLASS_REQUIRED]

    df["D_survival"] = pd.to_numeric(df["D_survival"], errors="coerce")
    df["N_matched"] = pd.to_numeric(df["N_matched"], errors="coerce")

    df = df.dropna(subset=["D_survival", "N_matched"])
    df["N_matched"] = df["N_matched"].astype(int)

    selected_parts = []

    for (cancer, db), sub in df.groupby(["cancer", "database"]):
        sub = sub.copy()
        sub = sub.sort_values("D_survival", ascending=False)

        if PILOT_MODE == "highD_only":
            pick = sub[sub["D_survival"] >= D_THRESHOLD].copy()

        elif PILOT_MODE == "topN_per_cancer_db":
            pick = sub.head(TOP_N_PER_CANCER_DB).copy()

        elif PILOT_MODE == "mixed":
            top_part = sub.head(TOP_N_PER_CANCER_DB).copy()

            remaining = sub.drop(index=top_part.index, errors="ignore")

            if remaining.shape[0] > 0:
                n_random = min(RANDOM_READY_N_PER_CANCER_DB, remaining.shape[0])
                random_part = remaining.sample(
                    n=n_random,
                    random_state=RANDOM_SEED
                ).copy()
                pick = pd.concat([top_part, random_part], axis=0)
            else:
                pick = top_part.copy()

        else:
            raise ValueError(f"Unknown PILOT_MODE: {PILOT_MODE}")

        pick = pick.drop_duplicates(subset=["cancer", "database", "bp_name"])
        selected_parts.append(pick)

        log_message(
            f"[SELECT] {cancer} | {db}: available={sub.shape[0]}, selected={pick.shape[0]}"
        )

    if len(selected_parts) == 0:
        raise ValueError("No pilot terms selected. Check filters and input result table.")

    selected = pd.concat(selected_parts, axis=0).reset_index(drop=True)

    selected["pilot_term_id"] = [
        f"T{i+1:05d}" for i in range(selected.shape[0])
    ]

    selected.to_csv(TABLE_DIR / "Table_00_selected_pilot_terms.csv", index=False)

    log_message(f"[OK] Total selected pilot terms: {selected.shape[0]}")

    return selected


# ============================================================
# 8. RANDOM BASELINE CORE
# ============================================================

def build_real_score_from_matched_genes(ge_z, matched_genes):
    genes = [g for g in matched_genes if g in ge_z.index]
    if len(genes) == 0:
        return None
    return ge_z.loc[genes].mean(axis=0, skipna=True)


def sample_random_genes(gene_universe, n, exclude_genes=None):
    gene_universe = list(gene_universe)

    if exclude_genes is not None and len(exclude_genes) > 0:
        exclude_set = set(exclude_genes)
        gene_universe = [g for g in gene_universe if g not in exclude_set]

    if len(gene_universe) < n:
        return None

    return list(rng.choice(gene_universe, size=n, replace=False))


def run_random_for_one_cancer(cancer_code, selected_terms):
    t0 = time.time()

    cancer_dir = find_cancer_dir(cancer_code)
    if cancer_dir is None:
        log_message(f"[SKIP] {cancer_code}: cancer folder not found.")
        return [], []

    ge_file = find_ge_file(cancer_dir)
    if ge_file is None:
        log_message(f"[SKIP] {cancer_code}: GE file not found.")
        return [], []

    clinical_file = find_survival_or_clinical_file(cancer_dir, cancer_code)
    if clinical_file is None:
        log_message(f"[SKIP] {cancer_code}: survival/clinical file not found.")
        return [], []

    try:
        ge = load_ge_matrix(ge_file)
        ge_z = zscore_by_gene(ge)
        survival_df, used_time_col, used_event_col = load_survival_from_any_clinical_file(clinical_file)
    except Exception as e:
        log_message(f"[SKIP] {cancer_code}: failed loading data. Error={e}")
        return [], []

    # Align to patients with survival
    common_patients = sorted(list(set(ge_z.columns).intersection(set(survival_df["patient_id"]))))
    if len(common_patients) < MIN_SURVIVAL_PATIENTS:
        log_message(
            f"[SKIP] {cancer_code}: insufficient common patients: {len(common_patients)}"
        )
        return [], []

    ge_z = ge_z[common_patients]
    survival_df = survival_df[survival_df["patient_id"].isin(common_patients)].copy()

    gene_universe = sorted([g for g in ge_z.index if isinstance(g, str) and g != ""])

    log_message(
        f"[RUN] {cancer_code}: terms={selected_terms.shape[0]}, "
        f"genes={len(gene_universe)}, patients={len(common_patients)}, "
        f"events={int(survival_df['OS_event'].sum())}"
    )

    all_random_records = []
    checkpoint_records = []

    for idx, row in selected_terms.iterrows():
        term_start = time.time()

        pilot_term_id = row["pilot_term_id"]
        db = row["database"]
        bp_name = row["bp_name"]
        n_matched = int(row["N_matched"])
        d_real_table = row.get("D_survival", np.nan)

        matched_genes = parse_matched_genes(row["matched_genes"])
        matched_genes = [g for g in matched_genes if g in ge_z.index]

        if len(matched_genes) < 1:
            continue

        # Recompute real D from matched genes for consistency
        real_score = build_real_score_from_matched_genes(ge_z, matched_genes)

        if real_score is None:
            continue

        d_real_recomputed = compute_logrank_D(real_score, survival_df)

        # Use effective matched count from actual available genes
        n_random = len(matched_genes)

        random_D_values = []

        for repeat_i in range(1, T_MAX + 1):
            exclude = matched_genes if EXCLUDE_REAL_GENES_FROM_RANDOM else None

            random_genes = sample_random_genes(
                gene_universe=gene_universe,
                n=n_random,
                exclude_genes=exclude
            )

            if random_genes is None:
                d_random = np.nan
            else:
                random_score = ge_z.loc[random_genes].mean(axis=0, skipna=True)
                d_random = compute_logrank_D(random_score, survival_df)

            random_D_values.append(d_random)

            if SAVE_ALL_RANDOM_REPEATS:
                all_random_records.append({
                    "cancer": cancer_code,
                    "database": db,
                    "bp_name": bp_name,
                    "pilot_term_id": pilot_term_id,
                    "repeat_i": repeat_i,
                    "N_matched_real": len(matched_genes),
                    "N_random": n_random,
                    "D_real_from_table": d_real_table,
                    "D_real_recomputed": d_real_recomputed,
                    "D_random": d_random
                })

        # Checkpoint summaries for this term
        random_arr = np.array(random_D_values, dtype=float)

        for T in T_CHECKPOINTS:
            vals = random_arr[:T]
            vals = vals[np.isfinite(vals)]

            if len(vals) == 0:
                median_random = np.nan
                mean_random = np.nan
                p95_random = np.nan
                max_random = np.nan
                frac_random_ge_real = np.nan
                delta_real_minus_random_median = np.nan
                structured_favored = np.nan
            else:
                median_random = float(np.nanmedian(vals))
                mean_random = float(np.nanmean(vals))
                p95_random = float(np.nanpercentile(vals, 95))
                max_random = float(np.nanmax(vals))

                if pd.isna(d_real_recomputed):
                    frac_random_ge_real = np.nan
                    delta_real_minus_random_median = np.nan
                    structured_favored = np.nan
                else:
                    frac_random_ge_real = float(np.mean(vals >= d_real_recomputed))
                    delta_real_minus_random_median = float(d_real_recomputed - median_random)
                    structured_favored = int(d_real_recomputed > median_random)

            checkpoint_records.append({
                "cancer": cancer_code,
                "database": db,
                "bp_name": bp_name,
                "pilot_term_id": pilot_term_id,
                "T": T,
                "N_matched_real": len(matched_genes),
                "D_real_from_table": d_real_table,
                "D_real_recomputed": d_real_recomputed,
                "random_mean_D": mean_random,
                "random_median_D": median_random,
                "random_p95_D": p95_random,
                "random_max_D": max_random,
                "delta_real_minus_random_median": delta_real_minus_random_median,
                "frac_random_ge_real": frac_random_ge_real,
                "structured_favored": structured_favored,
                "n_valid_random": int(len(vals))
            })

        if (idx + 1) % 10 == 0:
            elapsed = time.time() - t0
            log_message(
                f"[PROGRESS] {cancer_code}: processed {idx+1}/{selected_terms.shape[0]} "
                f"terms | elapsed={elapsed/60:.1f} min"
            )

    elapsed = time.time() - t0
    log_message(f"[DONE] {cancer_code}: elapsed={elapsed/60:.1f} min")

    return all_random_records, checkpoint_records


# ============================================================
# 9. SUMMARIES AND CONVERGENCE ANALYSIS
# ============================================================

def summarize_checkpoint_results(checkpoint_df):
    df = checkpoint_df.copy()

    # --------------------------------------------------------
    # Overall by T
    # --------------------------------------------------------
    overall = (
        df.groupby("T")
        .agg(
            n_terms=("pilot_term_id", "nunique"),
            n_rows=("pilot_term_id", "size"),
            median_random_median_D=("random_median_D", "median"),
            mean_random_median_D=("random_median_D", "mean"),
            median_random_p95_D=("random_p95_D", "median"),
            mean_random_p95_D=("random_p95_D", "mean"),
            median_delta_real_minus_random=("delta_real_minus_random_median", "median"),
            mean_delta_real_minus_random=("delta_real_minus_random_median", "mean"),
            structured_favored_fraction=("structured_favored", "mean"),
            median_frac_random_ge_real=("frac_random_ge_real", "median")
        )
        .reset_index()
    )

    overall = overall.sort_values("T")
    overall["change_median_random_median_D"] = overall["median_random_median_D"].diff().abs()
    overall["change_median_random_p95_D"] = overall["median_random_p95_D"].diff().abs()
    overall["change_structured_favored_fraction"] = overall["structured_favored_fraction"].diff().abs()

    overall.to_csv(
        TABLE_DIR / "Table_02_convergence_overall_by_T.csv",
        index=False
    )

    # --------------------------------------------------------
    # By cancer/database/T
    # --------------------------------------------------------
    by_cd = (
        df.groupby(["cancer", "database", "T"])
        .agg(
            n_terms=("pilot_term_id", "nunique"),
            median_random_median_D=("random_median_D", "median"),
            mean_random_median_D=("random_median_D", "mean"),
            median_random_p95_D=("random_p95_D", "median"),
            mean_random_p95_D=("random_p95_D", "mean"),
            median_delta_real_minus_random=("delta_real_minus_random_median", "median"),
            mean_delta_real_minus_random=("delta_real_minus_random_median", "mean"),
            structured_favored_fraction=("structured_favored", "mean"),
            median_frac_random_ge_real=("frac_random_ge_real", "median")
        )
        .reset_index()
        .sort_values(["cancer", "database", "T"])
    )

    by_cd["change_median_random_median_D"] = (
        by_cd.groupby(["cancer", "database"])["median_random_median_D"]
        .diff()
        .abs()
    )

    by_cd["change_median_random_p95_D"] = (
        by_cd.groupby(["cancer", "database"])["median_random_p95_D"]
        .diff()
        .abs()
    )

    by_cd["change_structured_favored_fraction"] = (
        by_cd.groupby(["cancer", "database"])["structured_favored_fraction"]
        .diff()
        .abs()
    )

    by_cd.to_csv(
        TABLE_DIR / "Table_03_convergence_by_cancer_database_T.csv",
        index=False
    )

    # --------------------------------------------------------
    # By cancer/T
    # --------------------------------------------------------
    by_cancer = (
        df.groupby(["cancer", "T"])
        .agg(
            n_terms=("pilot_term_id", "nunique"),
            median_random_median_D=("random_median_D", "median"),
            median_random_p95_D=("random_p95_D", "median"),
            median_delta_real_minus_random=("delta_real_minus_random_median", "median"),
            structured_favored_fraction=("structured_favored", "mean")
        )
        .reset_index()
        .sort_values(["cancer", "T"])
    )

    by_cancer["change_median_random_median_D"] = (
        by_cancer.groupby("cancer")["median_random_median_D"]
        .diff()
        .abs()
    )

    by_cancer.to_csv(
        TABLE_DIR / "Table_04_convergence_by_cancer_T.csv",
        index=False
    )

    # --------------------------------------------------------
    # By database/T
    # --------------------------------------------------------
    by_db = (
        df.groupby(["database", "T"])
        .agg(
            n_terms=("pilot_term_id", "nunique"),
            median_random_median_D=("random_median_D", "median"),
            median_random_p95_D=("random_p95_D", "median"),
            median_delta_real_minus_random=("delta_real_minus_random_median", "median"),
            structured_favored_fraction=("structured_favored", "mean")
        )
        .reset_index()
        .sort_values(["database", "T"])
    )

    by_db["change_median_random_median_D"] = (
        by_db.groupby("database")["median_random_median_D"]
        .diff()
        .abs()
    )

    by_db.to_csv(
        TABLE_DIR / "Table_05_convergence_by_database_T.csv",
        index=False
    )

    return overall, by_cd, by_cancer, by_db


def make_convergence_figures(overall, by_cd, by_cancer, by_db):
    # --------------------------------------------------------
    # Figure 1: Overall random median D convergence
    # --------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))
    plt.plot(
        overall["T"],
        overall["median_random_median_D"],
        marker="o",
        label="Median of term-level random median D"
    )
    plt.plot(
        overall["T"],
        overall["median_random_p95_D"],
        marker="o",
        label="Median of term-level random 95th percentile D"
    )
    plt.xlabel("Random repeat number T")
    plt.ylabel("Random-baseline D summary")
    plt.title("Overall random baseline convergence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Figure_RandomPilot_01_overall_convergence.png", dpi=300)
    plt.savefig(FIG_DIR / "Figure_RandomPilot_01_overall_convergence.pdf")
    plt.close()

    # --------------------------------------------------------
    # Figure 2: Overall checkpoint-to-checkpoint change
    # --------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))
    plt.plot(
        overall["T"],
        overall["change_median_random_median_D"],
        marker="o",
        label="Change in random median D"
    )
    plt.plot(
        overall["T"],
        overall["change_median_random_p95_D"],
        marker="o",
        label="Change in random 95th percentile D"
    )
    plt.axhline(0.02, linestyle="--", linewidth=1, label="0.02 reference")
    plt.axhline(0.05, linestyle=":", linewidth=1, label="0.05 reference")
    plt.xlabel("Random repeat number T")
    plt.ylabel("Absolute change from previous checkpoint")
    plt.title("Checkpoint-to-checkpoint convergence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Figure_RandomPilot_02_checkpoint_change.png", dpi=300)
    plt.savefig(FIG_DIR / "Figure_RandomPilot_02_checkpoint_change.pdf")
    plt.close()

    # --------------------------------------------------------
    # Figure 3: By cancer
    # --------------------------------------------------------
    plt.figure(figsize=(8.5, 5.5))
    for cancer, sub in by_cancer.groupby("cancer"):
        sub = sub.sort_values("T")
        plt.plot(
            sub["T"],
            sub["median_random_median_D"],
            marker="o",
            label=cancer
        )
    plt.xlabel("Random repeat number T")
    plt.ylabel("Median random median D")
    plt.title("Random baseline convergence by cancer")
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Figure_RandomPilot_03_by_cancer.png", dpi=300)
    plt.savefig(FIG_DIR / "Figure_RandomPilot_03_by_cancer.pdf")
    plt.close()

    # --------------------------------------------------------
    # Figure 4: By database
    # --------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))
    for db, sub in by_db.groupby("database"):
        sub = sub.sort_values("T")
        plt.plot(
            sub["T"],
            sub["median_random_median_D"],
            marker="o",
            label=db
        )
    plt.xlabel("Random repeat number T")
    plt.ylabel("Median random median D")
    plt.title("Random baseline convergence by database")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Figure_RandomPilot_04_by_database.png", dpi=300)
    plt.savefig(FIG_DIR / "Figure_RandomPilot_04_by_database.pdf")
    plt.close()

    # --------------------------------------------------------
    # Figure 5: Structured-favored fraction convergence
    # --------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))
    plt.plot(
        overall["T"],
        overall["structured_favored_fraction"],
        marker="o"
    )
    plt.xlabel("Random repeat number T")
    plt.ylabel("Fraction of terms with D_real > random median D")
    plt.title("Structured-favored fraction across random-repeat checkpoints")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Figure_RandomPilot_05_structured_favored_fraction.png", dpi=300)
    plt.savefig(FIG_DIR / "Figure_RandomPilot_05_structured_favored_fraction.pdf")
    plt.close()


# ============================================================
# 10. RECOMMEND FINAL T
# ============================================================

def recommend_final_T(overall):
    """
    Rule-of-thumb recommendation:
        Prefer the smallest T >= 300 where:
            change in median random median D <= 0.02
            change in median random p95 D <= 0.05

        If none satisfies, recommend 500.
    """
    df = overall.copy()
    df = df.sort_values("T")

    candidate = df[
        (df["T"] >= 300)
        & (df["change_median_random_median_D"] <= 0.02)
        & (df["change_median_random_p95_D"] <= 0.05)
    ]

    if candidate.shape[0] > 0:
        rec_T = int(candidate.iloc[0]["T"])
        reason = (
            "Smallest T >= 300 with checkpoint changes below "
            "0.02 for median random D and below 0.05 for random 95th percentile D."
        )
    else:
        rec_T = 500
        reason = (
            "No checkpoint satisfied the strict convergence rule. "
            "Use T=500 as the safest pilot-tested choice, or consider T=1000 if "
            "large instability remains."
        )

    summary = {
        "recommended_T": rec_T,
        "reason": reason,
        "T_checkpoints": T_CHECKPOINTS,
        "pilot_cancers": PILOT_CANCERS,
        "pilot_databases": PILOT_DATABASES,
        "pilot_mode": PILOT_MODE,
        "observation_class_required": OBSERVATION_CLASS_REQUIRED,
        "exclude_real_genes_from_random": EXCLUDE_REAL_GENES_FROM_RANDOM,
        "random_seed": RANDOM_SEED
    }

    with open(TABLE_DIR / "Recommendation_final_random_T.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    with open(OUT_DIR / "README_RandomPilot_Summary.txt", "w", encoding="utf-8") as f:
        f.write("AIDO-h-Biology II Random Baseline Pilot Summary\n")
        f.write("================================================\n\n")
        f.write(f"Recommended final random repeat number T: {rec_T}\n\n")
        f.write(f"Reason:\n{reason}\n\n")
        f.write("Pilot configuration:\n")
        f.write(f"  Pilot cancers: {PILOT_CANCERS}\n")
        f.write(f"  Pilot databases: {PILOT_DATABASES}\n")
        f.write(f"  T checkpoints: {T_CHECKPOINTS}\n")
        f.write(f"  Pilot mode: {PILOT_MODE}\n")
        f.write(f"  Observation class: {OBSERVATION_CLASS_REQUIRED}\n")
        f.write(f"  Exclude real genes from random: {EXCLUDE_REAL_GENES_FROM_RANDOM}\n")
        f.write(f"  Random seed: {RANDOM_SEED}\n\n")
        f.write("Key output tables:\n")
        f.write("  tables/Table_00_selected_pilot_terms.csv\n")
        f.write("  tables/Table_01_random_checkpoint_by_term.csv\n")
        f.write("  tables/Table_02_convergence_overall_by_T.csv\n")
        f.write("  tables/Table_03_convergence_by_cancer_database_T.csv\n")
        f.write("  tables/Table_04_convergence_by_cancer_T.csv\n")
        f.write("  tables/Table_05_convergence_by_database_T.csv\n\n")
        f.write("Key figures:\n")
        f.write("  figures/Figure_RandomPilot_01_overall_convergence.png\n")
        f.write("  figures/Figure_RandomPilot_02_checkpoint_change.png\n")
        f.write("  figures/Figure_RandomPilot_03_by_cancer.png\n")
        f.write("  figures/Figure_RandomPilot_04_by_database.png\n")
        f.write("  figures/Figure_RandomPilot_05_structured_favored_fraction.png\n")

    log_message(f"[RECOMMENDATION] Final T recommendation: {rec_T}")
    log_message(f"[RECOMMENDATION] {reason}")

    return summary


# ============================================================
# 11. MAIN
# ============================================================

def main():
    start_time = time.time()

    log_message("============================================================")
    log_message("AIDO-h-Biology II Random Baseline Pilot T=50-500")
    log_message("============================================================")
    log_message(f"BASE_DATA_DIR = {BASE_DATA_DIR}")
    log_message(f"MULTIDB_RESULT_FILE = {MULTIDB_RESULT_FILE}")
    log_message(f"OUT_DIR = {OUT_DIR}")
    log_message(f"PILOT_CANCERS = {PILOT_CANCERS}")
    log_message(f"PILOT_DATABASES = {PILOT_DATABASES}")
    log_message(f"T_CHECKPOINTS = {T_CHECKPOINTS}")
    log_message(f"T_MAX = {T_MAX}")
    log_message(f"PILOT_MODE = {PILOT_MODE}")
    log_message(f"OBSERVATION_CLASS_REQUIRED = {OBSERVATION_CLASS_REQUIRED}")
    log_message(f"RANDOM_SEED = {RANDOM_SEED}")
    log_message("============================================================")

    # Save config
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
        "MIN_GROUP_SIZE": MIN_GROUP_SIZE,
        "EXCLUDE_REAL_GENES_FROM_RANDOM": EXCLUDE_REAL_GENES_FROM_RANDOM,
        "SAVE_ALL_RANDOM_REPEATS": SAVE_ALL_RANDOM_REPEATS,
        "RANDOM_SEED": RANDOM_SEED
    }

    with open(OUT_DIR / "config_random_pilot.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    # Load and select terms
    main_df = load_main_multidb_results()
    selected = select_pilot_terms(main_df)

    all_random_records = []
    all_checkpoint_records = []

    for cancer_code in PILOT_CANCERS:
        sub_terms = selected[selected["cancer"] == cancer_code].copy()

        if sub_terms.shape[0] == 0:
            log_message(f"[SKIP] {cancer_code}: no selected terms.")
            continue

        random_records, checkpoint_records = run_random_for_one_cancer(
            cancer_code=cancer_code,
            selected_terms=sub_terms
        )

        all_random_records.extend(random_records)
        all_checkpoint_records.extend(checkpoint_records)

        # Save partial outputs after each cancer
        if len(all_checkpoint_records) > 0:
            partial_checkpoint_df = pd.DataFrame(all_checkpoint_records)
            partial_checkpoint_df.to_csv(
                TABLE_DIR / "PARTIAL_Table_01_random_checkpoint_by_term.csv",
                index=False
            )

        if SAVE_ALL_RANDOM_REPEATS and len(all_random_records) > 0:
            partial_random_df = pd.DataFrame(all_random_records)
            partial_random_df.to_csv(
                TABLE_DIR / "PARTIAL_Table_99_all_random_repeats.csv",
                index=False
            )

    if len(all_checkpoint_records) == 0:
        raise RuntimeError("No checkpoint records generated. Check input data and logs.")

    checkpoint_df = pd.DataFrame(all_checkpoint_records)

    checkpoint_df.to_csv(
        TABLE_DIR / "Table_01_random_checkpoint_by_term.csv",
        index=False
    )

    if SAVE_ALL_RANDOM_REPEATS and len(all_random_records) > 0:
        random_df = pd.DataFrame(all_random_records)
        random_df.to_csv(
            TABLE_DIR / "Table_99_all_random_repeats.csv",
            index=False
        )

    # Summaries
    overall, by_cd, by_cancer, by_db = summarize_checkpoint_results(checkpoint_df)

    # Figures
    make_convergence_figures(overall, by_cd, by_cancer, by_db)

    # Recommendation
    rec = recommend_final_T(overall)

    elapsed = time.time() - start_time

    log_message("============================================================")
    log_message("[FINISHED] Random pilot completed.")
    log_message(f"Elapsed time: {elapsed/60:.1f} minutes")
    log_message(f"Output folder: {OUT_DIR}")
    log_message(f"Recommended final T: {rec['recommended_T']}")
    log_message("============================================================")


if __name__ == "__main__":
    main()
