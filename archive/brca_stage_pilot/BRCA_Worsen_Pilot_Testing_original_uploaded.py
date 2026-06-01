# ============================================================
# AIDO-h Pilot Test FIXED VERSION
# BRCA GE + Hallmark / GO BP / Reactome
# AIDO-h Lite + Stage I/II vs III/IV discriminability
#
# Output:
#   D:/AIDO-Temp/AIDOh_BRCA_GE_Stage_FIXED/
#
# Main fixes:
#   1. Robust stage parser
#   2. Debug output for stage columns and Early/Late counts
#   3. Figure 2/3 no longer silently blank if D is invalid
#   4. AIDO-h Lite completeness + representation + mu-D/W-D
# ============================================================

import os
import re
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")


# ============================================================
# 0. CONFIG
# ============================================================

GSEA_DIR = Path(r"D:/AIDO-Data/GSEA")
BRCA_DIR = Path(r"D:/AIDO-Data/UCSC_XENA/Breast Cancer (BRCA)")

OUT_DIR = Path(r"D:/AIDO-Temp/AIDOh_BRCA_GE_Stage_FIXED")
OUT_DIR.mkdir(parents=True, exist_ok=True)

GE_FILE = BRCA_DIR / "GE.tsv"
STAGE_FILE = BRCA_DIR / "BRCA_stage_groups_from_survival.tsv"

GMT_FILES = {
    "Hallmark": GSEA_DIR / "h.all.v2026.1.Hs.symbols.gmt",
    "GO_BP": GSEA_DIR / "c5.go.bp.v2026.1.Hs.symbols.gmt",
    "Reactome": GSEA_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
}

MIN_MATCHED_GENES = 10

MAX_GENES_FOR_COHERENCE = 300
MAX_GENES_FOR_STABILITY = 150
BOOTSTRAP_N = 30

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


# ============================================================
# 1. BASIC UTILITIES
# ============================================================

def normalize_tcga_barcode(x):
    """
    Normalize TCGA sample barcode to patient-level barcode.
    Example:
        TCGA-XX-YYYY-01A -> TCGA-XX-YYYY
    """
    if pd.isna(x):
        return np.nan

    x = str(x).strip().upper()
    x = x.replace(".", "-")
    return x[:12]


def clean_gene_symbol(x):
    """
    Clean gene symbol.
    """
    if pd.isna(x):
        return np.nan

    x = str(x).strip()
    x = re.sub(r"\.\d+$", "", x)
    return x.upper()


def read_table_auto(path):
    """
    Robust TSV/CSV reader.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

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


def zscore_rows(mat):
    """
    Z-score each row across columns.
    Input: genes x samples.
    """
    arr = mat.astype(float).values

    mean = np.nanmean(arr, axis=1, keepdims=True)
    std = np.nanstd(arr, axis=1, ddof=1, keepdims=True)

    std[std == 0] = np.nan

    z = (arr - mean) / std
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)

    return pd.DataFrame(z, index=mat.index, columns=mat.columns)


def iqr(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) == 0:
        return np.nan

    return np.nanpercentile(x, 75) - np.nanpercentile(x, 25)


def neglog10p(p):
    if p is None or pd.isna(p):
        return np.nan

    p = max(float(p), 1e-300)
    return -np.log10(p)


# ============================================================
# 2. LOAD GE MATRIX
# ============================================================

def load_ge_matrix(path):
    """
    Load TCGA GE matrix.

    Expected common format:
        rows = genes
        columns = samples

    First column is treated as gene symbol.
    """
    print(f"[LOAD] GE file: {path}")

    raw = pd.read_csv(path, sep="\t", low_memory=False)

    first_col = raw.columns[0]
    raw = raw.rename(columns={first_col: "gene"})

    raw["gene"] = raw["gene"].map(clean_gene_symbol)
    raw = raw.dropna(subset=["gene"])
    raw = raw[raw["gene"] != ""]

    raw = raw.set_index("gene")
    raw = raw.apply(pd.to_numeric, errors="coerce")

    raw = raw.dropna(axis=0, how="all")
    raw = raw.dropna(axis=1, how="all")

    # Normalize sample barcodes to TCGA patient-level barcode
    raw.columns = [normalize_tcga_barcode(c) for c in raw.columns]

    # Keep TCGA-like columns
    keep_cols = [
        c for c in raw.columns
        if isinstance(c, str) and c.startswith("TCGA-")
    ]
    raw = raw.loc[:, keep_cols]

    # Collapse duplicated patient columns by mean
    raw = raw.T.groupby(level=0).mean().T

    # Collapse duplicated genes by mean
    raw = raw.groupby(raw.index).mean()

    print(f"[OK] GE matrix loaded: genes={raw.shape[0]}, patients={raw.shape[1]}")

    return raw


# ============================================================
# 3. LOAD STAGE GROUPS - FIXED
# ============================================================

def infer_stage_group_from_text(x):
    """
    Convert stage text or numeric code to Early / Late.

    Early = Stage I / II
    Late  = Stage III / IV

    Accepts:
        Stage I, Stage IA, I, IA, 1
        Stage II, Stage IIA, II, IIA, 2
        Stage III, IIIA, 3
        Stage IV, IVA, 4
        I/II, I-II, III/IV, III-IV
        Early, Late
    """
    if pd.isna(x):
        return np.nan

    s = str(x).strip().upper()

    if s == "":
        return np.nan

    # Normalize common strings
    s = s.replace("_", " ")
    s = s.replace("-", " ")
    s = s.replace(".", " ")
    s = s.replace("PATHOLOGIC", "")
    s = s.replace("CLINICAL", "")
    s = s.replace("AJCC", "")
    s = s.replace("STAGE", "")
    s = s.replace("STG", "")
    s = re.sub(r"\s+", " ", s).strip()

    # Common missing values
    if s in [
        "NA", "NAN", "NONE", "NULL", "NOT AVAILABLE", "NOT REPORTED",
        "UNKNOWN", "X", "[NOT AVAILABLE]", "[NOT REPORTED]"
    ]:
        return np.nan

    # Already grouped labels
    if s in [
        "EARLY", "LOW", "EARLY STAGE", "LOW STAGE",
        "I II", "I/II", "I-II", "STAGE I II"
    ]:
        return "Early"

    if s in [
        "LATE", "HIGH", "ADVANCED", "LATE STAGE", "HIGH STAGE",
        "III IV", "III/IV", "III-IV", "STAGE III IV"
    ]:
        return "Late"

    # Numeric encoding
    if s in ["0", "1", "2"]:
        return "Early"

    if s in ["3", "4"]:
        return "Late"

    # Important: check IV and III before II and I
    if re.search(r"\bIV[A-C]?\b", s):
        return "Late"

    if re.search(r"\bIII[A-C]?\b", s):
        return "Late"

    if re.search(r"\bII[A-C]?\b", s):
        return "Early"

    if re.search(r"\bI[A-C]?\b", s):
        return "Early"

    # Loose fallback
    if "III" in s or "IV" in s:
        return "Late"

    if "II" in s or s == "I":
        return "Early"

    return np.nan


def load_stage_groups(path):
    """
    Robust loader for BRCA stage group file.

    It tries to identify:
    - patient/sample barcode column
    - stage or group column

    Then converts the stage column into Early / Late.
    """
    print(f"[LOAD] Stage file: {path}")

    df = read_table_auto(path)
    df.columns = [str(c).strip() for c in df.columns]

    print("[DEBUG] Stage file columns:")
    for c in df.columns:
        print("   ", c)

    # -----------------------------
    # 1. Find patient/sample column
    # -----------------------------
    sample_candidates = [
        c for c in df.columns
        if c.lower() in [
            "sample",
            "patient",
            "patient_id",
            "barcode",
            "bcr_patient_barcode",
            "tcga_id",
            "submitter_id",
            "case_id",
            "sample_id"
        ]
        or "sample" in c.lower()
        or "patient" in c.lower()
        or "barcode" in c.lower()
        or "submitter" in c.lower()
        or "case" in c.lower()
    ]

    if len(sample_candidates) == 0:
        sample_col = df.columns[0]
    else:
        sample_col = sample_candidates[0]

    print(f"[DEBUG] Using patient column: {sample_col}")

    df["patient"] = df[sample_col].map(normalize_tcga_barcode)

    # -----------------------------
    # 2. Find stage/group column
    # -----------------------------
    priority_cols = []

    exact_priority = [
        "stage_group",
        "stage_binary",
        "early_late",
        "stage_category",
        "pathologic_stage",
        "clinical_stage",
        "ajcc_pathologic_stage",
        "ajcc_clinical_stage",
        "stage"
    ]

    for name in exact_priority:
        for c in df.columns:
            if c.lower() == name:
                priority_cols.append(c)

    # Contains stage
    for c in df.columns:
        if "stage" in c.lower() and c not in priority_cols:
            priority_cols.append(c)

    # Contains group
    for c in df.columns:
        if (
            "group" in c.lower()
            and c not in priority_cols
            and c != sample_col
            and "patient" not in c.lower()
            and "sample" not in c.lower()
        ):
            priority_cols.append(c)

    # Fallback: try all non-patient columns
    for c in df.columns:
        if c not in priority_cols and c not in [sample_col, "patient"]:
            priority_cols.append(c)

    # Try each candidate and choose the one that produces most Early/Late labels
    best_col = None
    best_valid = -1
    best_labels = None

    print("[DEBUG] Testing possible stage columns:")

    for c in priority_cols:
        try:
            labels = df[c].map(infer_stage_group_from_text)
            valid = int(labels.notna().sum())
            vc = labels.value_counts(dropna=True).to_dict()
            print(f"   {c}: valid={valid}, counts={vc}")

            if valid > best_valid:
                best_valid = valid
                best_col = c
                best_labels = labels

        except Exception as e:
            print(f"   {c}: failed with {e}")

    if best_col is None or best_valid == 0:
        preview_path = OUT_DIR / "DEBUG_stage_file_preview.csv"
        df.head(50).to_csv(preview_path, index=False)
        raise ValueError(
            "No usable stage labels were detected. "
            f"A preview was saved to: {preview_path}"
        )

    print(f"[DEBUG] Using stage column: {best_col}")

    df["stage_group"] = best_labels

    out = df[["patient", "stage_group"]].dropna()
    out = out.drop_duplicates("patient")

    print(f"[OK] Stage groups loaded: patients={out.shape[0]}")
    print(out["stage_group"].value_counts())

    if out["stage_group"].nunique() < 2:
        raise ValueError(
            "Only one stage group detected after parsing. "
            "Need both Early and Late."
        )

    out.to_csv(OUT_DIR / "DEBUG_stage_groups_parsed.csv", index=False)

    return out


# ============================================================
# 4. LOAD GMT FILES
# ============================================================

def load_gmt(path, database_name):
    """
    Load GMT file.

    Format:
        gene_set_name \t description \t gene1 \t gene2 ...
    """
    print(f"[LOAD] GMT: {database_name} | {path}")

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"GMT file not found: {path}")

    gene_sets = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            name = parts[0].strip()
            desc = parts[1].strip()

            genes = [
                clean_gene_symbol(g)
                for g in parts[2:]
                if str(g).strip() != ""
            ]

            genes = sorted(list(set([
                g for g in genes
                if isinstance(g, str) and g != ""
            ])))

            gene_sets.append({
                "database": database_name,
                "bp_id": name,
                "bp_name": name,
                "description": desc,
                "genes": genes,
                "original_gene_count": len(genes)
            })

    print(f"[OK] Loaded {len(gene_sets)} gene sets from {database_name}")

    return gene_sets


# ============================================================
# 5. AIDO-h DIAGNOSTICS
# ============================================================

def compute_coherence(z_sub, max_genes=300):
    """
    Mean pairwise gene-gene correlation.
    Input: genes x samples z-scored matrix.
    """
    n_genes = z_sub.shape[0]

    if n_genes < 2:
        return np.nan

    if n_genes > max_genes:
        selected = np.random.choice(z_sub.index, size=max_genes, replace=False)
        z_sub = z_sub.loc[selected]

    arr = z_sub.values.astype(float)

    try:
        corr = np.corrcoef(arr)
    except Exception:
        return np.nan

    if corr.ndim != 2:
        return np.nan

    iu = np.triu_indices_from(corr, k=1)
    vals = corr[iu]
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return np.nan

    return float(np.nanmean(vals))


def compute_pca1_stats(z_sub):
    """
    PCA1 variance explained and mean-score vs PC1 correlation.
    Input: genes x samples z-scored matrix.
    """
    if z_sub.shape[0] < 2 or z_sub.shape[1] < 3:
        return np.nan, np.nan

    X = z_sub.T.values.astype(float)  # samples x genes
    X = X - np.nanmean(X, axis=0, keepdims=True)
    X = np.nan_to_num(X, nan=0.0)

    try:
        U, S, Vt = np.linalg.svd(X, full_matrices=False)

        denom = np.sum(S ** 2)

        if denom <= 0:
            return np.nan, np.nan

        var_explained = (S ** 2) / denom
        pca1_var = float(var_explained[0])

        pc1_score = U[:, 0] * S[0]
        mean_score = z_sub.mean(axis=0).values

        if np.std(mean_score) == 0 or np.std(pc1_score) == 0:
            corr = np.nan
        else:
            corr = float(np.corrcoef(mean_score, pc1_score)[0, 1])

        # PCA sign is arbitrary, use absolute agreement
        return pca1_var, abs(corr)

    except Exception:
        return np.nan, np.nan


def compute_leave_one_gene_stability(z_sub, max_genes=150):
    """
    Correlation between full mean score and leave-one-gene-out mean scores.
    Returns median LOO correlation.
    """
    n_genes = z_sub.shape[0]

    if n_genes < 3:
        return np.nan

    if n_genes > max_genes:
        selected = np.random.choice(z_sub.index, size=max_genes, replace=False)
        z_sub = z_sub.loc[selected]
        n_genes = z_sub.shape[0]

    arr = z_sub.values.astype(float)
    full = np.nanmean(arr, axis=0)

    corrs = []
    total = np.nansum(arr, axis=0)

    for i in range(n_genes):
        loo = (total - arr[i, :]) / max(n_genes - 1, 1)

        if np.std(loo) == 0 or np.std(full) == 0:
            continue

        r = np.corrcoef(full, loo)[0, 1]

        if np.isfinite(r):
            corrs.append(r)

    if len(corrs) == 0:
        return np.nan

    return float(np.nanmedian(corrs))


def compute_bootstrap_stability(z_sub, n_boot=30):
    """
    Bootstrap gene-subset stability.
    Correlate full mean score with bootstrap mean scores.
    """
    n_genes = z_sub.shape[0]

    if n_genes < 5:
        return np.nan

    arr = z_sub.values.astype(float)
    full = np.nanmean(arr, axis=0)

    corrs = []
    sample_size = max(3, int(0.8 * n_genes))

    for _ in range(n_boot):
        idx = np.random.choice(np.arange(n_genes), size=sample_size, replace=True)
        boot = np.nanmean(arr[idx, :], axis=0)

        if np.std(boot) == 0 or np.std(full) == 0:
            continue

        r = np.corrcoef(full, boot)[0, 1]

        if np.isfinite(r):
            corrs.append(r)

    if len(corrs) == 0:
        return np.nan

    return float(np.nanmedian(corrs))


def compute_stage_D(score_series, stage_df):
    """
    Compute stage discriminability:
        Early vs Late

    Uses Mann-Whitney U test and AUROC.

    D_stage = -log10(p)
    """
    tmp = pd.DataFrame({
        "patient": score_series.index,
        "score": score_series.values
    })

    tmp = tmp.merge(stage_df, on="patient", how="inner")
    tmp = tmp.dropna(subset=["score", "stage_group"])

    # Normalize again just in case
    tmp["stage_group"] = tmp["stage_group"].map(infer_stage_group_from_text)
    tmp = tmp.dropna(subset=["stage_group"])

    early = tmp.loc[tmp["stage_group"] == "Early", "score"].astype(float).values
    late = tmp.loc[tmp["stage_group"] == "Late", "score"].astype(float).values

    n_early = len(early)
    n_late = len(late)

    if n_early < 5 or n_late < 5:
        return np.nan, np.nan, np.nan, n_early, n_late

    try:
        p = stats.mannwhitneyu(early, late, alternative="two-sided").pvalue
    except Exception:
        p = np.nan

    D = neglog10p(p)

    try:
        y = (tmp["stage_group"] == "Late").astype(int).values
        s = tmp["score"].astype(float).values
        auc = roc_auc_score(y, s)
        auc_abs = max(auc, 1 - auc)
    except Exception:
        auc_abs = np.nan

    return D, p, auc_abs, n_early, n_late


# ============================================================
# 6. MAIN AIDO-h PROCESSING
# ============================================================

def process_gene_sets_for_database(database_name, gene_sets, ge, stage_df):
    """
    Run AIDO-h diagnostics for one database.
    """
    print(f"\n[PROCESS] {database_name}")

    measured_genes = set(ge.index)
    all_records = []

    for idx, gs in enumerate(gene_sets, start=1):
        if idx % 500 == 0:
            print(f"  processed {idx}/{len(gene_sets)} gene sets...")

        genes = gs["genes"]
        matched = sorted(list(set(genes).intersection(measured_genes)))

        original_n = len(genes)
        matched_n = len(matched)
        coverage = matched_n / original_n if original_n > 0 else np.nan
        pass_min = matched_n >= MIN_MATCHED_GENES

        base = {
            "database": database_name,
            "bp_id": gs["bp_id"],
            "bp_name": gs["bp_name"],
            "description": gs.get("description", ""),
            "original_gene_count": original_n,
            "matched_gene_count": matched_n,
            "coverage_fraction": coverage,
            "pass_min_gene_count": pass_min,
        }

        if not pass_min:
            rec = dict(base)

            rec.update({
                "mu_variance": np.nan,
                "mu_IQR": np.nan,
                "W_median": np.nan,
                "W_IQR": np.nan,
                "coherence_mean": np.nan,
                "PCA1_variance": np.nan,
                "mean_PCA_corr": np.nan,
                "loo_stability": np.nan,
                "bootstrap_stability": np.nan,
                "D_stage_mu": np.nan,
                "p_stage_mu": np.nan,
                "AUROC_stage_mu": np.nan,
                "D_stage_W": np.nan,
                "p_stage_W": np.nan,
                "AUROC_stage_W": np.nan,
                "n_early": np.nan,
                "n_late": np.nan,
                "recommended_representation": "excluded",
                "h_status": "insufficient_coverage",
            })

            all_records.append(rec)
            continue

        sub = ge.loc[matched]
        z_sub = zscore_rows(sub)

        # mu and W
        mu = z_sub.mean(axis=0)
        W = z_sub.std(axis=0, ddof=1)

        mu_var = float(np.nanvar(mu.values, ddof=1))
        mu_iqr = float(iqr(mu.values))
        W_median = float(np.nanmedian(W.values))
        W_iqr = float(iqr(W.values))

        coherence = compute_coherence(z_sub, max_genes=MAX_GENES_FOR_COHERENCE)
        pca1_var, mean_pca_corr = compute_pca1_stats(z_sub)
        loo_stability = compute_leave_one_gene_stability(
            z_sub,
            max_genes=MAX_GENES_FOR_STABILITY
        )
        boot_stability = compute_bootstrap_stability(
            z_sub,
            n_boot=BOOTSTRAP_N
        )

        D_mu, p_mu, auc_mu, n_early, n_late = compute_stage_D(mu, stage_df)
        D_W, p_W, auc_W, _, _ = compute_stage_D(W, stage_df)

        rec = dict(base)

        rec.update({
            "mu_variance": mu_var,
            "mu_IQR": mu_iqr,
            "W_median": W_median,
            "W_IQR": W_iqr,
            "coherence_mean": coherence,
            "PCA1_variance": pca1_var,
            "mean_PCA_corr": mean_pca_corr,
            "loo_stability": loo_stability,
            "bootstrap_stability": boot_stability,
            "D_stage_mu": D_mu,
            "p_stage_mu": p_mu,
            "AUROC_stage_mu": auc_mu,
            "D_stage_W": D_W,
            "p_stage_W": p_W,
            "AUROC_stage_W": auc_W,
            "n_early": n_early,
            "n_late": n_late,
        })

        all_records.append(rec)

    df = pd.DataFrame(all_records)

    passed = df["pass_min_gene_count"] == True

    if passed.sum() > 0:
        mu_iqr_q05 = df.loc[passed, "mu_IQR"].quantile(0.05)
        W_iqr_q75 = df.loc[passed, "W_IQR"].quantile(0.75)
    else:
        mu_iqr_q05 = np.nan
        W_iqr_q75 = np.nan

    def assign_status(row):
        if not row["pass_min_gene_count"]:
            return "insufficient_coverage", "excluded"

        mu_iqr_val = row["mu_IQR"]
        W_iqr_val = row["W_IQR"]
        coh = row["coherence_mean"]
        pca = row["PCA1_variance"]
        mpc = row["mean_PCA_corr"]
        stable = row["bootstrap_stability"]

        # Low resolution
        if pd.notna(mu_iqr_val) and pd.notna(mu_iqr_q05) and mu_iqr_val <= mu_iqr_q05:
            if pd.notna(W_iqr_val) and pd.notna(W_iqr_q75) and W_iqr_val >= W_iqr_q75:
                return "heterogeneity_resolvable", "W"

            return "low_resolution", "mean_cautious"

        # Mean-ready
        if (
            pd.notna(pca) and pca >= 0.25 and
            pd.notna(mpc) and mpc >= 0.70 and
            pd.notna(stable) and stable >= 0.80
        ):
            return "ready", "mean"

        # W recommended
        if (
            pd.notna(W_iqr_val) and pd.notna(W_iqr_q75) and W_iqr_val >= W_iqr_q75 and
            (
                (pd.notna(coh) and coh < 0.10) or
                (pd.notna(pca) and pca < 0.20)
            )
        ):
            return "partial_W_recommended", "W"

        # PCA recommended
        if pd.notna(pca) and pca >= 0.25 and (pd.isna(mpc) or mpc < 0.70):
            return "partial_PCA_recommended", "PCA"

        # Unstable
        if pd.notna(stable) and stable < 0.60:
            return "unstable_aggregation", "review"

        return "partial_mean_cautious", "mean_cautious"

    status_rep = df.apply(assign_status, axis=1)
    df["h_status"] = [x[0] for x in status_rep]
    df["recommended_representation"] = [x[1] for x in status_rep]

    return df


# ============================================================
# 7. FIGURES
# ============================================================

def save_fig_effective_gene_distribution(df, out_dir):
    fig_path = out_dir / "Figure1_effective_gene_count_distribution.png"

    plt.figure(figsize=(8, 5))

    for db in df["database"].unique():
        sub = df[df["database"] == db]
        vals = sub["matched_gene_count"].dropna().astype(float)

        plt.hist(
            vals,
            bins=50,
            alpha=0.45,
            label=db
        )

    plt.axvline(
        MIN_MATCHED_GENES,
        linestyle="--",
        linewidth=2,
        label=f"minimum matched genes = {MIN_MATCHED_GENES}"
    )

    plt.xlabel("Matched/effective gene count after projection")
    plt.ylabel("Number of BP/gene sets")
    plt.title("AIDO-h representation completeness")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()

    return fig_path


def save_fig_coverage_vs_D(df, out_dir):
    fig_path = out_dir / "Figure2_coverage_fraction_vs_stage_D.png"

    sub = df[
        (df["pass_min_gene_count"] == True) &
        (df["coverage_fraction"].notna()) &
        (df["D_stage_mu"].notna())
    ].copy()

    plt.figure(figsize=(7, 5))

    if sub.empty:
        plt.text(
            0.5,
            0.5,
            "No valid stage D values.\nCheck stage label parsing and Early/Late counts.",
            ha="center",
            va="center",
            fontsize=12
        )
        plt.axis("off")
    else:
        for db in sub["database"].unique():
            sdb = sub[sub["database"] == db]
            plt.scatter(
                sdb["coverage_fraction"],
                sdb["D_stage_mu"],
                s=16,
                alpha=0.55,
                label=db
            )

        plt.xlabel("Coverage fraction")
        plt.ylabel("Stage D from mean BP score, -log10(p)")
        plt.title("Coverage fraction vs stage discriminability")
        plt.legend()

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()

    return fig_path


def save_fig_muD_vs_WD(df, out_dir):
    fig_path = out_dir / "Figure3_mu_D_vs_W_D.png"

    sub = df[
        (df["pass_min_gene_count"] == True) &
        (df["D_stage_mu"].notna()) &
        (df["D_stage_W"].notna())
    ].copy()

    plt.figure(figsize=(6, 6))

    if sub.empty:
        plt.text(
            0.5,
            0.5,
            "No valid mu-D / W-D values.\nCheck stage label parsing and Early/Late counts.",
            ha="center",
            va="center",
            fontsize=12
        )
        plt.axis("off")
    else:
        for db in sub["database"].unique():
            sdb = sub[sub["database"] == db]
            plt.scatter(
                sdb["D_stage_mu"],
                sdb["D_stage_W"],
                s=16,
                alpha=0.55,
                label=db
            )

        max_val = np.nanmax([
            sub["D_stage_mu"].max(),
            sub["D_stage_W"].max()
        ])

        if np.isfinite(max_val) and max_val > 0:
            plt.plot(
                [0, max_val],
                [0, max_val],
                linestyle="--",
                linewidth=1
            )

        plt.xlabel("Stage D from mean activity score, mu")
        plt.ylabel("Stage D from width / heterogeneity score, W")
        plt.title("Mean-score D vs width-score D")
        plt.legend()

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()

    return fig_path


def save_fig_h_status_composition(df, out_dir):
    fig_path = out_dir / "Figure4_AIDOh_status_composition.png"

    counts = df["h_status"].value_counts()

    plt.figure(figsize=(9, 5))
    plt.bar(counts.index.astype(str), counts.values)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Number of BP/gene sets")
    plt.title("AIDO-h status composition")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()

    return fig_path


def save_fig_database_pass_fail(df, out_dir):
    fig_path = out_dir / "Figure5_database_pass_fail_AIDOh.png"

    summary = (
        df.groupby("database")["pass_min_gene_count"]
        .agg(["sum", "count"])
        .reset_index()
    )

    summary["pass"] = summary["sum"]
    summary["fail"] = summary["count"] - summary["sum"]

    x = np.arange(len(summary))
    width = 0.35

    plt.figure(figsize=(7, 5))

    plt.bar(x - width / 2, summary["pass"], width, label="Pass >=10 matched genes")
    plt.bar(x + width / 2, summary["fail"], width, label="Excluded <10 matched genes")

    plt.xticks(x, summary["database"], rotation=30, ha="right")
    plt.ylabel("Number of gene sets")
    plt.title("AIDO-h minimum-gene rule by database")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()

    return fig_path


# ============================================================
# 8. SUMMARY REPORT
# ============================================================

def write_summary_report(df, out_dir):
    report_path = out_dir / "AIDOh_Summary_Report.txt"

    total = len(df)
    pass_n = int(df["pass_min_gene_count"].sum())
    fail_n = total - pass_n

    valid_mu_D = int(df["D_stage_mu"].notna().sum())
    valid_W_D = int(df["D_stage_W"].notna().sum())

    lines = []

    lines.append("AIDO-h Pilot Summary Report")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Output directory: {out_dir}")
    lines.append(f"Minimum matched gene threshold: {MIN_MATCHED_GENES}")
    lines.append("")
    lines.append("Overall AIDO-h completeness")
    lines.append("-" * 70)
    lines.append(f"Total BP/gene sets analyzed: {total}")
    lines.append(f"Passed matched gene cutoff >= {MIN_MATCHED_GENES}: {pass_n}")
    lines.append(f"Excluded due to matched genes < {MIN_MATCHED_GENES}: {fail_n}")

    if total > 0:
        lines.append(f"Exclusion fraction: {fail_n / total:.3f}")

    lines.append("")
    lines.append("Valid D results")
    lines.append("-" * 70)
    lines.append(f"Valid D_stage_mu values: {valid_mu_D}")
    lines.append(f"Valid D_stage_W values: {valid_W_D}")
    lines.append("")

    lines.append("By database")
    lines.append("-" * 70)

    for db, sub in df.groupby("database"):
        t = len(sub)
        p = int(sub["pass_min_gene_count"].sum())
        f = t - p
        valid_d = int(sub["D_stage_mu"].notna().sum())

        lines.append(
            f"{db}: total={t}, pass={p}, excluded={f}, "
            f"excluded_fraction={f/t:.3f}, valid_D_mu={valid_d}"
        )

    lines.append("")
    lines.append("AIDO-h status composition")
    lines.append("-" * 70)

    for k, v in df["h_status"].value_counts().items():
        lines.append(f"{k}: {v}")

    lines.append("")
    lines.append("Top 20 BP by mean-score stage D")
    lines.append("-" * 70)

    top_mu = (
        df[df["pass_min_gene_count"] == True]
        .dropna(subset=["D_stage_mu"])
        .sort_values("D_stage_mu", ascending=False)
        .head(20)
    )

    if top_mu.empty:
        lines.append("No valid D_stage_mu values. Check stage parsing.")
    else:
        for _, r in top_mu.iterrows():
            lines.append(
                f"[{r['database']}] {r['bp_name']} | "
                f"D_mu={r['D_stage_mu']:.3f}, "
                f"D_W={r['D_stage_W']:.3f}, "
                f"matched={int(r['matched_gene_count'])}, "
                f"coverage={r['coverage_fraction']:.3f}, "
                f"h={r['h_status']}, "
                f"rep={r['recommended_representation']}"
            )

    lines.append("")
    lines.append("Top 20 BP by width-score stage D")
    lines.append("-" * 70)

    top_w = (
        df[df["pass_min_gene_count"] == True]
        .dropna(subset=["D_stage_W"])
        .sort_values("D_stage_W", ascending=False)
        .head(20)
    )

    if top_w.empty:
        lines.append("No valid D_stage_W values. Check stage parsing.")
    else:
        for _, r in top_w.iterrows():
            lines.append(
                f"[{r['database']}] {r['bp_name']} | "
                f"D_W={r['D_stage_W']:.3f}, "
                f"D_mu={r['D_stage_mu']:.3f}, "
                f"matched={int(r['matched_gene_count'])}, "
                f"coverage={r['coverage_fraction']:.3f}, "
                f"h={r['h_status']}, "
                f"rep={r['recommended_representation']}"
            )

    lines.append("")
    lines.append("Interpretation guide")
    lines.append("-" * 70)
    lines.append("matched genes < 10: excluded as insufficient observation readiness.")
    lines.append("ready: mean BP score is relatively well-supported.")
    lines.append("partial_W_recommended: internal heterogeneity may be informative.")
    lines.append("low_resolution: current GE interface has limited inter-patient resolution.")
    lines.append("unstable_aggregation: mean score may be sensitive to gene membership.")
    lines.append("")
    lines.append("AIDO-h does not decide whether a BP is biologically important.")
    lines.append("It diagnoses whether the current molecular measurement and representation are suitable for observing it.")

    report_path.write_text("\n".join(lines), encoding="utf-8")

    return report_path


# ============================================================
# 9. MAIN
# ============================================================

def main():
    print("=" * 80)
    print("AIDO-h PILOT TEST FIXED VERSION")
    print("BRCA GE + Hallmark / GO BP / Reactome + Stage Early vs Late")
    print("=" * 80)

    # -----------------------------
    # Load GE
    # -----------------------------
    ge = load_ge_matrix(GE_FILE)

    # -----------------------------
    # Load stage
    # -----------------------------
    stage_df = load_stage_groups(STAGE_FILE)

    # -----------------------------
    # Align patients
    # -----------------------------
    common_patients = sorted(list(set(ge.columns).intersection(set(stage_df["patient"]))))

    print(f"[ALIGN] Common GE-stage patients: {len(common_patients)}")

    if len(common_patients) < 30:
        debug_ge_patients = OUT_DIR / "DEBUG_GE_patients.txt"
        debug_stage_patients = OUT_DIR / "DEBUG_stage_patients.txt"

        debug_ge_patients.write_text("\n".join(sorted(map(str, ge.columns))), encoding="utf-8")
        debug_stage_patients.write_text("\n".join(sorted(map(str, stage_df["patient"]))), encoding="utf-8")

        raise ValueError(
            "Too few common patients between GE and stage file. "
            f"GE patients saved to {debug_ge_patients}; "
            f"stage patients saved to {debug_stage_patients}"
        )

    ge = ge[common_patients]
    stage_df = stage_df[stage_df["patient"].isin(common_patients)].copy()

    print("[ALIGN] Stage counts after GE-stage matching:")
    print(stage_df["stage_group"].value_counts())

    # Save debug alignment
    stage_df.to_csv(OUT_DIR / "DEBUG_stage_groups_after_alignment.csv", index=False)

    # -----------------------------
    # Load and process GMT databases
    # -----------------------------
    all_results = []

    for db_name, gmt_path in GMT_FILES.items():
        gene_sets = load_gmt(gmt_path, db_name)

        df_db = process_gene_sets_for_database(
            database_name=db_name,
            gene_sets=gene_sets,
            ge=ge,
            stage_df=stage_df
        )

        df_db.to_csv(OUT_DIR / f"AIDOh_{db_name}_full_table.csv", index=False)
        all_results.append(df_db)

    df_all = pd.concat(all_results, ignore_index=True)

    # -----------------------------
    # Output tables
    # -----------------------------
    completeness_cols = [
        "database",
        "bp_id",
        "bp_name",
        "original_gene_count",
        "matched_gene_count",
        "coverage_fraction",
        "pass_min_gene_count",
        "h_status"
    ]

    representation_cols = [
        "database",
        "bp_id",
        "bp_name",
        "mu_variance",
        "mu_IQR",
        "W_median",
        "W_IQR",
        "coherence_mean",
        "PCA1_variance",
        "mean_PCA_corr",
        "loo_stability",
        "bootstrap_stability",
        "recommended_representation",
        "h_status"
    ]

    d_cols = [
        "database",
        "bp_id",
        "bp_name",
        "matched_gene_count",
        "coverage_fraction",
        "h_status",
        "recommended_representation",
        "D_stage_mu",
        "p_stage_mu",
        "AUROC_stage_mu",
        "D_stage_W",
        "p_stage_W",
        "AUROC_stage_W",
        "n_early",
        "n_late"
    ]

    df_all.to_csv(
        OUT_DIR / "AIDOh_BRCA_GE_Stage_FULL_TABLE.csv",
        index=False
    )

    df_all[completeness_cols].to_csv(
        OUT_DIR / "AIDOh_BP_Completeness_Table.csv",
        index=False
    )

    df_all[representation_cols].to_csv(
        OUT_DIR / "AIDOh_Representation_Table.csv",
        index=False
    )

    df_all[d_cols].to_csv(
        OUT_DIR / "AIDOh_D_Comparison_Table.csv",
        index=False
    )

    # Top rankings
    raw_top = (
        df_all
        .dropna(subset=["D_stage_mu"])
        .sort_values("D_stage_mu", ascending=False)
        .head(100)
        .copy()
    )

    filtered_top = (
        df_all[
            (df_all["pass_min_gene_count"] == True) &
            (~df_all["h_status"].isin(["insufficient_coverage", "unstable_aggregation"]))
        ]
        .dropna(subset=["D_stage_mu"])
        .sort_values("D_stage_mu", ascending=False)
        .head(100)
        .copy()
    )

    raw_top.to_csv(
        OUT_DIR / "AIDOh_Top100_Raw_MuD.csv",
        index=False
    )

    filtered_top.to_csv(
        OUT_DIR / "AIDOh_Top100_AIDOh_Filtered_MuD.csv",
        index=False
    )

    # Database summary
    db_summary = (
        df_all.groupby("database")
        .agg(
            total_gene_sets=("bp_id", "count"),
            passed_min_genes=("pass_min_gene_count", "sum"),
            median_original_genes=("original_gene_count", "median"),
            median_matched_genes=("matched_gene_count", "median"),
            median_coverage=("coverage_fraction", "median"),
            valid_D_mu=("D_stage_mu", lambda x: x.notna().sum()),
            valid_D_W=("D_stage_W", lambda x: x.notna().sum())
        )
        .reset_index()
    )

    db_summary["excluded_min_genes"] = (
        db_summary["total_gene_sets"] - db_summary["passed_min_genes"]
    )

    db_summary["excluded_fraction"] = (
        db_summary["excluded_min_genes"] / db_summary["total_gene_sets"]
    )

    db_summary.to_csv(
        OUT_DIR / "AIDOh_Database_Observation_Readiness_Summary.csv",
        index=False
    )

    # Excel workbook
    excel_path = OUT_DIR / "AIDOh_BRCA_GE_Stage_summary_tables.xlsx"

    try:
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            df_all.to_excel(writer, sheet_name="FULL", index=False)
            df_all[completeness_cols].to_excel(writer, sheet_name="Completeness", index=False)
            df_all[representation_cols].to_excel(writer, sheet_name="Representation", index=False)
            df_all[d_cols].to_excel(writer, sheet_name="D_Comparison", index=False)
            db_summary.to_excel(writer, sheet_name="Database_Summary", index=False)
            raw_top.to_excel(writer, sheet_name="Top100_Raw_MuD", index=False)
            filtered_top.to_excel(writer, sheet_name="Top100_hFiltered_MuD", index=False)

        print(f"[OK] Excel workbook saved: {excel_path}")

    except Exception as e:
        print(f"[WARN] Excel workbook not saved: {e}")

    # -----------------------------
    # Figures
    # -----------------------------
    fig1 = save_fig_effective_gene_distribution(df_all, OUT_DIR)
    fig2 = save_fig_coverage_vs_D(df_all, OUT_DIR)
    fig3 = save_fig_muD_vs_WD(df_all, OUT_DIR)
    fig4 = save_fig_h_status_composition(df_all, OUT_DIR)
    fig5 = save_fig_database_pass_fail(df_all, OUT_DIR)

    # -----------------------------
    # Summary report
    # -----------------------------
    report_path = write_summary_report(df_all, OUT_DIR)

    print("\n" + "=" * 80)
    print("[DONE] AIDO-h pilot completed.")
    print(f"Output directory: {OUT_DIR}")
    print("")
    print("Main outputs:")
    print("  FULL table:       ", OUT_DIR / "AIDOh_BRCA_GE_Stage_FULL_TABLE.csv")
    print("  Completeness:     ", OUT_DIR / "AIDOh_BP_Completeness_Table.csv")
    print("  Representation:   ", OUT_DIR / "AIDOh_Representation_Table.csv")
    print("  D comparison:     ", OUT_DIR / "AIDOh_D_Comparison_Table.csv")
    print("  DB summary:       ", OUT_DIR / "AIDOh_Database_Observation_Readiness_Summary.csv")
    print("  Excel workbook:   ", excel_path)
    print("  Summary report:   ", report_path)
    print("")
    print("Figures:")
    print("  Figure 1:", fig1)
    print("  Figure 2:", fig2)
    print("  Figure 3:", fig3)
    print("  Figure 4:", fig4)
    print("  Figure 5:", fig5)
    print("=" * 80)


if __name__ == "__main__":
    main()
