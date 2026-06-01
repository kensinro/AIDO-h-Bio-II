# ============================================================
# AIDO-h-Biology II | Main manuscript figures 2�C5
# Data-backed manuscript figures from latest experiment outputs
# ============================================================

from pathlib import Path
import os
import re
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# -----------------------------
# 0. User paths
# -----------------------------
MAIN_ROOT = Path(os.environ.get("AIDOH_MAIN_ROOT", r"D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB"))
RANDOM_ROOT = Path(os.environ.get("AIDOH_RANDOM_ROOT", r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300"))

MAIN_TABLE_DIR = MAIN_ROOT / "tables"
RANDOM_TABLE_DIR = RANDOM_ROOT / "tables"

OUT_DIR = Path(os.environ.get("AIDOH_FIGURE_OUT", r"D:/AIDO-Temp/AIDO-h-Biology-II-Manuscript-Figures"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# 1. Style settings
# -----------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 16,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

COL_READY = "#1f77b4"
COL_LOW = "#ff7f0e"
COL_NEAR = "#2ca02c"
COL_PURPLE = "#7b4ab8"
COL_GRAY = "#666666"
COL_LIGHTBLUE = "#aec7e8"

DATABASE_ORDER = ["Hallmark", "Reactome", "GO_BP"]
DATABASE_LABEL = {"GO_BP": "GO BP", "Hallmark": "Hallmark", "Reactome": "Reactome"}

# -----------------------------
# 2. Helper functions
# -----------------------------
def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return path

def clean_database_name(x):
    if pd.isna(x):
        return x
    x = str(x)
    if x.upper() in ["GO_BP", "GO BP"]:
        return "GO_BP"
    if x.lower() == "hallmark":
        return "Hallmark"
    if x.lower() == "reactome":
        return "Reactome"
    return x

def safe_col(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise KeyError(f"None of these columns found: {candidates}\nAvailable columns: {list(df.columns)}")
    return None

def save_figure(fig, name):
    png = OUT_DIR / f"{name}.png"
    pdf = OUT_DIR / f"{name}.pdf"
    svg = OUT_DIR / f"{name}.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    print(f"Saved: {png}")
    print(f"Saved: {pdf}")
    print(f"Saved: {svg}")

def comma_int(x):
    if pd.isna(x):
        return "NA"
    return f"{int(round(float(x))):,}"

def fmt3(x):
    if pd.isna(x):
        return "NA"
    return f"{float(x):.3f}"

def short_term_name(s, max_len=42):
    s = str(s)
    s = s.replace("HALLMARK_", "")
    s = s.replace("GOBP_", "")
    s = s.replace("REACTOME_", "")
    s = s.replace("_", " ")
    s = s.title()
    if len(s) > max_len:
        s = s[:max_len - 1] + "..."
    return s

def add_panel_label(ax, label):
    ax.text(
        -0.12, 1.08, label,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="left"
    )

def horizontal_bar_labels(ax, values, y_positions, color=COL_READY, dx_frac=0.01, fmt="int"):
    xmax = max(values) if len(values) else 1
    for v, y in zip(values, y_positions):
        if fmt == "int":
            label = comma_int(v)
        elif fmt == "float3":
            label = fmt3(v)
        else:
            label = str(v)
        ax.text(v + xmax * dx_frac, y, label, va="center", ha="left",
                fontsize=10, fontweight="bold", color=color)

# -----------------------------
# 3. Load main tables
# -----------------------------
table1_path = require_file(MAIN_TABLE_DIR / "Table_1_cancer_database_observation_readiness_summary.csv")
table2_path = require_file(MAIN_TABLE_DIR / "Table_2_database_overall_observation_readiness_summary.csv")
full_path = require_file(MAIN_TABLE_DIR / "Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv")

df_cd = pd.read_csv(table1_path)
df_db = pd.read_csv(table2_path)
df_full = pd.read_csv(full_path)

for df in [df_cd, df_db, df_full]:
    if "database" in df.columns:
        df["database"] = df["database"].map(clean_database_name)

# -----------------------------
# 4. Load random baseline tables
# -----------------------------
random_term_path = RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_term_level_mapped_results.csv"
random_db_path = RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_database_summary.csv"
random_cd_path = RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_cancer_database_summary.csv"
random_bin_path = RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_bin_level_baselines.csv"
random_top_path = RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_top_structured_favored_terms.csv"
random_high_not_path = RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_highD_but_not_above_random_q95.csv"

df_rand_term = pd.read_csv(require_file(random_term_path))
df_rand_db = pd.read_csv(require_file(random_db_path))
df_rand_cd = pd.read_csv(require_file(random_cd_path))
df_rand_bin = pd.read_csv(require_file(random_bin_path))

df_rand_top = pd.read_csv(random_top_path) if random_top_path.exists() else df_rand_term.copy()
df_high_not = pd.read_csv(random_high_not_path) if random_high_not_path.exists() else pd.DataFrame()

for df in [df_rand_term, df_rand_db, df_rand_cd, df_rand_bin, df_rand_top, df_high_not]:
    if isinstance(df, pd.DataFrame) and "database" in df.columns:
        df["database"] = df["database"].map(clean_database_name)

# ============================================================
# Figure 2. Observation-readiness landscape across databases
# ============================================================

def make_figure_2():
    d = df_db.copy()
    d["database"] = d["database"].map(clean_database_name)
    d = d.set_index("database").loc[DATABASE_ORDER].reset_index()

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))
    fig.suptitle(
        "Figure 2. Pan-cancer observation-readiness landscape across annotation databases",
        y=0.98, fontsize=17, fontweight="bold"
    )
    fig.text(0.5, 0.945, "Summary across 25 TCGA cancer labels",
             ha="center", va="center", fontsize=12, style="italic", color=COL_GRAY)

    # A. Total observations
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    order = ["GO_BP", "Reactome", "Hallmark"]
    da = d.set_index("database").loc[order].reset_index()
    y = np.arange(len(da))
    vals = da["n_bp_cancer_observations"].values
    labels = [f"{DATABASE_LABEL[x]}\n(n = {comma_int(n)})" for x, n in zip(da["database"], da["n_unique_bp"])]
    ax.barh(y, vals, color=COL_READY, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Number of cancer�Cdatabase�Cterm observations")
    ax.set_title("Total observations by database", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    horizontal_bar_labels(ax, vals, y, color=COL_READY, dx_frac=0.02, fmt="int")

    # B. Fraction observation-ready
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    db = d.set_index("database").loc[DATABASE_ORDER].reset_index()
    y = np.arange(len(db))
    vals = db["fraction_observation_ready"].values
    ax.barh(y, vals, color=COL_READY, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([DATABASE_LABEL[x] for x in db["database"]])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction observation-ready")
    ax.set_title("Fraction observation-ready", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    for v, yy in zip(vals, y):
        ax.text(v + 0.02, yy, f"{v:.3f}", va="center", ha="left",
                fontsize=10, fontweight="bold", color=COL_READY)

    # C. Observation-readiness class counts
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    dc = d.set_index("database").loc[["GO_BP", "Hallmark", "Reactome"]].reset_index()
    x = np.arange(len(dc))
    ready = dc["n_observation_ready"].values
    low = dc["n_low_resolution"].values
    near = dc["n_near_unobservable"].values
    ax.bar(x, ready, color=COL_READY, label="Observation-ready")
    ax.bar(x, low, bottom=ready, color=COL_LOW, label="Low-resolution")
    ax.bar(x, near, bottom=ready+low, color=COL_NEAR, label="Near-unobservable")
    ax.set_xticks(x)
    ax.set_xticklabels([DATABASE_LABEL[z] for z in dc["database"]])
    ax.set_ylabel("Number of cancer�Cdatabase�Cterm observations")
    ax.set_title("Observation-readiness class counts", fontweight="bold")
    ax.legend(frameon=True)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    for i in range(len(dc)):
        if ready[i] > 0:
            ax.text(x[i], ready[i] / 2, comma_int(ready[i]), ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
        if low[i] > 0:
            ax.text(x[i], ready[i] + low[i] / 2, comma_int(low[i]), ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
        if near[i] > 0:
            ax.text(x[i], ready[i] + low[i] + near[i] + max(ready+low+near)*0.015,
                    comma_int(near[i]), ha="center", va="bottom",
                    fontsize=9, color=COL_NEAR, fontweight="bold")

    # D. Effective gene-set representation
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    dd = d.set_index("database").loc[["Hallmark", "Reactome", "GO_BP"]].reset_index()
    y = np.arange(len(dd))
    h = 0.28
    ndef = dd["median_N_defined"].values
    nmatch = dd["median_N_matched"].values
    frac = dd["median_matched_fraction"].values

    ax.barh(y - h/2, ndef, height=h, color=COL_READY, label=r"Median $N_{\mathrm{defined}}$")
    ax.barh(y + h/2, nmatch, height=h, color=COL_LIGHTBLUE, label=r"Median $N_{\mathrm{matched}}$")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{DATABASE_LABEL[z]}\n(n = {comma_int(n)})"
                        for z, n in zip(dd["database"], dd["n_unique_bp"])])
    ax.invert_yaxis()
    ax.set_xlabel("Median genes per gene set")
    ax.set_title("Effective gene-set representation", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left", ncol=2)

    xmax = max(ndef) * 1.25
    ax.set_xlim(0, xmax)
    for i, (a, b, f) in enumerate(zip(ndef, nmatch, frac)):
        ax.text(a + xmax*0.01, i - h/2, f"{int(a)}", va="center", ha="left",
                fontsize=9, color=COL_READY, fontweight="bold")
        ax.text(b + xmax*0.01, i + h/2, f"{int(b)}", va="center", ha="left",
                fontsize=9, color=COL_LIGHTBLUE, fontweight="bold")
        ax.text(xmax*0.94, i, f"{f:.3f}", va="center", ha="center",
                fontsize=10, color="black", fontweight="bold")
    ax.text(xmax*0.94, -0.55, "Matched\nfraction", ha="center", va="bottom",
            fontsize=9, color="black")

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    save_figure(fig, "Figure_2_observation_readiness_landscape")
    plt.close(fig)

# ============================================================
# Figure 3. Survival discriminability after filtering
# ============================================================

def make_figure_3():
    d = df_full.copy()

    db_col = safe_col(d, ["database"])
    class_col = safe_col(d, ["observation_class", "readiness_class", "observation_readiness_class"], required=False)
    d_col = safe_col(d, ["D_survival", "Dsurvival", "D", "survival_D"], required=False)
    nmatch_col = safe_col(d, ["N_matched", "n_matched", "matched_gene_count", "Nmatched"], required=False)
    cancer_col = safe_col(d, ["cancer", "cancer_label"], required=False)
    term_col = safe_col(d, ["term", "bp", "gene_set", "gene_set_name", "term_name"], required=False)

    if d_col is None:
        raise KeyError("Cannot find survival discriminability column in full result table.")

    d[db_col] = d[db_col].map(clean_database_name)
    d[d_col] = pd.to_numeric(d[d_col], errors="coerce")
    if nmatch_col:
        d[nmatch_col] = pd.to_numeric(d[nmatch_col], errors="coerce")

    if class_col:
        obs_ready = d[d[class_col].astype(str).str.lower().isin(
            ["observation_ready", "observation-ready", "observation ready"]
        )].copy()
    elif nmatch_col:
        obs_ready = d[d[nmatch_col] >= 10].copy()
    else:
        obs_ready = d.copy()

    obs_ready_valid = obs_ready[np.isfinite(obs_ready[d_col])].copy()

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))
    fig.suptitle("Figure 3. Survival-discriminability landscape after observation-readiness filtering",
                 y=0.98, fontsize=17, fontweight="bold")
    fig.text(0.5, 0.945, r"Observation-ready terms; nominal high-D threshold: $D_{\mathrm{survival}} = 1.301$",
             ha="center", va="center", fontsize=12, style="italic", color=COL_GRAY)

    # A. D distributions by database
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    data = [obs_ready_valid.loc[obs_ready_valid[db_col] == db, d_col].dropna().values
            for db in DATABASE_ORDER]
    bp = ax.boxplot(data, labels=[DATABASE_LABEL[x] for x in DATABASE_ORDER],
                    patch_artist=True, showfliers=False)
    for box in bp["boxes"]:
        box.set(facecolor="#d9e8f5", edgecolor=COL_READY)
    for med in bp["medians"]:
        med.set(color=COL_LOW, linewidth=2)
    ax.axhline(1.301, color=COL_GRAY, linestyle="--", linewidth=1.4)
    ax.set_ylabel(r"$D_{\mathrm{survival}}$")
    ax.set_title("Distribution of survival discriminability", fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    # B. Max D by cancer and database, heatmap
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    if cancer_col is None:
        ax.text(0.5, 0.5, "Cancer column not found", ha="center", va="center")
    else:
        pivot = obs_ready_valid.pivot_table(index=cancer_col, columns=db_col, values=d_col, aggfunc="max")
        pivot = pivot[[c for c in DATABASE_ORDER if c in pivot.columns]]
        pivot = pivot.sort_index()
        im = ax.imshow(pivot.values, aspect="auto", cmap="Blues")
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels([DATABASE_LABEL[x] for x in pivot.columns], rotation=30, ha="right")
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels(pivot.index, fontsize=7)
        ax.set_title(r"Maximum $D_{\mathrm{survival}}$ by cancer", fontweight="bold")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label(r"Max $D_{\mathrm{survival}}$")

    # C. High-D count by cancer and database
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    if cancer_col is None:
        ax.text(0.5, 0.5, "Cancer column not found", ha="center", va="center")
    else:
        tmp = obs_ready_valid.copy()
        tmp["highD"] = tmp[d_col] >= 1.301
        count = tmp.pivot_table(index=cancer_col, columns=db_col, values="highD", aggfunc="sum", fill_value=0)
        count = count[[c for c in DATABASE_ORDER if c in count.columns]]
        count["total"] = count.sum(axis=1)
        count = count.sort_values("total", ascending=True).drop(columns="total")
        bottom = np.zeros(len(count))
        colors = {"Hallmark": "#4c78a8", "Reactome": "#72b7b2", "GO_BP": "#f58518"}
        for db in [c for c in DATABASE_ORDER if c in count.columns]:
            vals = count[db].values
            ax.barh(count.index, vals, left=bottom, label=DATABASE_LABEL[db], color=colors.get(db, COL_READY))
            bottom += vals
        ax.set_xlabel(r"Number of nominal high-D terms")
        ax.set_title(r"Nominal high-D terms by cancer", fontweight="bold")
        ax.legend(frameon=False)
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", labelsize=7)

    # D. Nmatched vs D survival scatter
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    if nmatch_col is None:
        ax.text(0.5, 0.5, "Matched-gene column not found", ha="center", va="center")
    else:
        sample = obs_ready_valid.dropna(subset=[nmatch_col, d_col]).copy()
        if len(sample) > 12000:
            sample = sample.sample(12000, random_state=1)
        color_map = {"Hallmark": "#4c78a8", "Reactome": "#72b7b2", "GO_BP": "#f58518"}
        for db in DATABASE_ORDER:
            sub = sample[sample[db_col] == db]
            ax.scatter(sub[nmatch_col], sub[d_col], s=8, alpha=0.35,
                       label=DATABASE_LABEL[db], color=color_map.get(db, COL_READY), edgecolors="none")
        ax.axhline(1.301, color=COL_GRAY, linestyle="--", linewidth=1.4)
        ax.axvline(10, color=COL_GRAY, linestyle="--", linewidth=1.4)
        ax.set_xlabel(r"Effective matched-gene count, $N_{\mathrm{matched}}$")
        ax.set_ylabel(r"$D_{\mathrm{survival}}$")
        ax.set_title("Representation versus survival discriminability", fontweight="bold")
        ax.legend(frameon=False, markerscale=2)
        ax.grid(alpha=0.25)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    save_figure(fig, "Figure_3_survival_discriminability_landscape")
    plt.close(fig)

# ============================================================
# Figure 4. Size-bin random baseline calibration
# ============================================================

def make_figure_4():
    term = df_rand_term.copy()
    dbsum = df_rand_db.copy()
    cdsum = df_rand_cd.copy()

    for df in [term, dbsum, cdsum]:
        if "database" in df.columns:
            df["database"] = df["database"].map(clean_database_name)

    db_col = safe_col(term, ["database"])
    cancer_col = safe_col(term, ["cancer", "cancer_label"], required=False)
    dreal_col = safe_col(term, ["D_real", "Dreal", "D_survival", "real_D", "D"], required=False)
    q95_col = safe_col(term, ["D_random_q95", "Drandom_q95", "random_q95", "D_random, q95", "D_random_q95_mapped"], required=False)
    sf_col = safe_col(term, ["structured_favored", "is_structured_favored", "above_random_q95"], required=False)
    dexcess_col = safe_col(term, ["D_excess_q95", "Dexcess_q95", "D_excess_over_q95"], required=False)

    if dreal_col is None:
        raise KeyError("Cannot find real D column in mapped random baseline table.")

    term[dreal_col] = pd.to_numeric(term[dreal_col], errors="coerce")
    if q95_col:
        term[q95_col] = pd.to_numeric(term[q95_col], errors="coerce")
    if dexcess_col:
        term[dexcess_col] = pd.to_numeric(term[dexcess_col], errors="coerce")

    if sf_col is None:
        if q95_col:
            term["structured_favored_calc"] = term[dreal_col] > term[q95_col]
            sf_col = "structured_favored_calc"
        elif dexcess_col:
            term["structured_favored_calc"] = term[dexcess_col] > 0
            sf_col = "structured_favored_calc"
        else:
            raise KeyError("Cannot infer structured-favored column.")

    term[sf_col] = term[sf_col].astype(str).str.lower().isin(["true", "1", "yes", "y"]) | (term[sf_col] == True)

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))
    fig.suptitle("Figure 4. Size-bin random-baseline calibration of observation-ready terms",
                 y=0.98, fontsize=17, fontweight="bold")
    fig.text(0.5, 0.945, r"Structured-favored criterion: $D_{\mathrm{real}} > D_{\mathrm{random},q95}$",
             ha="center", va="center", fontsize=12, style="italic", color=COL_GRAY)

    # A. Structured-favored fraction by database
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    frac = term.groupby(db_col)[sf_col].mean().reindex(DATABASE_ORDER)
    y = np.arange(len(frac))
    ax.barh(y, frac.values, color=COL_PURPLE, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([DATABASE_LABEL[x] for x in frac.index])
    ax.invert_yaxis()
    ax.set_xlim(0, max(0.05, np.nanmax(frac.values)*1.35))
    ax.set_xlabel("Structured-favored fraction")
    ax.set_title("Structured-favored fraction by database", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    for v, yy in zip(frac.values, y):
        ax.text(v + ax.get_xlim()[1]*0.02, yy, f"{v:.3f}", va="center",
                ha="left", fontsize=10, fontweight="bold", color=COL_PURPLE)

    # B. Dreal vs random q95
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    if q95_col is None:
        ax.text(0.5, 0.5, "Random q95 column not found", ha="center", va="center")
    else:
        plot = term.dropna(subset=[dreal_col, q95_col]).copy()
        if len(plot) > 15000:
            plot = plot.sample(15000, random_state=1)
        colors = np.where(plot[sf_col], COL_PURPLE, "#999999")
        ax.scatter(plot[q95_col], plot[dreal_col], s=10, alpha=0.35, c=colors, edgecolors="none")
        lim = max(plot[q95_col].max(), plot[dreal_col].max()) * 1.05
        ax.plot([0, lim], [0, lim], color=COL_GRAY, linestyle="--", linewidth=1.3)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel(r"Size-bin random q95, $D_{\mathrm{random},q95}$")
        ax.set_ylabel(r"Real term discriminability, $D_{\mathrm{real}}$")
        ax.set_title(r"$D_{\mathrm{real}}$ versus size-bin random q95", fontweight="bold")
        ax.grid(alpha=0.25)

    # C. Cancer-database heatmap of structured-favored fraction
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    if cancer_col is None:
        cancer_col = safe_col(cdsum, ["cancer", "cancer_label"], required=False)
    if cancer_col is None:
        ax.text(0.5, 0.5, "Cancer column not found", ha="center", va="center")
    else:
        heat = term.pivot_table(index=cancer_col, columns=db_col, values=sf_col, aggfunc="mean")
        heat = heat[[c for c in DATABASE_ORDER if c in heat.columns]]
        heat = heat.sort_index()
        im = ax.imshow(heat.values, aspect="auto", cmap="Purples", vmin=0, vmax=np.nanmax(heat.values))
        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels([DATABASE_LABEL[x] for x in heat.columns], rotation=30, ha="right")
        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels(heat.index, fontsize=7)
        ax.set_title("Structured-favored fraction by cancer and database", fontweight="bold")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("Structured-favored fraction")

    # D. High-D and structured-favored overlap
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    term["highD_calc"] = term[dreal_col] >= 1.301
    summary = term.groupby(db_col).agg(
        observation_ready=("highD_calc", "size"),
        highD=("highD_calc", "sum"),
        structured_favored=(sf_col, "sum"),
        both=("highD_calc", lambda x: 0)
    )
    both = term.groupby(db_col).apply(lambda x: ((x["highD_calc"]) & (x[sf_col])).sum())
    summary["both"] = both
    summary = summary.reindex(DATABASE_ORDER)

    x = np.arange(len(summary))
    w = 0.25
    ax.bar(x - w, summary["highD"], width=w, color=COL_READY, label="Nominal high-D")
    ax.bar(x, summary["structured_favored"], width=w, color=COL_PURPLE, label="Structured-favored")
    ax.bar(x + w, summary["both"], width=w, color=COL_LOW, label="Both")
    ax.set_xticks(x)
    ax.set_xticklabels([DATABASE_LABEL[z] for z in summary.index])
    ax.set_ylabel("Number of observation-ready terms")
    ax.set_title("Overlap of high-D and structured-favored terms", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    save_figure(fig, "Figure_4_size_bin_random_baseline_calibration")
    plt.close(fig)

# ============================================================
# Supplementary Figure S2 code block retained for reference; not run by this main-figure script
# ============================================================

def make_figure_5():
    term = df_rand_term.copy()
    for df in [term]:
        if "database" in df.columns:
            df["database"] = df["database"].map(clean_database_name)

    db_col = safe_col(term, ["database"])
    cancer_col = safe_col(term, ["cancer", "cancer_label"], required=False)
    term_col = safe_col(term, ["term", "bp", "gene_set", "gene_set_name", "term_name"], required=False)
    dreal_col = safe_col(term, ["D_real", "Dreal", "D_survival", "real_D", "D"], required=False)
    q95_col = safe_col(term, ["D_random_q95", "Drandom_q95", "random_q95", "D_random_q95_mapped"], required=False)
    dmedian_col = safe_col(term, ["D_random_median", "Drandom_median", "random_median"], required=False)
    dexcess_q95_col = safe_col(term, ["D_excess_q95", "Dexcess_q95", "D_excess_over_q95"], required=False)
    dexcess_med_col = safe_col(term, ["D_excess_median", "Dexcess_median", "D_excess_over_median"], required=False)
    sf_col = safe_col(term, ["structured_favored", "is_structured_favored", "above_random_q95"], required=False)

    if dreal_col is None:
        raise KeyError("Cannot find D_real / D_survival column in random mapped table.")
    if term_col is None:
        raise KeyError("Cannot find term name column in random mapped table.")

    term[dreal_col] = pd.to_numeric(term[dreal_col], errors="coerce")
    if q95_col:
        term[q95_col] = pd.to_numeric(term[q95_col], errors="coerce")
    if dmedian_col:
        term[dmedian_col] = pd.to_numeric(term[dmedian_col], errors="coerce")

    if dexcess_q95_col:
        term[dexcess_q95_col] = pd.to_numeric(term[dexcess_q95_col], errors="coerce")
    elif q95_col:
        term["D_excess_q95_calc"] = term[dreal_col] - term[q95_col]
        dexcess_q95_col = "D_excess_q95_calc"
    else:
        raise KeyError("Cannot calculate D excess over q95.")

    if dexcess_med_col:
        term[dexcess_med_col] = pd.to_numeric(term[dexcess_med_col], errors="coerce")
    elif dmedian_col:
        term["D_excess_median_calc"] = term[dreal_col] - term[dmedian_col]
        dexcess_med_col = "D_excess_median_calc"

    if sf_col is None:
        term["structured_favored_calc"] = term[dexcess_q95_col] > 0
        sf_col = "structured_favored_calc"
    else:
        term[sf_col] = term[sf_col].astype(str).str.lower().isin(["true", "1", "yes", "y"]) | (term[sf_col] == True)

    term["highD"] = term[dreal_col] >= 1.301

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5))
    fig.suptitle("Supplementary Figure S2. Term-level prioritization after random-baseline calibration",
                 y=0.98, fontsize=17, fontweight="bold")
    fig.text(0.5, 0.945, "Observation-ready structured terms ranked by survival discriminability and random-baseline excess",
             ha="center", va="center", fontsize=12, style="italic", color=COL_GRAY)

    # A. Top structured-favored terms by excess q95
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    top = term[term[sf_col]].dropna(subset=[dexcess_q95_col]).sort_values(dexcess_q95_col, ascending=False).head(15)
    if len(top) == 0:
        ax.text(0.5, 0.5, "No structured-favored terms found", ha="center", va="center")
    else:
        labels = []
        for _, r in top.iterrows():
            t = short_term_name(r[term_col], 35)
            if cancer_col:
                labels.append(f"{str(r[cancer_col])}: {t}")
            else:
                labels.append(t)
        y = np.arange(len(top))
        ax.barh(y, top[dexcess_q95_col].values, color=COL_PURPLE)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel(r"Excess over random q95")
        ax.set_title(r"Top structured-favored terms by $D_{\mathrm{excess},q95}$", fontweight="bold")
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        ax.set_axisbelow(True)

    # B. Top real D terms and whether they exceed q95
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    topD = term.dropna(subset=[dreal_col]).sort_values(dreal_col, ascending=False).head(15)
    labels = []
    for _, r in topD.iterrows():
        t = short_term_name(r[term_col], 34)
        if cancer_col:
            labels.append(f"{str(r[cancer_col])}: {t}")
        else:
            labels.append(t)
    y = np.arange(len(topD))
    colors = np.where(topD[sf_col], COL_PURPLE, COL_GRAY)
    ax.barh(y, topD[dreal_col].values, color=colors)
    ax.axvline(1.301, color=COL_GRAY, linestyle="--", linewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(r"$D_{\mathrm{real}}$")
    ax.set_title("Top terms by real survival discriminability", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(
        handles=[
            Patch(facecolor=COL_PURPLE, label="Above random q95"),
            Patch(facecolor=COL_GRAY, label="Not above random q95")
        ],
        frameon=False, loc="lower right"
    )

    # C. Database-level fractions: highD, structured favored, both
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    summary = []
    for db in DATABASE_ORDER:
        sub = term[term[db_col] == db]
        if len(sub) == 0:
            continue
        summary.append({
            "database": db,
            "High-D": sub["highD"].mean(),
            "Structured-favored": sub[sf_col].mean(),
            "Both": ((sub["highD"]) & (sub[sf_col])).mean()
        })
    s = pd.DataFrame(summary).set_index("database")
    x = np.arange(len(s))
    width = 0.25
    ax.bar(x - width, s["High-D"], width=width, color=COL_READY, label="High-D")
    ax.bar(x, s["Structured-favored"], width=width, color=COL_PURPLE, label="Structured-favored")
    ax.bar(x + width, s["Both"], width=width, color=COL_LOW, label="Both")
    ax.set_xticks(x)
    ax.set_xticklabels([DATABASE_LABEL[z] for z in s.index])
    ax.set_ylabel("Fraction of observation-ready terms")
    ax.set_title("Fractions of prioritized term classes", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    # D. Excess over q95 distribution
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    data = [term.loc[term[db_col] == db, dexcess_q95_col].dropna().values for db in DATABASE_ORDER]
    bp = ax.boxplot(data, labels=[DATABASE_LABEL[x] for x in DATABASE_ORDER],
                    patch_artist=True, showfliers=False)
    for box in bp["boxes"]:
        box.set(facecolor="#eadff5", edgecolor=COL_PURPLE)
    for med in bp["medians"]:
        med.set(color=COL_LOW, linewidth=2)
    ax.axhline(0, color=COL_GRAY, linestyle="--", linewidth=1.2)
    ax.set_ylabel(r"$D_{\mathrm{real}} - D_{\mathrm{random},q95}$")
    ax.set_title("Distribution of excess over random q95", fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    save_figure(fig, "Supplementary_Figure_S2_term_prioritization_after_random_calibration")
    plt.close(fig)

# -----------------------------
# Run all figures
# -----------------------------
if __name__ == "__main__":
    print("Main table directory:", MAIN_TABLE_DIR)
    print("Random table directory:", RANDOM_TABLE_DIR)
    print("Output directory:", OUT_DIR)

    make_figure_2()
    make_figure_3()
    make_figure_4()

    print("\nDone. Main manuscript Figures 2-4 saved to:")
    print(OUT_DIR)
