# ============================================================
# AIDO-h-Biology II | Fix and regenerate Figure 4 and Figure 5
# Data-backed manuscript figures from random baseline outputs
#
# Input:
#   D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300/tables
#
# Output:
#   D:/AIDO-Temp/AIDO-h-Biology-II-Manuscript-Figures
#
# Figures:
#   Figure 4. Size-bin random-baseline calibration
#   Figure 5. Term-level prioritization after random calibration
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ============================================================
# 0. Paths
# ============================================================

RANDOM_ROOT = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300")
RANDOM_TABLE_DIR = RANDOM_ROOT / "tables"

OUT_DIR = Path(r"D:/AIDO-Temp/AIDO-h-Biology-II-Manuscript-Figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. Style
# ============================================================

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
COL_LIGHT_GRAY = "#999999"
COL_TEAL = "#72b7b2"

DATABASE_ORDER = ["Hallmark", "Reactome", "GO_BP"]
DATABASE_LABEL = {
    "GO_BP": "GO BP",
    "Hallmark": "Hallmark",
    "Reactome": "Reactome",
}


# ============================================================
# 2. Helper functions
# ============================================================

def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return path


def clean_database_name(x):
    if pd.isna(x):
        return x
    x = str(x).strip()
    if x.upper() in ["GO_BP", "GO BP", "GOBP"]:
        return "GO_BP"
    if x.lower() == "hallmark":
        return "Hallmark"
    if x.lower() == "reactome":
        return "Reactome"
    return x


def safe_col(df: pd.DataFrame, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c

    if required:
        raise KeyError(
            f"None of these columns found: {candidates}\n"
            f"Available columns:\n{list(df.columns)}"
        )

    return None


def fuzzy_col(df: pd.DataFrame, include_any=None, include_all=None, exclude_any=None):
    """
    Flexible column finder.
    include_any: at least one of these substrings must appear.
    include_all: all of these substrings must appear.
    exclude_any: none of these substrings should appear.
    """
    include_any = include_any or []
    include_all = include_all or []
    exclude_any = exclude_any or []

    hits = []

    for c in df.columns:
        lc = c.lower()

        ok_any = True
        if include_any:
            ok_any = any(x.lower() in lc for x in include_any)

        ok_all = all(x.lower() in lc for x in include_all)
        ok_exclude = not any(x.lower() in lc for x in exclude_any)

        if ok_any and ok_all and ok_exclude:
            hits.append(c)

    return hits


def save_figure(fig, name: str):
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


def short_term_name(s, max_len=44):
    s = str(s)
    s = s.replace("HALLMARK_", "")
    s = s.replace("GOBP_", "")
    s = s.replace("GO_", "")
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


def to_bool_series(s):
    if s.dtype == bool:
        return s.fillna(False)

    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).astype(float) > 0

    return s.astype(str).str.lower().isin(["true", "1", "yes", "y", "t"])


def infer_structured_flag(df, sf_col, dreal_col, q95_col, dexcess_q95_col):
    """
    Return boolean structured-favored flag:
    structured_favored = D_real > random q95
    """
    if sf_col is not None:
        return to_bool_series(df[sf_col])

    if q95_col is not None:
        return (
            pd.to_numeric(df[dreal_col], errors="coerce")
            >
            pd.to_numeric(df[q95_col], errors="coerce")
        ).fillna(False)

    if dexcess_q95_col is not None:
        return (
            pd.to_numeric(df[dexcess_q95_col], errors="coerce") > 0
        ).fillna(False)

    raise KeyError("Cannot infer structured-favored status.")


def detect_random_mapped_columns(term: pd.DataFrame):
    """
    Robustly detect columns in Table_SizeBinRandom_T300_term_level_mapped_results.csv.
    Also infers random q95 if needed.
    """

    # Basic columns
    db_col = safe_col(term, ["database"])
    cancer_col = safe_col(term, ["cancer", "cancer_label"], required=False)

    term_col = safe_col(
        term,
        ["bp_name", "term", "bp", "gene_set", "gene_set_name", "term_name", "pathway", "pathway_name"],
        required=False
    )

    # Real D
    dreal_col = safe_col(
        term,
        [
            "D_real",
            "Dreal",
            "D_survival",
            "Dsurvival",
            "real_D",
            "D",
            "D_survival_real",
            "real_D_survival",
        ],
        required=False
    )

    # Random q95
    q95_col = safe_col(
        term,
        [
            "random_D_q95",
            "random_D_q95_by_bin",
            "D_random_q95",
            "D_random_q95_by_bin",
            "Drandom_q95",
            "random_q95",
            "q95_random_D",
            "q95_D_random",
            "D_random_95",
            "random_D_95",
            "D_random_p95",
            "random_D_p95",
            "random_q95_D",
            "random_D_q95_mapped",
            "D_random_q95_mapped",
            "bin_random_D_q95",
            "bin_D_random_q95",
        ],
        required=False
    )

    if q95_col is None:
        q95_candidates = fuzzy_col(
            term,
            include_any=["q95", "p95", "95"],
            include_all=["d"],
            exclude_any=[]
        )
        q95_candidates = [
            c for c in q95_candidates
            if ("random" in c.lower() or "rand" in c.lower())
        ]

        print("Detected q95 candidates:", q95_candidates)

        if len(q95_candidates) > 0:
            q95_col = q95_candidates[0]
            print("Using q95 column:", q95_col)

    # Random median
    dmedian_col = safe_col(
        term,
        [
            "random_D_median",
            "random_D_median_by_bin",
            "D_random_median",
            "D_random_median_by_bin",
            "Drandom_median",
            "random_median",
            "median_random_D",
            "D_random_med",
            "random_D_med",
            "bin_random_D_median",
        ],
        required=False
    )

    if dmedian_col is None:
        median_candidates = fuzzy_col(
            term,
            include_any=["median", "med"],
            include_all=["d"],
            exclude_any=[]
        )
        median_candidates = [
            c for c in median_candidates
            if ("random" in c.lower() or "rand" in c.lower())
        ]

        print("Detected random median candidates:", median_candidates)

        if len(median_candidates) > 0:
            dmedian_col = median_candidates[0]
            print("Using random median column:", dmedian_col)

    # D excess over q95
    dexcess_q95_col = safe_col(
        term,
        [
            "D_excess_vs_random_q95",
            "D_excess_q95",
            "Dexcess_q95",
            "D_excess_over_q95",
            "D_real_minus_random_q95",
            "excess_vs_random_q95",
            "D_excess_vs_q95",
            "D_minus_random_q95",
            "D_real_minus_q95",
        ],
        required=False
    )

    if dexcess_q95_col is None:
        excess_candidates = fuzzy_col(
            term,
            include_any=["excess", "minus"],
            include_all=[],
            exclude_any=[]
        )
        excess_candidates = [
            c for c in excess_candidates
            if ("q95" in c.lower() or "p95" in c.lower() or "95" in c.lower())
        ]

        print("Detected D-excess q95 candidates:", excess_candidates)

        if len(excess_candidates) > 0:
            dexcess_q95_col = excess_candidates[0]
            print("Using D-excess q95 column:", dexcess_q95_col)

    # D excess over random median
    dexcess_med_col = safe_col(
        term,
        [
            "D_excess_vs_random_median",
            "D_excess_median",
            "Dexcess_median",
            "D_excess_over_median",
            "D_real_minus_random_median",
            "excess_vs_random_median",
            "D_excess_vs_median",
            "D_minus_random_median",
            "D_real_minus_median",
        ],
        required=False
    )

    if dexcess_med_col is None:
        excess_med_candidates = fuzzy_col(
            term,
            include_any=["excess", "minus"],
            include_all=[],
            exclude_any=[]
        )
        excess_med_candidates = [
            c for c in excess_med_candidates
            if ("median" in c.lower() or "med" in c.lower())
        ]

        print("Detected D-excess median candidates:", excess_med_candidates)

        if len(excess_med_candidates) > 0:
            dexcess_med_col = excess_med_candidates[0]
            print("Using D-excess median column:", dexcess_med_col)

    # Structured flag
    sf_col = safe_col(
        term,
        [
            "structured_favored",
            "is_structured_favored",
            "above_random_q95",
            "above_q95",
            "real_above_random_q95",
        ],
        required=False
    )

    # High-D flag
    highD_col = safe_col(
        term,
        [
            "high_D_real",
            "highD_real",
            "high_D",
            "is_high_D",
            "nominal_high_D",
            "highD",
        ],
        required=False
    )

    # Validate required columns
    if dreal_col is None:
        raise KeyError(
            "Cannot find real D column in random mapped table.\n"
            f"Available columns:\n{list(term.columns)}"
        )

    if term_col is None:
        raise KeyError(
            "Cannot find term name column in random mapped table.\n"
            "Expected one of: bp_name, term, bp, gene_set, gene_set_name, term_name, pathway_name.\n"
            f"Available columns:\n{list(term.columns)}"
        )

    # Numeric conversion
    term[dreal_col] = pd.to_numeric(term[dreal_col], errors="coerce")

    if q95_col is not None:
        term[q95_col] = pd.to_numeric(term[q95_col], errors="coerce")

    if dmedian_col is not None:
        term[dmedian_col] = pd.to_numeric(term[dmedian_col], errors="coerce")

    if dexcess_q95_col is not None:
        term[dexcess_q95_col] = pd.to_numeric(term[dexcess_q95_col], errors="coerce")

    if dexcess_med_col is not None:
        term[dexcess_med_col] = pd.to_numeric(term[dexcess_med_col], errors="coerce")

    # Infer q95 from D_real - D_excess_q95
    if q95_col is None and dexcess_q95_col is not None:
        term["random_q95_inferred"] = (
            pd.to_numeric(term[dreal_col], errors="coerce")
            -
            pd.to_numeric(term[dexcess_q95_col], errors="coerce")
        )
        q95_col = "random_q95_inferred"
        print("Using inferred random q95 column: D_real - D_excess_vs_random_q95")

    # Infer D_excess_q95 from D_real - q95
    if dexcess_q95_col is None and q95_col is not None:
        term["D_excess_q95_inferred"] = (
            pd.to_numeric(term[dreal_col], errors="coerce")
            -
            pd.to_numeric(term[q95_col], errors="coerce")
        )
        dexcess_q95_col = "D_excess_q95_inferred"
        print("Using inferred D_excess_q95 column: D_real - random_q95")

    # Infer D_excess_median from D_real - random median
    if dexcess_med_col is None and dmedian_col is not None:
        term["D_excess_median_inferred"] = (
            pd.to_numeric(term[dreal_col], errors="coerce")
            -
            pd.to_numeric(term[dmedian_col], errors="coerce")
        )
        dexcess_med_col = "D_excess_median_inferred"
        print("Using inferred D_excess_median column: D_real - random_median")

    if q95_col is None:
        raise KeyError(
            "Could not find or infer random q95 column. "
            "Please inspect the table columns printed above."
        )

    if dexcess_q95_col is None:
        raise KeyError(
            "Could not find or infer D_excess_vs_random_q95 column. "
            "Please inspect the table columns printed above."
        )

    term["structured_flag"] = infer_structured_flag(
        term, sf_col, dreal_col, q95_col, dexcess_q95_col
    )

    if highD_col is not None:
        term["highD"] = to_bool_series(term[highD_col])
    else:
        term["highD"] = term[dreal_col] >= 1.301

    print("\nColumn mapping used:")
    print("database column       :", db_col)
    print("cancer column         :", cancer_col)
    print("term column           :", term_col)
    print("D_real column         :", dreal_col)
    print("random q95 column     :", q95_col)
    print("random median column  :", dmedian_col)
    print("D excess q95 column   :", dexcess_q95_col)
    print("D excess median column:", dexcess_med_col)
    print("structured flag column:", sf_col)
    print("high-D flag column    :", highD_col)

    return {
        "term": term,
        "db_col": db_col,
        "cancer_col": cancer_col,
        "term_col": term_col,
        "dreal_col": dreal_col,
        "q95_col": q95_col,
        "dmedian_col": dmedian_col,
        "dexcess_q95_col": dexcess_q95_col,
        "dexcess_med_col": dexcess_med_col,
    }


# ============================================================
# 3. Load random baseline tables
# ============================================================

term_path = require_file(RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_term_level_mapped_results.csv")
db_path = require_file(RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_database_summary.csv")
cd_path = require_file(RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_cancer_database_summary.csv")
bin_path = require_file(RANDOM_TABLE_DIR / "Table_SizeBinRandom_T300_bin_level_baselines.csv")

df_rand_term = pd.read_csv(term_path)
df_rand_db = pd.read_csv(db_path)
df_rand_cd = pd.read_csv(cd_path)
df_rand_bin = pd.read_csv(bin_path)

for df in [df_rand_term, df_rand_db, df_rand_cd, df_rand_bin]:
    if "database" in df.columns:
        df["database"] = df["database"].map(clean_database_name)

print("Loaded random mapped table:")
print(term_path)
print("Shape:", df_rand_term.shape)
print("Columns:")
for c in df_rand_term.columns:
    print(" -", c)


# ============================================================
# 4. Figure 4
# Size-bin random baseline calibration
# ============================================================

def make_figure_4():
    info = detect_random_mapped_columns(df_rand_term.copy())

    term = info["term"]
    db_col = info["db_col"]
    cancer_col = info["cancer_col"]
    dreal_col = info["dreal_col"]
    q95_col = info["q95_col"]

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))

    fig.suptitle(
        "Figure 4. Size-bin random-baseline calibration of observation-ready terms",
        y=0.98,
        fontsize=17,
        fontweight="bold"
    )

    fig.text(
        0.5, 0.945,
        r"Structured-favored criterion: $D_{\mathrm{real}} > D_{\mathrm{random},q95}$",
        ha="center",
        va="center",
        fontsize=12,
        style="italic",
        color=COL_GRAY
    )

    # --------------------------------------------------------
    # A. Structured-favored fraction by database
    # --------------------------------------------------------
    ax = axes[0, 0]
    add_panel_label(ax, "A")

    frac = term.groupby(db_col)["structured_flag"].mean().reindex(DATABASE_ORDER)

    y = np.arange(len(frac))

    ax.barh(y, frac.values, color=COL_PURPLE, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([DATABASE_LABEL[x] for x in frac.index])
    ax.invert_yaxis()

    xmax = max(0.05, np.nanmax(frac.values) * 1.35)
    ax.set_xlim(0, xmax)

    ax.set_xlabel("Structured-favored fraction")
    ax.set_title("Structured-favored fraction by database", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)

    for v, yy in zip(frac.values, y):
        ax.text(
            v + xmax * 0.02,
            yy,
            f"{v:.3f}",
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold",
            color=COL_PURPLE
        )

    # --------------------------------------------------------
    # B. Dreal versus random q95
    # --------------------------------------------------------
    ax = axes[0, 1]
    add_panel_label(ax, "B")

    plot = term.dropna(subset=[dreal_col, q95_col]).copy()

    if len(plot) > 15000:
        plot = plot.sample(15000, random_state=1)

    colors = np.where(plot["structured_flag"], COL_PURPLE, COL_LIGHT_GRAY)

    ax.scatter(
        plot[q95_col],
        plot[dreal_col],
        s=10,
        alpha=0.35,
        c=colors,
        edgecolors="none"
    )

    lim = max(plot[q95_col].max(), plot[dreal_col].max()) * 1.05

    ax.plot(
        [0, lim],
        [0, lim],
        color=COL_GRAY,
        linestyle="--",
        linewidth=1.3
    )

    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)

    ax.set_xlabel(r"Size-bin random q95, $D_{\mathrm{random},q95}$")
    ax.set_ylabel(r"Real term discriminability, $D_{\mathrm{real}}$")
    ax.set_title(r"$D_{\mathrm{real}}$ versus size-bin random q95", fontweight="bold")
    ax.grid(alpha=0.25)

    ax.legend(
        handles=[
            Patch(facecolor=COL_PURPLE, label="Structured-favored"),
            Patch(facecolor=COL_LIGHT_GRAY, label="Not above random q95")
        ],
        frameon=False,
        loc="lower right"
    )

    # --------------------------------------------------------
    # C. Structured-favored heatmap
    # --------------------------------------------------------
    ax = axes[1, 0]
    add_panel_label(ax, "C")

    if cancer_col is None:
        ax.text(0.5, 0.5, "Cancer column not found", ha="center", va="center")
    else:
        heat = term.pivot_table(
            index=cancer_col,
            columns=db_col,
            values="structured_flag",
            aggfunc="mean"
        )

        heat = heat[[c for c in DATABASE_ORDER if c in heat.columns]]
        heat = heat.sort_index()

        vmax = np.nanmax(heat.values)
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1

        im = ax.imshow(
            heat.values,
            aspect="auto",
            cmap="Purples",
            vmin=0,
            vmax=vmax
        )

        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels(
            [DATABASE_LABEL[x] for x in heat.columns],
            rotation=30,
            ha="right"
        )

        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels(heat.index, fontsize=7)

        ax.set_title("Structured-favored fraction by cancer and database", fontweight="bold")

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("Structured-favored fraction")

    # --------------------------------------------------------
    # D. Overlap of high-D and structured-favored
    # --------------------------------------------------------
    ax = axes[1, 1]
    add_panel_label(ax, "D")

    summary = term.groupby(db_col).agg(
        observation_ready=("highD", "size"),
        highD=("highD", "sum"),
        structured_favored=("structured_flag", "sum"),
    )

    both = term.groupby(db_col).apply(
        lambda x: ((x["highD"]) & (x["structured_flag"])).sum()
    )

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
# 5. Figure 5
# Term-level prioritization after random calibration
# ============================================================

def make_figure_5():
    info = detect_random_mapped_columns(df_rand_term.copy())

    term = info["term"]
    db_col = info["db_col"]
    cancer_col = info["cancer_col"]
    term_col = info["term_col"]
    dreal_col = info["dreal_col"]
    q95_col = info["q95_col"]
    dexcess_q95_col = info["dexcess_q95_col"]

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5))

    fig.suptitle(
        "Figure 5. Term-level prioritization after random-baseline calibration",
        y=0.98,
        fontsize=17,
        fontweight="bold"
    )

    fig.text(
        0.5, 0.945,
        "Observation-ready structured terms ranked by survival discriminability and random-baseline excess",
        ha="center",
        va="center",
        fontsize=12,
        style="italic",
        color=COL_GRAY
    )

    # --------------------------------------------------------
    # A. Top structured-favored terms by D_excess_q95
    # --------------------------------------------------------
    ax = axes[0, 0]
    add_panel_label(ax, "A")

    top = (
        term[term["structured_flag"]]
        .dropna(subset=[dexcess_q95_col])
        .sort_values(dexcess_q95_col, ascending=False)
        .head(15)
    )

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

    # --------------------------------------------------------
    # B. Top real-D terms and random-q95 status
    # --------------------------------------------------------
    ax = axes[0, 1]
    add_panel_label(ax, "B")

    topD = (
        term.dropna(subset=[dreal_col])
        .sort_values(dreal_col, ascending=False)
        .head(15)
    )

    labels = []

    for _, r in topD.iterrows():
        t = short_term_name(r[term_col], 34)
        if cancer_col:
            labels.append(f"{str(r[cancer_col])}: {t}")
        else:
            labels.append(t)

    y = np.arange(len(topD))
    colors = np.where(topD["structured_flag"], COL_PURPLE, COL_LIGHT_GRAY)

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
            Patch(facecolor=COL_LIGHT_GRAY, label="Not above random q95")
        ],
        frameon=False,
        loc="lower right"
    )

    # --------------------------------------------------------
    # C. Database-level fractions
    # --------------------------------------------------------
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
            "Structured-favored": sub["structured_flag"].mean(),
            "Both": ((sub["highD"]) & (sub["structured_flag"])).mean()
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

    # --------------------------------------------------------
    # D. D_real - D_random,q95 distribution
    # --------------------------------------------------------
    ax = axes[1, 1]
    add_panel_label(ax, "D")

    data = [
        term.loc[term[db_col] == db, dexcess_q95_col].dropna().values
        for db in DATABASE_ORDER
    ]

    bp = ax.boxplot(
        data,
        labels=[DATABASE_LABEL[x] for x in DATABASE_ORDER],
        patch_artist=True,
        showfliers=False
    )

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
    save_figure(fig, "Figure_5_term_prioritization_after_random_calibration")
    plt.close(fig)


# ============================================================
# 6. Run Figure 4 and Figure 5
# ============================================================

print("Random table directory:", RANDOM_TABLE_DIR)
print("Output directory:", OUT_DIR)

make_figure_4()
make_figure_5()

print("\nDone. Fixed Figure 4 and Figure 5 saved to:")
print(OUT_DIR)
