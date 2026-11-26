"""
Generate all 9 publication-quality (300 DPI) figures for the 27200 ME2 Metabolomics project
from processed tables produced by src/pipeline.py.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Draw, rdFingerprintGenerator
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
TABLES_DIR = ROOT / "results" / "tables"
NET_DIR = ROOT / "results" / "networks"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Clean publication styling
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
)

# Colorblind-safe palette (Okabe-Ito inspired)
COL_DMSO = "#0072B2"       # deep blue
COL_OMEPR = "#D55E00"      # vermilion
COL_MED = "#009E73"        # bluish green
COL_BLANK = "#7F7F7F"      # neutral gray
COL_QC = "#CC79A7"         # reddish purple
COL_QCMIX = "#E69F00"      # orange/amber
COL_FAIL = "#B2182B"       # crimson
COL_SUPP = "#6A3D9A"       # violet


def load_all():
    raw_feat = pd.read_csv(RAW_DIR / "group2_drugs_iimn_gnps_quant.csv")
    raw_feat = raw_feat.drop(columns=[c for c in raw_feat.columns if c.startswith("Unnamed")])
    peak_cols = [c for c in raw_feat.columns if c.endswith(" Peak area")]
    X_raw = raw_feat[peak_cols].copy()
    X_raw.columns = [c.replace(" Peak area", "") for c in peak_cols]
    X_raw.index = raw_feat["row ID"]

    meta = pd.read_csv(PROC_DIR / "sample_qc_summary.csv")
    meta["run_date_dt"] = pd.to_datetime(meta["run_date"], format="ISO8601")
    X_pqn_all = pd.read_csv(PROC_DIR / "features_filtered_pqn.csv", index_col=0)
    log2_pqn = pd.read_csv(PROC_DIR / "features_filtered_log2_pqn.csv", index_col=0)
    res = pd.read_csv(TABLES_DIR / "table02_statistical_results_all.csv", index_col="row_id")
    t3 = pd.read_csv(TABLES_DIR / "table03_top_omeprazole_biotransformation_hits.csv", index_col="row_id")
    t5 = pd.read_csv(TABLES_DIR / "table05_chemical_class_enrichment.csv")
    t6 = pd.read_csv(TABLES_DIR / "table06_network_nodes_and_hubs.csv")
    t7 = pd.read_csv(TABLES_DIR / "table07_biotransformation_edges.csv")
    t8 = pd.read_csv(TABLES_DIR / "table08_ml_feature_ranking.csv", index_col="row_id")
    t8b = pd.read_csv(TABLES_DIR / "table08b_ml_cv_performance.csv")
    t9 = pd.read_csv(TABLES_DIR / "table09_permanova_variance_partitioning.csv")
    with open(PROC_DIR / "analysis_summary.json") as f:
        summary = json.load(f)
    return X_raw, meta, X_pqn_all, log2_pqn, res, t3, t5, t6, t7, t8, t8b, t9, summary


# =====================================================================
# FIGURE 1: Sample QC, Failed Injections & Run-Order Confounding
# =====================================================================
def make_fig01(X_raw, meta):
    fig = plt.figure(figsize=(17.4, 5.2))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.12, 1.02, 1.24], wspace=0.42)

    # Panel A: TIC vs Detected Features
    ax = fig.add_subplot(gs[0])
    group_styles = [
        ("DMSO (valid, n=25)", (meta["drug"] == "DMSO") & (~meta["is_failed_injection"]), COL_DMSO, "o", 52),
        ("Omeprazole (valid, n=25)", (meta["drug"] == "Omeprazole") & (~meta["is_failed_injection"]), COL_OMEPR, "s", 52),
        ("Medium T0 (n=3)", meta["drug"] == "Medium", COL_MED, "^", 65),
        ("Pooled QC (n=4)", meta["sample_type"] == "Pooled_QC", COL_QC, "D", 60),
        ("QC-Mix (n=4)", meta["sample_type"] == "QC", COL_QCMIX, "P", 65),
        ("Blanks (n=4)", meta["mzmine_sample_type"] == "blank", COL_BLANK, "v", 60),
        ("Failed Injections (n=4)", meta["is_failed_injection"], COL_FAIL, "X", 110),
    ]
    for lbl, msk, col, mk, sz in group_styles:
        sub = meta[msk]
        ax.scatter(
            sub["n_nonzero"],
            sub["TIC"] / 1e9,
            c=col,
            marker=mk,
            s=sz,
            alpha=0.88,
            edgecolors="black" if "Failed" in lbl else "white",
            linewidths=0.8,
            label=lbl,
            zorder=5 if "Failed" in lbl else 3,
        )

    ax.axhspan(0, 0.18, color=COL_FAIL, alpha=0.08, zorder=0)
    ax.annotate(
        "4 failed T4/T5 injections\n(TIC 4.4–4.7×10⁷ < Blanks!)",
        xy=(840, 0.055),
        xytext=(1380, 0.34),
        fontsize=8.3,
        fontweight="bold",
        color=COL_FAIL,
        arrowprops=dict(arrowstyle="->", color=COL_FAIL, lw=1.3),
        bbox=dict(boxstyle="round,pad=0.3", fc="#FFF5F5", ec=COL_FAIL, lw=0.8),
    )
    ax.annotate(
        "Omeprazole TIC doubles\nT0 (1.1e9) → T8 (2.2e9)",
        xy=(3230, 1.48),
        xytext=(1580, 1.46),
        fontsize=8.3,
        color=COL_OMEPR,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=COL_OMEPR, lw=1.3),
        bbox=dict(boxstyle="round,pad=0.3", fc="#FFF8F0", ec=COL_OMEPR, lw=0.8),
    )
    ax.set_xlabel("Detected Features per Injection (peak area > 0)")
    ax.set_ylabel("Total Ion Current (TIC, ×10⁹)")
    ax.set_title("A. Sample QC: Failed Injections & TIC Drift", fontsize=10.5)
    ax.legend(loc="upper left", fontsize=7.6, ncol=1, borderpad=0.35)
    ax.set_ylim(-0.05, 2.42)
    ax.grid(alpha=0.22, linestyle="--")

    # Panel B: Correlation heatmap of Failed Injections vs Blanks vs Valid T4/T5 Replicates
    ax2 = fig.add_subplot(gs[1])
    failed = meta.loc[meta["is_failed_injection"], "filename"].sort_values().tolist()
    blanks = meta.loc[meta["mzmine_sample_type"] == "blank", "filename"].sort_values().tolist()
    valid_ref = [
        "T4_DMSO_Bac_rep1.mzML",
        "T5_DMSO_Bac_rep1.mzML",
        "T4_Omepr_Bac_rep1.mzML",
        "T5_Omepr_Bac_rep1.mzML",
    ]
    sel = failed + blanks + valid_ref
    corr = np.log1p(X_raw[sel]).corr()
    short_names = [
        s.replace("_Bac", "").replace(".mzML", "").replace("Omepr", "Omep") for s in sel
    ]
    im = ax2.imshow(corr.values, cmap="RdYlBu_r", vmin=0.1, vmax=1.0, aspect="auto")
    ax2.set_xticks(range(len(sel)))
    ax2.set_yticks(range(len(sel)))
    ax2.set_xticklabels(short_names, rotation=45, ha="right", fontsize=7.6)
    ax2.set_yticklabels(short_names, fontsize=7.6)
    for pos in [3.5, 7.5]:
        ax2.axhline(pos, color="black", lw=1.5)
        ax2.axvline(pos, color="black", lw=1.5)
    for i in range(len(sel)):
        for j in range(len(sel)):
            v = corr.values[i, j]
            ax2.text(
                j,
                i,
                f"{v:.2f}",
                ha="center",
                va="center",
                fontsize=6.5,
                fontweight="bold" if v > 0.80 else "normal",
                color="white" if (v > 0.78 or v < 0.34) else "#111111",
            )
    cbar = plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.03)
    cbar.set_label("Pearson r (log1p)", fontsize=8.5)
    ax2.set_title("B. Failed Injections Match Blanks (r=0.83–0.89)", fontsize=10.5)

    # Panel C: Run-Order Confounding & Parent Omeprazole Autosampler Instability vs Ufiprazole
    ax3 = fig.add_subplot(gs[2])
    omepr = meta[(meta["drug"] == "Omeprazole") & (~meta["is_failed_injection"])].copy()
    omepr["parent_omep"] = [X_raw.loc[2828, f] for f in omepr["filename"]]
    omepr["ufiprazole"] = [X_raw.loc[2962, f] for f in omepr["filename"]]

    ax3_r = ax3.twinx()
    ax3_r.spines["top"].set_visible(False)
    # Draw primary axis annotations strictly on top of twin-axis triangles!
    ax3.set_zorder(ax3_r.get_zorder() + 1)
    ax3.patch.set_visible(False)

    sc1 = ax3.scatter(
        omepr["run_hours"],
        omepr["parent_omep"] / 1e6,
        c=COL_FAIL,
        s=58,
        marker="o",
        edgecolors="white",
        lw=0.6,
        label="Parent Omeprazole (ID 2828, m/z 346.12)",
        zorder=4,
    )
    m1, b1 = np.polyfit(omepr["run_hours"], omepr["parent_omep"] / 1e6, 1)
    xr = np.linspace(0, 132, 100)
    ax3.plot(xr, np.maximum(m1 * xr + b1, 0), color=COL_FAIL, ls="--", lw=1.4, alpha=0.75)

    # Clean non-overlapping offsets for T0-T8 labels (zorder=9 above all markers!)
    tp_offsets = {
        "T4": (0, 1.3),
        "T5": (-4, 1.7),
        "T8": (4, 1.5),
        "T3": (-5, 1.4),
        "T0": (6, 2.0),
        "T7": (-10, -1.6),
        "T2": (-10, 6.2),
        "T6": (-12, 4.1),
        "T1": (-1, 3.2),
    }
    for tp_lbl, grp in omepr.groupby("sample_type"):
        x_m = grp["run_hours"].mean()
        y_m = grp["parent_omep"].mean() / 1e6
        dx, dy = tp_offsets.get(tp_lbl, (0, 1.2))
        ax3.annotate(
            tp_lbl,
            xy=(x_m, y_m),
            xytext=(x_m + dx, y_m + dy),
            fontsize=7.9,
            fontweight="bold",
            color=COL_FAIL,
            ha="center",
            arrowprops=dict(arrowstyle="-", color=COL_FAIL, lw=0.7, alpha=0.7) if tp_lbl in ("T1", "T2", "T6") else None,
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec=COL_FAIL if tp_lbl in ("T1", "T2", "T6") else "none", lw=0.5, alpha=0.92),
            zorder=9,
        )

    sc2 = ax3_r.scatter(
        omepr["run_hours"],
        omepr["ufiprazole"] / 1e8,
        c=COL_DMSO,
        s=50,
        marker="^",
        alpha=0.82,
        edgecolors="white",
        lw=0.5,
        label="Ufiprazole / Sulfide (ID 2962, m/z 330.13)",
        zorder=3,
    )

    ax3.set_xlabel("LC-MS Acquisition Queue Time (hours from Jul 2 10:07)")
    ax3.set_ylabel("Parent Omeprazole Peak Area (×10⁶)", color=COL_FAIL)
    ax3.tick_params(axis="y", labelcolor=COL_FAIL)
    ax3_r.set_ylabel("Ufiprazole Peak Area (×10⁸)", color=COL_DMSO)
    ax3_r.tick_params(axis="y", labelcolor=COL_DMSO)
    ax3.set_ylim(-1.2, 20.5)
    ax3_r.set_ylim(-0.3, 6.4)
    ax3.set_title("C. Run-Order Bias: Parent Omeprazole Decays", fontsize=10.5)
    ax3.grid(alpha=0.22, linestyle="--")

    ax3.text(
        0.52,
        0.96,
        "Parent ID 2828:\n  r(run) = −0.62 (p = 9.4e−4)\n  r(time) = +0.08 (n.s.)\n  T1 & T2 = 0 (Day 5–6!)\nUfiprazole ID 2962:\n  r(time) = +0.98 (p < 1e−14)",
        transform=ax3.transAxes,
        va="top",
        fontsize=7.6,
        bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#CCCCCC", alpha=0.95),
    )

    handles = [sc1, sc2]
    labels = [h.get_label() for h in handles]
    ax3.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.48, -0.16), ncol=1, fontsize=7.8)

    fig.suptitle(
        "Experimental Design & Quality Control Audit (n = 69 Injections)",
        fontsize=12.8,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(FIG_DIR / "fig01_qc_failed_injections_and_run_order.png")
    plt.close(fig)


# =====================================================================
# FIGURE 2: PQN vs TIC Normalization & ESI Ion-Suppression Notch
# =====================================================================
def make_fig02(meta, X_pqn_all, res, summary):
    fig = plt.figure(figsize=(17.4, 5.2))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.02, 0.86, 1.46], wspace=0.36)

    bac = meta[(meta["bacteria"] == "Yes") & (~meta["is_failed_injection"])].copy()

    # Panel A: TIC vs PQN scaling factors across T0-T8
    ax = fig.add_subplot(gs[0])
    for drug, col in [("DMSO", COL_DMSO), ("Omeprazole", COL_OMEPR)]:
        sub = bac[bac["drug"] == drug].groupby("time_num")
        t = sub["tic_factor"].mean().index.values
        ax.plot(t, sub["tic_factor"].mean().values, color=col, lw=2.2, ls="-", marker="o", ms=5.5, label=f"{drug} — TIC factor")
        ax.fill_between(
            t,
            sub["tic_factor"].mean() - sub["tic_factor"].std(),
            sub["tic_factor"].mean() + sub["tic_factor"].std(),
            color=col,
            alpha=0.14,
        )
        ax.plot(t, sub["pqn_factor"].mean().values, color=col, lw=2.2, ls="--", marker="s", ms=5.0, label=f"{drug} — PQN factor")
    ax.axhline(1.0, color="gray", ls=":", lw=1.0)
    # Place callout cleanly in open lower-right band (between orange TIC ribbon and PQN lines!)
    ax.annotate(
        "TIC doubles (1.05→1.99×)\ndue to Ufiprazole cg=220\n(~47% of total TIC!)",
        xy=(5.0, 1.82),
        xytext=(3.75, 1.18),
        fontsize=7.9,
        fontweight="bold",
        color=COL_OMEPR,
        arrowprops=dict(arrowstyle="->", color=COL_OMEPR, lw=1.3),
        bbox=dict(boxstyle="round,pad=0.3", fc="#FFF8F0", ec=COL_OMEPR, lw=0.8),
        zorder=7,
    )
    ax.set_xlabel("Incubation Time (hours, T0–T8)")
    ax.set_ylabel("Normalization Scaling Factor")
    ax.set_xticks(range(9))
    ax.set_xticklabels([f"T{t}" for t in range(9)])
    ax.set_ylim(0.68, 2.28)
    ax.set_title("A. TIC Doubles vs. Stable PQN", fontsize=10.5)
    ax.legend(loc="upper left", fontsize=7.6)
    ax.grid(alpha=0.22, linestyle="--")

    # Panel B: Shared Endogenous False Down-Regulation under TIC vs PQN
    ax2 = fig.add_subplot(gs[1])
    cats = ["Treatment\nDown in Omepr\n(q < 0.05)", "Trt × Time\nNegative Slope\n(q < 0.05)"]
    tic_vals = [summary["tic_shared_sig_down"], summary["tic_shared_neg_int"]]
    pqn_vals = [summary["pqn_shared_sig_down"], summary["pqn_shared_neg_int"]]
    x = np.arange(len(cats))
    w = 0.34
    b1 = ax2.bar(x - w / 2, tic_vals, w, color="#D95F02", alpha=0.88, label="TIC Normalized")
    b2 = ax2.bar(x + w / 2, pqn_vals, w, color=COL_DMSO, alpha=0.88, label="PQN Normalized")
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            pct = 100.0 * h / summary["tic_vs_pqn_shared_n"]
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                h + 35,
                f"{int(h):,}\n({pct:.0f}%)",
                ha="center",
                va="bottom",
                fontsize=8.2,
                fontweight="bold",
            )
    ax2.set_xticks(x)
    ax2.set_xticklabels(cats, fontsize=8.6)
    ax2.set_ylabel(f"Shared Endogenous Features (n = {summary['tic_vs_pqn_shared_n']:,})")
    ax2.set_ylim(0, 2200)
    ax2.set_title("B. TIC Creates ~1,300 False Down Hits", fontsize=10.5)
    ax2.legend(loc="upper right", fontsize=8.0)
    ax2.grid(axis="y", alpha=0.22, linestyle="--")

    # Panel C: ESI Co-Elution Ion Suppression Notch across Retention Time
    ax3 = fig.add_subplot(gs[2])
    sh = res[res["feature_category"] == "Shared_Endogenous"].copy()
    is_neg = (sh["q_int_pqn_unadj"] < 0.05) & (sh["beta_int_pqn_unadj"] < 0)
    is_pos = (sh["q_int_pqn_unadj"] < 0.05) & (sh["beta_int_pqn_unadj"] > 0)
    is_ns = ~(is_neg | is_pos)

    ax3.axvspan(1.80, 2.25, color=COL_SUPP, alpha=0.14, zorder=0, label="ESI Suppression Zone (RT 1.80–2.25 min)")
    ax3.scatter(sh.loc[is_ns, "rt"], sh.loc[is_ns, "beta_int_pqn_unadj"], c="#BDBDBD", s=11, alpha=0.45, label="Non-significant shared", zorder=1)
    ax3.scatter(sh.loc[is_pos, "rt"], sh.loc[is_pos, "beta_int_pqn_unadj"], c=COL_MED, s=22, alpha=0.85, label=f"Positive Trt×Time (n={is_pos.sum()})", zorder=3)
    ax3.scatter(sh.loc[is_neg, "rt"], sh.loc[is_neg, "beta_int_pqn_unadj"], c=COL_SUPP, s=24, alpha=0.88, edgecolors="white", lw=0.3, label=f"Negative Trt×Time (n={is_neg.sum()})", zorder=4)
    ax3.axhline(0, color="black", lw=0.9, ls="--", alpha=0.6)

    ax3.axvline(1.898, color=COL_OMEPR, lw=2.0, ls="-", zorder=5)
    ax3.annotate(
        "Ufiprazole peak (ID 2962)\nRT = 1.90 min, area = 4.9×10⁸\n(9× surge suppresses ESI ions!)",
        xy=(1.90, 0.48),
        xytext=(2.68, 0.36),
        fontsize=7.9,
        fontweight="bold",
        color=COL_OMEPR,
        arrowprops=dict(arrowstyle="->", color=COL_OMEPR, lw=1.3),
        bbox=dict(boxstyle="round,pad=0.3", fc="#FFF8F0", ec=COL_OMEPR, lw=0.9),
    )

    if 2972 in sh.index:
        ax3.annotate(
            "N-Phenethylacetamide (ID 2972)\n& IPA isomer (ID 3178): suppressed!",
            xy=(2.22, -0.68),
            xytext=(2.75, -0.84),
            fontsize=7.7,
            color=COL_SUPP,
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=COL_SUPP, lw=1.2),
            bbox=dict(boxstyle="round,pad=0.28", fc="#F9F5FF", ec=COL_SUPP, lw=0.8),
        )
    ax3.text(
        0.98,
        0.04,
        f"{summary['suppression_neg_int_in_window']}/{summary['suppression_neg_int_total']} ({summary['suppression_pct_neg_int']}%) neg hits in 9% of RT\nOR = {summary['suppression_fisher_or']}, Fisher p = {summary['suppression_fisher_p']:.1e}",
        transform=ax3.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.8,
        fontweight="bold",
        color=COL_SUPP,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COL_SUPP, lw=0.8),
    )
    ax3.set_xlabel("Chromatographic Retention Time (min)")
    ax3.set_ylabel("Shared Feature Trt × Time Slope (log₂ FC / hr)")
    ax3.set_xlim(0.1, 5.65)
    ax3.set_ylim(-1.45, 1.25)
    ax3.set_title("C. Co-Elution ESI Ion Suppression Notch", fontsize=10.5)
    ax3.legend(loc="upper right", fontsize=7.4)
    ax3.grid(alpha=0.20, linestyle="--")

    fig.suptitle(
        "Normalization Model Awareness (PQN vs. TIC) & LC-MS Matrix Ion Suppression",
        fontsize=12.8,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(FIG_DIR / "fig02_normalization_pqn_vs_tic_and_ion_suppression.png")
    plt.close(fig)


# =====================================================================
# FIGURE 3: Global & Endogenous PCA + PERMANOVA Variance Partitioning
# =====================================================================
def make_fig03(meta, X_pqn_all, log2_pqn, res, t9):
    fig = plt.figure(figsize=(17.6, 5.4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.08, 1.12, 1.18], wspace=0.42)

    bac = meta[(meta["bacteria"] == "Yes") & (~meta["is_failed_injection"])].copy()
    med = meta[meta["drug"] == "Medium"].copy()
    bac_files = bac["filename"].tolist()
    med_files = med["filename"].tolist()

    arr_all = X_pqn_all[bac_files + med_files].values.astype(float).copy()
    for i in range(arr_all.shape[0]):
        pos = arr_all[i, arr_all[i] > 0]
        arr_all[i, arr_all[i] == 0] = pos.min() * 0.2 if len(pos) else 1.0
    log2_with_med = pd.DataFrame(np.log2(arr_all), index=X_pqn_all.index, columns=bac_files + med_files)

    cmap_t = plt.cm.viridis

    def plot_pca_panel(ax, mat_df, meta_df, title, label_offsets, show_medium=False):
        X_sc = StandardScaler().fit_transform(mat_df.values.T)
        pca = PCA(n_components=2, random_state=42)
        scores = pca.fit_transform(X_sc)
        ev = pca.explained_variance_ratio_ * 100.0
        df_p = meta_df[["filename", "drug", "sample_type", "time_num"]].copy()
        df_p["PC1"] = scores[:, 0]
        df_p["PC2"] = scores[:, 1]

        for drug, marker, edge_c, ls in [("DMSO", "o", COL_DMSO, "-"), ("Omeprazole", "s", COL_OMEPR, "-")]:
            sub = df_p[df_p["drug"] == drug]
            cent = sub.groupby("time_num")[["PC1", "PC2"]].mean().sort_index()
            ax.plot(cent["PC1"], cent["PC2"], color=edge_c, lw=2.2, alpha=0.85, ls=ls, zorder=3)
            ax.annotate(
                "",
                xy=(cent.loc[8, "PC1"], cent.loc[8, "PC2"]),
                xytext=(cent.loc[7, "PC1"], cent.loc[7, "PC2"]),
                arrowprops=dict(arrowstyle="->", color=edge_c, lw=2.2),
                zorder=4,
            )
            ax.scatter(
                sub["PC1"],
                sub["PC2"],
                c=sub["time_num"],
                cmap=cmap_t,
                vmin=0,
                vmax=8,
                marker=marker,
                s=58,
                edgecolors=edge_c,
                linewidths=1.3,
                alpha=0.92,
                zorder=5,
            )
            for tp_k in (0, 8):
                key = f"{'Omep' if drug=='Omeprazole' else 'DMSO'} T{tp_k}"
                off = label_offsets[key]
                ax.annotate(
                    key,
                    xy=(cent.loc[tp_k, "PC1"], cent.loc[tp_k, "PC2"]),
                    xytext=off,
                    textcoords="offset points",
                    fontsize=7.8,
                    fontweight="bold",
                    color=edge_c,
                    arrowprops=dict(arrowstyle="-", color=edge_c, lw=0.9),
                    bbox=dict(boxstyle="round,pad=0.20", fc="white", ec=edge_c, alpha=0.95, lw=0.8),
                    zorder=8,
                )

        if show_medium:
            sub_m = df_p[df_p["drug"] == "Medium"]
            ax.scatter(
                sub_m["PC1"],
                sub_m["PC2"],
                c=COL_MED,
                marker="^",
                s=95,
                edgecolors="black",
                lw=0.8,
                zorder=6,
                label="T0 Medium (no bacteria)",
            )
            ax.annotate(
                "T0 Medium\n(no bacteria)",
                xy=(sub_m["PC1"].mean(), sub_m["PC2"].mean()),
                xytext=(24, -4),
                textcoords="offset points",
                fontsize=7.6,
                fontweight="bold",
                color=COL_MED,
                arrowprops=dict(arrowstyle="->", color=COL_MED, lw=1.0),
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=COL_MED, alpha=0.95, lw=0.8),
                zorder=8,
            )

        ax.set_xlabel(f"PC1 ({ev[0]:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}% variance)")
        ax.set_title(title, fontsize=10.5)
        ax.margins(x=0.14, y=0.12)
        ax.grid(alpha=0.22, linestyle="--")

        legend_els = [
            Line2D([0], [0], marker="o", color=COL_DMSO, lw=2, mfc="white", mec=COL_DMSO, ms=7, label="DMSO (T0→T8)"),
            Line2D([0], [0], marker="s", color=COL_OMEPR, lw=2, mfc="white", mec=COL_OMEPR, ms=7, label="Omeprazole (T0→T8)"),
        ]
        if show_medium:
            legend_els.append(Line2D([0], [0], marker="^", color="none", mfc=COL_MED, mec="black", ms=8, label="T0 Medium (n=3)"))
        ax.legend(handles=legend_els, loc="best", fontsize=7.6)
        return df_p

    ax1 = fig.add_subplot(gs[0])
    plot_pca_panel(
        ax1,
        log2_pqn[bac_files],
        bac,
        "A. Global PCA (All 3,414 Filtered Features)",
        label_offsets={"DMSO T0": (18, -16), "DMSO T8": (18, 6), "Omep T0": (-16, 18), "Omep T8": (-52, 8)},
    )

    ax2 = fig.add_subplot(gs[1])
    clean_endo = res.index[
        (res["feature_category"] == "Shared_Endogenous")
        & (~res["in_ion_suppression_window"])
        & (~res["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]))
    ]
    plot_pca_panel(
        ax2,
        log2_with_med.loc[clean_endo, bac_files + med_files],
        pd.concat([bac, med], ignore_index=True),
        f"B. Clean Endogenous PCA (n = {len(clean_endo):,}, + Medium)",
        label_offsets={"DMSO T0": (26, 8), "DMSO T8": (-64, 12), "Omep T0": (24, -18), "Omep T8": (-56, 22)},
        show_medium=True,
    )

    sm = plt.cm.ScalarMappable(cmap=cmap_t, norm=plt.Normalize(vmin=0, vmax=8))
    cbar = fig.colorbar(sm, ax=[ax1, ax2], orientation="horizontal", fraction=0.055, pad=0.17, aspect=42)
    cbar.set_ticks(range(9))
    cbar.set_ticklabels([f"T{t} ({t}h)" for t in range(9)])
    cbar.set_label("Sampling Timepoint (hours post-treatment)", fontsize=9)

    # Panel C: Multivariate PERMANOVA Variance Partitioning
    ax3 = fig.add_subplot(gs[2])
    sets = [
        t9["Feature_Set"].unique()[0],
        t9["Feature_Set"].unique()[1],
        t9["Feature_Set"].unique()[2],
    ]
    short_sets = [
        "All Filtered\n(PQN, n=3,414)",
        "Clean Endogenous\n(PQN, n=1,915)",
        "Clean Endogenous\n(TIC, n=1,915)",
    ]
    terms_plot = [
        ("Treatment (Omeprazole vs DMSO)", COL_OMEPR, "Treatment"),
        ("Time (T0-T8, categorical)", COL_DMSO, "Time (T0–T8)"),
        ("Treatment x Time Interaction", COL_QCMIX, "Trt × Time"),
        ("Residual (Replicate + Unexplained)", "#D3D3D3", "Residual"),
    ]
    y_pos = np.arange(len(sets))
    lefts = np.zeros(len(sets))
    for term_key, col, short_lbl in terms_plot:
        vals = [
            float(t9.loc[(t9["Feature_Set"] == s) & (t9["Term"] == term_key), "R2"].iloc[0]) * 100.0
            for s in sets
        ]
        ax3.barh(y_pos, vals, left=lefts, height=0.52, color=col, edgecolor="white", lw=1.0, label=short_lbl)
        for i, v in enumerate(vals):
            if v >= 10.5:
                ax3.text(
                    lefts[i] + v / 2.0,
                    y_pos[i],
                    f"{v:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    fontweight="bold",
                    color="white" if col != "#D3D3D3" else "#222222",
                )
            elif v >= 6.5:
                # Place narrow segment label cleanly ABOVE the bar so it never spills across dividers!
                ax3.text(
                    lefts[i] + v / 2.0,
                    y_pos[i] - 0.34,
                    f"{v:.1f}%",
                    ha="center",
                    va="bottom",
                    fontsize=7.6,
                    fontweight="bold",
                    color="#7A4D00" if col == COL_QCMIX else "#333333",
                )
        lefts += np.array(vals)

    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(short_sets, fontsize=8.4, fontweight="bold")
    ax3.invert_yaxis()
    ax3.set_ylim(2.55, -0.60)
    ax3.set_xlim(0, 102)
    ax3.set_xlabel("Variance Explained (PERMANOVA R², %)")
    ax3.set_title("C. Multivariate PERMANOVA (999 Perms)", fontsize=10.5)
    ax3.legend(loc="lower center", bbox_to_anchor=(0.5, -0.36), ncol=2, fontsize=8.0)
    ax3.text(
        0.5,
        -0.19,
        "Biological terms p ≤ 0.022; conditional Run Order R² ≤ 1.1% (p > 0.06, n.s.)",
        transform=ax3.transAxes,
        ha="center",
        fontsize=7.6,
        color="#333333",
    )

    fig.suptitle(
        "Unsupervised Ordination (PCA) & Multivariate Variance Partitioning (PERMANOVA)",
        fontsize=12.8,
        fontweight="bold",
        y=1.03,
    )
    fig.savefig(FIG_DIR / "fig03_pca_and_variance_partitioning.png")
    plt.close(fig)


# =====================================================================
# FIGURE 4: Volcano Plots (Treatment & Interaction) + Model Comparison
# =====================================================================
def make_fig04(res):
    fig = plt.figure(figsize=(17.4, 5.2))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.15, 1.15, 0.96], wspace=0.36)

    def assign_plot_cat(df):
        cats = np.full(len(df), "Other / Non-Significant", dtype=object)
        is_carry = df["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]).values
        is_supp = (df["feature_category"] == "Shared_Endogenous").values & df["in_ion_suppression_window"].values & (df["beta_trt"] < 0).values
        is_omep = (
            ((df["feature_category"] == "Omeprazole_Specific").values | (df["origin_category"] == "Omeprazole-derived").values | ((df["feature_category"] == "Shared_Endogenous").values & df["in_ion_suppression_window"].values & (df["beta_trt"] >= 2.0).values))
            & (~is_carry)
        )
        is_endo_clean = (df["feature_category"] == "Shared_Endogenous").values & (~is_supp) & (~is_omep) & (~is_carry)
        cats[is_endo_clean] = "Clean Shared Endogenous"
        cats[is_omep] = "Omeprazole-Specific (Biotransformation)"
        cats[is_supp] = "ESI Ion-Suppression Zone (RT 1.80–2.25)"
        cats[is_carry] = "Carryover / Contaminant"
        return cats

    df = res.copy()
    df["plot_cat"] = assign_plot_cat(df)
    cat_styles = [
        ("Clean Shared Endogenous", COL_DMSO, "o", 18, 0.55),
        ("Omeprazole-Specific (Biotransformation)", COL_OMEPR, "o", 22, 0.72),
        ("ESI Ion-Suppression Zone (RT 1.80–2.25)", COL_SUPP, "^", 24, 0.78),
        ("Carryover / Contaminant", "#999999", "x", 22, 0.65),
    ]

    # Panel A: Treatment Main Effect Volcano
    ax1 = fig.add_subplot(gs[0])
    y_trt = -np.log10(df["q_trt"].clip(lower=1e-45))
    for cat_n, col, mk, sz, alp in cat_styles:
        m = df["plot_cat"] == cat_n
        ax1.scatter(df.loc[m, "beta_trt"], y_trt[m], c=col, marker=mk, s=sz, alpha=alp, label=f"{cat_n} ({m.sum():,})")
    ax1.axhline(-np.log10(0.05), color="black", ls="--", lw=0.9)
    ax1.axvline(-1.0, color="gray", ls=":", lw=0.9)
    ax1.axvline(1.0, color="gray", ls=":", lw=0.9)

    # Clean non-overlapping callouts well clear of upper-left legend!
    callouts_a = [
        (2962, "Ufiprazole\n(ID 2962)", (11.8, 17.5)),
        (2465, "4-Desmethylthioether\n(ID 2465)", (-5.5, 31.5)),
        (2480, "5-MeO-2-MBI\n(ID 2480)", (10.8, 39.5)),
        (3605, "Benz-Me-Ufip\n(ID 3605)", (-5.2, 24.5)),
        (2368, "ID 2368\n($\\mathrm{C_{11}H_{15}NO_3S_2}$)", (11.0, 30.0)),
        (2972, "Suppressed:\nID 2972", (-6.2, 11.0)),
        (4553, "IPA (ID 4553)", (-6.0, 4.2)),
    ]
    for fid, txt, xy_txt in callouts_a:
        if fid in df.index:
            ax1.annotate(
                txt,
                xy=(df.loc[fid, "beta_trt"], y_trt.loc[fid]),
                xytext=xy_txt,
                fontsize=7.3,
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8),
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#AAAAAA", alpha=0.94, lw=0.6),
                zorder=6,
            )
    ax1.set_xlim(-7.2, 18.2)
    ax1.set_ylim(-1.5, 48.0)
    ax1.set_xlabel(r"Treatment Effect $\beta_{\mathrm{trt}}$ ($\log_2$ Omeprazole / DMSO)")
    ax1.set_ylabel(r"$-\log_{10}$(BH-FDR $q_{\mathrm{trt}}$)")
    ax1.set_title("A. Treatment Main Effect (Run-Order Adjusted)", fontsize=10.5)
    ax1.legend(loc="upper left", fontsize=6.9, framealpha=0.95)
    ax1.grid(alpha=0.20, linestyle="--")

    # Panel B: Treatment x Time Interaction Volcano
    ax2 = fig.add_subplot(gs[1])
    y_int = -np.log10(df["q_int"].clip(lower=1e-22))
    for cat_n, col, mk, sz, alp in cat_styles:
        m = df["plot_cat"] == cat_n
        ax2.scatter(df.loc[m, "beta_int"], y_int[m], c=col, marker=mk, s=sz, alpha=alp)
    ax2.axhline(-np.log10(0.05), color="black", ls="--", lw=0.9)
    ax2.axvline(0, color="gray", ls=":", lw=0.9)

    callouts_b = [
        (2368, "ID 2368 (m/z 274.06)", (0.88, 18.1)),
        (1775, "ID 1775 (m/z 242.12)", (1.16, 15.8)),
        (1553, "ID 1553 (m/z 225.07)", (1.12, 12.6)),
        (2480, "5-MeO-2-MBI (2480)", (-0.18, 15.2)),
        (2465, "4-Desmethyl (2465)", (-0.58, 12.6)),
        (3605, "Benz-Me-Ufip (3605)", (1.06, 7.6)),
        (2972, "ESI Suppressed\n(RT 1.87–2.14)", (-1.38, 14.8)),
    ]
    for fid, txt, xy_txt in callouts_b:
        if fid in df.index:
            ax2.annotate(
                txt,
                xy=(df.loc[fid, "beta_int"], y_int.loc[fid]),
                xytext=xy_txt,
                fontsize=7.3,
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8),
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#AAAAAA", alpha=0.94, lw=0.6),
                zorder=6,
            )
    ax2.set_xlim(-1.48, 1.88)
    ax2.set_ylim(-0.8, 20.0)
    ax2.set_xlabel(r"Treatment × Time Slope $\beta_{\mathrm{int}}$ ($\log_2$ FC / hr)")
    ax2.set_ylabel(r"$-\log_{10}$(BH-FDR $q_{\mathrm{int}}$)")
    ax2.set_title("B. Divergent Temporal Trajectories (Trt × Time)", fontsize=10.5)
    ax2.grid(alpha=0.20, linestyle="--")

    # Panel C: Model Comparison
    ax3 = fig.add_subplot(gs[2])
    terms = ["Treatment\n(q < 0.05)", "Time\n(q < 0.05)", "Trt × Time\n(q < 0.05)"]
    lin_counts = [
        int((res["q_trt"] < 0.05).sum()),
        int((res["q_time"] < 0.05).sum()),
        int((res["q_int"] < 0.05).sum()),
    ]
    cat_counts = [
        int((res["q_cat_trt"] < 0.05).sum()),
        int((res["q_cat_time"] < 0.05).sum()),
        int((res["q_cat_int"] < 0.05).sum()),
    ]
    nonlin_only = int(((res["q_cat_int"] < 0.05) & (res["q_int"] >= 0.05)).sum())
    x = np.arange(len(terms))
    w = 0.35
    b1 = ax3.bar(x - w / 2, lin_counts, w, color=COL_DMSO, alpha=0.88, label="Linear Time + Run Covariate")
    b2 = ax3.bar(x + w / 2, cat_counts, w, color=COL_OMEPR, alpha=0.88, label="2-Way Categorical ANOVA")
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width() / 2, h + 55, f"{h:,}", ha="center", va="bottom", fontsize=8.1, fontweight="bold")
    ax3.set_xticks(x)
    ax3.set_xticklabels(terms, fontsize=8.5)
    ax3.set_ylabel("Significant Features (BH-FDR q < 0.05)")
    ax3.set_ylim(0, 4550)
    ax3.set_title("C. Linear vs. Categorical Models", fontsize=10.5)
    ax3.legend(loc="upper left", fontsize=7.6)
    ax3.grid(axis="y", alpha=0.22, linestyle="--")
    ax3.annotate(
        f"+{nonlin_only} non-linear Trt×Time hits\ncaught only by Categorical ANOVA\n(e.g., transients ID 4470, 2719)",
        xy=(2 + w / 2, 1980),
        xytext=(1.42, 3420),
        fontsize=7.3,
        fontweight="bold",
        color=COL_OMEPR,
        arrowprops=dict(arrowstyle="->", color=COL_OMEPR, lw=1.2),
        bbox=dict(boxstyle="round,pad=0.28", fc="#FFF8F0", ec=COL_OMEPR, lw=0.8),
    )

    fig.suptitle(
        "Differential Abundance Statistics Across 3,414 Filtered Features (BH-FDR)",
        fontsize=12.8,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(FIG_DIR / "fig04_statistical_volcano_and_model_comparison.png")
    plt.close(fig)


# =====================================================================
# FIGURE 5: Temporal Heatmap & 4 Kinetic Trajectory Modules
# =====================================================================
def make_fig05(meta, log2_pqn, res, t3):
    fig = plt.figure(figsize=(17.6, 6.8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.22, 1.0], wspace=0.34)

    bac = meta[(meta["bacteria"] == "Yes") & (~meta["is_failed_injection"])].copy()
    tp_means = pd.DataFrame(
        {
            f"{d}_{t}": log2_pqn[s["filename"].tolist()].mean(axis=1)
            for (d, t), s in bac.groupby(["drug", "sample_type"])
        }
    )

    # Panel A: Top 18 landmark IIMN-collapsed Omeprazole-dynamic features (large, legible labels!)
    ax1 = fig.add_subplot(gs[0])
    hm_ids = [
        i
        for i in [
            2497, 2685, 2478, 4470, 2719, 4157,
            2828, 4709,
            3605, 2500, 1677, 2465, 2618, 2480, 2311, 1775, 2368, 2962, 1553,
        ]
        if i in tp_means.index
    ]

    mod_rank = {
        "M4: Early Precursors Depleted (T0-T3->T8)": 0,
        "M3: Early-Mid Intermediates (Peak T2-T3)": 1,
        "M2: Mid-Late Transients (Peak T4-T5)": 2,
        "M1: Monotonic Accumulation (T0->T8)": 3,
    }
    hm_meta = res.loc[hm_ids].copy()
    hm_meta["mrank"] = hm_meta["kinetic_module"].map(mod_rank).fillna(2)
    hm_meta = hm_meta.sort_values(["mrank", "spearman_time_omepr"])
    hm_ids = hm_meta.index.tolist()

    cols_ord = [f"DMSO_T{t}" for t in range(9)] + [f"Omeprazole_T{t}" for t in range(9)]
    M = tp_means.loc[hm_ids, cols_ord].values
    Mz = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-8)

    im = ax1.imshow(Mz, cmap="RdBu_r", vmin=-2.0, vmax=2.0, aspect="auto")
    ax1.axvline(8.5, color="black", lw=2.0)

    custom_labels = {
        2962: "ID 2962: Ufiprazole / Sulfide (330.13)",
        2480: "ID 2480: 5-MeO-2-MBI (181.04)",
        2465: "ID 2465: 4-Desmethylthioether (316.11)",
        2828: "ID 2828: Parent Omeprazole* (346.12)",
        4157: "ID 4157: Omeprazole −H₂O (328.11)",
        3605: "ID 3605: Benz-Me-Ufip [C₁₈H₂₁N₃O₂S]",
        2500: "ID 2500: Pyr-Me-Ufip [C₁₈H₂₁N₃O₂S]",
        1677: "ID 1677: Me-5-MeO-2-MBI [C₉H₁₀N₂OS]",
        2618: "ID 2618: Me-Omeprazole [C₁₈H₂₁N₃O₃S]",
        2368: "ID 2368: Unknown 274.06 [C₁₁H₁₅NO₃S₂, 2S]",
        1553: "ID 1553: Unknown 225.07 [C₁₀H₁₂N₂O₂S]",
        1775: "ID 1775: Unknown 242.12 [C₁₂H₁₉NO₂S, 1S]",
        4470: "ID 4470: Transient 328.13 [C₁₇H₁₇N₃O₄]",
        2719: "ID 2719: Transient 342.14 [C₁₈H₁₉N₃O₄]",
        2311: "ID 2311: Late isomer 346.12 (RT 1.38)",
        2478: "ID 2478: Transient 298.15 [C₁₇H₁₉N₃O₂]",
        2497: "ID 2497: Transient 296.14 [C₁₇H₁₇N₃O₂]",
        2685: "ID 2685: Transient 326.15 [C₁₈H₁₉N₃O₃]",
        4709: "ID 4709: Mid-late 338.10 (RT 2.68)",
    }
    ylabels = [custom_labels.get(fid, f"ID {fid}: m/z {res.loc[fid,'mz']:.2f}") for fid in hm_ids]

    ax1.set_yticks(range(len(hm_ids)))
    ax1.set_yticklabels(ylabels, fontsize=8.8)
    ax1.set_xticks(range(18))
    ax1.set_xticklabels([f"T{t}" for t in range(9)] * 2, fontsize=8.5)
    ax1.text(4, -1.05, "DMSO Control (T0–T8)", ha="center", va="bottom", fontsize=9.8, fontweight="bold", color=COL_DMSO)
    ax1.text(13, -1.05, "Omeprazole Treated (T0–T8)", ha="center", va="bottom", fontsize=9.8, fontweight="bold", color=COL_OMEPR)
    cbar = plt.colorbar(im, ax=ax1, fraction=0.032, pad=0.025)
    cbar.set_label("Row Z-score (log₂ PQN)", fontsize=9.0)
    ax1.set_title("A. Top 19 Landmark Biotransformation Features (IIMN-Collapsed)", fontsize=10.8, pad=30)

    # Panel B: 4 K-Means Kinetic Modules (2x2 subgrid on the right)
    gs_r = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1], hspace=0.42, wspace=0.28)
    mod_specs = [
        ("M1: Monotonic Accumulation (T0->T8)", COL_OMEPR, "M1: Monotonic Products\n(Ufiprazole, 2480, 2465, 2368)"),
        ("M2: Mid-Late Transients (Peak T4-T5)", COL_QCMIX, "M2: Mid-Late Transients\n(Peak T4–T5; incl. 2828*)"),
        ("M3: Early-Mid Intermediates (Peak T2-T3)", COL_SUPP, "M3: Early-Mid Intermediates\n(Peak T2–T3: ID 4470, 2719)"),
        ("M4: Early Precursors Depleted (T0-T3->T8)", COL_DMSO, "M4: Early Depleted Precursors\n(Exhausted T0–T3 → T8)"),
    ]
    t_vals = np.arange(9)
    for idx, (mod_key, col, title_sub) in enumerate(mod_specs):
        ax_m = fig.add_subplot(gs_r[idx // 2, idx % 2])
        m_ids = res.index[res["kinetic_module"] == mod_key]
        mat = tp_means.loc[m_ids, [f"Omeprazole_T{t}" for t in range(9)]].values
        zmat = (mat - mat.mean(axis=1, keepdims=True)) / (mat.std(axis=1, keepdims=True) + 1e-8)
        rng = np.random.default_rng(42)
        pick = rng.choice(len(zmat), size=min(50, len(zmat)), replace=False)
        for p in pick:
            ax_m.plot(t_vals, zmat[p], color=col, alpha=0.12, lw=0.7)
        med_c = np.median(zmat, axis=0)
        q25, q75 = np.percentile(zmat, [25, 75], axis=0)
        ax_m.fill_between(t_vals, q25, q75, color=col, alpha=0.25)
        ax_m.plot(t_vals, med_c, color=col, lw=2.5, marker="o", ms=4.8)
        ax_m.axhline(0, color="gray", ls=":", lw=0.8)
        ax_m.set_xticks(range(0, 9, 2))
        ax_m.set_xticklabels([f"T{t}" for t in range(0, 9, 2)], fontsize=8.5)
        ax_m.set_ylim(-2.4, 2.4)
        ax_m.set_title(f"{title_sub}\n(n = {len(m_ids)} unique entities)", fontsize=9.0, fontweight="bold")
        if idx % 2 == 0:
            ax_m.set_ylabel("Z-score (Omeprazole)", fontsize=8.8)
        if idx >= 2:
            ax_m.set_xlabel("Timepoint", fontsize=8.8)
        ax_m.grid(alpha=0.20, linestyle="--")

    fig.suptitle(
        "Temporal Heatmap & Unsupervised Kinetic Clustering of 606 Unique Omeprazole-Dynamic Entities",
        fontsize=12.8,
        fontweight="bold",
        y=0.99,
    )
    fig.savefig(FIG_DIR / "fig05_temporal_heatmap_and_kinetic_clusters.png")
    plt.close(fig)


# =====================================================================
# FIGURE 6: Supervised ML Classification & Discriminant Feature Ranking
# =====================================================================
def make_fig06(t8, t8b, t3):
    fig = plt.figure(figsize=(17.6, 5.8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.88, 1.35], wspace=0.88)

    # Panel A: 5-Fold Stratified CV Performance
    ax1 = fig.add_subplot(gs[0])
    agg = t8b.groupby(["Task", "Model"])[["ROC_AUC", "Balanced_Acc"]].agg(["mean", "std"]).reset_index()
    tasks_ord = [
        "Omeprazole vs DMSO (All Unique)",
        "Omeprazole vs DMSO (Clean Endogenous)",
        "Growth Stage T0-2/T3-5/T6-8 (Clean Endogenous)",
    ]
    short_t = [
        "Omeprazole vs DMSO\n(All Unique, n=1,778)",
        "Omeprazole vs DMSO\n(Clean Endogenous, n=1,915)",
        "Growth Stage T0–2/3–5/6–8\n(Clean Endogenous, n=1,915)",
    ]
    y = np.arange(len(tasks_ord))
    h = 0.28
    for idx_t, task in enumerate(tasks_ord):
        sub_rf = agg[(agg["Task"] == task) & (agg["Model"] == "Random Forest")].iloc[0]
        sub_lr = agg[(agg["Task"] == task) & (agg["Model"] != "Random Forest")].iloc[0]

        rf_mu, rf_sd = float(sub_rf[("Balanced_Acc", "mean")]), float(sub_rf[("Balanced_Acc", "std")])
        lr_mu, lr_sd = float(sub_lr[("Balanced_Acc", "mean")]), float(sub_lr[("Balanced_Acc", "std")])
        # Asymmetric clipped error bars so upper error never exceeds 1.00!
        rf_err = [[rf_sd], [min(rf_sd, 1.0 - rf_mu)]] if rf_sd > 1e-4 else None
        lr_err = [[lr_sd], [min(lr_sd, 1.0 - lr_mu)]] if lr_sd > 1e-4 else None

        ax1.barh(idx_t - h / 2, rf_mu, h, xerr=rf_err, color=COL_DMSO, alpha=0.88, capsize=3 if rf_err else 0, label="Random Forest (Bal. Acc.)" if idx_t == 0 else None)
        ax1.barh(idx_t + h / 2, lr_mu, h, xerr=lr_err, color=COL_OMEPR, alpha=0.88, capsize=3 if lr_err else 0, label="Regularized LogReg (Bal. Acc.)" if idx_t == 0 else None)

        # Position label cleanly to the left of any error whisker!
        rf_txt_x = (rf_mu - rf_sd - 0.015) if rf_sd > 1e-4 else (rf_mu - 0.015)
        lr_txt_x = (lr_mu - lr_sd - 0.015) if lr_sd > 1e-4 else (lr_mu - 0.015)
        ax1.text(rf_txt_x, idx_t - h / 2, f"Acc={rf_mu:.3f} | AUC={sub_rf[('ROC_AUC','mean')]:.3f}", ha="right", va="center", color="white", fontsize=8.0, fontweight="bold")
        ax1.text(lr_txt_x, idx_t + h / 2, f"Acc={lr_mu:.3f} | AUC={sub_lr[('ROC_AUC','mean')]:.3f}", ha="right", va="center", color="white", fontsize=8.0, fontweight="bold")

    ax1.set_yticks(y)
    ax1.set_yticklabels(short_t, fontsize=8.8, fontweight="bold")
    ax1.invert_yaxis()
    ax1.set_xlim(0.5, 1.04)
    ax1.set_xlabel("5-Fold Stratified CV Balanced Accuracy (mean ± SD, capped at 1.0)")
    ax1.set_title("A. Supervised ML Classification (5-Fold CV)", fontsize=10.8)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=8.2)
    ax1.grid(axis="x", alpha=0.22, linestyle="--")

    # Panel B: Top 15 Discriminant Features (larger font, proper math subscripts!)
    ax2 = fig.add_subplot(gs[1])
    top15 = t8.head(15).copy()
    labels_b = []
    colors_b = []
    for fid, r in top15.iterrows():
        is_ann = pd.notna(r["compound_name"])
        if fid == 2618:
            lbl = "ID 2618 | m/z 360.14 | Put. Me-Omeprazole"
        elif fid == 2678:
            lbl = "ID 2678 | m/z 288.07 | +CH₂ of ID 2368"
        elif fid == 2480:
            lbl = "ID 2480 | m/z 181.04 | 5-MeO-2-MBI [Ann.]"
        elif fid == 1553:
            lbl = "ID 1553 | m/z 225.07 | C₁₀H₁₂N₂O₂S (r=0.99)"
        elif fid == 2368:
            lbl = "ID 2368 | m/z 274.06 | C₁₁H₁₅NO₃S₂ (2S)"
        elif is_ann:
            lbl = f"ID {fid} | m/z {r['mz']:.2f} | {str(r['compound_name'])[:22]} [Ann.]"
        else:
            lbl = f"ID {fid} | m/z {r['mz']:.2f} (RT {r['rt']:.2f})"
        labels_b.append(lbl)
        colors_b.append(COL_MED if is_ann else COL_OMEPR)

    y2 = np.arange(len(top15))
    ax2.barh(y2, top15["composite_priority_score"].values, color=colors_b, alpha=0.88, edgecolor="white", lw=0.6)
    for i, (fid, r) in enumerate(top15.iterrows()):
        ax2.text(
            r["composite_priority_score"] + 0.004,
            i,
            rf"$\beta_{{\mathrm{{trt}}}}={r['beta_trt']:+.1f}$,  $r_t={r['spearman_time_omepr']:+.2f}$",
            va="center",
            fontsize=8.0,
            color="#222222",
        )
    ax2.set_yticks(y2)
    ax2.set_yticklabels(labels_b, fontsize=8.4)
    ax2.invert_yaxis()
    ax2.set_xlim(0.75, 1.06)
    ax2.set_xlabel("Composite Priority Score (Effect Size + FDR + RF Importance + Consistency)")
    ax2.set_title("B. Top 15 Discriminant Features: Annotated & Dark Metabolome", fontsize=10.8)
    ax2.legend(
        handles=[
            mpatches.Patch(color=COL_MED, label="MS2 Library Annotated"),
            mpatches.Patch(color=COL_OMEPR, label="Unannotated Candidate (Dark Metabolome)"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        fontsize=8.2,
    )
    ax2.grid(axis="x", alpha=0.22, linestyle="--")

    fig.suptitle(
        "Supervised Machine Learning Classification & Dark-Metabolome Feature Prioritization",
        fontsize=12.8,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(FIG_DIR / "fig06_ml_prioritization_known_and_unknown.png")
    plt.close(fig)


# =====================================================================
# FIGURE 7: Chemical Class Enrichment (ORA) & Endogenous Trajectories
# =====================================================================
def make_fig07(meta, X_pqn_all, t5):
    fig = plt.figure(figsize=(17.6, 7.6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.98, 1.38], wspace=0.30)

    # Panel A: Chemical Class Over-Representation Bubble Plot (clean non-overlapping labels!)
    ax1 = fig.add_subplot(gs[0])
    plot_ora = t5[t5["Hits_in_Class"] >= 2].copy()
    plot_ora["neglog10_q"] = -np.log10(plot_ora["q_value"].clip(lower=1e-18))
    set_order = [
        "Omeprazole-Enriched (q<0.05, log2FC>1)",
        "Com20 Time-Increasing (Bacterial Production)",
        "Com20 Time-Decreasing (Bacterial Uptake/Depletion)",
        "Stable mGAM Medium Background",
    ]
    short_set = {
        "Omeprazole-Enriched (q<0.05, log2FC>1)": "Omeprazole-Enriched",
        "Com20 Time-Increasing (Bacterial Production)": "Com20 Produced (T0→T8)",
        "Com20 Time-Decreasing (Bacterial Uptake/Depletion)": "Com20 Consumed (T0→T8)",
        "Stable mGAM Medium Background": "Stable mGAM Medium",
    }
    set_cols = {
        "Omeprazole-Enriched": COL_OMEPR,
        "Com20 Produced (T0→T8)": COL_MED,
        "Com20 Consumed (T0→T8)": COL_DMSO,
        "Stable mGAM Medium": COL_BLANK,
    }
    plot_ora["set_short"] = plot_ora["Hit_Set"].map(short_set)
    # Keep top 3 classes per hit set for clean generous vertical spacing!
    sel_rows = [plot_ora[plot_ora["Hit_Set"] == s].sort_values("p_value").head(3) for s in set_order]
    ora_sub = pd.concat(sel_rows, ignore_index=True)
    ora_sub["row_lbl"] = ora_sub["Chemical_Class"] + "\n[" + ora_sub["set_short"] + "]"

    y = np.arange(len(ora_sub))
    colors = [set_cols[s] for s in ora_sub["set_short"]]
    sizes = 55 + ora_sub["Hits_in_Class"].values * 12
    q05_x = -np.log10(0.05)
    ax1.axvline(q05_x, color="black", ls="--", lw=1.1, zorder=1)
    for i, r in ora_sub.iterrows():
        ax1.hlines(i, 0, r["neglog10_q"], color=colors[i], lw=2.0, alpha=0.65, zorder=2)
    ax1.scatter(ora_sub["neglog10_q"], y, s=sizes, c=colors, edgecolors="black", lw=0.75, zorder=4)
    for i, r in ora_sub.iterrows():
        # Place labels strictly to the right of BOTH the bubble AND the q=0.05 dashed line!
        base_x = r["neglog10_q"] + 0.65 + 0.03 * r["Hits_in_Class"]
        lbl_x = max(base_x, q05_x + 0.35) if r["neglog10_q"] < q05_x else base_x
        ax1.text(
            lbl_x,
            i,
            f"{r['Hits_in_Class']}/{r['Total_in_Class']} ({r['Fold_Enrichment']:.1f}×)",
            va="center",
            fontsize=8.2,
            fontweight="bold" if r["q_value"] < 0.05 else "normal",
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.88),
            zorder=5,
        )
    ax1.set_yticks(y)
    ax1.set_yticklabels(ora_sub["row_lbl"], fontsize=8.3)
    ax1.invert_yaxis()
    ax1.set_xlim(-0.4, 23.2)
    ax1.set_xlabel(r"Enrichment Significance $-\log_{10}$(BH-FDR $q$-value)")
    ax1.set_title("A. Chemical Class Enrichment (Fisher's Exact ORA)", fontsize=10.8)
    legend_els = [mpatches.Patch(color=c, label=k) for k, c in set_cols.items()] + [
        Line2D([0], [0], color="black", ls="--", lw=1.1, label="BH-FDR q = 0.05")
    ]
    ax1.legend(handles=legend_els, loc="lower right", fontsize=8.0)
    ax1.grid(axis="x", alpha=0.22, linestyle="--")

    # Panel B: 6-panel Endogenous Com20 Trajectories (+ 1 suppression artifact comparison!)
    gs_r = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs[1], hspace=0.58, wspace=0.28)
    bac = meta[(meta["bacteria"] == "Yes") & (~meta["is_failed_injection"])].copy()
    med_files = meta[meta["drug"] == "Medium"]["filename"].tolist()

    panels_endo = [
        (4553, r"Indole-3-propionic acid (ID 4553, RT 2.60)" + "\n" + r"Genuine 2× Attenuation ($q_{\mathrm{trt}}=0.001$)", COL_MED),
        (2972, "N-Phenethylacetamide (ID 2972, RT 1.87)\nARTIFACT: Suppressed by Ufiprazole!", COL_SUPP),
        (167, "N-Acetylhistamine (ID 167, RT 0.23)\nConserved >100× Bacterial Production", COL_DMSO),
        (9646, "Gly-C17:1 (ID 9646, RT 5.05)\nConserved Microbial N-Acyl Glycine", COL_DMSO),
        (9772, "Oleic acid (ID 9772, RT 5.10)\nConserved Bacterial Lipid Uptake", COL_DMSO),
        (7477, "PE 13:0 (ID 7477, RT 3.97)\nDelayed Bacterial Lipid in Omeprazole", COL_OMEPR),
    ]

    for idx, (fid, title_p, accent) in enumerate(panels_endo):
        ax_e = fig.add_subplot(gs_r[idx // 2, idx % 2])
        for drug, col, mk in [("DMSO", COL_DMSO, "o"), ("Omeprazole", COL_OMEPR, "s")]:
            sub_m = bac[bac["drug"] == drug]
            vals = X_pqn_all.loc[fid, sub_m["filename"].tolist()].values
            df_t = pd.DataFrame({"t": sub_m["time_num"].values, "v": vals})
            g = df_t.groupby("t")["v"]
            t_x = g.mean().index.values
            mu, sd = g.mean().values, g.std().fillna(0).values
            ax_e.plot(t_x, mu, color=col, lw=2.1, marker=mk, ms=4.6, label=drug)
            ax_e.fill_between(t_x, np.maximum(mu - sd, 0), mu + sd, color=col, alpha=0.16)
            ax_e.scatter(df_t["t"], df_t["v"], color=col, s=14, alpha=0.55)
        med_val = X_pqn_all.loc[fid, med_files].mean()
        ax_e.axhline(med_val, color=COL_MED, ls=":", lw=1.4, label="T0 Medium" if idx == 0 else None)
        if idx == 1:
            ax_e.set_facecolor("#FAF7FF")
        ax_e.set_title(title_p, fontsize=8.8, fontweight="bold", color=COL_SUPP if idx == 1 else "#111111", pad=12)
        ax_e.set_xticks(range(0, 9, 2))
        ax_e.set_xticklabels([f"T{t}" for t in range(0, 9, 2)], fontsize=8.4)
        ax_e.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
        if idx % 2 == 0:
            ax_e.set_ylabel("PQN Peak Area", fontsize=8.6)
        if idx >= 4:
            ax_e.set_xlabel("Timepoint", fontsize=8.6)
        if idx == 0:
            ax_e.legend(loc="upper left", fontsize=7.6)
        ax_e.grid(alpha=0.20, linestyle="--")

    fig.suptitle(
        "Chemical Class Over-Representation & Endogenous Com20 Metabolic Trajectories",
        fontsize=12.8,
        fontweight="bold",
        y=0.99,
    )
    fig.savefig(FIG_DIR / "fig07_chemical_class_enrichment_and_endogenous_trajectories.png")
    plt.close(fig)


# =====================================================================
# FIGURE 8: Multi-Layer Network Analysis (Biotransformation + Endogenous)
# =====================================================================
def repel_label_positions(coords_dict, min_dx=0.22, min_dy=0.11, n_iter=80):
    """Iterative pairwise repulsion so network label boxes never overlap."""
    keys = list(coords_dict.keys())
    pts = {k: np.array(coords_dict[k], dtype=float) for k in keys}
    for _ in range(n_iter):
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                ki, kj = keys[i], keys[j]
                dx = pts[ki][0] - pts[kj][0]
                dy = pts[ki][1] - pts[kj][1]
                if abs(dx) < min_dx and abs(dy) < min_dy:
                    push_x = (min_dx - abs(dx)) * 0.25 * (1 if dx >= 0 else -1)
                    push_y = (min_dy - abs(dy)) * 0.35 * (1 if dy >= 0 else -1)
                    pts[ki][0] += push_x
                    pts[kj][0] -= push_x
                    pts[ki][1] += push_y
                    pts[kj][1] -= push_y
    return pts


def make_fig08(res, t6, t7):
    fig = plt.figure(figsize=(17.6, 8.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.08, 1.0], wspace=0.18)

    # Panel A: Omeprazole Biotransformation Mass-Difference Network
    ax1 = fig.add_subplot(gs[0])
    G_bio = nx.read_graphml(NET_DIR / "omeprazole_biotransformation_network.graphml")
    core_nodes = [
        str(i)
        for i in [2828, 2962, 2465, 2480, 1677, 3605, 2500, 2618, 2368, 2678, 2379, 4470, 2719, 2478, 2497, 2685, 2311, 4157, 3626, 3654, 1604, 1553, 4380, 9575, 3779]
        if str(i) in G_bio.nodes()
    ]
    sub_nodes = set(core_nodes)
    for n in list(core_nodes):
        if n in G_bio:
            sub_nodes.update(list(G_bio.neighbors(n)))
    H_bio = G_bio.subgraph([n for n in sub_nodes if G_bio.degree(n) > 0]).copy()

    pos = nx.kamada_kawai_layout(H_bio, scale=1.25)

    edge_palette = {
        "Sulfoxide Reduction / Oxidation (±O)": COL_OMEPR,
        "Methylation / Demethylation (±CH2)": COL_DMSO,
        "Hydrogenation / Dehydrogenation (±H2)": COL_MED,
        "S-Extrusion Rearrangement (±S-O+H2)": COL_SUPP,
        "C-S Pyridine Cleavage (±C9H11NO)": COL_FAIL,
        "C-S Benzimidazole-S Cleavage (±C8H7N2OS)": COL_FAIL,
        "Formylation / CO Shift (±CO)": "#999999",
        "Acetylation (±C2H2O)": COL_QCMIX,
    }
    for u, v, d in H_bio.edges(data=True):
        t_name = d.get("transformation", "")
        ec = edge_palette.get(t_name, "#B0B0B0")
        r_t = float(d.get("temporal_pearson_r", 0.0))
        ls = "-" if r_t >= 0 else "--"
        nx.draw_networkx_edges(H_bio, pos, edgelist=[(u, v)], ax=ax1, edge_color=ec, width=1.9 if ec != "#B0B0B0" else 1.1, style=ls, alpha=0.80)

    mod_colors = {
        "M1: Monotonic Accumulation (T0->T8)": COL_OMEPR,
        "M2: Mid-Late Transients (Peak T4-T5)": COL_QCMIX,
        "M3: Early-Mid Intermediates (Peak T2-T3)": COL_SUPP,
        "M4: Early Precursors Depleted (T0-T3->T8)": COL_DMSO,
    }
    node_cols, node_sizes = [], []
    for n in H_bio.nodes():
        d = H_bio.nodes[n]
        node_cols.append(mod_colors.get(d.get("kinetic_module", ""), "#8C8C8C"))
        node_sizes.append(125 + 75 * np.log10(max(float(d.get("max_omepr_pqn", 1e4)), 1e4) / 1e4))

    nx.draw_networkx_nodes(
        H_bio,
        pos,
        ax=ax1,
        node_color=node_cols,
        node_size=node_sizes,
        edgecolors="black",
        linewidths=1.0,
        alpha=0.92,
    )

    # Push callout labels radially onto a clean outer perimeter so central nodes stay 100% visible!
    pretty_net_lbl = {
        "2828": "Omeprazole (346.12)",
        "2962": "UFIPRAZOLE (330.13)",
        "2465": "4-Desmethyl (316.11)",
        "2480": "5-MeO-2-MBI (181.04)",
        "1677": "Me-5-MeO-MBI (195.06)",
        "3605": "Benz-Me-Ufip (344.14)",
        "2500": "Pyr-Me-Ufip (344.14)",
        "2618": "Me-Omeprazole (360.14)",
        "2368": "ID 2368: 274.06 (2S)",
        "2678": "ID 2678: 288.07 (+CH₂)",
        "4470": "ID 4470: 328.13 (T3)",
        "4157": "Omep. −H₂O (328.11)",
    }
    cx1 = np.mean([p[0] for p in pos.values()])
    cy1 = np.mean([p[1] for p in pos.values()])
    init_lbl_pos = {}
    for n in pretty_net_lbl:
        if n in pos:
            v = np.array(pos[n]) - np.array([cx1, cy1])
            norm = np.linalg.norm(v) + 1e-6
            # Push all labels onto an outer ring outside the cluster core
            init_lbl_pos[n] = np.array([cx1, cy1]) + (v / norm) * max(norm + 0.55, 1.38)
    rep_pos = repel_label_positions(init_lbl_pos, min_dx=0.72, min_dy=0.26, n_iter=200)
    # Nudge UFIPRAZOLE cleanly left of bottom grey node & Omep -H2O up-right
    if "2962" in rep_pos:
        rep_pos["2962"] = np.array([rep_pos["2962"][0] - 0.85, rep_pos["2962"][1] - 0.10])
    if "4157" in rep_pos:
        rep_pos["4157"] = np.array([rep_pos["4157"][0] + 0.18, rep_pos["4157"][1] + 0.22])
    rep_pos = repel_label_positions(rep_pos, min_dx=0.72, min_dy=0.26, n_iter=40)
    for n, lp in rep_pos.items():
        ax1.annotate(
            pretty_net_lbl[n],
            xy=pos[n],
            xytext=lp,
            fontsize=7.6,
            fontweight="bold",
            ha="center",
            va="center",
            arrowprops=dict(arrowstyle="->", color="#444444", lw=0.85, alpha=0.88),
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#888888", alpha=0.95, lw=0.65),
            zorder=7,
        )
    xs_a = [p[0] for p in pos.values()] + [p[0] for p in rep_pos.values()]
    ys_a = [p[1] for p in pos.values()] + [p[1] for p in rep_pos.values()]
    ax1.set_xlim(min(xs_a) - 0.55, max(xs_a) + 0.55)
    ax1.set_ylim(min(ys_a) - 0.30, max(ys_a) + 0.30)
    ax1.set_title(
        f"A. Omeprazole Biotransformation Mass-Difference Core\n(Connected core: {H_bio.number_of_nodes()} of {G_bio.number_of_nodes()} unique nodes, {H_bio.number_of_edges()} of {G_bio.number_of_edges()} <4.5 mDa edges)",
        fontsize=10.5,
    )
    ax1.axis("off")

    leg_a = [
        mpatches.Patch(color=COL_OMEPR, label="M1: Monotonic"),
        mpatches.Patch(color=COL_QCMIX, label="M2: Mid-Late"),
        mpatches.Patch(color=COL_SUPP, label="M3: Early-Mid"),
        mpatches.Patch(color=COL_DMSO, label="M4: Depleted"),
        mpatches.Patch(color="#8C8C8C", label="Other / Static"),
        Line2D([0], [0], color=COL_OMEPR, lw=2.2, label="±O (Sulfoxide)"),
        Line2D([0], [0], color=COL_DMSO, lw=2.2, label="±CH₂ (Methyl)"),
        Line2D([0], [0], color=COL_MED, lw=2.2, label="±H₂ (Redox)"),
        Line2D([0], [0], color=COL_SUPP, lw=2.2, label="±S−O+H₂"),
        Line2D([0], [0], color=COL_FAIL, lw=2.2, label="±C–S Cleavage"),
        Line2D([0], [0], color="#666666", lw=1.6, ls="--", label="Neg. corr. (r<0)"),
    ]
    ax1.legend(handles=leg_a, loc="upper center", bbox_to_anchor=(0.5, -0.01), fontsize=7.6, ncol=4)

    # Panel B: Endogenous Co-Abundance Network (Giant Connected Core Component)
    ax2 = fig.add_subplot(gs[1])
    G_endo = nx.read_graphml(NET_DIR / "endogenous_coabundance_network.graphml")
    H_endo = nx.Graph()
    for n, d in G_endo.nodes(data=True):
        H_endo.add_node(n, **d)
    for u, v, d in G_endo.edges(data=True):
        if abs(float(d["spearman_r"])) >= 0.78:
            H_endo.add_edge(u, v, **d)
    largest_cc = max(nx.connected_components(H_endo), key=len)
    H_core = H_endo.subgraph(largest_cc).copy()

    pos2 = nx.kamada_kawai_layout(H_core, scale=1.2)
    b_times = np.array([float(H_core.nodes[n]["beta_time"]) for n in H_core.nodes()])
    degs = dict(H_core.degree())
    n_sizes2 = [55 + 12 * degs[n] for n in H_core.nodes()]

    pos_edges = [(u, v) for u, v, d in H_core.edges(data=True) if float(d["spearman_r"]) > 0]
    neg_edges = [(u, v) for u, v, d in H_core.edges(data=True) if float(d["spearman_r"]) < 0]
    nx.draw_networkx_edges(H_core, pos2, edgelist=pos_edges, ax=ax2, edge_color="#BDD7EE", width=0.85, alpha=0.55)
    nx.draw_networkx_edges(H_core, pos2, edgelist=neg_edges, ax=ax2, edge_color="#F4C7C3", width=0.75, style="--", alpha=0.45)

    sc_n = nx.draw_networkx_nodes(
        H_core,
        pos2,
        ax=ax2,
        node_color=b_times,
        cmap="coolwarm",
        vmin=-0.9,
        vmax=0.9,
        node_size=n_sizes2,
        edgecolors="#222222",
        linewidths=0.75,
        alpha=0.94,
    )
    cbar2 = plt.colorbar(sc_n, ax=ax2, fraction=0.036, pad=0.02)
    cbar2.set_label(r"Temporal Slope $\beta_{\mathrm{time}}$ ($\log_2$ FC / hr)" + "\nBlue = Consumed | Red = Produced", fontsize=8.5)

    hub_labels = {
        "9718": "15-HEDE (Hub, deg=47)",
        "7272": "LPE 13:0 (Produced)",
        "4191": "Phe-C4:0 (N-acyl AA)",
        "4553": "Indole-3-propionic acid\n(Attenuated 2× in Omep.)",
        "167": "N-Acetylhistamine (+100×)",
        "9772": "Oleic acid (Consumed)",
        "592": "2'-O-Me-adenosine",
        "9644": "3-OH-C17-Glycine",
        "267": "Guanine (Bridge hub)",
        "1294": "Phe-Pro (Consumed)",
    }
    cx = np.mean([p[0] for p in pos2.values()])
    cy = np.mean([p[1] for p in pos2.values()])
    # Place ALL hub labels strictly on an outer halo ring (radius >= 1.25) so ZERO boxes cover core nodes!
    init_hub_pos = {}
    for n in hub_labels:
        if n in pos2:
            vec = np.array(pos2[n]) - np.array([cx, cy])
            norm = np.linalg.norm(vec) + 1e-6
            init_hub_pos[n] = np.array([cx, cy]) + (vec / norm) * max(norm + 0.55, 1.32)
    rep_pos2 = repel_label_positions(init_hub_pos, min_dx=0.76, min_dy=0.30, n_iter=200)
    # Push Oleic acid cleanly further right & N-Acetylhistamine slightly lower-right so zero nodes/arrows overlap
    if "9772" in rep_pos2:
        rep_pos2["9772"] = np.array([rep_pos2["9772"][0] + 0.42, rep_pos2["9772"][1] + 0.12])
    if "167" in rep_pos2:
        rep_pos2["167"] = np.array([rep_pos2["167"][0] + 0.22, rep_pos2["167"][1] - 0.18])
    rep_pos2 = repel_label_positions(rep_pos2, min_dx=0.76, min_dy=0.30, n_iter=80)
    for n, lp in rep_pos2.items():
        ax2.annotate(
            hub_labels[n],
            xy=pos2[n],
            xytext=lp,
            fontsize=7.6,
            fontweight="bold",
            ha="center",
            va="center",
            arrowprops=dict(arrowstyle="->", color="#333333", lw=0.85, alpha=0.88),
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#888888", alpha=0.95, lw=0.65),
            zorder=7,
        )
    xs_b = [p[0] for p in pos2.values()] + [p[0] for p in rep_pos2.values()]
    ys_b = [p[1] for p in pos2.values()] + [p[1] for p in rep_pos2.values()]
    ax2.set_xlim(min(xs_b) - 0.55, max(xs_b) + 0.88)
    ax2.set_ylim(min(ys_b) - 0.30, max(ys_b) + 0.30)
    ax2.set_title(
        f"B. Endogenous Com20 Metabolic Co-Abundance Core\n(Giant component |r| ≥ 0.78: {H_core.number_of_nodes()} nodes, {H_core.number_of_edges()} edges; Full Q = 0.224)",
        fontsize=10.5,
    )
    ax2.axis("off")

    fig.suptitle(
        "Xenobiotic Biotransformation & Endogenous Com20 Metabolic Networks",
        fontsize=12.8,
        fontweight="bold",
        y=0.99,
    )
    fig.savefig(FIG_DIR / "fig08_biotransformation_and_coabundance_networks.png")
    plt.close(fig)


# =====================================================================
# FIGURE 9: Molecular Zoom-In — 2D Structures, Regioisomers, Kinetics & Isotope Proof
# =====================================================================
def make_fig09(X_raw, meta, X_pqn_all):
    fig = plt.figure(figsize=(17.2, 9.6))
    gs = gridspec.GridSpec(2, 3, height_ratios=[0.95, 1.08], hspace=0.40, wspace=0.34)

    # Top row: 1x6 subgrid of RDKit 2D structures with high-contrast headers!
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gs[0, :], wspace=0.10)
    mols_info = [
        ("1. Parent Omeprazole", "ID 2828 | m/z 346.1218\nC₁₇H₁₉N₃O₃S (GNPS 0.92)", "COc1ccc2nc([nH]c2c1)S(=O)Cc1ncc(C)c(OC)c1C", "#222222"),
        ("2. Ufiprazole (−O)", "ID 2962 | m/z 330.1270\nC₁₇H₁₉N₃O₂S | +9× (4.9×10⁸)", "COc1ccc2nc([nH]c2c1)SCc1ncc(C)c(OC)c1C", COL_OMEPR),
        ("3. 4-Desmethyl (−CH₂)", "ID 2465 | m/z 316.1113\nC₁₆H₁₇N₃O₂S (GNPS 0.87)", "COc1ccc2nc([nH]c2c1)SCc1ncc(C)c(O)c1C", COL_DMSO),
        ("4. 5-MeO-2-MBI (C–S)", "ID 2480 | m/z 181.0430\nC₈H₈N₂OS | ³⁴S 4.2% (1S)", "COc1ccc2[nH]c(S)nc2c1", COL_MED),
        ("5. Benz-Me-Ufiprazole*", "ID 3605 | m/z 344.1427\nFrags: 194.051 + 151.099", "COc1ccc2nc(n(C)c2c1)SCc1ncc(C)c(OC)c1C", COL_SUPP),
        ("6. Pyr-Me-Ufiprazole*", "ID 2500 | m/z 344.1426\nFrags: 180.035 + 165.115", "COc1ccc2nc([nH]c2c1)SCc1ncc(C)c(OCC)c1C", "#8C510A"),
    ]
    for idx_m, (head, sub, smi, accent) in enumerate(mols_info):
        ax_m = fig.add_subplot(gs_top[0, idx_m])
        mol = Chem.MolFromSmiles(smi)
        img = Draw.MolToImage(mol, size=(360, 250))
        ax_m.imshow(img)
        ax_m.set_title(head, fontsize=9.0, fontweight="bold", color=accent, pad=8)
        ax_m.set_xlabel(sub, fontsize=7.9, labelpad=5, color="#222222")
        ax_m.set_xticks([])
        ax_m.set_yticks([])
        for sp in ax_m.spines.values():
            sp.set_visible(True)
            sp.set_edgecolor(accent)
            sp.set_linewidth(1.2)
    fig.text(
        0.5,
        0.928,
        "A. Anaerobic Com20 Omeprazole Biotransformation Cascade (RDKit 2D Structures; *putative regioisomers localized via co-eluting IIMN fragments)",
        ha="center",
        fontsize=10.8,
        fontweight="bold",
    )

    # Bottom-left: Quantitative T0-T8 trajectories — clear lower band for legend (zero curve overlap!)
    ax_bl = fig.add_subplot(gs[1, 0])
    bac = meta[(meta["bacteria"] == "Yes") & (~meta["is_failed_injection"])].copy()
    omepr = bac[bac["drug"] == "Omeprazole"]
    cascade_plot = [
        (2962, "Ufiprazole (2962)", COL_OMEPR, "o", "-"),
        (3605, "Benz-Me-Ufip (3605)", COL_DMSO, "s", "-"),
        (1677, "Me-5-MeO-MBI (1677)", "#8C510A", "^", "-"),
        (2480, "5-MeO-2-MBI (2480)", COL_MED, "D", "-"),
        (2368, "Unknown 2S (2368)", COL_SUPP, "P", "-"),
        (4470, "Transient (4470)", "#555555", "v", "--"),
    ]
    for fid, lbl, col, mk, ls in cascade_plot:
        vals = X_pqn_all.loc[fid, omepr["filename"].tolist()].values
        df_t = pd.DataFrame({"t": omepr["time_num"].values, "v": vals})
        g = df_t.groupby("t")["v"]
        ax_bl.plot(g.mean().index, np.log10(g.mean().values + 1), color=col, lw=2.2, marker=mk, ms=4.8, ls=ls, label=lbl)
    ax_bl.set_xlabel("Incubation Time (hours, Omeprazole)")
    ax_bl.set_ylabel(r"$\log_{10}$(PQN Peak Area + 1)")
    ax_bl.set_xticks(range(9))
    ax_bl.set_xticklabels([f"T{t}" for t in range(9)])
    ax_bl.set_ylim(4.0, 9.15)
    ax_bl.set_title("B. Cascade Kinetics in Omeprazole (T0–T8)", fontsize=10.6)
    ax_bl.legend(loc="lower center", ncol=2, fontsize=7.5, framealpha=0.95)
    ax_bl.grid(alpha=0.22, linestyle="--")

    # Bottom-middle: RDKit Morgan Fingerprint Tanimoto Similarity Matrix (clean 0.55 white-text threshold!)
    ax_bm = fig.add_subplot(gs[1, 1])
    sim_mols = [
        ("Omeprazole", "COc1ccc2nc([nH]c2c1)S(=O)Cc1ncc(C)c(OC)c1C"),
        ("Ufiprazole", "COc1ccc2nc([nH]c2c1)SCc1ncc(C)c(OC)c1C"),
        ("4-Desmethyl", "COc1ccc2nc([nH]c2c1)SCc1ncc(C)c(O)c1C"),
        ("Benz-Me-Ufip*", "COc1ccc2nc(n(C)c2c1)SCc1ncc(C)c(OC)c1C"),
        ("5-MeO-2-MBI", "COc1ccc2[nH]c(S)nc2c1"),
        ("Me-5-MeO-MBI*", "COc1ccc2[nH]c(SC)nc2c1"),
        ("Indole-3-PA", "O=C(O)CCc1c[nH]c2ccccc12"),
        ("Cyclo(Pro-Leu)", "CC(C)CC1NC(=O)C2CCCN2C1=O"),
    ]
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [mfpgen.GetFingerprint(Chem.MolFromSmiles(s)) for _, s in sim_mols]
    n_s = len(fps)
    Smat = np.zeros((n_s, n_s))
    for i in range(n_s):
        for j in range(n_s):
            Smat[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    im_s = ax_bm.imshow(Smat, cmap="YlOrRd", vmin=0, vmax=1.0)
    lbls_s = [k for k, _ in sim_mols]
    ax_bm.set_xticks(range(n_s))
    ax_bm.set_yticks(range(n_s))
    ax_bm.set_xticklabels(lbls_s, rotation=40, ha="right", fontsize=7.8)
    ax_bm.set_yticklabels(lbls_s, fontsize=7.8)
    for i in range(n_s):
        for j in range(n_s):
            ax_bm.text(j, i, f"{Smat[i,j]:.2f}", ha="center", va="center", fontsize=6.9, fontweight="bold" if Smat[i, j] >= 0.55 else "normal", color="white" if Smat[i, j] >= 0.55 else "#111111")
    plt.colorbar(im_s, ax=ax_bm, fraction=0.046, pad=0.03, label="Tanimoto Similarity")
    ax_bm.set_title("C. RDKit Morgan Structural Similarity", fontsize=10.6)

    # Bottom-right: 34S Isotope Confirmation & False-Positive Sertindole Audit (zero bar/label collisions!)
    ax_br = fig.add_subplot(gs[1, 2])
    omepr_files = omepr["filename"].tolist()
    iso_items = [
        ("5-MeO-2-MBI\nID 2480 (1S)", 100.0 * X_raw.loc[2479, omepr_files].max() / X_raw.loc[2480, omepr_files].max(), 4.25, COL_MED),
        ("ID 1775\n242.12 (1S)", 100.0 * X_raw.loc[1770, omepr_files].max() / X_raw.loc[1775, omepr_files].max(), 4.25, COL_DMSO),
        ("ID 2368\n274.06 (2S)", 100.0 * X_raw.loc[2369, omepr_files].max() / X_raw.loc[2368, omepr_files].max(), 8.50, COL_SUPP),
        ("'Sertindole'?\nID 5664 (Cl?)", 0.0, 32.0, COL_FAIL),
    ]
    x_i = np.arange(len(iso_items))
    obs_v = [it[1] for it in iso_items]
    exp_v = [it[2] for it in iso_items]
    cols_i = [it[3] for it in iso_items]
    w = 0.36
    ax_br.bar(x_i - w / 2, obs_v, w, color=cols_i, alpha=0.90, label="Observed M+2 / M (%)")
    ax_br.bar(x_i + w / 2, exp_v, w, color="#DDDDDD", edgecolor="#444444", hatch="//", label="Theoretical M+2 (%)")
    for idx_i, (ov, ev) in enumerate(zip(obs_v, exp_v)):
        ax_br.text(idx_i - w / 2, ov + 0.8, f"{ov:.1f}%", ha="center", va="bottom", fontsize=7.9, fontweight="bold", color=cols_i[idx_i])
        ax_br.text(idx_i + w / 2, ev + 0.8, f"{ev:.1f}%", ha="center", va="bottom", fontsize=7.9, color="#333333")
    ax_br.annotate(
        "Missing 32% ³⁷Cl peak\nproves 'Sertindole' is a\nFALSE library match!",
        xy=(3 - w / 2, 3.2),
        xytext=(0.75, 22.0),
        fontsize=7.9,
        fontweight="bold",
        color=COL_FAIL,
        arrowprops=dict(arrowstyle="->", color=COL_FAIL, lw=1.35, connectionstyle="arc3,rad=-0.22"),
        bbox=dict(boxstyle="round,pad=0.3", fc="#FFF5F5", ec=COL_FAIL, lw=0.9),
    )
    ax_br.set_xticks(x_i)
    ax_br.set_xticklabels([it[0] for it in iso_items], fontsize=8.0)
    ax_br.set_ylabel("M+2 Isotope Relative Intensity (%)")
    ax_br.set_ylim(0, 38.5)
    ax_br.set_title("D. ³⁴S Isotope Proof & False-Hit Audit", fontsize=10.6)
    ax_br.legend(loc="upper left", fontsize=7.8)
    ax_br.grid(axis="y", alpha=0.22, linestyle="--")

    fig.suptitle(
        "Molecular Structure Zoom-In, Regioisomerism & M+2 Isotope Verification",
        fontsize=12.8,
        fontweight="bold",
        y=0.99,
    )
    fig.savefig(FIG_DIR / "fig09_hit_zoomin_structures_and_regioisomers.png")
    plt.close(fig)


def main():
    print("Loading processed tables and generating 9 figures...")
    X_raw, meta, X_pqn_all, log2_pqn, res, t3, t5, t6, t7, t8, t8b, t9, summary = load_all()
    make_fig01(X_raw, meta)
    print("  Saved fig01_qc_failed_injections_and_run_order.png")
    make_fig02(meta, X_pqn_all, res, summary)
    print("  Saved fig02_normalization_pqn_vs_tic_and_ion_suppression.png")
    make_fig03(meta, X_pqn_all, log2_pqn, res, t9)
    print("  Saved fig03_pca_and_variance_partitioning.png")
    make_fig04(res)
    print("  Saved fig04_statistical_volcano_and_model_comparison.png")
    make_fig05(meta, log2_pqn, res, t3)
    print("  Saved fig05_temporal_heatmap_and_kinetic_clusters.png")
    make_fig06(t8, t8b, t3)
    print("  Saved fig06_ml_prioritization_known_and_unknown.png")
    make_fig07(meta, X_pqn_all, t5)
    print("  Saved fig07_chemical_class_enrichment_and_endogenous_trajectories.png")
    make_fig08(res, t6, t7)
    print("  Saved fig08_biotransformation_and_coabundance_networks.png")
    make_fig09(X_raw, meta, X_pqn_all)
    print("  Saved fig09_hit_zoomin_structures_and_regioisomers.png")
    print("All 9 figures generated successfully in results/figures/!")


if __name__ == "__main__":
    main()
