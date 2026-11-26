"""
27200 Data-driven Bioengineering — ME2 Group Project Pipeline
Untargeted LC-MS/MS Metabolomics: Omeprazole vs DMSO in Com20 Synthetic Gut Community (T0-T8)

Reproducible end-to-end analysis covering:
  Step 1: Experimental design, failed-injection QC, and run-order confounder audit
  Step 2: Feature filtering (replicate prevalence + blank ratio), PQN vs TIC normalization,
          imputation, and log2 transformation + ESI co-elution ion-suppression diagnosis
  Step 3: Statistical modeling (run-order-adjusted linear model, two-way categorical ANOVA,
          within-treatment temporal tests, BH-FDR correction, multivariate PERMANOVA)
  Step 4: Unsupervised ML (global & clean-endogenous PCA, temporal k-means kinetic modules)
          and Supervised ML feature prioritization (Stratified CV Random Forest & L1-LogReg)
  Step 5: Annotation curation, carryover & false-positive audit, RDKit chemical class
          ontology, and Fisher's exact Over-Representation Analysis (ORA)
  Step 6: Multi-layer Network Analysis (IIMN ion-identity collapse, ChemProp-style exact-mass
          biotransformation network, and endogenous metabolic co-abundance network + hubs/modules)
  Step 7: Molecular zoom-in (exact-mass formula & 34S/37Cl isotope verification, regioisomer
          fragmentation evidence, and RDKit Morgan-fingerprint Tanimoto structural similarity)
"""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
TABLES_DIR = ROOT / "results" / "tables"
NET_DIR = ROOT / "results" / "networks"
FIG_DIR = ROOT / "results" / "figures"

for d in (PROC_DIR, TABLES_DIR, NET_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Constants
RNG_SEED = 42
SUPPRESSION_RT_MIN = 1.80
SUPPRESSION_RT_MAX = 2.25
MASS = {
    "C": 12.0,
    "H": 1.007825032,
    "N": 14.003074005,
    "O": 15.994914620,
    "S": 31.972071174,
}
PROTON = 1.007276467


# ---------------------------------------------------------------------------
# Annotation curation & RDKit chemical class ontology
# ---------------------------------------------------------------------------
CARRYOVER_DRUG_PATTERNS = [
    r"\berythromycin\b",
    r"\bclarithromycin\b",
    r"\bmetronidazole\b",
    r"\bket\b",
    r"\bketoconazole\b",
    r"\bsertraline\b",
    r"\bsimvastatin\b",
    r"\blovastatin\b",
    r"\bloratadine\b",
    r"\benclomiphene\b",
    r"\bclomifene",
    r"\bmontelukast\b",
    r"\btelmisartan\b",
    r"\blansoprazole\b",
]
QCMIX_STD_PATTERNS = [
    r"\bsulfamethizole\b",
    r"\bsulfachloropyridazine\b",
    r"\bsulfadimethoxine\b",
    r"\bsulfisomidine\b",
    r"\btrimethoprim\b",
    r"\bamitriptyline\b",
]
CONTAMINANT_PATTERNS = [
    r"tris\(2-butoxyethyl\) phosphate",
    r"cyclotetrasiloxane",
    r"lauryl diethanolamide",
    r"benzododecinium",
    r"dicyclohexylurea",
    r"triphenylphosphine oxide",
    r"2-mercaptobenzothiazole",
    r"2,4,6-tri-tert-butylphenol",
    r"diethylamino hydroxybenzoyl",
    r"4-hydroxy-1-\(2-hydroxyethyl\)-2,2,6,6-tetramethylpiperidine",
    r"2,2'-\[\(4-methylphenyl\)imino\]diethanol",
    r"triallyloxy-1,3,5-triazine",
]
OMEPRAZOLE_PATTERNS = [
    r"\bomeprazole\b",
    r"\bufiprazole\b",
    r"37052-78-1",
    r"2-mercaptobenzimidazole",
    r"benzodiazol-2-yl",
]


def classify_annotation(row: pd.Series) -> tuple[str, str]:
    """Return (chemical_class, origin_category) using curated rules + RDKit SMARTS."""
    name = str(row.get("compound_name", "") or "")
    nl = name.lower()
    method = str(row.get("method", "") or "")
    smiles = row.get("smiles", None)
    fid = int(row.name)

    # Suspect false-positive spectral library hits on Omeprazole-only features
    if fid in {5664, 2658}:
        return ("Omeprazole-Specific (Suspect Lib Match)", "Omeprazole-derived")

    for pat in OMEPRAZOLE_PATTERNS:
        if re.search(pat, nl):
            return ("Omeprazole & Benzimidazole Derivatives", "Omeprazole-derived")

    for pat in CARRYOVER_DRUG_PATTERNS:
        if re.search(pat, nl):
            return ("Study Carryover Drugs (Other 11 Drugs)", "Technical/Carryover")

    for pat in QCMIX_STD_PATTERNS:
        if re.search(pat, nl):
            return ("QC-Mix Standards (Carryover)", "Technical/Carryover")

    for pat in CONTAMINANT_PATTERNS:
        if re.search(pat, nl):
            return ("Plasticizers & Background Contaminants", "Technical/Contaminant")

    # Cyclic dipeptides / diketopiperazines
    if (
        re.search(r"cyclo[\-\(]", nl)
        or "piperazine-2,5-dione" in nl
        or "2,5-piperazinedione" in nl
        or "pyrrolo[1,2-a]pyrazine-1,4-dione" in nl
        or "phenylalanine, prolyl-" in nl
    ):
        return ("Cyclic Dipeptides (DKPs)", "Endogenous/Medium")

    # Bile acids
    if "cholic acid" in nl or "choldienic acid" in nl or "lithocholic" in nl:
        return ("Bile Acids", "Endogenous/Medium")

    # N-Acyl amino acids / microbial amino acid conjugates
    if (
        re.search(r"^(gly|leu|phe)-c\d+", nl)
        or re.search(r"candidate .+-c\d+:\d+", nl)
        or "heptadecanoyl)glycine" in nl
        or "oleoyl-phenylalanine" in nl
        or "oleyl sarcosine" in nl
        or "propanoyl]norleucine" in nl
    ):
        return ("N-Acyl Amino Acids", "Endogenous/Medium")

    # Indoles & tryptophan derivatives
    if (
        "indole" in nl
        or "tryptophan" in nl
        or "tryptophol" in nl
        or "tryptamine" in nl
        or "abrine" in nl
        or "serotonin" in nl
        or "ficz" in nl
        or "indolizine" in nl
        or "86704-55-4" in nl
    ):
        return ("Indoles & Tryptophan Derivatives", "Endogenous/Medium")

    # Nucleosides, purines & biogenic amines
    if (
        "adenosine" in nl
        or "guanine" in nl
        or "xanthine" in nl
        or "histamine" in nl
        or "2365-40-4" in nl
        or "pyridox" in nl
    ):
        return ("Nucleosides, Purines & Amines", "Endogenous/Medium")

    # Amino acids & linear peptides
    if (
        any(
            k in nl
            for k in [
                "arginine",
                "histidine",
                "tyrosine",
                "phenylalanine",
                "phe-leu",
                "phe-pro",
                "phe-phe",
                "phe-trp",
                "leu-val",
                "val-leu",
                "leu-gly",
                "tyramine",
                "phenethylacetamide",
                "amino-2-phenylacetate",
            ]
        )
    ):
        return ("Amino Acids & Linear Peptides", "Endogenous/Medium")

    # Isoflavones, polyphenols & plant-derived medium components
    if any(
        k in nl
        for k in [
            "daidzin",
            "daidzein",
            "genistein",
            "glycitein",
            "naringenin",
            "tricin",
            "vanillin",
            "schisandrin",
            "578-74-5",
            "trimethoxyphenol",
            "4-hydroxy-3-methoxyphenyl",
        ]
    ):
        return ("Plant Polyphenols & Isoflavones (mGAM)", "Endogenous/Medium")

    # Lipids, phospholipids, sphingolipids, fatty acids & oxylipins
    if method == "matchedlipid" or any(
        k in nl
        for k in [
            "oleic acid",
            "linoleic acid",
            "palmitic acid",
            "stearic acid",
            "myristoleic",
            "eicosenoic",
            "octadecadiynoic",
            "octadecenamide",
            "palmitamide",
            "docosenamide",
            "ethanolamide",
            "epome",
            "epoxy-9z",
            "ekode",
            "hode",
            "hede",
            "hydroxydecanoic",
            "hydroxydodecanoic",
            "hydroxynonenal",
            "octadecadienoic",
            "octadecadien-1-ol",
            "monoelaidin",
            "glycerol",
            "sphinganine",
            "cer ",
            "lpe ",
            "pe ",
            "nae ",
            "spb ",
            "lactone",
            "oxiran-2-yl",
            "nerolidol",
            "camphor",
        ]
    ):
        return ("Fatty Acids, Oxylipins & Lipids", "Endogenous/Medium")

    # RDKit SMARTS fallback if SMILES present
    if isinstance(smiles, str) and smiles.strip():
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            if mol.HasSubstructMatch(Chem.MolFromSmarts("C1(=O)NCC(=O)NC1")):
                return ("Cyclic Dipeptides (DKPs)", "Endogenous/Medium")
            if mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")):
                return ("Indoles & Tryptophan Derivatives", "Endogenous/Medium")
            if mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccc2[nH]cnc2c1")):
                return ("Omeprazole & Benzimidazole Derivatives", "Omeprazole-derived")
            if Descriptors.MolWt(mol) > 180 and mol.HasSubstructMatch(
                Chem.MolFromSmarts("CCCCCCCC")
            ):
                return ("Fatty Acids, Oxylipins & Lipids", "Endogenous/Medium")

    return ("Other Annotated Metabolites", "Endogenous/Medium")


def solve_exact_formulas(mz_obs: float, s_range=(0, 2), ppm_tol: float = 5.0) -> list[dict]:
    """Find neutral molecular formulas (C, H, N, O, S) matching [M+H]+ within ppm_tol."""
    target_neutral = mz_obs - PROTON
    hits = []
    max_c = min(28, int(target_neutral // 12))
    for s in range(s_range[0], s_range[1] + 1):
        for n in range(0, 7):
            for o in range(0, 8):
                rem = target_neutral - (s * MASS["S"] + n * MASS["N"] + o * MASS["O"])
                if rem < 0:
                    continue
                for c in range(4, max_c + 1):
                    h_rem = rem - c * MASS["C"]
                    if h_rem < 0:
                        break
                    h = int(round(h_rem / MASS["H"]))
                    if h < 1 or h > 2 * c + n + 2:
                        continue
                    if (h + n) % 2 != 0:  # Nitrogen rule
                        continue
                    rdb = c - h / 2.0 + n / 2.0 + 1.0
                    if rdb < 0 or rdb > 18:
                        continue
                    m_calc = (
                        c * MASS["C"]
                        + h * MASS["H"]
                        + n * MASS["N"]
                        + o * MASS["O"]
                        + s * MASS["S"]
                        + PROTON
                    )
                    ppm = (mz_obs - m_calc) / m_calc * 1e6
                    if abs(ppm) <= ppm_tol:
                        formula = (
                            f"C{c}H{h}"
                            + (f"N{n}" if n else "")
                            + (f"O{o}" if o else "")
                            + (f"S{s}" if s else "")
                        )
                        hits.append(
                            {
                                "formula": formula,
                                "ppm": round(ppm, 2),
                                "abs_ppm": abs(ppm),
                                "rdb": rdb,
                                "S": s,
                            }
                        )
    hits.sort(key=lambda d: (d["abs_ppm"], d["S"]))
    return hits


def impute_and_log2(df: pd.DataFrame, frac_min: float = 0.2) -> pd.DataFrame:
    """Per-feature left-censored imputation (1/5 of minimum positive value) followed by log2."""
    arr = df.values.astype(float).copy()
    for i in range(arr.shape[0]):
        pos = arr[i, arr[i] > 0]
        fill = pos.min() * frac_min if len(pos) > 0 else 1.0
        arr[i, arr[i] == 0] = fill
    return pd.DataFrame(np.log2(arr), index=df.index, columns=df.columns)


def ols_vectorized(Xmat: np.ndarray, Ymat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized OLS returning (beta, tstat, pvals) of shape (k_predictors, n_features)."""
    n, k = Xmat.shape
    df_res = n - k
    beta = np.linalg.lstsq(Xmat, Ymat, rcond=None)[0]
    resid = Ymat - Xmat @ beta
    s2 = (resid ** 2).sum(axis=0) / df_res
    xtx_inv_diag = np.diag(np.linalg.inv(Xmat.T @ Xmat))
    se = np.sqrt(np.maximum(np.outer(xtx_inv_diag, s2), 1e-30))
    tstat = beta / se
    pvals = 2.0 * stats.t.sf(np.abs(tstat), df=df_res)
    return beta, tstat, pvals


def permanova_test(
    Ymat: np.ndarray,
    design_df: pd.DataFrame,
    n_perm: int = 999,
    seed: int = RNG_SEED,
) -> pd.DataFrame:
    """
    Euclidean-distance multivariate PERMANOVA (exact equivalence to centered multivariate
    variance partitioning) with permutation p-values for Treatment, Time, Interaction, and Run Date.
    """
    rng = np.random.default_rng(seed)
    n = Ymat.shape[0]
    Yc = Ymat - Ymat.mean(axis=0, keepdims=True)
    sst = float((Yc ** 2).sum())

    D_trt = pd.get_dummies(design_df["ATTRIBUTE_Drug"], drop_first=True, dtype=float).values
    D_time = pd.get_dummies(design_df["ATRRIBUTE_SampleType"], drop_first=True, dtype=float).values
    D_int = D_trt * D_time
    D_run = design_df[["run_hours"]].values
    D_run = (D_run - D_run.mean()) / D_run.std()

    # Precompute orthonormal bases Q via QR for ultra-fast vectorized permutations
    def make_Q(X):
        Xa = np.column_stack([np.ones(n), X])
        Q, _ = np.linalg.qr(Xa)
        return Q

    Q_trt = make_Q(D_trt)
    Q_time = make_Q(D_time)
    Q_run = make_Q(D_run)
    Q_add = make_Q(np.column_stack([D_trt, D_time]))
    Q_full = make_Q(np.column_stack([D_trt, D_time, D_int]))
    Q_full_run = make_Q(np.column_stack([D_trt, D_time, D_int, D_run]))

    def ssr_q(Q, Ym):
        proj = Q.T @ Ym
        return float((proj ** 2).sum())

    # Marginal / Type-II breakdown
    ss_trt = ssr_q(Q_trt, Yc)
    ss_time = ssr_q(Q_time, Yc)
    ss_run = ssr_q(Q_run, Yc)
    ss_add = ssr_q(Q_add, Yc)
    ss_full = ssr_q(Q_full, Yc)
    ss_int = ss_full - ss_add
    ss_res = sst - ss_full
    df_res = n - (1 + D_trt.shape[1] + D_time.shape[1] + D_int.shape[1])

    ss_full_plus_run = ssr_q(Q_full_run, Yc)
    ss_run_cond = ss_full_plus_run - ss_full

    # Type II partial sums of squares for Treatment and Time
    ss_trt_part = ss_add - ss_time
    ss_time_part = ss_add - ss_trt

    F_obs = {
        "Treatment (Omeprazole vs DMSO)": (ss_trt_part / 1.0) / (ss_res / df_res),
        "Time (T0-T8, categorical)": (ss_time_part / 8.0) / (ss_res / df_res),
        "Treatment x Time Interaction": (ss_int / 8.0) / (ss_res / df_res),
        "Run Order (marginal)": (ss_run / 1.0) / ((sst - ss_run) / (n - 2)),
        "Run Order (conditional on Trt*Time)": (ss_run_cond / 1.0)
        / ((sst - ss_full_plus_run) / (df_res - 1)),
    }
    R2_obs = {
        "Treatment (Omeprazole vs DMSO)": ss_trt_part / sst,
        "Time (T0-T8, categorical)": ss_time_part / sst,
        "Treatment x Time Interaction": ss_int / sst,
        "Run Order (marginal)": ss_run / sst,
        "Run Order (conditional on Trt*Time)": ss_run_cond / sst,
    }
    df_map = {
        "Treatment (Omeprazole vs DMSO)": 1,
        "Time (T0-T8, categorical)": 8,
        "Treatment x Time Interaction": 8,
        "Run Order (marginal)": 1,
        "Run Order (conditional on Trt*Time)": 1,
    }

    # Fast QR permutation tests (row permutations of Yc)
    ge_counts = {k: 1 for k in F_obs}
    for _ in range(n_perm):
        perm = rng.permutation(n)
        Yp = Yc[perm]

        p_trt = ssr_q(Q_trt, Yp)
        p_time = ssr_q(Q_time, Yp)
        p_run = ssr_q(Q_run, Yp)
        p_add = ssr_q(Q_add, Yp)
        p_full = ssr_q(Q_full, Yp)
        p_full_run = ssr_q(Q_full_run, Yp)
        p_res = sst - p_full

        f_perm = {
            "Treatment (Omeprazole vs DMSO)": ((p_add - p_time) / 1.0) / (p_res / df_res),
            "Time (T0-T8, categorical)": ((p_add - p_trt) / 8.0) / (p_res / df_res),
            "Treatment x Time Interaction": ((p_full - p_add) / 8.0) / (p_res / df_res),
            "Run Order (marginal)": (p_run / 1.0) / ((sst - p_run) / (n - 2)),
            "Run Order (conditional on Trt*Time)": ((p_full_run - p_full) / 1.0)
            / ((sst - p_full_run) / (df_res - 1)),
        }
        for k in F_obs:
            if f_perm[k] >= F_obs[k]:
                ge_counts[k] += 1

    rows = []
    for k in F_obs:
        rows.append(
            {
                "Term": k,
                "df": df_map[k],
                "R2": R2_obs[k],
                "pseudo_F": F_obs[k],
                "p_perm": ge_counts[k] / (n_perm + 1),
            }
        )
    rows.append(
        {
            "Term": "Residual (Replicate + Unexplained)",
            "df": df_res,
            "R2": ss_res / sst,
            "pseudo_F": np.nan,
            "p_perm": np.nan,
        }
    )
    return pd.DataFrame(rows)


def main():
    print("=" * 72)
    print("27200 ME2 METABOLOMICS PIPELINE — STARTING FULL ANALYSIS")
    print("=" * 72)

    # ------------------------------------------------------------------
    # STEP 1: Load data & perform sample-level QC + run-order audit
    # ------------------------------------------------------------------
    meta = pd.read_csv(RAW_DIR / "group2_drugs_metadata.csv")
    feat = pd.read_csv(RAW_DIR / "group2_drugs_iimn_gnps_quant.csv")
    ann = pd.read_csv(RAW_DIR / "group2_drugs_annotations.csv")

    unnamed = [c for c in feat.columns if c.startswith("Unnamed")]
    feat = feat.drop(columns=unnamed)

    peak_cols = [c for c in feat.columns if c.endswith(" Peak area")]
    peak_files = [c.replace(" Peak area", "") for c in peak_cols]
    assert set(peak_files) == set(meta["filename"]), "Sample filename mismatch!"

    X_raw = feat[peak_cols].copy()
    X_raw.columns = peak_files
    X_raw.index = feat["row ID"]
    feat_meta = feat.drop(columns=peak_cols).set_index("row ID")

    # Standardize metadata
    meta["sample_type"] = meta["ATRRIBUTE_SampleType"]
    meta["drug"] = meta["ATTRIBUTE_Drug"]
    meta["replicate"] = meta["ATTRIBUTE_Replicate"]
    meta["bacteria"] = meta["ATTRIBUTE_Bacteria"]
    meta["time_num"] = meta["sample_type"].str.extract(r"^T(\d+)$").astype(float)
    meta["run_date_dt"] = pd.to_datetime(meta["run_date"], format="ISO8601")
    meta["run_hours"] = (
        meta["run_date_dt"] - meta["run_date_dt"].min()
    ).dt.total_seconds() / 3600.0
    meta["TIC"] = [X_raw[f].sum() for f in meta["filename"]]
    meta["n_nonzero"] = [(X_raw[f] > 0).sum() for f in meta["filename"]]
    meta["frac_zero"] = [(X_raw[f] == 0).mean() for f in meta["filename"]]
    meta["median_nonzero"] = [
        X_raw[f][X_raw[f] > 0].median() for f in meta["filename"]
    ]

    # Automated failed-injection detection: bacterial sample whose TIC < 2 * max(Blank TIC)
    blank_files = meta.loc[meta["mzmine_sample_type"] == "blank", "filename"].tolist()
    max_blank_tic = meta.loc[meta["filename"].isin(blank_files), "TIC"].max()
    meta["is_failed_injection"] = (meta["bacteria"] == "Yes") & (
        meta["TIC"] < 2.0 * max_blank_tic
    )
    failed_files = meta.loc[meta["is_failed_injection"], "filename"].tolist()
    print(f"[Step 1] Total injections: {len(meta)} | Features: {len(X_raw)}")
    print(
        f"[Step 1] Detected {len(failed_files)} failed/blank-like bacterial injections "
        f"(TIC < 2x max Blank = {2*max_blank_tic:.2e}): {failed_files}"
    )

    # Verify failed-injection correlation to blanks vs true replicates
    log1p_all = np.log1p(X_raw)
    sample_corr = log1p_all.corr()
    meta["mean_corr_to_blanks"] = [
        sample_corr.loc[f, blank_files].mean() for f in meta["filename"]
    ]

    pooled_qc_files = meta.loc[meta["sample_type"] == "Pooled_QC", "filename"].tolist()
    qc_mix_files = meta.loc[meta["sample_type"] == "QC", "filename"].tolist()
    medium_files = meta.loc[meta["drug"] == "Medium", "filename"].tolist()
    bac_meta = meta.loc[
        (meta["bacteria"] == "Yes") & (~meta["is_failed_injection"])
    ].copy()
    bac_files = bac_meta["filename"].tolist()
    dmso_files = bac_meta.loc[bac_meta["drug"] == "DMSO", "filename"].tolist()
    omepr_files = bac_meta.loc[bac_meta["drug"] == "Omeprazole", "filename"].tolist()

    # Check run-order confounding with timepoint
    r_run_dmso_design, p_run_dmso_design = stats.spearmanr(
        bac_meta.loc[bac_meta["drug"] == "DMSO", "time_num"],
        bac_meta.loc[bac_meta["drug"] == "DMSO", "run_hours"],
    )
    r_run_omepr_design, p_run_omepr_design = stats.spearmanr(
        bac_meta.loc[bac_meta["drug"] == "Omeprazole", "time_num"],
        bac_meta.loc[bac_meta["drug"] == "Omeprazole", "run_hours"],
    )
    print(
        f"[Step 1] Run-order vs Timepoint confounding: DMSO r={r_run_dmso_design:.3f} (p={p_run_dmso_design:.4f}), "
        f"Omeprazole r={r_run_omepr_design:.3f} (p={p_run_omepr_design:.4f})"
    )

    # ------------------------------------------------------------------
    # STEP 2: Feature filtering & PQN vs TIC normalization
    # ------------------------------------------------------------------
    rep_counts = pd.DataFrame(
        {
            f"{d}_{t}": (X_raw[s["filename"].tolist()] > 0).sum(axis=1)
            for (d, t), s in bac_meta.groupby(["drug", "sample_type"])
        }
    )
    tp_raw_means = pd.DataFrame(
        {
            f"{d}_{t}": X_raw[s["filename"].tolist()].mean(axis=1)
            for (d, t), s in bac_meta.groupby(["drug", "sample_type"])
        }
    )
    mean_blank = X_raw[blank_files].mean(axis=1)

    pass_prevalence = rep_counts.max(axis=1) >= 2
    pass_blank = tp_raw_means.max(axis=1) > 3.0 * mean_blank
    keep_mask = pass_prevalence & pass_blank

    # Document Pooled QC stats (diagnostic, NOT strict presence filter!)
    pqc_det = (X_raw[pooled_qc_files] > 0).sum(axis=1)
    pqc_mean = X_raw[pooled_qc_files].mean(axis=1)
    pqc_rsd = X_raw[pooled_qc_files].std(axis=1, ddof=1) / pqc_mean.replace(0, np.nan)

    print(
        f"[Step 2] Feature filtering: {len(X_raw)} raw -> "
        f"{pass_prevalence.sum()} (det>=2 in >=1 group) -> "
        f"{keep_mask.sum()} (max group mean > 3x blank)"
    )

    Xf_raw = X_raw.loc[keep_mask, bac_files + medium_files].copy()

    # Probabilistic Quotient Normalization (PQN)
    # Reference spectrum = median across DMSO bacterial samples for ubiquitous shared features (>=80% det in DMSO)
    shared_ref_mask = (Xf_raw[dmso_files] > 0).mean(axis=1) >= 0.80
    ref_spectrum = Xf_raw.loc[shared_ref_mask, dmso_files].median(axis=1)
    ref_spectrum = ref_spectrum[ref_spectrum > 0]
    print(f"[Step 2] PQN reference spectrum built from {len(ref_spectrum)} ubiquitous DMSO features")

    pqn_factors = {}
    tic_factors = {}
    med_tic_bac = Xf_raw[bac_files].sum(axis=0).median()
    for col in bac_files + medium_files:
        ratios = Xf_raw.loc[ref_spectrum.index, col] / ref_spectrum
        ratios = ratios[ratios > 0]
        pqn_factors[col] = float(ratios.median())
        tic_factors[col] = float(Xf_raw[col].sum() / med_tic_bac)

    meta["pqn_factor"] = meta["filename"].map(pqn_factors)
    meta["tic_factor"] = meta["filename"].map(tic_factors)

    # Save sample QC summary
    meta.to_csv(PROC_DIR / "sample_qc_summary.csv", index=False)
    meta.to_csv(TABLES_DIR / "table01_sample_qc_and_design.csv", index=False)

    X_pqn_all = Xf_raw.div(pd.Series(pqn_factors), axis=1)
    X_pqn = X_pqn_all[bac_files].copy()
    X_med_pqn = X_pqn_all[medium_files].copy()
    X_tic = Xf_raw[bac_files].div(pd.Series({c: tic_factors[c] for c in bac_files}), axis=1)

    # Impute & log2 transform (bacterial samples only for modeling; also build combined with medium for PCA)
    log2_pqn = impute_and_log2(X_pqn)
    log2_tic = impute_and_log2(X_tic)
    log2_pqn_with_med = impute_and_log2(X_pqn_all)

    X_pqn_all.to_csv(PROC_DIR / "features_filtered_pqn.csv")
    log2_pqn.to_csv(PROC_DIR / "features_filtered_log2_pqn.csv")

    # ------------------------------------------------------------------
    # STEP 5A (early): Curate annotations & RDKit descriptors
    # ------------------------------------------------------------------
    ann_sorted = ann.sort_values(["method", "score"], ascending=[False, False])
    ann_best = ann_sorted.drop_duplicates("id").set_index("id")
    ann_kept = ann_best.loc[ann_best.index.intersection(X_pqn.index)].copy()

    classes, origins, mws, logps = [], [], [], []
    for fid, row in ann_kept.iterrows():
        cls, orig = classify_annotation(row)
        classes.append(cls)
        origins.append(orig)
        smi = row.get("smiles", None)
        if isinstance(smi, str) and smi.strip():
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                mws.append(round(Descriptors.MolWt(m), 3))
                logps.append(round(Descriptors.MolLogP(m), 3))
            else:
                mws.append(np.nan)
                logps.append(np.nan)
        else:
            mws.append(np.nan)
            logps.append(np.nan)

    ann_kept["chemical_class"] = classes
    ann_kept["origin_category"] = origins
    ann_kept["rdkit_mw"] = mws
    ann_kept["rdkit_logp"] = logps
    ann_kept["n_library_rows"] = ann.groupby("id").size().reindex(ann_kept.index)
    ann_kept.to_csv(PROC_DIR / "feature_annotations_curated.csv")
    print(f"[Step 5A] Curated {len(ann_kept)} annotated features across {ann_kept['chemical_class'].nunique()} chemical classes")

    # ------------------------------------------------------------------
    # STEP 3: Statistical Modeling (Linear + Categorical ANOVA + Within-Trt + TIC check)
    # ------------------------------------------------------------------
    design = bac_meta.set_index("filename").loc[bac_files].copy()
    design["trt"] = (design["drug"] == "Omeprazole").astype(float)
    design["t_c"] = design["time_num"] - design["time_num"].mean()
    design["trt_t"] = design["trt"] * design["t_c"]
    design["run_z"] = (design["run_hours"] - design["run_hours"].mean()) / design["run_hours"].std()

    Y = log2_pqn[bac_files].values.T  # (50, 3414)

    # Model A: Run-order-adjusted linear model: y ~ 1 + trt + time + trt:time + run_z
    X_lin = np.column_stack(
        [np.ones(len(design)), design["trt"], design["t_c"], design["trt_t"], design["run_z"]]
    )
    beta_lin, t_lin, p_lin = ols_vectorized(X_lin, Y)

    # Also fit unadjusted PQN vs TIC linear models on SHARED features for direct TIC-bias quantification
    X_lin_no_run = np.column_stack(
        [np.ones(len(design)), design["trt"], design["t_c"], design["trt_t"]]
    )
    beta_tic, _, p_tic = ols_vectorized(X_lin_no_run, log2_tic[bac_files].values.T)
    beta_pqn_nr, _, p_pqn_nr = ols_vectorized(X_lin_no_run, Y)

    # Model B: Two-way Categorical ANOVA: y ~ Treatment * factor(Time) (captures non-linear transients!)
    D_trt = pd.get_dummies(design["drug"], drop_first=True, dtype=float).values
    D_time = pd.get_dummies(design["sample_type"], drop_first=True, dtype=float).values
    D_int = D_trt * D_time
    X_full = np.column_stack([np.ones(len(design)), D_trt, D_time, D_int])
    X_noint = np.column_stack([np.ones(len(design)), D_trt, D_time])
    X_notrt = np.column_stack([np.ones(len(design)), D_time])
    X_notime = np.column_stack([np.ones(len(design)), D_trt])

    def calc_rss(Xm):
        b = np.linalg.lstsq(Xm, Y, rcond=None)[0]
        return ((Y - Xm @ b) ** 2).sum(axis=0)

    rss_full = calc_rss(X_full)
    rss_noint = calc_rss(X_noint)
    rss_notrt = calc_rss(X_notrt)
    rss_notime = calc_rss(X_notime)
    df_res_cat = len(design) - X_full.shape[1]  # 32
    ms_res_cat = np.maximum(rss_full / df_res_cat, 1e-12)

    F_cat_int = ((rss_noint - rss_full) / 8.0) / ms_res_cat
    F_cat_trt = ((rss_notrt - rss_noint) / 1.0) / ms_res_cat
    F_cat_time = ((rss_notime - rss_noint) / 8.0) / ms_res_cat

    p_cat_int = stats.f.sf(F_cat_int, 8, df_res_cat)
    p_cat_trt = stats.f.sf(F_cat_trt, 1, df_res_cat)
    p_cat_time = stats.f.sf(F_cat_time, 8, df_res_cat)

    # Model C: Within-Omeprazole and Within-DMSO temporal Spearman correlations & 1-way ANOVA
    omepr_t = design.loc[omepr_files, "time_num"].values
    omepr_run = design.loc[omepr_files, "run_hours"].values
    dmso_t = design.loc[dmso_files, "time_num"].values

    r_time_omepr, p_time_omepr, r_run_omepr, r_time_dmso, p_time_dmso = [], [], [], [], []
    Y_omepr = log2_pqn[omepr_files].values
    Y_dmso = log2_pqn[dmso_files].values
    for i in range(len(log2_pqn)):
        rt_o, pt_o = stats.spearmanr(omepr_t, Y_omepr[i])
        rr_o, _ = stats.spearmanr(omepr_run, Y_omepr[i])
        rt_d, pt_d = stats.spearmanr(dmso_t, Y_dmso[i])
        r_time_omepr.append(0.0 if np.isnan(rt_o) else rt_o)
        p_time_omepr.append(1.0 if np.isnan(pt_o) else pt_o)
        r_run_omepr.append(0.0 if np.isnan(rr_o) else rr_o)
        r_time_dmso.append(0.0 if np.isnan(rt_d) else rt_d)
        p_time_dmso.append(1.0 if np.isnan(pt_d) else pt_d)

    # Per-timepoint PQN means & replicate consistency (median within-group CV)
    tp_pqn_means = pd.DataFrame(
        {
            f"{d}_{t}": X_pqn[s["filename"].tolist()].mean(axis=1)
            for (d, t), s in bac_meta.groupby(["drug", "sample_type"])
        }
    )
    tp_log2_means = pd.DataFrame(
        {
            f"{d}_{t}": log2_pqn[s["filename"].tolist()].mean(axis=1)
            for (d, t), s in bac_meta.groupby(["drug", "sample_type"])
        }
    )
    # Within-group SD on log2 scale (lower = higher replicate consistency)
    tp_log2_sds = pd.DataFrame(
        {
            f"{d}_{t}": log2_pqn[s["filename"].tolist()].std(axis=1, ddof=1)
            for (d, t), s in bac_meta.groupby(["drug", "sample_type"])
        }
    )

    det_dmso = (Xf_raw[dmso_files] > 0).sum(axis=1)
    det_omepr = (Xf_raw[omepr_files] > 0).sum(axis=1)
    det_med = (Xf_raw[medium_files] > 0).sum(axis=1)

    feature_group = np.where(
        (det_dmso >= 5) & (det_omepr >= 5),
        "Shared_Endogenous",
        np.where(
            (det_omepr >= 5) & (det_dmso < 5),
            "Omeprazole_Specific",
            np.where(
                (det_dmso >= 5) & (det_omepr < 5),
                "DMSO_Specific",
                "Sparse_Group_Specific",
            ),
        ),
    )

    res = pd.DataFrame(
        {
            "row_id": log2_pqn.index,
            "mz": feat_meta.loc[log2_pqn.index, "row m/z"].values,
            "rt": feat_meta.loc[log2_pqn.index, "row retention time"].values,
            "correlation_group_id": feat_meta.loc[log2_pqn.index, "correlation group ID"].values,
            "annotation_network_number": feat_meta.loc[
                log2_pqn.index, "annotation network number"
            ].values,
            "best_ion": feat_meta.loc[log2_pqn.index, "best ion"].values,
            "neutral_M_mass": feat_meta.loc[log2_pqn.index, "neutral M mass"].values,
            "feature_category": feature_group,
            "det_dmso": det_dmso.values,
            "det_omepr": det_omepr.values,
            "det_medium": det_med.values,
            "det_pooled_qc": pqc_det.loc[log2_pqn.index].values,
            "pooled_qc_rsd": pqc_rsd.loc[log2_pqn.index].values,
            "mean_blank_raw": mean_blank.loc[log2_pqn.index].values,
            "mean_medium_pqn": X_med_pqn.mean(axis=1).values,
            "mean_dmso_pqn": X_pqn[dmso_files].mean(axis=1).values,
            "mean_omepr_pqn": X_pqn[omepr_files].mean(axis=1).values,
            "max_omepr_pqn": X_pqn[omepr_files].max(axis=1).values,
            "max_dmso_pqn": X_pqn[dmso_files].max(axis=1).values,
            "omepr_T0_pqn": tp_pqn_means["Omeprazole_T0"].values,
            "omepr_T4_pqn": tp_pqn_means["Omeprazole_T4"].values,
            "omepr_T8_pqn": tp_pqn_means["Omeprazole_T8"].values,
            "dmso_T0_pqn": tp_pqn_means["DMSO_T0"].values,
            "dmso_T8_pqn": tp_pqn_means["DMSO_T8"].values,
            "peak_tp_omepr": tp_pqn_means[[f"Omeprazole_T{t}" for t in range(9)]]
            .idxmax(axis=1)
            .str.replace("Omeprazole_", "")
            .values,
            "median_rep_log2_sd": tp_log2_sds.median(axis=1).values,
            # Linear model (run-date adjusted)
            "beta_trt": beta_lin[1],
            "p_trt": p_lin[1],
            "q_trt": multipletests(p_lin[1], method="fdr_bh")[1],
            "beta_time": beta_lin[2],
            "p_time": p_lin[2],
            "q_time": multipletests(p_lin[2], method="fdr_bh")[1],
            "beta_int": beta_lin[3],
            "p_int": p_lin[3],
            "q_int": multipletests(p_lin[3], method="fdr_bh")[1],
            "beta_run": beta_lin[4],
            "p_run": p_lin[4],
            "q_run": multipletests(p_lin[4], method="fdr_bh")[1],
            # Categorical 2-way ANOVA
            "F_cat_trt": F_cat_trt,
            "p_cat_trt": p_cat_trt,
            "q_cat_trt": multipletests(p_cat_trt, method="fdr_bh")[1],
            "F_cat_time": F_cat_time,
            "p_cat_time": p_cat_time,
            "q_cat_time": multipletests(p_cat_time, method="fdr_bh")[1],
            "F_cat_int": F_cat_int,
            "p_cat_int": p_cat_int,
            "q_cat_int": multipletests(p_cat_int, method="fdr_bh")[1],
            # Within-treatment temporal trends
            "spearman_time_omepr": r_time_omepr,
            "q_spearman_time_omepr": multipletests(p_time_omepr, method="fdr_bh")[1],
            "spearman_run_omepr": r_run_omepr,
            "spearman_time_dmso": r_time_dmso,
            "q_spearman_time_dmso": multipletests(p_time_dmso, method="fdr_bh")[1],
            # TIC comparison (unadjusted linear model)
            "beta_trt_tic": beta_tic[1],
            "q_trt_tic": multipletests(p_tic[1], method="fdr_bh")[1],
            "beta_int_tic": beta_tic[3],
            "q_int_tic": multipletests(p_tic[3], method="fdr_bh")[1],
            "beta_trt_pqn_unadj": beta_pqn_nr[1],
            "q_trt_pqn_unadj": multipletests(p_pqn_nr[1], method="fdr_bh")[1],
            "beta_int_pqn_unadj": beta_pqn_nr[3],
            "q_int_pqn_unadj": multipletests(p_pqn_nr[3], method="fdr_bh")[1],
        }
    ).set_index("row_id")

    # Flag ESI co-elution ion suppression window
    res["in_ion_suppression_window"] = (res["rt"] >= SUPPRESSION_RT_MIN) & (
        res["rt"] <= SUPPRESSION_RT_MAX
    )

    # Join curated annotations
    for col in ["compound_name", "score", "mol_formula", "smiles", "adduct", "method", "chemical_class", "origin_category"]:
        res[col] = res.index.map(ann_kept[col])

    # Mark IIMN representative feature per correlation_group_id (highest mean peak area in group,
    # preferring [M+H]+ / annotated ion when within 20% of top intensity)
    res["iimn_is_representative"] = True
    for cg, grp in res[res["correlation_group_id"].notna()].groupby("correlation_group_id"):
        if len(grp) > 1:
            # Prefer annotated [M+H]+ or highest overall intensity
            score_rank = grp[["mean_omepr_pqn", "mean_dmso_pqn"]].max(axis=1).copy()
            # Boost annotated features slightly so parent [M+H]+ beats minor fragment/isotope
            is_ann = grp["compound_name"].notna()
            is_mh = grp["best_ion"] == "[M+H]+"
            # Keep highest-intensity ion unless annotated [M+H]+ is >30% of max
            max_i = score_rank.max()
            rep_id = score_rank.idxmax()
            ann_mh = grp[is_ann & (score_rank >= 0.3 * max_i)]
            if len(ann_mh) > 0:
                rep_id = score_rank.loc[ann_mh.index].idxmax()
            elif (is_mh & (score_rank >= 0.5 * max_i)).any():
                rep_id = score_rank.loc[grp[is_mh & (score_rank >= 0.5 * max_i)].index].idxmax()
            res.loc[grp.index, "iimn_is_representative"] = False
            res.loc[rep_id, "iimn_is_representative"] = True

    n_unique_entities = int(res["iimn_is_representative"].sum())
    print(
        f"[Step 3] Statistical summary across {len(res)} features ({n_unique_entities} IIMN-collapsed unique entities):"
    )
    print(
        f"  Linear (run-adj): q_trt<0.05 = {(res['q_trt']<0.05).sum()}, "
        f"q_time<0.05 = {(res['q_time']<0.05).sum()}, "
        f"q_int<0.05 = {(res['q_int']<0.05).sum()}, q_run<0.05 = {(res['q_run']<0.05).sum()}"
    )
    print(
        f"  Categorical ANOVA: q_trt<0.05 = {(res['q_cat_trt']<0.05).sum()}, "
        f"q_time<0.05 = {(res['q_cat_time']<0.05).sum()}, "
        f"q_int<0.05 = {(res['q_cat_int']<0.05).sum()}"
    )

    # Quantify TIC vs PQN bias and Ion Suppression Fisher's exact test
    sh = res[res["feature_category"] == "Shared_Endogenous"]
    tic_fake_down = int(((sh["q_trt_tic"] < 0.05) & (sh["beta_trt_tic"] < 0)).sum())
    pqn_down = int(((sh["q_trt_pqn_unadj"] < 0.05) & (sh["beta_trt_pqn_unadj"] < 0)).sum())
    tic_neg_int = int(((sh["q_int_tic"] < 0.05) & (sh["beta_int_tic"] < 0)).sum())
    pqn_neg_int = int(((sh["q_int_pqn_unadj"] < 0.05) & (sh["beta_int_pqn_unadj"] < 0)).sum())

    neg_int_sh = (sh["q_int_pqn_unadj"] < 0.05) & (sh["beta_int_pqn_unadj"] < 0)
    in_win_sh = sh["in_ion_suppression_window"]
    a_supp = int((neg_int_sh & in_win_sh).sum())
    b_supp = int((neg_int_sh & (~in_win_sh)).sum())
    c_supp = int(((~neg_int_sh) & in_win_sh).sum())
    d_supp = int(((~neg_int_sh) & (~in_win_sh)).sum())
    or_supp, p_supp = stats.fisher_exact([[a_supp, b_supp], [c_supp, d_supp]], alternative="greater")
    print(
        f"[Step 2/3] TIC bias on {len(sh)} shared features: {tic_fake_down} down (TIC) vs {pqn_down} down (PQN); "
        f"neg interaction: {tic_neg_int} (TIC) vs {pqn_neg_int} (PQN)"
    )
    print(
        f"[Step 2/3] ESI ion suppression (RT {SUPPRESSION_RT_MIN}-{SUPPRESSION_RT_MAX} min): "
        f"{a_supp}/{a_supp+b_supp} ({100*a_supp/(a_supp+b_supp):.1f}%) of neg shared interactions in {100*(a_supp+c_supp)/len(sh):.1f}% of features "
        f"(OR={or_supp:.1f}, Fisher p={p_supp:.2e})"
    )

    # ------------------------------------------------------------------
    # STEP 3B: Multivariate PERMANOVA (All vs Clean Endogenous)
    # ------------------------------------------------------------------
    clean_endo_ids = res.index[
        (res["feature_category"] == "Shared_Endogenous")
        & (~res["in_ion_suppression_window"])
        & (~res["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]))
    ]
    perm_all = permanova_test(log2_pqn[bac_files].values.T, design, n_perm=999)
    perm_all["Feature_Set"] = f"All Filtered Features (n={len(log2_pqn)})"
    perm_endo = permanova_test(log2_pqn.loc[clean_endo_ids, bac_files].values.T, design, n_perm=999)
    perm_endo["Feature_Set"] = f"Clean Endogenous Outside Suppression (n={len(clean_endo_ids)})"
    perm_tic = permanova_test(log2_tic.loc[clean_endo_ids, bac_files].values.T, design, n_perm=999)
    perm_tic["Feature_Set"] = f"Clean Endogenous under TIC Normalization (n={len(clean_endo_ids)})"
    perm_df = pd.concat([perm_all, perm_endo, perm_tic], ignore_index=True)
    perm_df.to_csv(TABLES_DIR / "table09_permanova_variance_partitioning.csv", index=False)
    print("[Step 3B] PERMANOVA complete:")
    print(perm_df[["Feature_Set", "Term", "R2", "pseudo_F", "p_perm"]].to_string(index=False))

    # ------------------------------------------------------------------
    # STEP 4: Kinetic Trajectory Clustering & Supervised ML Prioritization
    # ------------------------------------------------------------------
    # 4A. Temporal kinetic clustering of Omeprazole-enriched dynamic features (IIMN-collapsed!)
    omepr_dyn = res[
        (res["iimn_is_representative"])
        & (res["beta_trt"] > 1.5)
        & (res["q_cat_trt"] < 0.01)
        & ((res["q_cat_int"] < 0.01) | (res["q_spearman_time_omepr"] < 0.01))
        & (~res["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]))
    ].copy()
    traj_cols = [f"Omeprazole_T{t}" for t in range(9)]
    traj_mat = tp_log2_means.loc[omepr_dyn.index, traj_cols].values
    traj_z = (traj_mat - traj_mat.mean(axis=1, keepdims=True)) / (
        traj_mat.std(axis=1, keepdims=True) + 1e-8
    )
    km = KMeans(n_clusters=4, random_state=RNG_SEED, n_init=50)
    raw_labels = km.fit_predict(traj_z)
    # Classify cluster semantics accurately from each centroid's empirical shape
    t_vec = np.arange(9)
    centroids = np.array([traj_z[raw_labels == c].mean(axis=0) for c in range(4)])
    cluster_names = {}
    for c in range(4):
        cen = centroids[c]
        if cen[8] > 0.5 and cen[0] < -0.5:
            cluster_names[c] = "M1: Monotonic Accumulation (T0->T8)"
        elif cen[0] > 0.3 and cen[8] < -0.5:
            cluster_names[c] = "M4: Early Precursors Depleted (T0-T3->T8)"
        elif np.argmax(cen) in (2, 3):
            cluster_names[c] = "M3: Early-Mid Intermediates (Peak T2-T3)"
        else:
            cluster_names[c] = "M2: Mid-Late Transients (Peak T4-T5)"
    res["kinetic_module"] = pd.Series(pd.NA, index=res.index, dtype="object")
    res.loc[omepr_dyn.index, "kinetic_module"] = [cluster_names[c] for c in raw_labels]
    print(f"[Step 4A] Clustered {len(omepr_dyn)} unique Omeprazole-dynamic entities into 4 kinetic modules:")
    print(res["kinetic_module"].value_counts().to_string())

    # 4B. Supervised ML: Stratified CV Random Forest & L1-Logistic Regression
    # Task 1: Omeprazole vs DMSO classification (using IIMN-representative non-carryover features)
    ml_feats = res.index[
        (res["iimn_is_representative"])
        & (~res["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]))
    ]
    X_ml = log2_pqn.loc[ml_feats, bac_files].values.T
    y_trt = (design["drug"] == "Omeprazole").astype(int).values
    # Task 2: Early (T0-T2) vs Mid (T3-T5) vs Late (T6-T8) metabolic stage within clean endogenous features
    X_endo_ml = log2_pqn.loc[clean_endo_ids, bac_files].values.T
    y_stage = pd.cut(design["time_num"], bins=[-1, 2, 5, 9], labels=[0, 1, 2]).astype(int).values
    # Task 3: Omeprazole vs DMSO strictly on Clean Endogenous features outside suppression zone!
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG_SEED)

    ml_metrics = []
    rf_importances_trt = np.zeros(len(ml_feats))
    l1_coefs_trt = np.zeros(len(ml_feats))
    for fold, (tr, te) in enumerate(skf.split(X_ml, y_trt)):
        sc = StandardScaler()
        Xtr, Xte = sc.fit_transform(X_ml[tr]), sc.transform(X_ml[te])
        rf = RandomForestClassifier(n_estimators=400, max_features="sqrt", random_state=RNG_SEED + fold)
        rf.fit(Xtr, y_trt[tr])
        prob_rf = rf.predict_proba(Xte)[:, 1]
        rf_importances_trt += rf.feature_importances_ / 5.0

        lr = LogisticRegression(l1_ratio=1.0, C=0.5, solver="liblinear", max_iter=2000, random_state=RNG_SEED + fold)
        lr.fit(Xtr, y_trt[tr])
        prob_lr = lr.predict_proba(Xte)[:, 1]
        l1_coefs_trt += np.abs(lr.coef_[0]) / 5.0

        ml_metrics.append({"Task": "Omeprazole vs DMSO (All Unique)", "Model": "Random Forest", "Fold": fold, "ROC_AUC": roc_auc_score(y_trt[te], prob_rf), "Balanced_Acc": balanced_accuracy_score(y_trt[te], (prob_rf >= 0.5).astype(int))})
        ml_metrics.append({"Task": "Omeprazole vs DMSO (All Unique)", "Model": "L1 Logistic Regression", "Fold": fold, "ROC_AUC": roc_auc_score(y_trt[te], prob_lr), "Balanced_Acc": balanced_accuracy_score(y_trt[te], (prob_lr >= 0.5).astype(int))})

    # Clean endogenous treatment classification & temporal stage classification
    rf_importances_endo_trt = np.zeros(len(clean_endo_ids))
    for fold, (tr, te) in enumerate(skf.split(X_endo_ml, y_trt)):
        sc = StandardScaler()
        Xtr, Xte = sc.fit_transform(X_endo_ml[tr]), sc.transform(X_endo_ml[te])
        rf = RandomForestClassifier(n_estimators=400, max_features="sqrt", random_state=RNG_SEED + fold)
        rf.fit(Xtr, y_trt[tr])
        prob_rf = rf.predict_proba(Xte)[:, 1]
        rf_importances_endo_trt += rf.feature_importances_ / 5.0
        lr = LogisticRegression(l1_ratio=1.0, C=0.5, solver="liblinear", max_iter=2000, random_state=RNG_SEED + fold)
        lr.fit(Xtr, y_trt[tr])
        prob_lr = lr.predict_proba(Xte)[:, 1]
        ml_metrics.append({"Task": "Omeprazole vs DMSO (Clean Endogenous)", "Model": "Random Forest", "Fold": fold, "ROC_AUC": roc_auc_score(y_trt[te], prob_rf), "Balanced_Acc": balanced_accuracy_score(y_trt[te], (prob_rf >= 0.5).astype(int))})
        ml_metrics.append({"Task": "Omeprazole vs DMSO (Clean Endogenous)", "Model": "L1 Logistic Regression", "Fold": fold, "ROC_AUC": roc_auc_score(y_trt[te], prob_lr), "Balanced_Acc": balanced_accuracy_score(y_trt[te], (prob_lr >= 0.5).astype(int))})

    for fold, (tr, te) in enumerate(skf.split(X_endo_ml, y_stage)):
        sc = StandardScaler()
        Xtr, Xte = sc.fit_transform(X_endo_ml[tr]), sc.transform(X_endo_ml[te])
        rf = RandomForestClassifier(n_estimators=400, max_features="sqrt", random_state=RNG_SEED + fold)
        rf.fit(Xtr, y_stage[tr])
        prob_rf = rf.predict_proba(Xte)
        pred_rf = rf.predict(Xte)
        lr = LogisticRegression(C=0.5, solver="lbfgs", max_iter=2000, random_state=RNG_SEED + fold)
        lr.fit(Xtr, y_stage[tr])
        prob_lr = lr.predict_proba(Xte)
        pred_lr = lr.predict(Xte)
        ml_metrics.append({"Task": "Growth Stage T0-2/T3-5/T6-8 (Clean Endogenous)", "Model": "Random Forest", "Fold": fold, "ROC_AUC": roc_auc_score(y_stage[te], prob_rf, multi_class="ovr"), "Balanced_Acc": balanced_accuracy_score(y_stage[te], pred_rf)})
        ml_metrics.append({"Task": "Growth Stage T0-2/T3-5/T6-8 (Clean Endogenous)", "Model": "L2 Logistic Regression", "Fold": fold, "ROC_AUC": roc_auc_score(y_stage[te], prob_lr, multi_class="ovr"), "Balanced_Acc": balanced_accuracy_score(y_stage[te], pred_lr)})

    ml_metrics_df = pd.DataFrame(ml_metrics)
    ml_metrics_df.to_csv(TABLES_DIR / "table08b_ml_cv_performance.csv", index=False)

    res["rf_importance_trt"] = 0.0
    res.loc[ml_feats, "rf_importance_trt"] = rf_importances_trt
    res["l1_coef_trt"] = 0.0
    res.loc[ml_feats, "l1_coef_trt"] = l1_coefs_trt
    res["rf_importance_endo_trt"] = 0.0
    res.loc[clean_endo_ids, "rf_importance_endo_trt"] = rf_importances_endo_trt

    # Composite prioritization score for known & unknown discriminant features (effect size + significance + ML + consistency)
    cand = res.loc[ml_feats].copy()
    cand["neglog10_q"] = -np.log10(np.minimum(cand["q_cat_trt"], cand["q_cat_int"]).clip(lower=1e-50))
    cand["composite_priority_score"] = (
        0.30 * (cand["beta_trt"].abs().rank(pct=True))
        + 0.25 * (cand["neglog10_q"].rank(pct=True))
        + 0.25 * (cand["rf_importance_trt"].rank(pct=True))
        + 0.20 * ((1.0 / (cand["median_rep_log2_sd"] + 0.1)).rank(pct=True))
    )
    cand_sorted = cand.sort_values("composite_priority_score", ascending=False)
    cand_sorted.head(100).to_csv(TABLES_DIR / "table08_ml_feature_ranking.csv")
    print("[Step 4B] Supervised ML CV summary:")
    print(ml_metrics_df.groupby(["Task", "Model"])[["ROC_AUC", "Balanced_Acc"]].mean().round(3).to_string())

    # ------------------------------------------------------------------
    # STEP 5B: Chemical Class Over-Representation Analysis (ORA)
    # ------------------------------------------------------------------
    ann_res = res[res["compound_name"].notna()].copy()
    # Background for biological ORA: exclude study carryover & QC-mix standards (228 biological/drug features)
    bio_ann = ann_res[~ann_res["origin_category"].isin(["Technical/Carryover"])].copy()

    ora_rows = []
    # Test sets:
    # (1) Omeprazole-enriched (q_cat_trt < 0.05 & beta_trt > 1.0)
    # (2) Endogenous Time-Increasing (shared, q_time < 0.05 & beta_time > 0.15)
    # (3) Endogenous Time-Decreasing (shared, q_time < 0.05 & beta_time < -0.15)
    # (4) Stable mGAM Medium-Derived (det_medium >= 2, |beta_time| < 0.08, mean_medium > 1e5)
    test_sets = {
        "Omeprazole-Enriched (q<0.05, log2FC>1)": (bio_ann["q_cat_trt"] < 0.05) & (bio_ann["beta_trt"] > 1.0),
        "Com20 Time-Increasing (Bacterial Production)": (bio_ann["feature_category"] == "Shared_Endogenous")
        & (bio_ann["q_time"] < 0.05)
        & (bio_ann["beta_time"] > 0.15),
        "Com20 Time-Decreasing (Bacterial Uptake/Depletion)": (bio_ann["feature_category"] == "Shared_Endogenous")
        & (bio_ann["q_time"] < 0.05)
        & (bio_ann["beta_time"] < -0.15),
        "Stable mGAM Medium Background": (bio_ann["det_medium"] >= 2)
        & (bio_ann["beta_time"].abs() < 0.08)
        & (bio_ann["mean_medium_pqn"] > 1e5),
    }

    for set_name, hit_mask in test_sets.items():
        n_total = len(bio_ann)
        n_hits = int(hit_mask.sum())
        p_list, temp = [], []
        for cls_name, cls_grp in bio_ann.groupby("chemical_class"):
            in_cls = bio_ann["chemical_class"] == cls_name
            k = int((hit_mask & in_cls).sum())
            M = int(in_cls.sum())
            if M < 3:
                continue
            a = k
            b = n_hits - k
            c = M - k
            d = n_total - n_hits - M + k
            odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            fold_enr = (k / max(n_hits, 1)) / (M / n_total)
            p_list.append(pval)
            temp.append(
                {
                    "Hit_Set": set_name,
                    "Chemical_Class": cls_name,
                    "Hits_in_Class": k,
                    "Total_in_Class": M,
                    "Total_Hits": n_hits,
                    "Background_Total": n_total,
                    "Fold_Enrichment": round(fold_enr, 2),
                    "Odds_Ratio": round(odds, 2) if np.isfinite(odds) else np.inf,
                    "p_value": pval,
                    "Example_Compounds": "; ".join(
                        bio_ann.loc[hit_mask & in_cls, "compound_name"].head(4).tolist()
                    ),
                }
            )
        if p_list:
            qvals = multipletests(p_list, method="fdr_bh")[1]
            for r_dict, qv in zip(temp, qvals):
                r_dict["q_value"] = qv
                ora_rows.append(r_dict)

    ora_df = pd.DataFrame(ora_rows).sort_values(["Hit_Set", "p_value"])
    ora_df.to_csv(TABLES_DIR / "table05_chemical_class_enrichment.csv", index=False)
    print("[Step 5B] Significant Chemical Class Enrichments (q < 0.10):")
    print(
        ora_df[ora_df["q_value"] < 0.10][
            ["Hit_Set", "Chemical_Class", "Hits_in_Class", "Total_in_Class", "Fold_Enrichment", "q_value"]
        ].to_string(index=False)
    )

    # ------------------------------------------------------------------
    # STEP 6: Multi-Layer Network Analysis
    #   6A: ChemProp-Style Biotransformation Mass-Difference Network
    #   6B: Endogenous Metabolic Co-Abundance Network (Hubs & Modules)
    # ------------------------------------------------------------------
    # 6A: Select top 45 IIMN-collapsed Omeprazole-specific features (by peak area & dynamics) + key curated hits
    curated_key_ids = [
        2828, 2962, 2465, 2480, 1677, 3605, 2500, 2618, 2368, 1553, 1775, 2027,
        4470, 2719, 2311, 2626, 3654, 4157, 4380, 5664, 2658, 2478, 2497, 2685,
        3626, 6545, 5096, 1604, 2379, 2678,
    ]
    top_omepr_pool = res[
        (res["iimn_is_representative"])
        & (res["feature_category"].isin(["Omeprazole_Specific", "Shared_Endogenous"]))
        & (res["beta_trt"] > 2.5)
        & (~res["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]))
    ].sort_values("max_omepr_pqn", ascending=False).head(45).index.tolist()
    net_ids = sorted(set([i for i in curated_key_ids if i in res.index] + top_omepr_pool))

    TRANSFORMATIONS = [
        ("Sulfoxide Reduction / Oxidation (±O)", 15.994915),
        ("Methylation / Demethylation (±CH2)", 14.015650),
        ("Hydrogenation / Dehydrogenation (±H2)", 2.015650),
        ("Hydration / Dehydration (±H2O)", 18.010565),
        ("Acetylation (±C2H2O)", 42.010565),
        ("Hydroxyethyl / Ethoxy Shift (±C2H4O)", 44.026215),
        ("Formylation / CO Shift (±CO)", 27.994915),
        ("S-Extrusion Rearrangement (±S-O+H2)", 17.992806),
        ("C-S Pyridine Cleavage (±C9H11NO)", 149.084064),
        ("C-S Benzimidazole-S Cleavage (±C8H7N2OS)", 179.027909),
        ("Conjugation +C6H7NS (Pakkir Shah et al.)", 125.029920),
        ("Conjugation +C6H6OS (Pakkir Shah et al.)", 126.013936),
        ("Thioglycolate Disulfide Shift (±C2H2O2S)", 89.977550),
    ]

    G_bio = nx.Graph()
    for fid in net_ids:
        r = res.loc[fid]
        label = (
            str(r["compound_name"])[:28]
            if pd.notna(r["compound_name"]) and fid not in {5664, 2658}
            else f"ID {fid} ({r['mz']:.2f})"
        )
        G_bio.add_node(
            int(fid),
            label=label,
            mz=float(r["mz"]),
            rt=float(r["rt"]),
            max_omepr_pqn=float(r["max_omepr_pqn"]),
            spearman_time_omepr=float(r["spearman_time_omepr"]),
            peak_tp=str(r["peak_tp_omepr"]),
            kinetic_module=str(r["kinetic_module"]) if pd.notna(r["kinetic_module"]) else "Other",
            annotated=bool(pd.notna(r["compound_name"]) and fid not in {5664, 2658}),
        )

    bio_edges = []
    # Compute temporal trajectories in Omeprazole across T0-T8 for ChemProp2-like correlation
    omepr_tp_mat = tp_log2_means.loc[net_ids, [f"Omeprazole_T{t}" for t in range(9)]]
    for i_idx, u in enumerate(net_ids):
        for v in net_ids[i_idx + 1 :]:
            dmz = abs(res.loc[u, "mz"] - res.loc[v, "mz"])
            for tname, tmass in TRANSFORMATIONS:
                err_da = abs(dmz - tmass)
                if err_da <= 0.0045:  # strict high-resolution Orbitrap match (<4.5 mDa)
                    corr_t, _ = stats.pearsonr(omepr_tp_mat.loc[u].values, omepr_tp_mat.loc[v].values)
                    G_bio.add_edge(
                        int(u),
                        int(v),
                        transformation=tname,
                        delta_mz=round(float(res.loc[v, "mz"] - res.loc[u, "mz"]), 4),
                        mass_error_mDa=round(err_da * 1000.0, 2),
                        temporal_pearson_r=round(float(corr_t), 3),
                    )
                    bio_edges.append(
                        {
                            "source_id": int(u),
                            "source_mz": round(float(res.loc[u, "mz"]), 4),
                            "source_name": res.loc[u, "compound_name"] if pd.notna(res.loc[u, "compound_name"]) else f"Unknown_{u}",
                            "target_id": int(v),
                            "target_mz": round(float(res.loc[v, "mz"]), 4),
                            "target_name": res.loc[v, "compound_name"] if pd.notna(res.loc[v, "compound_name"]) else f"Unknown_{v}",
                            "transformation": tname,
                            "delta_mz": round(float(res.loc[v, "mz"] - res.loc[u, "mz"]), 4),
                            "mass_error_mDa": round(err_da * 1000.0, 2),
                            "temporal_pearson_r": round(float(corr_t), 3),
                        }
                    )
                    break

    bio_edges_df = pd.DataFrame(bio_edges).sort_values("mass_error_mDa")
    bio_edges_df.to_csv(TABLES_DIR / "table07_biotransformation_edges.csv", index=False)
    nx.write_graphml(G_bio, NET_DIR / "omeprazole_biotransformation_network.graphml")
    print(
        f"[Step 6A] Omeprazole Biotransformation Network: {G_bio.number_of_nodes()} nodes, "
        f"{G_bio.number_of_edges()} exact-mass (<4.5 mDa) transformation edges"
    )

    # 6B: Endogenous Annotated Co-Abundance Network (shared biological metabolites outside suppression zone)
    endo_ann_ids = res.index[
        (res["compound_name"].notna())
        & (res["feature_category"] == "Shared_Endogenous")
        & (~res["in_ion_suppression_window"])
        & (res["origin_category"] == "Endogenous/Medium")
        & (res["iimn_is_representative"])
    ]
    Y_ea = log2_pqn.loc[endo_ann_ids, bac_files].values
    corr_mat, p_corr_mat = stats.spearmanr(Y_ea, axis=1)
    n_ea = len(endo_ann_ids)
    triu_i, triu_j = np.triu_indices(n_ea, k=1)
    q_corr = multipletests(p_corr_mat[triu_i, triu_j], method="fdr_bh")[1]

    G_endo = nx.Graph()
    for fid in endo_ann_ids:
        r = res.loc[fid]
        G_endo.add_node(
            int(fid),
            label=str(r["compound_name"])[:35],
            chemical_class=str(r["chemical_class"]),
            beta_time=float(r["beta_time"]),
            beta_trt=float(r["beta_trt"]),
            beta_int=float(r["beta_int"]),
            mz=float(r["mz"]),
            rt=float(r["rt"]),
        )
    for idx_e, (ii, jj) in enumerate(zip(triu_i, triu_j)):
        r_val = float(corr_mat[ii, jj])
        if abs(r_val) >= 0.72 and q_corr[idx_e] < 0.001:
            G_endo.add_edge(
                int(endo_ann_ids[ii]),
                int(endo_ann_ids[jj]),
                spearman_r=round(r_val, 3),
                sign="positive" if r_val > 0 else "negative",
                weight=abs(r_val),
            )
    # Remove isolated nodes for cleaner module/hub metrics (keep track of connected core)
    G_endo_conn = G_endo.subgraph([n for n, deg in G_endo.degree() if deg > 0]).copy()
    communities = list(nx.community.greedy_modularity_communities(G_endo_conn, weight="weight"))
    modularity_score = nx.community.modularity(G_endo_conn, communities, weight="weight")
    node_to_mod = {}
    for m_idx, comm in enumerate(communities, 1):
        for n_id in comm:
            node_to_mod[n_id] = m_idx
            G_endo_conn.nodes[n_id]["module"] = m_idx

    deg_cent = nx.degree_centrality(G_endo_conn)
    bet_cent = nx.betweenness_centrality(G_endo_conn, weight="weight")
    endo_node_rows = []
    for n_id in G_endo_conn.nodes():
        r = res.loc[n_id]
        endo_node_rows.append(
            {
                "row_id": n_id,
                "compound_name": r["compound_name"],
                "chemical_class": r["chemical_class"],
                "module": node_to_mod[n_id],
                "degree": G_endo_conn.degree(n_id),
                "degree_centrality": round(deg_cent[n_id], 4),
                "betweenness_centrality": round(bet_cent[n_id], 4),
                "beta_time": round(float(r["beta_time"]), 3),
                "beta_trt": round(float(r["beta_trt"]), 3),
                "q_time": r["q_time"],
                "q_trt": r["q_trt"],
            }
        )
    endo_nodes_df = pd.DataFrame(endo_node_rows).sort_values(["degree", "betweenness_centrality"], ascending=False)
    endo_nodes_df.to_csv(TABLES_DIR / "table06_network_nodes_and_hubs.csv", index=False)
    nx.write_graphml(G_endo_conn, NET_DIR / "endogenous_coabundance_network.graphml")
    print(
        f"[Step 6B] Endogenous Co-Abundance Network: {G_endo_conn.number_of_nodes()} connected nodes, "
        f"{G_endo_conn.number_of_edges()} edges (|r|>=0.72, q<0.001), {len(communities)} modules "
        f"(modularity Q={modularity_score:.3f})"
    )

    # ------------------------------------------------------------------
    # STEP 7: Curated Hit Tables & Exact-Mass / Isotope / Structural Zoom-in
    # ------------------------------------------------------------------
    # Top 30 unique Omeprazole biotransformation hits (annotated + top unknowns)
    top_bio_candidates = res[
        (res["iimn_is_representative"])
        & (~res["origin_category"].isin(["Technical/Carryover", "Technical/Contaminant"]))
        & (
            ((res["feature_category"] == "Omeprazole_Specific") & (res["max_omepr_pqn"] > 1e5))
            | (res.index.isin([2962, 2828, 2480, 2465, 3605, 2500, 1677, 2368, 1553, 1775, 2027, 4470, 2719, 2311, 2626, 2618, 3654, 5664, 2658]))
        )
    ].copy()

    # Add exact-mass top formula and 34S / M+2 isotope evidence
    top_formulas, ppm_errs, iso_notes = [], [], []
    for fid in top_bio_candidates.index:
        mz_v = float(res.loc[fid, "mz"])
        rt_v = float(res.loc[fid, "rt"])
        mx_v = float(X_raw.loc[fid, omepr_files].max())
        # Check M+1.9958 (34S) within strict RT (<0.015 min), m/z (<0.004 Da), intensity ratio (2-15%), and co-abundance (r>0.7)
        co_m2 = feat_meta[
            ((feat_meta["row retention time"] - rt_v).abs() < 0.015)
            & ((feat_meta["row m/z"] - (mz_v + 1.9958)).abs() < 0.004)
        ]
        iso_found = False
        if len(co_m2) > 0 and mx_v > 0:
            for i_id in co_m2.index:
                ratio_pct = float(X_raw.loc[i_id, omepr_files].max() / mx_v * 100.0)
                r_iso, _ = stats.pearsonr(X_raw.loc[fid, omepr_files].values, X_raw.loc[i_id, omepr_files].values)
                if 2.0 <= ratio_pct <= 20.0 and r_iso > 0.70:
                    iso_notes.append(f"34S M+2 confirmed (ID {i_id}, {ratio_pct:.1f}% ~ {ratio_pct/4.25:.1f}S, r={r_iso:.2f})")
                    s_rng = (1, 2)
                    iso_found = True
                    break
        if not iso_found:
            iso_notes.append("No M+2 (low signal or 0-1 S)")
            s_rng = (0, 2)
        # Prefer known annotation formula if verified omeprazole derivative
        if pd.notna(res.loc[fid, "mol_formula"]) and fid not in {5664, 2658}:
            top_formulas.append(str(res.loc[fid, "mol_formula"]))
            ppm_errs.append(0.0)
        else:
            f_hits = solve_exact_formulas(mz_v, s_range=s_rng, ppm_tol=5.0)
            if f_hits:
                top_formulas.append(f_hits[0]["formula"])
                ppm_errs.append(f_hits[0]["ppm"])
            else:
                top_formulas.append("n/a")
                ppm_errs.append(np.nan)

    top_bio_candidates["putative_formula"] = top_formulas
    top_bio_candidates["formula_ppm_error"] = ppm_errs
    top_bio_candidates["isotope_evidence"] = iso_notes
    top_bio_out = top_bio_candidates.sort_values("max_omepr_pqn", ascending=False).head(35)
    top_bio_out.to_csv(TABLES_DIR / "table03_top_omeprazole_biotransformation_hits.csv")

    # Top endogenous hits: (a) strongest Time-changing Com20 bacterial metabolites + (b) Clean Trt/Int divergent outside suppression
    endo_clean_ann = res[
        (res["compound_name"].notna())
        & (res["feature_category"] == "Shared_Endogenous")
        & (~res["in_ion_suppression_window"])
        & (res["origin_category"] == "Endogenous/Medium")
    ].copy()
    endo_clean_ann["min_q_endo"] = np.minimum(endo_clean_ann["q_time"], np.minimum(endo_clean_ann["q_trt"], endo_clean_ann["q_int"]))
    endo_clean_ann.sort_values("min_q_endo").to_csv(TABLES_DIR / "table04_top_endogenous_hits.csv")

    # Save complete master statistical table (all 3,414 features)
    res.to_csv(TABLES_DIR / "table02_statistical_results_all.csv")

    # Save key numerical summary JSON for figures and notebooks
    summary_stats = {
        "n_injections_total": int(len(meta)),
        "n_bacterial_raw": 54,
        "n_failed_injections": len(failed_files),
        "failed_injection_files": failed_files,
        "n_bacterial_clean": len(bac_files),
        "n_medium_controls": len(medium_files),
        "n_blanks": len(blank_files),
        "n_pooled_qcs": len(pooled_qc_files),
        "n_qc_mix": len(qc_mix_files),
        "n_features_raw": int(len(X_raw)),
        "n_features_filtered": int(len(res)),
        "n_unique_iimn_entities": n_unique_entities,
        "n_annotated_raw": int(ann["id"].nunique()),
        "n_annotated_filtered": int(len(ann_kept)),
        "n_pqn_ref_features": int(len(ref_spectrum)),
        "run_order_spearman_dmso": round(float(r_run_dmso_design), 3),
        "run_order_p_dmso": float(p_run_dmso_design),
        "run_order_spearman_omepr": round(float(r_run_omepr_design), 3),
        "run_order_p_omepr": float(p_run_omepr_design),
        "tic_vs_pqn_shared_n": int(len(sh)),
        "tic_shared_sig_down": tic_fake_down,
        "pqn_shared_sig_down": pqn_down,
        "tic_shared_neg_int": tic_neg_int,
        "pqn_shared_neg_int": pqn_neg_int,
        "suppression_neg_int_in_window": a_supp,
        "suppression_neg_int_total": a_supp + b_supp,
        "suppression_pct_neg_int": round(100.0 * a_supp / (a_supp + b_supp), 1),
        "suppression_pct_features": round(100.0 * (a_supp + c_supp) / len(sh), 1),
        "suppression_fisher_or": round(float(or_supp), 1),
        "suppression_fisher_p": float(p_supp),
        "n_clean_endogenous": int(len(clean_endo_ids)),
        "bio_network_nodes": G_bio.number_of_nodes(),
        "bio_network_edges": G_bio.number_of_edges(),
        "endo_network_nodes": G_endo_conn.number_of_nodes(),
        "endo_network_edges": G_endo_conn.number_of_edges(),
        "endo_network_modules": len(communities),
        "endo_network_modularity": round(float(modularity_score), 3),
    }
    with open(PROC_DIR / "analysis_summary.json", "w") as f:
        json.dump(summary_stats, f, indent=2)

    print("=" * 72)
    print("PIPELINE FINISHED SUCCESSFULLY! Outputs saved to data/processed/ & results/")
    print("=" * 72)


if __name__ == "__main__":
    main()
