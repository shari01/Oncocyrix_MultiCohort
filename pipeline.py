# =========================
# 0. IMPORTS
# =========================
from typing import List, Optional
import argparse
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

# Use non-interactive backend (no Tkinter issues)
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path
import logging
import sys
from scipy import io as spio
from scipy import sparse as sp_sparse
import gzip
from typing import Dict, Optional

# Optional: mygene for Ensembl → gene symbol mapping
try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    MYGENE_AVAILABLE = False

# Optional: CellTypist for cell type annotation (supervised ML classifier)
try:
    import celltypist
    CELLTYPIST_AVAILABLE = True
except ImportError:
    CELLTYPIST_AVAILABLE = False

# Optional: BBKNN / other integration methods via scanpy.external
try:
    import scanpy.external as sce
    SC_EXT_AVAILABLE = True
except ImportError:
    SC_EXT_AVAILABLE = False

# Optional: rpy2 for DESeq2
try:
    import rpy2.robjects as ro
    from rpy2.rinterface_lib import callbacks
    RPY2_AVAILABLE = True
except ImportError:
    ro = None
    callbacks = None
    RPY2_AVAILABLE = False

# Optional gseapy for multi-database pathway enrichment (Enrichr)
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False


# =========================
# 1. CONFIG
# =========================

SINGLE_MODE = False      # True = one dataset only
MULTI_MODE = True        # True = many samples (e.g. pre/post, multi-patient)

# --- SINGLE MODE CONFIG ---
SINGLE_10X_DIR = Path(
    r"D:\AyassBio_Workspace_Downloads\SCANPY-SINGLECELL_PIPELINE\Cervical_cancer_sc_test\SINGLE_CELL_10X\ADC_(CANCER)\ADC_6"
)
SINGLE_SAMPLE_LABEL = "X_10"
SINGLE_GROUP_LABEL = "LUNG_CANCER"   # e.g. "CASE" / "CONTROL" or "TUMOR" / "NORMAL"

# --- MULTI MODE CONFIG ---
MULTI_BASE_DIR = Path(
    r"C:\Users\shery\Downloads\oncocyrix-multicohort\GSE208653_RAW"
)
# meta file for sample→group / CAR-T mapping
MULTI_META_FILENAME = "metadata.xlsx"

# SHORTER to avoid Excel 259-char path issues
OUTPUT_FOLDER_NAME = "SC_ANALYSIS_RESULTS"

# --- PATHWAY CLUSTERING CONTROL ---
DO_PATHWAY_CLUSTERING = True

# --- CAR-T ANALYSIS CONTROL ---
# If True: compute CAR-T / T-cell state signatures (TStemCM, TPEX, TEX, etc.)
# If False: skip all CAR-T scoring and plots (pipeline behaves as a generic scRNA-seq pipeline).
DO_CART_SCORING = False


# --- CLI OVERRIDES (optional) ---

parser = argparse.ArgumentParser(
    description="10x-only single-cell pipeline with per-sample + combined analysis (CAR-T aware)."
)
parser.add_argument("--mode", choices=["single", "multi"], help="Run single or multi mode.")
parser.add_argument("--single-10x-dir", type=str, help="Path to single 10x folder (single mode).")
parser.add_argument("--multi-base-dir", type=str, help="Base dir for multi mode (group/sample structure).")
parser.add_argument("--out-name", type=str, help="Custom combined output dir name.")
parser.add_argument("--no-pathway-clustering", action="store_true", help="Disable pathway clustering (cluster/celltype DEG enrichment + pseudobulk GSEA).")
parser.add_argument("--do-cart", action="store_true",help="Enable CAR-T signature scoring and CAR-T state analysis.")
                  
                   

args, _ = parser.parse_known_args()

if args.do_cart:
    DO_CART_SCORING = True

if args.mode is not None:
    if args.mode == "single":
        SINGLE_MODE, MULTI_MODE = True, False
    else:
        SINGLE_MODE, MULTI_MODE = False, True

if args.single_10x_dir is not None:
    SINGLE_10X_DIR = Path(args.single_10x_dir)

if args.multi_base_dir is not None:
    MULTI_BASE_DIR = Path(args.multi_base_dir)

if args.out_name is not None:
    OUTPUT_FOLDER_NAME = args.out_name

if args.no_pathway_clustering:
    DO_PATHWAY_CLUSTERING = False

if SINGLE_MODE == MULTI_MODE:
    raise ValueError("Set exactly one of SINGLE_MODE or MULTI_MODE to True.")

if SINGLE_MODE:
    IN_DIR = SINGLE_10X_DIR
    DATA_PATH = SINGLE_10X_DIR
else:
    IN_DIR = MULTI_BASE_DIR
    DATA_PATH = MULTI_BASE_DIR

COMBINED_OUT_DIR = IN_DIR / OUTPUT_FOLDER_NAME
COMBINED_OUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = COMBINED_OUT_DIR / "pipeline.log"

# ============= LOGGING =============
for h in logging.root.handlers[:]:
    logging.root.removeHandler(h)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s\t[%(levelname)s]\t%(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)

logger = logging.getLogger()  # root logger

logger.info(f"Logging to file: {LOG_FILE}")
logger.info(f"SINGLE_MODE: {SINGLE_MODE}")
logger.info(f"MULTI_MODE: {MULTI_MODE}")
logger.info(f"IN_DIR: {IN_DIR}")
logger.info(f"DATA_PATH: {DATA_PATH}")
logger.info(f"COMBINED_OUT_DIR: {COMBINED_OUT_DIR}")
logger.info(f"mygene available: {MYGENE_AVAILABLE}")
logger.info(f"celltypist available: {CELLTYPIST_AVAILABLE}")
logger.info(f"scanpy.external available: {SC_EXT_AVAILABLE}")
logger.info(f"rpy2 available: {RPY2_AVAILABLE}")
logger.info(f"gseapy available: {GSEAPY_AVAILABLE}")
logger.info(f"DO_PATHWAY_CLUSTERING: {DO_PATHWAY_CLUSTERING}")
logger.info(f"DO_CART_SCORING: {DO_CART_SCORING}")

# Redirect R console output to Python logger
if RPY2_AVAILABLE and callbacks is not None:
    def _r_console_write(x):
        x = x.strip()
        if x:
            logger.info(f"[R] {x}")
    callbacks.consolewrite_print = _r_console_write
    callbacks.consolewrite_warnerror = _r_console_write

# Global plotting aesthetics
sc.settings.verbosity = 3
sc.settings.figdir = COMBINED_OUT_DIR
sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings.set_figure_params(
    dpi=140,
    dpi_save=300,
    facecolor="white",
    format="png",
    figsize=(10, 8),
)

plt.rcParams["figure.figsize"] = (10, 8)
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 12

logger.info(f"Scanpy version: {sc.__version__}")


# =========================
# 1b. CAR-T SIGNATURES & SCORING
# =========================

# These signatures are intentionally a bit redundant/overlapping.
# Idea: capture pre-CAR T cell states (naive / Tcm / Tem), stem-cell memory,
# exhaustion, effector/differentiation, and context (TME, proliferation, etc.)
CART_SIGNATURES = {
    # --- Pre-CAR / baseline T cell states ---

    # Classic naive CD4/CD8 (CCR7+ SELL+ IL7R+)
    "Naive_T": [
        "CCR7", "SELL", "IL7R", "TCF7", "LEF1", "MAL",
    ],

    # Central memory (Tcm) – IL7R+ CCR7+ co-stimulatory
    "Central_Memory_Tcm": [
        "CCR7", "SELL", "IL7R", "CD27", "CD28", "ICOS",
    ],

    # Effector memory (Tem) – CCR7- but still IL7R+, more effector skew
    "Effector_Memory_Tem": [
        "IL7R", "CXCR3", "CXCR4", "GZMK", "CCL5",
    ],

    # Stem-cell memory / central memory–like (good prognosis)
    "Memory_TStemCM": [
        "TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "IL2RB",
    ],

    # Tissue-resident memory (Trm-like) – more relevant in tissues/TME
    "Resident_Memory_Trm": [
        "ITGAE",  # CD103
        "ITGA1",
        "CXCR6",
        "CD69",
        "NR4A1",
    ],

    # Memory-like CAR T associated with durable responses
    # (less-differentiated, TCF7/LEF1/BCL2/BACH2/KLF2+)
    "Memory_like_CAR_T_good": [
        "TCF7", "LEF1", "CCR7", "SELL", "IL7R",
        "CD27", "CD28",
        "BCL2",
        "BACH2",
        "KLF2",
        "CXCR3",
    ],

    # --- Exhaustion / activation continuum ---

    # Precursor exhausted (TPEX) – TCF7+ TOX+ PD1+ with memory-like features
    "TPEX": [
        "TCF7",
        "TOX",
        "PDCD1",
        "CXCR5",
        "SLAMF6",
        "BCL6",
    ],

    # Terminally exhausted (TEX) – multiple inhibitory receptors + exhaustion TFs
    "TEX": [
        "TOX", "TOX2",
        "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1",
        "EOMES",
        "BATF",
        "IRF4",
        "NR4A1", "NR4A2",
        "PRDM1",  # BLIMP1
    ],

    # General T cell activation (CD25, HLA, NFAT target genes)
    "Activation_Tcell": [
        "IL2RA",      # CD25
        "TNFRSF4",    # OX40
        "TNFRSF18",   # GITR
        "HLA-DRA", "HLA-DRB1", "HLA-DPA1",
        "CD40LG",
        "CD69",
        "CD38",
    ],

    # Costimulatory / good “fitness” signal
    "Costimulatory_Signaling": [
        "CD28", "ICOS", "TNFRSF9", "TNFSF9", "CD27", "CD40LG",
    ],

    # Inhibitory checkpoint module
    "Inhibitory_Checkpoints": [
        "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "BTLA",
        "CD244", "CD160",
    ],

    # Th1 / inflammatory cytokine module (linked to CAR-T effector function & CRS)
    "Th1_Inflammatory_Cytokines": [
        "IFNG",
        "TNF",
        "IL2",
        "CXCL9",
        "CXCL10",
        "CCL3", "CCL4",
        "STAT1",
        "TBX21",  # T-bet
    ],

    # Effector / cytotoxic
    "Effector_TEFF": [
        "GZMB", "GZMH", "GZMK", "PRF1", "IFNG", "GNLY", "NKG7",
    ],

    # Terminal differentiation (KLRG1+ CX3CR1+)
    "Terminal_diff": [
        "KLRG1", "CX3CR1", "GZMB", "PRF1",
    ],

    # NK-like dysfunctional CAR T (CD8→NK-like transition, poor outcomes)
    "NK_like_Dysfunctional_CAR_T": [
        # NK receptors
        "KLRB1",
        "KLRC1", "KLRC2", "KLRC3",
        "KLRD1",
        # cytotoxic granules
        "GNLY",
        "NKG7",
        # NK/T-cell surface markers & dysfunction genes
        "FCGR3A",   # CD16
        "LAYN",
        "CD9",
        "PHLDA1",
        "TNFRSF9",
        "RGS16",
        "SRGAP3",
        # exhaustion-associated TFs regulating NK-like fate
        "SOX4",
        "ID3",
    ],

    # Treg-like CAR T / suppressive CD4+ T cells (linked with poor response)
    "Treg_like_CAR_T": [
        "FOXP3",
        "IKZF2",    # Helios
        "IL2RA",    # CD25
        "CTLA4",
        "TIGIT",
        "LAG3",
        "CCR8",
        "TNFRSF18", # GITR
    ],

    # --- Metabolism / proliferation / TME context ---

    # Metabolic fitness (PGC1a axis)
    "Metabolic_fitness": [
        "PPARGC1A",  # PGC1A
        "PGC1A",     # some datasets use this alias
        "TFAM",
        "MYC",
        "SLC2A1",    # GLUT1
    ],

    # Hypoxia / glycolytic stress (often in exhausted CAR T in TME)
    "Hypoxia_Glycolysis_Stress": [
        "HIF1A",
        "SLC2A1",
        "LDHA",
        "ENO1",
        "GAPDH",
        "PGK1",
        "VEGFA",
    ],

    # Proliferation (Ki67+ cycling CAR T cells)
    "Proliferation": [
        "MKI67", "TOP2A", "PCNA", "CCNB1", "CDC20", "BIRC5",
    ],

    # More refined G2/M cell-cycle module often seen in cycling CAR T clusters
    "Cell_cycle_G2M": [
        "MKI67",
        "TOP2A",
        "CCNB1",
        "CDC6",
        "MCM2", "MCM4",
        "UBE2C",
        "AURKB",
    ],

    # Type I IFN / ISG response – sometimes enriched in dysfunctional CAR T
    "Interferon_ISG_Response": [
        "ISG15",
        "IFIT1", "IFIT3",
        "MX1",
        "OAS1",
        "IRF7",
        "STAT1",
    ],

    # Immunosuppressive / TME (myeloid + Treg + stromal)
    "Immunosuppressive_TME": [
        "IDO1", "ARG1", "IL10", "TGFB1", "VEGFA",
        "FOXP3",    # Treg
        "S100A8", "S100A9",  # MDSC / neutrophil-like
        "CSF3R",
    ],
}

# ============================================================
# === CAR-T SIGNATURE SCORING FUNCTIONS (DIRECTLY EMBEDDED) ===
# ============================================================


def compute_cart_scores(adata):
    """
    Compute CAR-T signature scores for every signature in CART_SIGNATURES.
    Adds columns to adata.obs named CART_<SignatureName>_score.
    """
    print("[CART] Computing CAR-T signature scores ...")

    # Map uppercase → real var name (case-insensitive match)
    var_map = {g.upper(): g for g in adata.var_names}

    for sig_name, genes in CART_SIGNATURES.items():
        matched = []
        for g in genes:
            up = g.upper()
            if up in var_map:
                matched.append(var_map[up])

        if len(matched) < 3:
            print(f"[CART] Signature '{sig_name}' → too few genes in dataset, skipping.")
            continue

        score_name = f"CART_{sig_name}_score"
        print(f"[CART] Scoring {score_name} (n={len(matched)})")

        sc.tl.score_genes(
            adata,
            gene_list=matched,
            score_name=score_name,
            use_raw=True
        )

    print("[CART] Finished computing CAR-T signatures.")


def plot_cart_state_proportions_by_patient_phase(
    adata,
    out_dir: Path,
    patient_col="patient_id",
    phase_col="cart_phase",
    state_col="CART_State_v2",
    min_cells=20,
):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    for col in [patient_col, phase_col, state_col]:
        if col not in adata.obs.columns:
            logger.info(f"[CART-PLOT] Missing {col} -> skip.")
            return

    df = adata.obs[[patient_col, phase_col, state_col]].astype(str).copy()
    # filter tiny groups
    keep = df.value_counts([patient_col, phase_col]).reset_index(name="n")
    keep = keep[keep["n"] >= min_cells][[patient_col, phase_col]]
    df = df.merge(keep, on=[patient_col, phase_col], how="inner")

    tab = df.value_counts([patient_col, phase_col, state_col]).reset_index(name="n_cells")
    tot = tab.groupby([patient_col, phase_col])["n_cells"].sum().reset_index(name="total")
    tab = tab.merge(tot, on=[patient_col, phase_col])
    tab["fraction"] = tab["n_cells"] / tab["total"]

    tab.to_csv(out_dir / "CART_state_fraction_patient_phase.csv", index=False)

    # plot (one panel per patient)
    patients = sorted(tab[patient_col].unique())
    for p in patients:
        sub = tab[tab[patient_col] == p].copy()
        pivot = sub.pivot(index=phase_col, columns=state_col, values="fraction").fillna(0)

        ax = pivot.plot(kind="bar", stacked=True, figsize=(12, 5))
        ax.set_ylabel("Fraction of cells")
        ax.set_title(f"{p}: CAR-T state composition by phase")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(out_dir / f"{p}_CART_state_composition_stacked.png", dpi=300)
        plt.close()


def _safe_overlap(adata: ad.AnnData, genes: list[str]) -> list[str]:
    """
    Case-insensitive intersection of a gene list with adata.var_names.
    Returns the original var_names entry (not uppercased).
    """
    var_upper = pd.Index([g.upper() for g in adata.var_names])
    wanted = [g.upper() for g in genes]
    mask = var_upper.isin(wanted)
    if not mask.any():
        return []
    return list(adata.var_names[mask])


def refine_cart_states(
    adata: ad.AnnData,
    celltype_col: str | None = "celltype",
    out_col: str = "CART_State_v2",
) -> None:
    """
    Build a more detailed CAR-T classifier to reduce 'Other'.

    Uses:
      - CART_*_score columns from compute_cart_scores
      - Optional celltype column to flag non-T cells
    Writes:
      adata.obs[out_col]
    """

    logger.info(f"[CART-REFINE] Building refined CAR-T states into obs['{out_col}'].")

    # 1) Pull all scores safely (fallback to zeros if missing)
    def _get_score(name: str):
        col = f"CART_{name}_score"
        if col in adata.obs.columns:
            return adata.obs[col]
        return pd.Series(0.0, index=adata.obs_names)

    s_naive   = _get_score("Naive_T")
    s_tcm     = _get_score("Central_Memory_Tcm")
    s_tem     = _get_score("Effector_Memory_Tem")
    s_mem     = _get_score("Memory_TStemCM")
    s_tpex    = _get_score("TPEX")
    s_tex     = _get_score("TEX")
    s_teff    = _get_score("Effector_TEFF")
    s_term    = _get_score("Terminal_diff")
    s_prolif  = _get_score("Proliferation")

    # 2) Compute thresholds (slightly softer than 75th percentile for some)
    def q(series, qv):
        return float(series.quantile(qv)) if len(series) > 0 else 0.0

    thr = {
        "naive":   q(s_naive, 0.70),
        "tcm":     q(s_tcm,   0.70),
        "tem":     q(s_tem,   0.70),
        "mem":     q(s_mem,   0.75),
        "tpex":    q(s_tpex,  0.75),
        "tex":     q(s_tex,   0.75),
        "teff":    q(s_teff,  0.70),
        "term":    q(s_term,  0.75),
        "prolif":  q(s_prolif, 0.80),
    }

    logger.info(f"[CART-REFINE] Thresholds: {thr}")

    # 3) Optional: detect non-T cells from celltype
    is_t_like = pd.Series(True, index=adata.obs_names)
    if celltype_col is not None and celltype_col in adata.obs.columns:
        ct = adata.obs[celltype_col].astype(str).str.upper()
        t_keywords = ["T", "CD4", "CD8", "NK", "CAR"]
        is_t_like = ct.apply(
            lambda x: any(kw in x for kw in t_keywords)
        )

    labels: list[str] = []

    for i in range(adata.n_obs):
        if not is_t_like.iloc[i]:
            labels.append("Non_T_or_other_immune")
            continue

        v = {
            "naive":   s_naive.iloc[i],
            "tcm":     s_tcm.iloc[i],
            "tem":     s_tem.iloc[i],
            "mem":     s_mem.iloc[i],
            "tpex":    s_tpex.iloc[i],
            "tex":     s_tex.iloc[i],
            "teff":    s_teff.iloc[i],
            "term":    s_term.iloc[i],
            "prolif":  s_prolif.iloc[i],
        }

        state = "Unclassified_T"

        # Priority order (from most specific / clinically relevant)
        if v["tex"]   >= thr["tex"]:
            state = "TEX_terminal"
        elif v["tpex"] >= thr["tpex"] and v["tex"] < thr["tex"]:
            state = "TPEX_like"
        elif v["mem"] >= thr["mem"] and v["tex"] < thr["tex"] and v["term"] < thr["term"]:
            state = "TStemCM_like"
        elif v["term"] >= thr["term"]:
            state = "Terminal_diff"
        elif v["teff"] >= thr["teff"]:
            state = "Effector_TEFF"
        elif v["prolif"] >= thr["prolif"]:
            state = "Proliferating_T"
        elif (v["naive"] >= thr["naive"]) or (v["tcm"] >= thr["tcm"]):
            state = "Naive_or_Tcm"
        elif v["tem"] >= thr["tem"]:
            state = "Tem_like"
        else:
            # low-signal T cell
            state = "Low_signal_T"

        labels.append(state)

    adata.obs[out_col] = pd.Categorical(labels)
    logger.info(f"[CART-REFINE] Added refined CAR-T labels: {adata.obs[out_col].value_counts().to_dict()}")


def export_cart_signature_gene_evidence(
    adata: ad.AnnData,
    cart_dir: Path,
    analysis_name: str,
    state_col: str = "CART_State_v2",
    use_layer: str | None = None,   # None -> use adata.raw if exists else adata.X
    expr_threshold: float = 0.0,
):
    """
    Writes:
      A) <analysis>_CART_signature_gene_evidence_overall.csv
         signature × gene: mean_expr + %cells expressing

      B) <analysis>_CART_signature_gene_evidence_by_state.csv
         signature × state × gene: mean_expr + %cells expressing
    """
    cart_dir = Path(cart_dir)
    cart_dir.mkdir(parents=True, exist_ok=True)

    # --- pick matrix ---
    if use_layer is not None and use_layer in adata.layers:
        X = adata.layers[use_layer]
        gene_names = adata.var_names
    elif adata.raw is not None:
        X = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        X = adata.X
        gene_names = adata.var_names

    # case-insensitive mapping
    var_map = {str(g).upper(): str(g) for g in gene_names}
    gene_index = pd.Index(gene_names)

    if state_col in adata.obs.columns:
        states = adata.obs[state_col].astype(str)
    else:
        states = pd.Series(["ALL"] * adata.n_obs, index=adata.obs_names)

    rows_overall = []
    rows_by_state = []

    def _to_1d(arr):
        if sp_sparse.issparse(arr):
            return arr.toarray().ravel()
        return np.array(arr).ravel()

    for sig_name, genes in CART_SIGNATURES.items():
        matched = []
        for g in genes:
            up = str(g).upper()
            if up in var_map:
                matched.append(var_map[up])

        if not matched:
            continue

        idx = gene_index.get_indexer(matched)
        idx = idx[idx >= 0]
        if len(idx) == 0:
            continue

        Xsub = X[:, idx]
        if sp_sparse.issparse(Xsub):
            Xsub = Xsub.tocsr()

        # overall
        for j, gname in enumerate(matched):
            data = _to_1d(Xsub[:, j])
            rows_overall.append({
                "analysis": analysis_name,
                "signature": sig_name,
                "gene": gname,
                "mean_expr": float(data.mean()),
                "pct_cells_expressing": float((data > expr_threshold).mean() * 100.0),
            })

        # by-state
        for st in sorted(states.unique()):
            mask = (states == st).values
            if mask.sum() < 5:
                continue
            Xst = Xsub[mask, :]
            for j, gname in enumerate(matched):
                data = _to_1d(Xst[:, j])
                rows_by_state.append({
                    "analysis": analysis_name,
                    "signature": sig_name,
                    "state": st,
                    "gene": gname,
                    "mean_expr": float(data.mean()),
                    "pct_cells_expressing": float((data > expr_threshold).mean() * 100.0),
                    "n_cells_state": int(mask.sum()),
                })

    df_overall = pd.DataFrame(rows_overall)
    df_by_state = pd.DataFrame(rows_by_state)

    out1 = cart_dir / f"{analysis_name}_CART_signature_gene_evidence_overall.csv"
    out2 = cart_dir / f"{analysis_name}_CART_signature_gene_evidence_by_state.csv"

    df_overall.to_csv(out1, index=False)
    df_by_state.to_csv(out2, index=False)

    logger.info(f"[CART] Wrote signature gene evidence overall: {out1}")
    logger.info(f"[CART] Wrote signature gene evidence by state: {out2}")


def export_cart_state_markers(
    adata: ad.AnnData,
    cart_dir: Path,
    analysis_name: str,
    state_col: str = "CART_State_v2",
    n_genes: int = 100,
):
    """
    Writes:
      <analysis>_CART_state_markers_ranked.csv  (Scanpy wilcoxon DE between CART states)
    """
    if state_col not in adata.obs.columns:
        logger.info(f"[CART] No '{state_col}' found -> skipping CART state marker DE.")
        return

    cart_dir = Path(cart_dir)
    cart_dir.mkdir(parents=True, exist_ok=True)

    if not pd.api.types.is_categorical_dtype(adata.obs[state_col]):
        adata.obs[state_col] = adata.obs[state_col].astype("category")

    sc.tl.rank_genes_groups(
        adata,
        groupby=state_col,
        method="wilcoxon",
        n_genes=n_genes,
    )

    df = sc.get.rank_genes_groups_df(adata, None)
    out = cart_dir / f"{analysis_name}_CART_state_markers_ranked.csv"
    df.to_csv(out, index=False)
    logger.info(f"[CART] Wrote CART state marker genes: {out}")


def export_cart_patient_phase_gene_summary(
    adata: ad.AnnData,
    cart_dir: Path,
    analysis_name: str,
    patient_col: str = "patient_id",
    phase_col: str = "cart_phase",          # meta column like Pre/Post
    fallback_phase_col: str = "group",      # if cart_phase missing, use group
    state_col: str = "CART_State_v2",
    use_layer: str | None = None,           # None -> adata.raw if exists else adata.X
    expr_threshold: float = 0.0,
    min_cells: int = 10,
):
    # ---- DEBUG: sanity checks ----
    print("\n[DEBUG export_cart_patient_phase_gene_summary]")
    print("n_obs:", adata.n_obs)
    print("Available obs columns:",
          [c for c in [patient_col, phase_col, fallback_phase_col, state_col]
           if c in adata.obs.columns])

    if patient_col in adata.obs:
        print("patient_id unique (first 10):",
              adata.obs[patient_col].astype(str).unique()[:10])
    else:
        print("MISSING patient_col:", patient_col)

    if phase_col in adata.obs:
        print("cart_phase unique:",
              adata.obs[phase_col].astype(str).unique())
    else:
        print("cart_phase missing, fallback group unique:",
              adata.obs[fallback_phase_col].astype(str).unique()
              if fallback_phase_col in adata.obs else "MISSING")

    if state_col in adata.obs:
        print("CART_State_v2 unique:",
              adata.obs[state_col].astype(str).unique())
    else:
        print("MISSING state_col:", state_col)
    # ---- END DEBUG ----

    """
    Writes:
      1) <analysis>_CART_patient_phase_gene_summary.csv
         patient × phase × gene

      2) <analysis>_CART_patient_phase_state_gene_summary.csv
         patient × phase × CART_State_v2 × gene   <-- THIS is your main goal

      3) <analysis>_CART_patient_phase_state_gene_delta_post_minus_pre.csv (optional)
         patient × state × gene, Post-Pre deltas (only if phase labels include Pre & Post)
    """
    cart_dir = Path(cart_dir)
    cart_dir.mkdir(parents=True, exist_ok=True)

    # --- pick matrix ---
    if use_layer is not None and use_layer in adata.layers:
        X = adata.layers[use_layer]
        gene_names = adata.var_names
    elif adata.raw is not None:
        X = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        X = adata.X
        gene_names = adata.var_names

    if patient_col not in adata.obs.columns:
        logger.warning(f"[CART-PATIENT] '{patient_col}' not in obs -> skipping patient-phase export.")
        return

    if phase_col in adata.obs.columns:
        phase_use = phase_col
    elif fallback_phase_col in adata.obs.columns:
        phase_use = fallback_phase_col
        logger.warning(f"[CART-PATIENT] '{phase_col}' missing; using '{fallback_phase_col}' as phase.")
    else:
        logger.warning(f"[CART-PATIENT] Neither '{phase_col}' nor '{fallback_phase_col}' found -> skipping.")
        return

    # Build CART gene list (all signature genes, deduped, matched)
    var_map = {str(g).upper(): str(g) for g in gene_names}
    all_cart_genes = []
    for _, genes in CART_SIGNATURES.items():
        for g in genes:
            up = str(g).upper()
            if up in var_map:
                all_cart_genes.append(var_map[up])
    all_cart_genes = sorted(set(all_cart_genes))
    if not all_cart_genes:
        logger.warning("[CART-PATIENT] No CART signature genes matched -> skipping.")
        return

    gene_index = pd.Index(gene_names)
    idx = gene_index.get_indexer(all_cart_genes)
    idx = idx[idx >= 0]
    if len(idx) == 0:
        logger.warning("[CART-PATIENT] CART gene indices empty after matching -> skipping.")
        return

    Xsub = X[:, idx]
    if sp_sparse.issparse(Xsub):
        Xsub = Xsub.tocsr()

    patients = adata.obs[patient_col].astype(str)
    phases = adata.obs[phase_use].astype(str)

    if state_col in adata.obs.columns:
        states = adata.obs[state_col].astype(str)
    else:
        states = pd.Series(["ALL"] * adata.n_obs, index=adata.obs_names)

    def _to_1d(arr):
        if sp_sparse.issparse(arr):
            return arr.toarray().ravel()
        return np.array(arr).ravel()

    # -------- 1) patient × phase × gene --------
    rows_pp = []
    grouped_pp = pd.DataFrame({
        "patient": patients.values,
        "phase": phases.values,
    }).groupby(["patient", "phase"]).indices

    for (p, ph), obs_idx in grouped_pp.items():
        if len(obs_idx) < min_cells:
            continue
        Xg = Xsub[obs_idx, :]
        for j, gname in enumerate(all_cart_genes):
            arr = _to_1d(Xg[:, j])
            rows_pp.append({
                "analysis": analysis_name,
                "patient_id": p,
                "phase": ph,
                "gene": gname,
                "mean_expr": float(arr.mean()),
                "pct_cells_expressing": float((arr > expr_threshold).mean() * 100.0),
                "n_cells": int(len(obs_idx)),
            })

    df_pp = pd.DataFrame(rows_pp)
    out_pp = cart_dir / f"{analysis_name}_CART_patient_phase_gene_summary.csv"
    df_pp.to_csv(out_pp, index=False)
    logger.info(f"[CART-PATIENT] Wrote: {out_pp}")

    # -------- 2) patient × phase × state × gene --------
    rows_pps = []
    grouped_pps = pd.DataFrame({
        "patient": patients.values,
        "phase": phases.values,
        "state": states.values,
    }).groupby(["patient", "phase", "state"]).indices

    for (p, ph, st), obs_idx in grouped_pps.items():
        if len(obs_idx) < min_cells:
            continue
        Xg = Xsub[obs_idx, :]
        for j, gname in enumerate(all_cart_genes):
            arr = _to_1d(Xg[:, j])
            rows_pps.append({
                "analysis": analysis_name,
                "patient_id": p,
                "phase": ph,
                "CART_State_v2": st,
                "gene": gname,
                "mean_expr": float(arr.mean()),
                "pct_cells_expressing": float((arr > expr_threshold).mean() * 100.0),
                "n_cells": int(len(obs_idx)),
            })

    df_pps = pd.DataFrame(rows_pps)
    out_pps = cart_dir / f"{analysis_name}_CART_patient_phase_state_gene_summary.csv"
    df_pps.to_csv(out_pps, index=False)
    logger.info(f"[CART-PATIENT] Wrote: {out_pps}")

    # -------- 3) Post - Pre delta (only if phases map to Pre & Post) --------
    if df_pps.empty:
        return

    tmp = df_pps.copy()
    tmp["phase_norm"] = tmp["phase"].str.strip().str.lower()

    pre_names = {"pre", "baseline", "before"}
    post_names = {"post", "after"}

    def map_phase(x: str):
        if x in pre_names:
            return "Pre"
        if x in post_names:
            return "Post"
        return x

    tmp["phase_bucket"] = tmp["phase_norm"].map(map_phase)

    if {"Pre", "Post"}.issubset(set(tmp["phase_bucket"].unique())):
        wide_mean = tmp.pivot_table(
            index=["patient_id", "CART_State_v2", "gene"],
            columns="phase_bucket",
            values="mean_expr",
            aggfunc="mean",
        )
        wide_pct = tmp.pivot_table(
            index=["patient_id", "CART_State_v2", "gene"],
            columns="phase_bucket",
            values="pct_cells_expressing",
            aggfunc="mean",
        )

        common_idx = wide_mean.dropna(subset=["Pre", "Post"]).index
        if len(common_idx) == 0:
            logger.info("[CART-PATIENT] Pre/Post detected but no common rows -> delta skipped.")
            return

        wide_mean = wide_mean.loc[common_idx]
        wide_pct = wide_pct.loc[common_idx]

        df_delta = pd.DataFrame({
            "patient_id": [i[0] for i in common_idx],
            "CART_State_v2": [i[1] for i in common_idx],
            "gene": [i[2] for i in common_idx],
            "mean_expr_Pre": wide_mean["Pre"].values,
            "mean_expr_Post": wide_mean["Post"].values,
            "mean_expr_delta_Post_minus_Pre": (wide_mean["Post"] - wide_mean["Pre"]).values,
            "pct_expr_Pre": wide_pct["Pre"].values,
            "pct_expr_Post": wide_pct["Post"].values,
            "pct_expr_delta_Post_minus_Pre": (wide_pct["Post"] - wide_pct["Pre"]).values,
        })

        out_delta = cart_dir / f"{analysis_name}_CART_patient_phase_state_gene_delta_post_minus_pre.csv"
        df_delta.to_csv(out_delta, index=False)
        logger.info(f"[CART-PATIENT] Wrote: {out_delta}")
    else:
        logger.info("[CART-PATIENT] Could not detect both Pre and Post phase labels; delta file not written.")

# =========================
# 2. LOADER FUNCTIONS (10x ONLY)
# =========================

def _find_first_matching(dir_path: Path, patterns) -> Path | None:
    for pat in patterns:
        matches = sorted(dir_path.glob(pat))
        if matches:
            return matches[0]
    return None


def _mmread_auto(path: Path):
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return spio.mmread(f)
    else:
        return spio.mmread(str(path))


def load_10x_feature_barcode_matrix(tenx_dir: Path) -> ad.AnnData:
    """
    Single 10x folder:
      matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv/genes.tsv[.gz]
    """
    tenx_dir = Path(tenx_dir)
    if not tenx_dir.exists():
        raise FileNotFoundError(f"10X folder not found: {tenx_dir}")

    logger.info(f"Loading 10X feature-barcode matrix from: {tenx_dir}")

    matrix_path = _find_first_matching(
        tenx_dir,
        ["matrix.mtx", "matrix.mtx.gz", "*.matrix.mtx", "*.matrix.mtx.gz", "*.mtx", "*.mtx.gz"],
    )
    if matrix_path is None:
        raise FileNotFoundError(
            f"No matrix.mtx[.gz] file found in {tenx_dir}."
        )

    barcodes_path = _find_first_matching(
        tenx_dir,
        [
            "barcodes.tsv",
            "barcodes.tsv.gz",
            "*barcodes.tsv",
            "*barcodes.tsv.gz",
            "barcode.tsv",
            "barcode.tsv.gz",
            "*barcode.tsv",
            "*barcode.tsv.gz",
        ],
    )
    if barcodes_path is None:
        raise FileNotFoundError(
            f"No barcodes.tsv[.gz] file found in {tenx_dir}."
        )

    features_path = _find_first_matching(
        tenx_dir,
        [
            "features.tsv",
            "features.tsv.gz",
            "*features.tsv",
            "*features.tsv.gz",
            "genes.tsv",
            "genes.tsv.gz",
            "*genes.tsv",
            "*genes.tsv.gz",
        ],
    )
    if features_path is None:
        raise FileNotFoundError(
            f"No features.tsv[.gz] or genes.tsv[.gz] file found in {tenx_dir}."
        )

    logger.info(f"matrix:   {matrix_path.name}")
    logger.info(f"barcodes: {barcodes_path.name}")
    logger.info(f"features: {features_path.name}")

    M = _mmread_auto(matrix_path)
    if not sp_sparse.issparse(M):
        M = sp_sparse.coo_matrix(M)
    M = M.tocsr()
    X = M.T  # cells x genes

    barcodes_df = pd.read_csv(barcodes_path, sep="\t", header=None, compression="infer")
    barcodes = barcodes_df.iloc[:, 0].astype(str).values

    feat_df = pd.read_csv(features_path, sep="\t", header=None, compression="infer")
    ncols = feat_df.shape[1]
    colnames = []
    if ncols >= 1:
        colnames.append("feature_id")
    if ncols >= 2:
        colnames.append("feature_name")
    if ncols >= 3:
        colnames.append("feature_type")
    while len(colnames) < ncols:
        colnames.append(f"extra_{len(colnames)}")
    feat_df.columns = colnames

    gene_ids = feat_df["feature_id"].astype(str).values
    if "feature_name" in feat_df.columns:
        gene_names = feat_df["feature_name"].astype(str).values
    else:
        gene_names = gene_ids

    adata_ = ad.AnnData(X=X)
    adata_.obs_names = barcodes
    adata_.obs["barcode"] = barcodes
    adata_.var_names = gene_names
    adata_.var["feature_id"] = gene_ids
    adata_.var["gene_symbol"] = gene_names
    if "feature_type" in feat_df.columns:
        adata_.var["feature_type"] = feat_df["feature_type"].astype(str).values
    adata_.var_names_make_unique()
    logger.info(f"Loaded AnnData from 10x: {adata_}")
    return adata_


# =========================
# 3. GENE NAME HANDLING
# =========================

def looks_like_ensembl(ids, prefix="ENSG"):
    ids = list(ids)
    if len(ids) == 0:
        return False
    sample = ids[: min(50, len(ids))]
    flags = [str(x).upper().startswith(prefix) for x in sample]
    return np.mean(flags) > 0.6


def update_gene_names(adata_in: ad.AnnData) -> ad.AnnData:
    candidate_cols = [
        "gene_symbol", "symbol", "GeneSymbol", "SYMBOL",
        "gene_name", "GeneName", "name"
    ]
    for col in candidate_cols:
        if col in adata_in.var.columns:
            col_vals = adata_in.var[col].astype(str)
            non_na_ratio = (col_vals != "nan").mean()
            if non_na_ratio > 0.5:
                logger.info(f"Using adata.var['{col}'] as gene names.")
                new_names = []
                for old, new in zip(adata_in.var_names, col_vals):
                    if new != "nan" and new is not None and len(new) > 0:
                        new_names.append(new)
                    else:
                        new_names.append(old)
                adata_in.var_names = new_names
                adata_in.var_names_make_unique()
                return adata_in

    if MYGENE_AVAILABLE and looks_like_ensembl(adata_in.var_names, prefix="ENSG"):
        logger.info("var_names look like Ensembl – mapping via mygene.")
        mg = mygene.MyGeneInfo()
        gene_ids = adata_in.var_names.tolist()
        chunk_size = 1000
        all_results = []

        for i in range(0, len(gene_ids), chunk_size):
            try:
                chunk = gene_ids[i:i+chunk_size]
                res = mg.querymany(
                    chunk,
                    scopes="ensembl.gene",
                    fields="symbol,name",
                    species="human",
                    as_dataframe=True,
                    df_index=True,
                )
                all_results.append(res)
            except Exception as e:
                logger.warning(f"mygene chunk {i} failed: {e}")
                continue

        if len(all_results) > 0:
            mapping_df = pd.concat(all_results, axis=0)
            mapping_df = mapping_df[~mapping_df.index.duplicated(keep="first")]
            symbols = []
            names = []
            for gid in adata_in.var_names:
                if gid in mapping_df.index:
                    row = mapping_df.loc[gid]
                    symbols.append(row.get("symbol", np.nan))
                    names.append(row.get("name", np.nan))
                else:
                    symbols.append(np.nan)
                    names.append(np.nan)
            adata_in.var["mapped_symbol"] = pd.Series(symbols, index=adata_in.var.index, dtype="string")
            adata_in.var["mapped_name"] = pd.Series(names, index=adata_in.var.index, dtype="string")
            mapped_sym = adata_in.var["mapped_symbol"]
            valid_mask = mapped_sym.notna() & (mapped_sym != "") & (mapped_sym != "nan")
            mapped_count = int(valid_mask.sum())
            logger.info(f"Mapped Ensembl → symbol for ~{mapped_count} genes.")
            new_names = []
            for old, sym, ok in zip(adata_in.var_names, mapped_sym, valid_mask):
                if ok:
                    new_names.append(str(sym))
                else:
                    new_names.append(old)
            adata_in.var_names = new_names
            adata_in.var_names_make_unique()
        else:
            logger.warning("No successful mygene mapping. Keeping original names.")
    else:
        if not MYGENE_AVAILABLE:
            logger.info("mygene not installed; skipping Ensembl mapping.")
        else:
            logger.info("var_names do not look like Ensembl; skipping mapping.")
    return adata_in


# =========================
# 4. SIMPLE QC FILTER (per-sample pre-combine)
# =========================

def simple_qc_filter(
    adata_in: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_mt_pct: float = 15.0,
) -> ad.AnnData:
    """
    Basic QC filter: per-sample or pre-combine (no plots).
    """
    adata = adata_in.copy()
    logger.info(f"[simple_qc_filter] Input shape: {adata.n_obs} cells x {adata.n_vars} genes")

    if "counts" not in adata.layers:
        if sp_sparse.issparse(adata.X):
            adata.layers["counts"] = adata.X.copy()
        else:
            adata.layers["counts"] = np.array(adata.X)

    adata = update_gene_names(adata)
    if getattr(adata, "raw", None) is not None:
        adata.raw = None

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    adata = adata[
        (adata.obs["n_genes_by_counts"] > min_genes)
        & (adata.obs["n_genes_by_counts"] < max_genes)
        & (adata.obs["pct_counts_mt"] < max_mt_pct)
    ].copy()

    logger.info(f"[simple_qc_filter] Output shape after QC: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


# =========================
# 5. PSEUDOBULK HELPERS
# =========================

def compute_pseudobulk_matrices(
    adata_in: ad.AnnData,
    sample_col: str = "sample",
    group_col: str = "group",
    out_dir: Path | None = None,
):
    """
    Pseudobulk count aggregation at sample level.
    - Writes counts + design into out_dir (matrices subfolder in pipeline).
    - Returns can_do_deg flag for group comparisons.
    """
    logger.info("[PSEUDOBULK] compute_pseudobulk_matrices called.")

    if out_dir is None:
        out_dir = COMBINED_OUT_DIR
    out_dir = Path(out_dir)

    if sample_col not in adata_in.obs.columns:
        logger.info(f"[PSEUDOBULK] No '{sample_col}' in obs → skip.")
        return None, None, None, False

    if "counts" in adata_in.layers:
        X = adata_in.layers["counts"]
        logger.info("[PSEUDOBULK] Using adata.layers['counts'].")
    else:
        logger.info("[PSEUDOBULK] Using adata.X (may be normalized).")
        X = adata_in.X

    if hasattr(X, "toarray"):
        mat = X.toarray()
    else:
        mat = np.array(X)

    counts_df = pd.DataFrame(
        mat,
        index=adata_in.obs_names,
        columns=adata_in.var_names,
    )

    sample_ids = adata_in.obs[sample_col].astype(str)
    pb_sample = counts_df.groupby(sample_ids).sum()
    pb_sample_file = out_dir / "pseudobulk_counts_by_sample.csv"
    pb_sample.to_csv(pb_sample_file)
    logger.info(f"[PSEUDOBULK] sample-level counts: {pb_sample_file}")

    pb_group = None
    design = None
    can_do_deg = False

    if group_col in adata_in.obs.columns:
        design = (
            adata_in.obs[[sample_col, group_col]]
            .drop_duplicates()
            .set_index(sample_col)
        )
        design = design.loc[pb_sample.index]
        design_file = out_dir / "pseudobulk_design_sample_group.csv"
        design.to_csv(design_file)
        logger.info(f"[PSEUDOBULK] design table: {design_file}")

        group_ids = adata_in.obs[group_col].astype(str)
        pb_group = counts_df.groupby(group_ids).sum()
        pb_group_file = out_dir / "pseudobulk_counts_by_group.csv"
        pb_group.to_csv(pb_group_file)
        logger.info(f"[PSEUDOBULK] group-level counts: {pb_group_file}")

        group_sizes = design[group_col].value_counts()
        logger.info(f"[PSEUDOBULK] samples per group: {group_sizes.to_dict()}")
        if group_sizes.min() >= 2:
            can_do_deg = True
            logger.info("[PSEUDOBULK] Each group has >=2 samples → OK for DESeq2.")
        else:
            logger.info("[PSEUDOBULK] At least one group has <2 samples → DO NOT run DESeq2.")
    else:
        logger.info(f"[PSEUDOBULK] no '{group_col}' in obs → only sample-level counts.")

    return pb_sample, design, pb_group, can_do_deg


def compute_pseudobulk_by_sample_celltype(
    adata_in: ad.AnnData,
    sample_col: str = "sample",
    group_col: str = "group",
    celltype_col: str = "celltype",
    min_cells: int = 20,
    out_dir: Path | None = None,
):
    """
    Pseudobulk at sample × celltype level.
    Writes pseudobulk_counts_by_sample_celltype.csv and design in out_dir.
    """
    logger.info("[PB-CT] compute_pseudobulk_by_sample_celltype called.")

    if out_dir is None:
        out_dir = COMBINED_OUT_DIR
    out_dir = Path(out_dir)

    if sample_col not in adata_in.obs.columns or group_col not in adata_in.obs.columns:
        logger.info("[PB-CT] sample or group column missing → skipping cell-type pseudobulk.")
        return None, None

    if celltype_col not in adata_in.obs.columns:
        logger.info(f"[PB-CT] '{celltype_col}' not in obs → skipping cell-type pseudobulk.")
        return None, None

    if "counts" in adata_in.layers:
        X = adata_in.layers["counts"]
        logger.info("[PB-CT] Using adata.layers['counts'].")
    else:
        X = adata_in.X
        logger.info("[PB-CT] Using adata.X (may be normalized).")

    if hasattr(X, "toarray"):
        mat = X.toarray()
    else:
        mat = np.array(X)

    counts_df = pd.DataFrame(mat, index=adata_in.obs_names, columns=adata_in.var_names)

    comb_key = (
        adata_in.obs[sample_col].astype(str)
        + "|"
        + adata_in.obs[celltype_col].astype(str)
    )
    adata_in.obs["sample_celltype"] = comb_key

    ct_sizes = comb_key.value_counts()
    keep_keys = ct_sizes[ct_sizes >= min_cells].index
    logger.info(f"[PB-CT] keeping {len(keep_keys)} sample×celltype combos with ≥{min_cells} cells.")
    mask = adata_in.obs["sample_celltype"].isin(keep_keys)

    counts_df = counts_df.loc[mask, :].copy()
    comb_key = adata_in.obs.loc[mask, "sample_celltype"]

    pb_ct = counts_df.groupby(comb_key).sum()

    meta = adata_in.obs.loc[mask, [sample_col, group_col, celltype_col, "sample_celltype"]]
    design_ct = meta.drop_duplicates().set_index("sample_celltype")
    pb_ct = pb_ct.loc[design_ct.index]

    pb_ct_file = out_dir / "pseudobulk_counts_by_sample_celltype.csv"
    design_ct_file = out_dir / "pseudobulk_design_sample_celltype.csv"
    pb_ct.to_csv(pb_ct_file)
    design_ct.to_csv(design_ct_file)
    logger.info(f"[PB-CT] cell-type pseudobulk counts: {pb_ct_file}")
    logger.info(f"[PB-CT] cell-type design table: {design_ct_file}")

    return pb_ct, design_ct


# =========================
# 5b. DESeq2 (R) WITHOUT PATHWAYS
# =========================

def run_deseq2_and_pathway_analysis(matrix_dir: Path, deseq2_out_dir: Path):
    """
    Run DESeq2 on pseudobulk counts via rpy2.
    Uses:
      - matrix_dir/pseudobulk_counts_by_sample.csv
      - matrix_dir/pseudobulk_design_sample_group.csv
    Writes all DESeq2 outputs to deseq2_out_dir.
    """
    matrix_dir = Path(matrix_dir)
    deseq2_out_dir = Path(deseq2_out_dir)
    deseq2_out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"[RPY2] run_deseq2_and_pathway_analysis called. matrix_dir={matrix_dir}, out={deseq2_out_dir}")

    if not RPY2_AVAILABLE:
        logger.error("[RPY2] rpy2 not available. Install rpy2 to run DESeq2 automatically.")
        return

    counts_file = matrix_dir / "pseudobulk_counts_by_sample.csv"
    design_file = matrix_dir / "pseudobulk_design_sample_group.csv"
    if not counts_file.exists() or not design_file.exists():
        logger.error(f"[RPY2] Missing pseudobulk files: {counts_file} or {design_file}")
        return

    # Prepare R working directory + paths
    r_workdir = str(deseq2_out_dir)
    r_workdir_escaped = r_workdir.replace("\\", "/").replace('"', '\\"')
    counts_path_escaped = str(counts_file).replace("\\", "/").replace('"', '\\"')
    design_path_escaped = str(design_file).replace("\\", "/").replace('"', '\\"')

    logger.info(f"[RPY2] Using R working directory: {r_workdir_escaped}")

    r_code = r"""
setwd("{wd}")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

message("=== DESeq2 pseudobulk (all group pairs) ===")

counts_path <- "{counts_file}"
design_path <- "{design_file}"

if (!file.exists(counts_path) || !file.exists(design_path)) {
  stop("Missing counts or design file for DESeq2.")
}

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
design <- read.csv(design_path, row.names = 1, check.names = FALSE)

common <- intersect(rownames(counts), rownames(design))
if (length(common) < 2) {
  stop("Not enough overlapping samples between counts and design.")
}

counts <- counts[common, , drop = FALSE]
design <- design[common, , drop = FALSE]

# transpose → genes x samples
counts_t <- t(counts)

if (ncol(counts_t) < 2) {
  stop("Need at least 2 samples for DESeq2.")
}

if (ncol(counts_t) < 4) {
  message("Warning: <4 samples in pseudobulk, DESeq2 may be unstable.")
}

design$group <- as.factor(design$group)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_t),
  colData   = design,
  design    = ~ group
)

dds <- DESeq(dds)

handle_contrast <- function(dds, contrast_vec, label) {
  message("=== Contrast: ", label, " ===")

  res <- results(dds, contrast = contrast_vec)

  # LFC shrink
  res_shrunk <- tryCatch({
    lfcShrink(dds, contrast = contrast_vec, res = res, type = "apeglm")
  }, error = function(e) {
    message("lfcShrink (apeglm) failed for ", label, ": ", e$message,
            " – returning unshrunk results.")
    res
  })

  res_df <- as.data.frame(res_shrunk)
  res_df$gene <- rownames(res_df)
  res_df$comparison <- label

  # classify up/down
  res_df$regulation <- "no_change"
  res_df$regulation[!is.na(res_df$log2FoldChange) &
                    res_df$log2FoldChange > 0 &
                    res_df$padj < 0.05] <- "up"
  res_df$regulation[!is.na(res_df$log2FoldChange) &
                    res_df$log2FoldChange < 0 &
                    res_df$padj < 0.05] <- "down"

  base <- paste0("deseq2_", label)
  write.csv(res_df, file = paste0(base, "_results.csv"), row.names = FALSE)
  write.csv(res_df[res_df$regulation == "up",   ], file = paste0(base, "_upregulated_genes.csv"),   row.names = FALSE)
  write.csv(res_df[res_df$regulation == "down", ], file = paste0(base, "_downregulated_genes.csv"), row.names = FALSE)

  # volcano
  png(paste0(base, "_volcano.png"), width = 7, height = 6, units = "in", res = 300)
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("down" = "blue", "no_change" = "grey70", "up" = "red")) +
    theme_bw() +
    ggtitle(paste("DESeq2 pseudobulk DEG (volcano)", label))
  dev.off()

  # MA plot
  png(paste0(base, "_MAplot.png"), width = 7, height = 6, units = "in", res = 300)
  plotMA(res_shrunk, ylim = c(-5, 5))
  dev.off()

  # top-50 heatmap if enough DEGs
  sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
  sig <- sig[order(sig$padj), ]
  top_genes <- head(sig$gene, 50)

  if (length(top_genes) > 1) {
    mat <- counts_t[top_genes, , drop = FALSE]
    mat <- t(scale(t(log2(mat + 1))))
    png(paste0(base, "_top50_heatmap.png"), width = 8, height = 10, units = "in", res = 300)
    pheatmap(mat, show_rownames = TRUE, main = paste("Top 50 DEGs:", label))
    dev.off()
  } else {
    message("Not enough significant DEGs for heatmap: ", label)
  }
}

groups <- levels(design$group)
if (is.null(groups)) {
  groups <- sort(unique(design$group))
}

if (length(groups) < 2) {
  stop("Need at least 2 groups for DESeq2 contrasts.")
}

message("Groups detected for DESeq2: ", paste(groups, collapse = ", "))

# Default contrast: second vs first
if (length(groups) >= 2) {
  default_label <- paste0(groups[2], "_vs_", groups[1])
  handle_contrast(dds, c("group", groups[2], groups[1]), default_label)
}

# All pairwise
if (length(groups) >= 2) {
  for (i in seq_len(length(groups) - 1)) {
    for (j in seq((i + 1), length(groups))) {
      gA <- groups[i]
      gB <- groups[j]
      label <- paste0(gB, "_vs_", gA)
      if (exists("default_label") && label == default_label) {
        next
      }
      handle_contrast(dds, c("group", gB, gA), label)
    }
  }
}

message("DESeq2 pseudobulk DE finished.")
"""

    r_code = r_code.replace("{wd}", r_workdir_escaped)
    r_code = r_code.replace("{counts_file}", counts_path_escaped)
    r_code = r_code.replace("{design_file}", design_path_escaped)

    try:
        logger.info("[RPY2] Running DESeq2 pseudobulk in R via rpy2 (all group pairs)...")
        ro.r(r_code)
        logger.info("[RPY2] Finished DESeq2 pseudobulk via rpy2.")
    except Exception as e:
        logger.error(f"[RPY2] Error while running DESeq2: {e}")


# =========================
# 6. GROUP COMPARISON HELPERS (single-cell DE)
# =========================

def compute_de_between_groups_per_cluster(
    adata: ad.AnnData,
    group_col: str = "group",
    cluster_col: str = "leiden",
    out_dir: Path | None = None,
):
    """
    For each cluster, run DE for ALL group pairs (focus vs ref).
    Output: one CSV with all comparisons per cluster.
    """
    if out_dir is None:
        out_dir = COMBINED_OUT_DIR / "06_intercluster_analysis_deg" / "intercluster_group_DE"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if group_col not in adata.obs.columns:
        logger.info(f"[DE-CLUSTER-GROUP] '{group_col}' not in obs → skipping.")
        return

    if cluster_col not in adata.obs.columns:
        logger.info(f"[DE-CLUSTER-GROUP] '{cluster_col}' not in obs → skipping.")
        return

    groups = sorted(adata.obs[group_col].astype(str).unique().tolist())
    if len(groups) < 2:
        logger.info(
            f"[DE-CLUSTER-GROUP] Need at least 2 groups in '{group_col}', found: {groups}"
        )
        return

    logger.info(
        f"[DE-CLUSTER-GROUP] Comparing ALL group pairs within each cluster ({cluster_col}). "
        f"Groups: {groups}"
    )

    all_rows = []
    clusters = sorted(adata.obs[cluster_col].astype(str).unique().tolist())

    for cl in clusters:
        sub = adata[adata.obs[cluster_col].astype(str) == cl].copy()
        sub_group_counts = sub.obs[group_col].astype(str).value_counts()

        if sub_group_counts.size < 2:
            logger.info(
                f"[DE-CLUSTER-GROUP] Cluster {cl}: only one group present "
                f"{sub_group_counts.to_dict()} → skip."
            )
            continue

        for i in range(len(groups) - 1):
            for j in range(i + 1, len(groups)):
                group_ref = groups[i]
                group_focus = groups[j]

                if group_ref not in sub_group_counts.index or group_focus not in sub_group_counts.index:
                    continue

                logger.info(
                    f"[DE-CLUSTER-GROUP] Cluster {cl}: {group_focus} vs {group_ref}, "
                    f"counts={sub_group_counts.to_dict()}"
                )

                try:
                    sc.tl.rank_genes_groups(
                        sub,
                        groupby=group_col,
                        groups=[group_focus],
                        reference=group_ref,
                        method="wilcoxon",
                        n_genes=sub.n_vars,
                    )

                    try:
                        df = sc.get.rank_genes_groups_df(sub, group_focus)
                    except Exception:
                        rg = sub.uns["rank_genes_groups"]
                        names = rg["names"][group_focus]
                        scores = rg["scores"][group_focus]
                        pvals_adj = rg["pvals_adj"][group_focus]
                        logfc = None
                        if "logfoldchanges" in rg:
                            logfc = rg["logfoldchanges"][group_focus]
                        df_dict = {
                            "names": names,
                            "scores": scores,
                            "pvals_adj": pvals_adj,
                        }
                        if logfc is not None:
                            df_dict["logfoldchanges"] = logfc
                        df = pd.DataFrame(df_dict)

                    df["cluster"] = cl
                    df["group_focus"] = group_focus
                    df["group_ref"] = group_ref
                    df["comparison"] = f"{group_focus}_vs_{group_ref}"

                    # Up / down classification + direction_simple
                    if "logfoldchanges" in df.columns and "pvals_adj" in df.columns:
                        df["regulation"] = "no_change"
                        up_mask = (df["logfoldchanges"] > 1.0) & (df["pvals_adj"] < 0.05)
                        down_mask = (df["logfoldchanges"] < -1.0) & (df["pvals_adj"] < 0.05)
                        df.loc[up_mask, "regulation"] = "up"
                        df.loc[down_mask, "regulation"] = "down"

                        df["direction_simple"] = np.where(
                            df["logfoldchanges"] > 0,
                            "upregulated",
                            np.where(
                                df["logfoldchanges"] < 0,
                                "downregulated",
                                "no_expression_change",
                            ),
                        )

                    all_rows.append(df)

                except Exception as e:
                    logger.warning(
                        f"[DE-CLUSTER-GROUP] DE failed for cluster {cl}, "
                        f"{group_focus} vs {group_ref}: {e}"
                    )
                    continue

    if all_rows:
        de_all = pd.concat(all_rows, axis=0, ignore_index=True)
        out_file = out_dir / f"DE_{cluster_col}_ALL_GROUP_PAIRS.csv"
        de_all.to_csv(out_file, index=False)
        logger.info(
            f"[DE-CLUSTER-GROUP] Wrote group-wise DE per cluster (all pairs): {out_file}"
        )
    else:
        logger.info("[DE-CLUSTER-GROUP] No DE tables were generated.")


def plot_groupwise_celltype_proportions(
    adata: ad.AnnData,
    group_col: str = "group",
    celltype_col: str | None = None,
    out_dir: Path | None = None,
):
    if out_dir is None:
        out_dir = COMBINED_OUT_DIR
    out_dir = Path(out_dir)

    if group_col not in adata.obs.columns:
        logger.info(f"[CT-PROP] '{group_col}' not in obs → skip.")
        return

    if celltype_col is None or celltype_col not in adata.obs.columns:
        logger.info("[CT-PROP] celltype column missing → skip.")
        return

    df = (
        adata.obs[[group_col, celltype_col]]
        .astype(str)
        .value_counts()
        .reset_index(name="n_cells")
    )
    total_per_group = df.groupby(group_col)["n_cells"].sum().rename("total_cells")
    df = df.merge(total_per_group, on=group_col, how="left")
    df["fraction"] = df["n_cells"] / df["total_cells"]

    prop_file = out_dir / "celltype_proportions_by_group.csv"
    df.to_csv(prop_file, index=False)
    logger.info(f"[CT-PROP] cell type proportions: {prop_file}")

    pivot = df.pivot(
        index=celltype_col,
        columns=group_col,
        values="fraction"
    ).fillna(0.0)

    fig, ax = plt.subplots(figsize=(12, 5))
    pivot.plot(kind="bar", ax=ax)
    ax.set_ylabel("fraction of cells")
    ax.set_title("Cell type proportions by group")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_dir / "celltype_proportions_by_group.png", dpi=300)
    plt.close(fig)


def plot_group_specific_umaps(
    adata: ad.AnnData,
    group_col: str = "group",
    color_col: str = "leiden",
    out_dir: Path | None = None,
):
    if out_dir is None:
        out_dir = COMBINED_OUT_DIR
    out_dir = Path(out_dir)

    if group_col not in adata.obs.columns:
        logger.info(f"[UMAP-BY-GROUP] '{group_col}' not in obs → skip.")
        return

    sc.settings.figdir = out_dir

    groups = sorted(adata.obs[group_col].astype(str).unique().tolist())
    for g in groups:
        sub = adata[adata.obs[group_col].astype(str) == g, :].copy()
        logger.info(f"[UMAP-BY-GROUP] Plotting group={g}, color={color_col}")
        try:
            sc.pl.umap(
                sub,
                color=[color_col],
                size=15,
                show=False,
                save=f"_group_{g}_{color_col}.png",
            )
        except Exception as e:
            logger.warning(f"[UMAP-BY-GROUP] Failed to plot group {g}: {e}")


def compute_de_by_celltype(
    adata: ad.AnnData,
    celltype_col: str,
    group_col: str = "group",
    deg_root_dir: Path | None = None,
    deg_pairs: Optional[List[tuple[str, str]]] = None,   # [(focus, ref), ...]
    deg_all_pairs: bool = False,                         # if True: all pairwise
    pval_cutoff: float = 0.05,
):

    """
    Single-cell DE per cell type (group comparisons).
    Folder structure:

      06_groupwise_deg/
          celltype_specific_deg/
              Tcell/
                  Tcell_CASE_vs_CONTROL.csv
              Bcell/
                  Bcell_CASE_vs_CONTROL.csv
              ...

    Format:

      names, logfoldchanges, pvals_adj, celltype, comparison, logFC_str
    """
    if deg_root_dir is None:
        deg_root_dir = COMBINED_OUT_DIR / "06_groupwise_deg"

    if celltype_col not in adata.obs.columns or group_col not in adata.obs.columns:
        logger.info(f"[DE-CELLTYPE] Missing '{celltype_col}' or '{group_col}' → skip.")
        return

    groups = sorted(adata.obs[group_col].astype(str).unique().tolist())
    if len(groups) < 2:
        logger.info(f"[DE-CELLTYPE] Need >=2 groups in '{group_col}', found: {groups}")
        return

    # only now create DEG folder
    celltype_deg_root = deg_root_dir / "celltype_specific_deg"
    celltype_deg_root.mkdir(parents=True, exist_ok=True)
    celltypes = sorted(adata.obs[celltype_col].astype(str).unique())
    # Build comparisons
    if deg_pairs is not None and len(deg_pairs) > 0:
        pairs = [(str(f), str(r)) for (f, r) in deg_pairs]
    elif deg_all_pairs:
        pairs = []
        for i in range(len(groups) - 1):
            for j in range(i + 1, len(groups)):
                pairs.append((groups[j], groups[i]))  # later vs earlier
    else:
        pairs = [(groups[1], groups[0])]  # default second vs first

    logger.info(f"[DE-CELLTYPE] DE comparisons (focus vs ref): {pairs}")

    celltypes = adata.obs[celltype_col].astype(str).unique().tolist()
    for ct in celltypes:
        sub = adata[adata.obs[celltype_col].astype(str) == ct, :].copy()
        tab = sub.obs[group_col].value_counts()

        # if any group has <2 cells, we skip DEG for this cell type
        if (tab < 2).any() or tab.size < 2:
            logger.info(
                f"[DE-CELLTYPE] celltype={ct}: not enough cells per group {tab.to_dict()} → skip."
            )
            continue

        for (group_focus, group_ref) in pairs:
            tab = sub.obs[group_col].astype(str).value_counts()

            if group_focus not in tab.index or group_ref not in tab.index:
                logger.info(f"[DE-CELLTYPE] celltype={ct}: missing groups for {group_focus} vs {group_ref} -> skip.")
                continue

            if tab[group_focus] < 2 or tab[group_ref] < 2:
                logger.info(f"[DE-CELLTYPE] celltype={ct}: not enough cells for {group_focus} vs {group_ref} -> {tab.to_dict()} skip.")
                continue

            logger.info(f"[DE-CELLTYPE] celltype={ct}: running {group_focus} vs {group_ref} ...")
            try:
                sc.tl.rank_genes_groups(
                    sub,
                    groupby=group_col,
                    groups=[group_focus],
                    reference=group_ref,
                    method="wilcoxon",
                    n_genes=sub.n_vars,
                )

                try:
                    df = sc.get.rank_genes_groups_df(sub, group_focus)
                except Exception:
                    rg = sub.uns["rank_genes_groups"]
                    names = rg["names"][group_focus]
                    scores = rg["scores"][group_focus] if "scores" in rg else None
                    pvals_adj = rg["pvals_adj"][group_focus] if "pvals_adj" in rg else None

                    logfc = None
                    if "logfoldchanges" in rg:
                        logfc = rg["logfoldchanges"][group_focus]

                    df_dict = {
                        "names": names,
                        "scores": scores if scores is not None else np.nan,
                        "pvals_adj": pvals_adj if pvals_adj is not None else np.nan,
                    }
                    if logfc is not None:
                        df_dict["logfoldchanges"] = logfc
                    df = pd.DataFrame(df_dict)

                keep_cols = ["names", "logfoldchanges", "pvals_adj"]
                missing = [c for c in keep_cols if c not in df.columns]
                if missing:
                    logger.warning(f"[DE-CELLTYPE] Missing columns {missing} for celltype={ct}, {group_focus} vs {group_ref}; skipping.")
                    continue

                df = df[keep_cols].copy()
                df["celltype"] = ct
                df["comparison"] = f"{group_focus}_vs_{group_ref}"
                df["logFC_str"] = df["logfoldchanges"].map(lambda x: f"{x:+.2f}" if pd.notnull(x) else "")

                # regulation
                df["regulation"] = "no_change"
                df.loc[(df["logfoldchanges"] > 0) & (df["pvals_adj"] < pval_cutoff), "regulation"] = "up"
                df.loc[(df["logfoldchanges"] < 0) & (df["pvals_adj"] < pval_cutoff), "regulation"] = "down"

                df_up = df[df["regulation"] == "up"].copy()
                df_down = df[df["regulation"] == "down"].copy()

                ct_safe = str(ct).replace(" ", "_").replace("/", "_")
                ct_dir = celltype_deg_root / ct_safe
                ct_dir.mkdir(parents=True, exist_ok=True)

                base = f"{ct_safe}_{group_focus}_vs_{group_ref}"
                df.to_csv(ct_dir / f"{base}_ALL.csv", index=False)
                df_up.to_csv(ct_dir / f"{base}_UP.csv", index=False)
                df_down.to_csv(ct_dir / f"{base}_DOWN.csv", index=False)

                logger.info(f"[DE-CELLTYPE] Wrote: {ct_dir / f'{base}_ALL.csv'} (+UP/+DOWN)")

            except Exception as e:
                logger.warning(f"[DE-CELLTYPE] Failed for celltype={ct} {group_focus} vs {group_ref}: {e}")
                continue


# =========================
# 7. CELLTYPE MARKERS (UPDATED)
# =========================

def compute_celltype_markers(
    adata: ad.AnnData,
    celltype_col: str,
    out_dir: Path,
    analysis_name: str,
    n_markers_per_type: int = 50,
    reference_dir: Path | None = None,
):
    """
    Celltype-specific markers: for each cell type, DE vs all other cell types.

    Writes:
      - global rankplot / heatmap / dotplot
      - celltype_marker_genes_{celltype_col}_ALL.csv
      - celltype_marker_genes_{celltype_col}_{celltype}.csv per cell type
      - per-celltype dotplots and rankplots in sc_dot_plot_vis / sc_rank_plot_vis
      - copy of ALL markers in reference_dir (if provided)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"[CT-MARKERS] compute_celltype_markers called for column '{celltype_col}'.")

    if celltype_col not in adata.obs.columns:
        logger.info(f"[CT-MARKERS] '{celltype_col}' not in obs → skip.")
        return

    sc.settings.figdir = out_dir

    if not pd.api.types.is_categorical_dtype(adata.obs[celltype_col]):
        adata.obs[celltype_col] = adata.obs[celltype_col].astype("category")

    celltypes = adata.obs[celltype_col].cat.categories.tolist()
    logger.info(f"[CT-MARKERS] Found {len(celltypes)} cell types: {celltypes}")

    heatmap_height = max(8.0, 0.6 * len(celltypes))

    # Subfolders for extra visualizations
    dotplot_dir = out_dir / "sc_dot_plot_vis"
    rankplot_dir = out_dir / "sc_rank_plot_vis"
    dotplot_dir.mkdir(parents=True, exist_ok=True)
    rankplot_dir.mkdir(parents=True, exist_ok=True)

    try:
        logger.info(f"[CT-MARKERS] Running rank_genes_groups by '{celltype_col}' (wilcoxon)...")
        sc.tl.rank_genes_groups(
            adata,
            groupby=celltype_col,
            method="wilcoxon",
            n_genes=n_markers_per_type,
        )

        # Global combined plots
        sc.pl.rank_genes_groups(
            adata,
            n_genes=20,
            sharey=False,
            show=False,
            save=f"_{analysis_name}_celltype_markers_rankplot.png",
        )

        sc.pl.rank_genes_groups_heatmap(
            adata,
            n_genes=10,
            show=False,
            save=f"_{analysis_name}_celltype_markers_heatmap.png",
            figsize=(10, heatmap_height),
        )
        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=10,
            show=False,
            save=f"_{analysis_name}_celltype_markers_dotplot.png",
        )

        # Extract markers_all
        try:
            markers_all = sc.get.rank_genes_groups_df(adata, None)
        except Exception:
            rg = adata.uns["rank_genes_groups"]
            groups = rg["names"].dtype.names
            rows = []
            for g in groups:
                names = rg["names"][g]
                scores = rg["scores"][g] if "scores" in rg else None
                pvals_adj = rg["pvals_adj"][g] if "pvals_adj" in rg else None

                logfc = None
                if "logfoldchanges" in rg:
                    logfc = rg["logfoldchanges"][g]

                for idx, gene in enumerate(names):
                    row = {
                        "group": g,
                        "names": gene,
                        "scores": float(scores[idx]) if scores is not None else np.nan,
                        "pvals_adj": float(pvals_adj[idx]) if pvals_adj is not None else np.nan,
                        "rank": idx + 1,
                    }
                    if logfc is not None:
                        row["logfoldchanges"] = float(logfc[idx])
                    rows.append(row)

            markers_all = pd.DataFrame(rows)


        if "logfoldchanges" in markers_all.columns:
            markers_all["direction_simple"] = np.where(
                markers_all["logfoldchanges"] > 0,
                "upregulated",
                np.where(
                    markers_all["logfoldchanges"] < 0,
                    "downregulated",
                    "no_expression_change",
                ),
            )

        all_file = out_dir / f"celltype_marker_genes_{celltype_col}_ALL.csv"
        markers_all.to_csv(all_file, index=False)
        logger.info(f"[CT-MARKERS] Wrote all celltype markers: {all_file}")

        if reference_dir is not None:
            reference_dir = Path(reference_dir)
            reference_dir.mkdir(parents=True, exist_ok=True)
            ref_file = reference_dir / f"{analysis_name}_celltype_markers_{celltype_col}_ALL.csv"
            markers_all.to_csv(ref_file, index=False)
            logger.info(f"[CT-MARKERS] Copied celltype markers to reference dir: {ref_file}")

        # Per-celltype marker tables + per-celltype dotplots & rankplots
        for ct in celltypes:
           # celltypes = sorted(adata.obs[celltype_col].astype(str).unique())

            sub = markers_all[markers_all["group"] == ct].copy()
            ct_safe = str(ct).replace(" ", "_").replace("/", "_")
            ct_file = out_dir / f"celltype_marker_genes_{celltype_col}_{ct_safe}.csv"
            sub.to_csv(ct_file, index=False)
            logger.info(f"[CT-MARKERS] Wrote markers for celltype '{ct}': {ct_file}")

            # Top markers for this celltype
            if "pvals_adj" in sub.columns:
                sub_sorted = sub.sort_values("pvals_adj")
            else:
                sub_sorted = sub.sort_values("rank")
            top_genes = sub_sorted["names"].head(min(n_markers_per_type, len(sub_sorted))).tolist()
            if not top_genes:
                continue

            # Dotplot: show expression of top markers across all cell types
            sc.settings.figdir = dotplot_dir
            try:
                sc.pl.dotplot(
                    adata,
                    var_names=top_genes,
                    groupby=celltype_col,
                    show=False,
                    save=f"_{analysis_name}_dotplot_{ct_safe}.png",
                )
            except Exception as e:
                logger.warning(f"[CT-MARKERS] Dotplot failed for {ct}: {e}")

            # Rank plot: only this celltype
            sc.settings.figdir = rankplot_dir
            try:
                sc.pl.rank_genes_groups(
                    adata,
                    groups=[ct],
                    n_genes=20,
                    sharey=False,
                    show=False,
                    save=f"_{analysis_name}_rankplot_{ct_safe}.png",
                )
            except Exception as e:
                logger.warning(f"[CT-MARKERS] Rankplot (single celltype) failed for {ct}: {e}")

        # Restore figdir to out_dir
        sc.settings.figdir = out_dir

    except Exception as e:
        logger.warning(f"[CT-MARKERS] Failed to compute celltype markers: {e}")


# =========================
# 8. PATHWAY ENRICHMENT HELPERS (gseapy / Enrichr)
# =========================

def deduplicate_pathways_semantic(
    df: pd.DataFrame,
    combined_dir: Path,
    prefix: str,
    sim_threshold: float = 0.9,
) -> pd.DataFrame:
    """
    Optional semantic deduplication using SentenceTransformer (MiniLM) + FAISS.
    Falls back to simple string-based dedup if libraries are not available.
    Writes a log file describing which pathways were removed and why.
    """
    combined_dir = Path(combined_dir)
    combined_dir.mkdir(parents=True, exist_ok=True)
    log_file = combined_dir / f"{prefix}_pathway_dedup_log.txt"

    # Reset index so we can refer to integer positions
    df = df.reset_index(drop=True)

    # Priority: lower Adjusted P-value first, then higher Combined Score
    def _score(row):
        p = row.get("Adjusted P-value", np.nan)
        cs = row.get("Combined Score", 0.0)
        if pd.isna(p):
            p = 1.0
        return (p, -cs)

    order = sorted(range(len(df)), key=lambda i: _score(df.iloc[i]))

    log_lines = []
    log_lines.append("=== Pathway deduplication (semantic) ===")
    log_lines.append(f"Original rows: {len(df)}")
    log_lines.append(f"Similarity threshold: {sim_threshold}")
    log_lines.append("Priority: min(Adjusted P-value), then max(Combined Score)")
    log_lines.append("")

    # Try semantic dedup
    try:
        from sentence_transformers import SentenceTransformer
        import faiss

        # Friendly MiniLM model
        model = SentenceTransformer("all-MiniLM-L6-v2")

        names = df["Pathways"].astype(str).tolist()
        embeddings = model.encode(names, convert_to_numpy=True, show_progress_bar=False)
        embeddings = embeddings.astype("float32")

        faiss.normalize_L2(embeddings)
        dim = embeddings.shape[1]
        index = faiss.IndexFlatIP(dim)

        keep_mask = np.zeros(len(df), dtype=bool)

        kept_indices = []
        for idx in order:
            if not kept_indices:
                index.add(embeddings[idx:idx+1])
                kept_indices.append(idx)
                keep_mask[idx] = True
                log_lines.append(f"KEEP idx={idx}: {df.loc[idx, 'Biological_Database']} :: {df.loc[idx, 'Pathways']}")
                continue

            D, I = index.search(embeddings[idx:idx+1], k=1)
            sim = float(D[0][0])
            dup_idx = int(I[0][0])

            if sim >= sim_threshold:
                log_lines.append(
                    f"REMOVE idx={idx} (sim={sim:.3f} vs idx={dup_idx}) → "
                    f"{df.loc[idx, 'Biological_Database']} :: {df.loc[idx, 'Pathways']}"
                )
            else:
                index.add(embeddings[idx:idx+1])
                kept_indices.append(idx)
                keep_mask[idx] = True
                log_lines.append(f"KEEP idx={idx}: {df.loc[idx, 'Biological_Database']} :: {df.loc[idx, 'Pathways']}")

        new_df = df.loc[keep_mask].reset_index(drop=True)
        log_lines.append("")
        log_lines.append(f"Final rows after semantic dedup: {len(new_df)}")

    except Exception as e:
        # Fallback: simple dedup by (Biological_Database, Pathways) with best Adjusted P-value
        log_lines.append("")
        log_lines.append(f"Semantic dedup skipped (reason: {e}).")
        log_lines.append("Using simple exact-string dedup instead.")
        df["_score_adj"] = df["Adjusted P-value"].fillna(1.0)
        df["_score_comb"] = -df.get("Combined Score", 0.0).fillna(0.0)
        df = df.sort_values(["Biological_Database", "Pathways", "_score_adj", "_score_comb"])
        new_df = df.drop_duplicates(subset=["Biological_Database", "Pathways"], keep="first").copy()
        new_df = new_df.drop(columns=["_score_adj", "_score_comb"])
        new_df = new_df.reset_index(drop=True)
        log_lines.append(f"Final rows after string dedup: {len(new_df)}")

    log_file.write_text("\n".join(log_lines), encoding="utf-8")
    logger.info(f"[ENRICHR-DEDUP] Wrote pathway deduplication log: {log_file}")
    return new_df


def run_enrichr_multidb(
    gene_list,
    out_dir: Path,
    prefix: str,
    pval_cutoff: float = 0.05,
):
    """
    Run Enrichr-based enrichment on a list of gene symbols across multiple databases.

    Output columns (per file):

        Biological_Database    Pathways    Overlap    P-value    Adjusted P-value
        Odds Ratio    Combined Score    Genes

    GO results are stored in separate folders: GO_BP, GO_MF, GO_CC, KEGG, Reactome, WikiPathways.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not GSEAPY_AVAILABLE:
        logger.warning(f"[ENRICHR] gseapy not available; skipping enrichment for {prefix}.")
        return

    genes = [g for g in set(map(str, gene_list)) if g not in (None, "", "nan")]
    if len(genes) < 5:
        logger.info(f"[ENRICHR] Not enough genes for enrichment ({len(genes)}) for {prefix}. Skipping.")
        return

    # Use current Enrichr library names that are typically available
    gene_sets = [
        "GO_Biological_Process_2021",
        "GO_Molecular_Function_2021",
        "GO_Cellular_Component_2021",
        "KEGG_2021_Human",
        "Reactome_2022",
        # Use 2019 WikiPathways Human, more widely available than 2023
        "WikiPathways_2019_Human",
    ]

    # (friendly_name, folder_name)
    gs_meta = {
        "GO_Biological_Process_2021": (
            "Gene Ontology – Biological Process (GO BP)",
            "GO_BP",
        ),
        "GO_Molecular_Function_2021": (
            "Gene Ontology – Molecular Function (GO MF)",
            "GO_MF",
        ),
        "GO_Cellular_Component_2021": (
            "Gene Ontology – Cellular Component (GO CC)",
            "GO_CC",
        ),
        "KEGG_2021_Human": (
            "KEGG Pathway Database",
            "KEGG",
        ),
        "Reactome_2022": (
            "Reactome Pathway Database",
            "Reactome",
        ),
        "WikiPathways_2019_Human": (
            "WikiPathways Database",
            "WikiPathways",
        ),
    }

    logger.info(f"[ENRICHR] Running enrichment for {prefix} on {len(genes)} genes...")

    base_pathway_dir = out_dir / "pathways"
    base_pathway_dir.mkdir(parents=True, exist_ok=True)
    combined_dir = base_pathway_dir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)

    combined_results = []

    def _plot_top_bar(df: pd.DataFrame, out_png: Path, top_n: int = 20):
        if "Combined Score" not in df.columns:
            return
        df_plot = df.sort_values("Combined Score", ascending=False).head(top_n)
        if df_plot.empty:
            return

        # More colorful & readable barplot
        plt.figure(figsize=(12, max(5.0, 0.5 * len(df_plot))))
        colors = plt.cm.tab20(np.linspace(0, 1, len(df_plot)))
        plt.barh(df_plot["Pathways"], df_plot["Combined Score"], color=colors)
        plt.gca().invert_yaxis()
        plt.xlabel("Combined Score", fontsize=13)
        plt.ylabel("Pathways", fontsize=12)
        plt.title(out_png.stem.replace("_", " "), fontsize=14)
        plt.xticks(fontsize=11)
        plt.yticks(fontsize=9)
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
        plt.close()

    for gs in gene_sets:
        friendly_name, db_folder = gs_meta.get(gs, (gs, "Other"))
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gs,
                outdir=None,
                cutoff=pval_cutoff,
            )
        except Exception as e:
            logger.warning(
                f"[ENRICHR] Enrichment failed for {prefix}, gene_set={gs}: {e}"
            )
            continue

        if enr is None or getattr(enr, "results", None) is None or enr.results.empty:
            logger.info(f"[ENRICHR] No significant terms for {prefix}, gene_set={gs}.")
            continue

        df_res = enr.results.copy()

        # Remove legacy P-values
        for col in ["Old P-value", "Old Adjusted P-value"]:
            if col in df_res.columns:
                df_res = df_res.drop(columns=[col])

        # Enrichr outputs: ['Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes', 'Gene_set', ...]
        if "Gene_set" not in df_res.columns:
            df_res["Gene_set"] = gs

        if "Genes" not in df_res.columns:
            df_res["Genes"] = ""

        # Rename + set human-readable Biological_Database and Pathways
        df_res["Biological_Database"] = friendly_name
        if "Term" in df_res.columns:
            df_res = df_res.rename(columns={"Term": "Pathways"})

        # Drop the original Gene_set column, we now use Biological_Database
        if "Gene_set" in df_res.columns:
            df_res = df_res.drop(columns=["Gene_set"])

        # Reorder columns to the exact layout you want
        desired_cols = [
            "Biological_Database",
            "Pathways",
            "Overlap",
            "P-value",
            "Adjusted P-value",
            "Odds Ratio",
            "Combined Score",
            "Genes",
        ]
        cols = [c for c in desired_cols if c in df_res.columns] + [
            c for c in df_res.columns if c not in desired_cols
        ]
        df_res = df_res[cols]

        # Write per-database result into its own folder
        db_dir = base_pathway_dir / db_folder
        db_dir.mkdir(parents=True, exist_ok=True)

        out_file = db_dir / f"{prefix}_{db_folder}_enrichment.csv"
        df_res.to_csv(out_file, index=False)
        logger.info(f"[ENRICHR] Saved {out_file} with {df_res.shape[0]} rows.")

        barplot_file = db_dir / f"{prefix}_{db_folder}_top20_barplot.png"
        _plot_top_bar(df_res, barplot_file, top_n=20)

        combined_results.append(df_res)

    if combined_results:
        combined_df = pd.concat(combined_results, axis=0, ignore_index=True)

        # Semantic + exact dedup
        combined_dedup = deduplicate_pathways_semantic(
            combined_df,
            combined_dir=combined_dir,
            prefix=prefix,
        )

        combined_file_raw = combined_dir / f"{prefix}_combined_pathways_RAW.csv"
        combined_df.to_csv(combined_file_raw, index=False)

        combined_file_clean = combined_dir / f"{prefix}_combined_pathways_DEDUP.csv"
        combined_dedup.to_csv(combined_file_clean, index=False)

        logger.info(
            f"[ENRICHR] Saved combined pathways for {prefix}: "
            f"{combined_file_raw.name} (raw), {combined_file_clean.name} (deduplicated)"
        )
    else:
        logger.info(f"[ENRICHR] No enrichment results to combine for {prefix}.")


def run_cluster_marker_enrichment(
    markers_all: pd.DataFrame,
    out_dir: Path,
    analysis_name: str,
    pval_col: str = "pvals_adj",
    pval_cutoff: float = 0.05,
    top_n: int = 200,
    cluster_celltype_map: Optional[Dict[str, str]] = None,
):
    """
    For each Leiden cluster (markers_all['group']),
    run multi-DB Enrichr enrichment on top marker genes.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if "group" not in markers_all.columns or "names" not in markers_all.columns:
        logger.warning("[CLUSTER-ENRICH] markers_all must have 'group' and 'names' columns.")
        return

    use_pval = pval_col in markers_all.columns

    for cl in sorted(markers_all["group"].unique()):
        sub = markers_all[markers_all["group"] == cl].copy()

        if use_pval:
            sub = sub.sort_values(pval_col).dropna(subset=[pval_col])
            sig = sub[sub[pval_col] < pval_cutoff]
            if sig.empty:
                logger.info(f"[CLUSTER-ENRICH] No significant markers for cluster {cl}. Skipping.")
                continue
            gene_list = sig["names"].head(top_n).tolist()
        else:
            if "rank" not in sub.columns:
                logger.warning(f"[CLUSTER-ENRICH] No '{pval_col}' or 'rank' column for cluster {cl}. Skipping.")
                continue
            sub = sub.sort_values("rank")
            gene_list = sub["names"].head(top_n).tolist()

        if not gene_list:
            logger.info(f"[CLUSTER-ENRICH] Empty gene list for cluster {cl}. Skipping.")
            continue

        cl_str = str(cl)
        if cluster_celltype_map is not None and cl_str in cluster_celltype_map:
            cl_label = cluster_celltype_map[cl_str]   # e.g. "1_Tcell"
        else:
            cl_label = cl_str

        prefix = f"{analysis_name}_cluster_{cl_label}"
        run_enrichr_multidb(
            gene_list=gene_list,
            out_dir=out_dir,
            prefix=prefix,
            pval_cutoff=pval_cutoff,
        )


def run_group_cluster_deg_enrichment_from_file(
    deg_table_path: Path,
    out_dir: Path,
    analysis_name: str,
    pval_col: str = "pvals_adj",
    logfc_col: str = "logfoldchanges",
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.25,
):
    """
    Use the table DE_<cluster_col>_ALL_GROUP_PAIRS.csv (from
    compute_de_between_groups_per_cluster) and run Enrichr (gseapy)
    for EACH (cluster, comparison) pair.

    Prefix example:
        combined_all_samples_cluster0_CASE_vs_CONTROL
    """
    deg_table_path = Path(deg_table_path)
    if not deg_table_path.exists():
        logger.info(f"[GRP-CL-DEG-ENRICH] DEG table not found: {deg_table_path}")
        return

    df = pd.read_csv(deg_table_path)
    required_cols = {"names", "cluster", "comparison", pval_col, logfc_col}
    missing = required_cols - set(df.columns)
    if missing:
        logger.warning(
            f"[GRP-CL-DEG-ENRICH] Missing columns {missing} in {deg_table_path.name}; skipping."
        )
        return

    wrote_any = False
    for (cl, comp), sub in df.groupby(["cluster", "comparison"]):
        sub = sub.dropna(subset=[pval_col, logfc_col])
        sub = sub[sub[pval_col] < pval_cutoff]
        sub = sub[sub[logfc_col].abs() >= logfc_cutoff]

        if sub.empty:
            continue

        if not wrote_any:
            out_dir = Path(out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            wrote_any = True

        prefix = f"{analysis_name}_cluster{cl}_{comp}"
        gene_list = sub["names"].tolist()
        logger.info(
            f"[GRP-CL-DEG-ENRICH] Running enrichment for cluster={cl}, comparison={comp} "
            f"({len(gene_list)} genes, prefix={prefix})"
        )
        run_enrichr_multidb(
            gene_list=gene_list,
            out_dir=out_dir,
            prefix=prefix,
            pval_cutoff=pval_cutoff,
        )


def run_gseapy_on_deseq2_results(
    deseq2_dir: Path,
    pathway_out_dir: Path,
    pval_col: str = "padj",
    logfc_col: str = "log2FoldChange",
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.0,
):
    """
    Scan DESeq2 *_results.csv files from pseudobulk, pick significant DEGs,
    and run multi-database enrichment via gseapy (Enrichr).
    Also writes a filtered DEG file with 'direction_simple'.
    """
    deseq2_dir = Path(deseq2_dir)
    pathway_out_dir = Path(pathway_out_dir)
    pathway_out_dir.mkdir(parents=True, exist_ok=True)

    if not GSEAPY_AVAILABLE:
        logger.warning("[DESEQ2-GSEAPY] gseapy not available; skipping pseudobulk pathway enrichment.")
        return

    files = list(deseq2_dir.glob("deseq2_*_results.csv"))
    if not files:
        logger.info("[DESEQ2-GSEAPY] No deseq2_*_results.csv files found.")
        return

    for f in files:
        df = pd.read_csv(f)
        if pval_col not in df.columns or logfc_col not in df.columns or "gene" not in df.columns:
            logger.warning(f"[DESEQ2-GSEAPY] {f.name} missing required columns; skipping.")
            continue

        df2 = df.dropna(subset=[pval_col, logfc_col]).copy()
        df2 = df2[df2[pval_col] < pval_cutoff]
        df2 = df2[abs(df2[logfc_col]) > logfc_cutoff]

        if df2.empty:
            logger.info(f"[DESEQ2-GSEAPY] No significant DEGs in {f.name}; skipping.")
            continue

        df2["direction_simple"] = np.where(
            df2[logfc_col] > 0,
            "upregulated",
            np.where(df2[logfc_col] < 0, "downregulated", "no_expression_change"),
        )

        filtered_out = pathway_out_dir / f"{f.stem}_sigDEGs_with_direction.csv"
        df2.to_csv(filtered_out, index=False)

        gene_list = df2["gene"].astype(str).tolist()
        prefix = f.stem
        logger.info(f"[DESEQ2-GSEAPY] Running enrichment for {prefix} on {len(gene_list)} genes...")
        run_enrichr_multidb(
            gene_list=gene_list,
            out_dir=pathway_out_dir,
            prefix=prefix,
            pval_cutoff=0.05,
        )


# =========================
# 8b. CELLTYPE DEG + MARKER + PATHWAY SUMMARY
# =========================

def summarize_celltype_degs_markers_pathways(
    out_dir: Path,
    analysis_name: str,
    deg_dir: Path,
    celltype_dir: Path,
    ct_deg_pathway_dir: Path | None = None,
    top_n_genes: int = 30,
    top_n_pathways: int = 10,
):
    """
    Summary across:
      - celltype-specific DEGs (group A vs group B)
      - whether each DEG is also a cell-type marker
      - a simple score (2 = DEG+marker, 1 = DEG only)
      - top pathways per celltype comparison
    Outputs in 09_reference_summary.
    """
    out_dir = Path(out_dir)
    deg_dir = Path(deg_dir)
    celltype_dir = Path(celltype_dir)
    ref_dir = out_dir / "09_reference_summary"
    ref_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"[SUMMARY-CT-DEG] Building celltype DEG + marker + pathway summary for {analysis_name}")

    markers_all_path = celltype_dir / "celltype_marker_genes_celltype_ALL.csv"
    marker_pairs = set()
    if markers_all_path.exists():
        try:
            markers_all = pd.read_csv(markers_all_path)
            if {"group", "names"}.issubset(markers_all.columns):
                for _, row in markers_all.iterrows():
                    ct = str(row["group"])
                    g = str(row["names"])
                    marker_pairs.add((ct, g))
                logger.info(f"[SUMMARY-CT-DEG] Loaded {len(markers_all)} marker rows from {markers_all_path.name}")
            else:
                logger.warning(
                    f"[SUMMARY-CT-DEG] Marker file {markers_all_path.name} missing 'group' or 'names' columns."
                )
        except Exception as e:
            logger.warning(f"[SUMMARY-CT-DEG] Failed to read {markers_all_path}: {e}")
    else:
        logger.info(
            f"[SUMMARY-CT-DEG] No marker ALL file found at {markers_all_path}; "
            "DEGs will still be summarised but without marker scores."
        )

    pathways_combined_dir = None
    if ct_deg_pathway_dir is not None:
        ct_deg_pathway_dir = Path(ct_deg_pathway_dir)
        candidate = ct_deg_pathway_dir / "pathways" / "combined"
        if candidate.exists():
            pathways_combined_dir = candidate
            logger.info(f"[SUMMARY-CT-DEG] Using combined pathway dir: {pathways_combined_dir}")

    # new layout: deg_dir/celltype_specific_deg/<CT>/<CT_CASE_vs_CONTROL>.csv
    ct_deg_files = sorted(deg_dir.glob("*/*.csv"))
    if not ct_deg_files:
        logger.info("[SUMMARY-CT-DEG] No celltype-specific DEG CSVs found.")
        return

    genelevel_rows = []
    overview_rows = []
    summary_txt_lines = [
        f"=== Cell-type DEGs + markers + pathways summary for {analysis_name} ===",
        "",
        "Score legend:",
        "  2 = gene is DEG AND cell-type marker",
        "  1 = gene is DEG only",
        "",
    ]

    for f in ct_deg_files:
        try:
            df = pd.read_csv(f)
        except Exception as e:
            logger.warning(f"[SUMMARY-CT-DEG] Failed to read {f.name}: {e}")
            continue

        if not {"names", "celltype", "comparison", "logfoldchanges", "pvals_adj"}.issubset(df.columns):
            logger.warning(
                f"[SUMMARY-CT-DEG] File {f.name} missing required columns; skipping."
            )
            continue

        ct_name = str(df["celltype"].iloc[0])
        comparison_label = str(df["comparison"].iloc[0])

        pcol = "pvals_adj"
        lfc_col = "logfoldchanges"

        df_deg = df.dropna(subset=[pcol, lfc_col]).copy()
        df_deg = df_deg[df_deg[pcol] < 0.05]
        df_deg = df_deg[df_deg[lfc_col].abs() > 0.25]

        if df_deg.empty:
            summary_txt_lines.append(
                f"Example: {ct_name} DEGs ({comparison_label}) → 0 DEGs passing filters."
            )
            summary_txt_lines.append("")
            continue

        is_marker_list = []
        score_list = []
        for _, row in df_deg.iterrows():
            gene = str(row["names"])
            is_marker = int((ct_name, gene) in marker_pairs)
            score = 2 if is_marker == 1 else 1
            is_marker_list.append(is_marker)
            score_list.append(score)

        df_deg["is_marker_for_celltype"] = is_marker_list
        df_deg["celltype_DEG_marker_score"] = score_list

        genelevel_rows.append(
            df_deg.assign(
                celltype=ct_name,
                comparison=comparison_label,
            )[["celltype", "comparison", "names", lfc_col, pcol,
               "is_marker_for_celltype", "celltype_DEG_marker_score"]].rename(
                columns={
                    "names": "gene",
                    lfc_col: "logFC",
                    pcol: "pval_adj_or_p",
                }
            )
        )

        n_deg = df_deg.shape[0]
        n_markers = int(df_deg["is_marker_for_celltype"].sum())

        df_markers_only = df_deg[df_deg["is_marker_for_celltype"] == 1].copy()
        df_markers_only = df_markers_only.sort_values(pcol)
        top_marker_genes = df_markers_only["names"].head(min(10, len(df_markers_only))).tolist()
        top_marker_genes_str = "; ".join(map(str, top_marker_genes)) if top_marker_genes else ""

        top_pathway_str_list = []
        if pathways_combined_dir is not None:
            prefix = f.stem
            # Prefer deduplicated, then raw, then legacy filename
            candidate_files = [
                pathways_combined_dir / f"{prefix}_combined_pathways_DEDUP.csv",
                pathways_combined_dir / f"{prefix}_combined_pathways_RAW.csv",
                pathways_combined_dir / f"{prefix}_combined_pathways.csv",
            ]
            comb_file = None
            for cf in candidate_files:
                if cf.exists():
                    comb_file = cf
                    break

            if comb_file is not None:
                try:
                    pdf = pd.read_csv(comb_file)
                    if not pdf.empty:
                        if "Combined Score" in pdf.columns:
                            pdf = pdf.sort_values("Combined Score", ascending=False)
                        elif "Adjusted P-value" in pdf.columns:
                            pdf = pdf.sort_values("Adjusted P-value", ascending=True)
                        pdf_top = pdf.head(top_n_pathways)
                        for _, prow in pdf_top.iterrows():
                            gs = str(prow.get("Biological_Database", "NA"))
                            term = str(prow.get("Pathways", "NA"))
                            adjp = float(prow.get("Adjusted P-value", np.nan)) if "Adjusted P-value" in pdf.columns else np.nan
                            if np.isnan(adjp):
                                top_pathway_str_list.append(f"{gs}:: {term}")
                            else:
                                top_pathway_str_list.append(f"{gs}:: {term} (adj_p={adjp:.2e})")
                except Exception as e:
                    logger.warning(f"[SUMMARY-CT-DEG] Failed to load pathways for {prefix}: {e}")

        top_pathways_str = "; ".join(top_pathway_str_list) if top_pathway_str_list else ""

        overview_rows.append(
            {
                "celltype": ct_name,
                "comparison": comparison_label,
                "n_deg_filtered": n_deg,
                "n_deg_markers": n_markers,
                "top_deg_marker_genes": top_marker_genes_str,
                "top_pathways": top_pathways_str,
            }
        )

        summary_txt_lines.append(
            f"Example: {ct_name} DEGs ({comparison_label}) "
            f"→ {n_deg} DEGs passing filters (|logFC|>0.25, p<0.05)"
        )
        summary_txt_lines.append(
            f"  marker DEGs (also cell-type markers): {n_markers}"
        )
        if top_marker_genes:
            summary_txt_lines.append("  Top marker DEGs:")
            for g in top_marker_genes[: min(5, len(top_marker_genes))]:
                summary_txt_lines.append(f"    - {g} (score=2)")
        else:
            summary_txt_lines.append("  Top marker DEGs: none (no DEG overlaps with markers).")

        if top_pathway_str_list:
            summary_txt_lines.append("  Top pathways:")
            for pw in top_pathway_str_list[: min(5, len(top_pathway_str_list))]:
                summary_txt_lines.append(f"    - {pw}")
        else:
            summary_txt_lines.append("  Top pathways: none (no enriched terms found).")

        summary_txt_lines.append("")

    if genelevel_rows:
        genelevel_df = pd.concat(genelevel_rows, axis=0, ignore_index=True)
        genelevel_file = ref_dir / "celltype_DEG_marker_genelevel_summary.csv"
        genelevel_df.to_csv(genelevel_file, index=False)
        logger.info(f"[SUMMARY-CT-DEG] Wrote gene-level DEG+marker summary: {genelevel_file}")

    if overview_rows:
        overview_df = pd.DataFrame(overview_rows)
        overview_file = ref_dir / "celltype_DEG_marker_pathway_overview.csv"
        overview_df.to_csv(overview_file, index=False)
        logger.info(f"[SUMMARY-CT-DEG] Wrote celltype DEG+marker+pathway overview: {overview_file}")

    summary_txt_file = ref_dir / "celltype_DEG_marker_pathway_summary.txt"
    summary_txt_file.write_text("\n".join(summary_txt_lines), encoding="utf-8")
    logger.info(f"[SUMMARY-CT-DEG] Wrote text summary: {summary_txt_file}")



# =========================
# 10. MAIN SCANPY PIPELINE (UPDATED WITH CAR-T)
# =========================

def run_scanpy_pipeline(
    adata: ad.AnnData,
    out_dir: Path,
    analysis_name: str,
    do_pseudobulk: bool = False,
    batch_key: str | None = None,
    integration_method: str | None = None,  # "bbknn" or None
    do_groupwise_de: bool = False,
    group_col: str = "group",
    cluster_col: str = "leiden",
    do_dpt: bool = True,  # optional pseudotime
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Structured subfolders for this analysis
    summary_dir = out_dir / "00_analysis_summary"
    qc_dir = out_dir / "01_qc_and_filtering"
    hvg_dir = out_dir / "02_highly_variable_genes"
    dimred_dir = out_dir / "03_dimensionality_reduction_and_embeddings"
    clustering_dir = out_dir / "04_clustering_and_cell_states"

    # root for cell-type things
    celltype_root_dir = out_dir / "05_celltype_analysis"
    celltype_anno_dir = celltype_root_dir / "celltype_annotation"
    celltype_markers_dir = celltype_root_dir / "celltype_specific_markers"

    #car-t
    cart_dir           = out_dir / "05_CART_analysis"

    # DEG / pathway / pseudobulk dirs
    deg_root_dir = out_dir / "06_groupwise_deg"
    pathway_root_dir = out_dir / "07_pathway_enrichment"
    pseudobulk_root_dir = out_dir / "08_pseudobulk"
    reference_dir = out_dir / "09_reference_summary"

    for d in [
        summary_dir,
        qc_dir,
        hvg_dir,
        dimred_dir,
        clustering_dir,
        celltype_root_dir,
        cart_dir,       # 🔴 make sure this is in the list
        celltype_anno_dir,
        celltype_markers_dir,
        reference_dir,
    ]:
        d.mkdir(parents=True, exist_ok=True)

    logger.info(f"=== Running Scanpy pipeline: {analysis_name} ===")
    logger.info(f"Output root folder: {out_dir}")

    initial_cells = adata.n_obs
    initial_genes = adata.n_vars

    if "counts" not in adata.layers:
        if sp_sparse.issparse(adata.X):
            adata.layers["counts"] = adata.X.copy()
        else:
            adata.layers["counts"] = np.array(adata.X)

    adata = update_gene_names(adata)
    if getattr(adata, "raw", None) is not None:
        adata.raw = None

    sc.pp.filter_genes(adata, min_cells=3)
    genes_after_min_cells = adata.n_vars

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # ========= QC plots =========
    groupby_col = "sample" if "sample" in adata.obs.columns else None
    sc.settings.figdir = qc_dir

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        groupby=groupby_col,
        multi_panel=True,
        show=False,
        save=f"_{analysis_name}_qc_violin.png",
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        show=False,
        save=f"_{analysis_name}_qc_total_vs_mito.png",
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        show=False,
        save=f"_{analysis_name}_qc_total_vs_genes.png",
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].hist(adata.obs["n_genes_by_counts"], bins=60)
    axes[0].set_title("n_genes_by_counts")
    axes[1].hist(adata.obs["total_counts"], bins=60)
    axes[1].set_title("total_counts")
    axes[2].hist(adata.obs["pct_counts_mt"], bins=60)
    axes[2].set_title("pct_counts_mt")
    plt.tight_layout()
    plt.savefig(qc_dir / f"{analysis_name}_qc_metric_histograms.png", dpi=300)
    plt.close()

    # ========= QC filtering =========
    MIN_GENES = 200
    MAX_GENES = 6000
    MAX_MT_PCT = 15.0
    adata = adata[
        (adata.obs["n_genes_by_counts"] > MIN_GENES)
        & (adata.obs["n_genes_by_counts"] < MAX_GENES)
        & (adata.obs["pct_counts_mt"] < MAX_MT_PCT)
    ].copy()

    # ========= Normalization =========
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # ========= HVG selection =========
    n_cells = adata.n_obs
    if n_cells < 4000:
        N_HVG = 2000
    elif n_cells > 200000:
        N_HVG = 4000
    else:
        N_HVG = 4000

    if "sample" in adata.obs.columns:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=N_HVG,
            batch_key="sample",
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=N_HVG,
        )

    sc.settings.figdir = hvg_dir
    sc.pl.highly_variable_genes(adata, show=False, save=f"_{analysis_name}_highly_variable_genes_plot.png")
    hvg_table = adata.var.copy()
    hvg_table.to_csv(hvg_dir / f"{analysis_name}_highly_variable_genes_table.csv")
    hvg_count = int(adata.var["highly_variable"].sum())

    # ========= Dimensionality reduction =========
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", use_highly_variable=True)

    sc.settings.figdir = dimred_dir
    sc.pl.pca_variance_ratio(adata, log=True, show=False, save=f"_{analysis_name}_PCA_variance_explained.png")

    if integration_method == "bbknn" and batch_key is not None and batch_key in adata.obs.columns:
        if SC_EXT_AVAILABLE:
            logger.info(f"[INTEGRATION] Using BBKNN with batch_key='{batch_key}'.")
            sce.pp.bbknn(adata, batch_key=batch_key)
        else:
            logger.warning("[INTEGRATION] scanpy.external (bbknn) not available. Falling back to standard neighbors.")
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    else:
        logger.info("[INTEGRATION] Using standard neighbors (no explicit integration).")
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

    sc.tl.umap(adata)

    if adata.n_obs <= 50000:
        sc.tl.tsne(adata, n_pcs=30, use_rep="X_pca")
    else:
        logger.info(f"[TSNE] Skipping t-SNE because n_obs={adata.n_obs} > 50000.")

    # embeddings coloured by sample / group
    color_cols = []
    if "sample" in adata.obs.columns:
        color_cols.append("sample")
    if "group" in adata.obs.columns:
        color_cols.append("group")

    if color_cols:
        sc.pl.umap(
            adata,
            color=color_cols,
            wspace=0.4,
            size=15,
            show=False,
            save=f"_{analysis_name}_UMAP_samples_groups.png",
        )
        if "X_tsne" in adata.obsm.keys():
            sc.pl.tsne(
                adata,
                color=color_cols,
                wspace=0.4,
                size=15,
                show=False,
                save=f"_{analysis_name}_TSNE_samples_groups.png",
            )
    else:
        sc.pl.umap(adata, size=15, show=False, save=f"_{analysis_name}_UMAP.png")
        if "X_tsne" in adata.obsm.keys():
            sc.pl.tsne(adata, size=15, show=False, save=f"_{analysis_name}_TSNE.png")

    for qc_col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
        if qc_col in adata.obs.columns:
            sc.pl.umap(
                adata,
                color=[qc_col],
                size=15,
                show=False,
                save=f"_{analysis_name}_UMAP_{qc_col}.png",
            )
            if "X_tsne" in adata.obsm.keys():
                sc.pl.tsne(
                    adata,
                    color=[qc_col],
                    size=15,
                    show=False,
                    save=f"_{analysis_name}_TSNE_{qc_col}.png",
                )

    # ========= Clustering =========
    sc.settings.figdir = clustering_dir
    sc.tl.leiden(adata, resolution=0.5)
    clusters = sorted(adata.obs["leiden"].unique().tolist(), key=lambda x: int(x))
    cluster_sizes = adata.obs["leiden"].value_counts().sort_index()
    cluster_sizes.to_csv(
        clustering_dir / f"{analysis_name}_cluster_cell_counts_leiden.csv",
        header=["n_cells"],
    )

    sc.pl.umap(
        adata,
        color=["leiden"],
        legend_loc="on data",
        size=15,
        show=False,
        save=f"_{analysis_name}_UMAP_leiden.png",
    )
    if "X_tsne" in adata.obsm.keys():
        sc.pl.tsne(
            adata,
            color=["leiden"],
            legend_loc="on data",
            size=15,
            show=False,
            save=f"_{analysis_name}_TSNE_leiden.png",
        )
    sc.pl.pca(
        adata,
        color=["leiden"],
        show=False,
        save=f"_{analysis_name}_PCA_leiden.png",
    )

    # ========= DPT (optional pseudotime) =========
    if do_dpt:
        try:
            logger.info(f"[{analysis_name}] Computing diffusion map + DPT (trajectory inference)...")
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
            sc.settings.figdir = dimred_dir
            sc.pl.umap(
                adata,
                color=["dpt_pseudotime"],
                show=False,
                save=f"_{analysis_name}_UMAP_dpt_pseudotime.png",
            )
            if "X_tsne" in adata.obsm.keys():
                sc.pl.tsne(
                    adata,
                    color=["dpt_pseudotime"],
                    show=False,
                    save=f"_{analysis_name}_TSNE_dpt_pseudotime.png",
                )
            sc.pl.diffmap(
                adata,
                color=["dpt_pseudotime"],
                show=False,
                save=f"_{analysis_name}_DIFFMAP_dpt_pseudotime.png",
            )
        except Exception as e:
            logger.warning(f"[{analysis_name}] DPT computation/plotting failed: {e}")

    # ========= Cell type detection / prediction =========
    celltype_col_raw = None
    celltype_source = None
    standard_celltype_col = "celltype"

    candidate_celltype_cols = [
        "cell_type",
        "celltype",
        "CellType",
        "cell_types",
        "celltype_major",
        "cell_type_major",
        "cell_identity",
        "cell_ontology_class",
        "celltype_new",
        "celltype_celltypist",
    ]
    for col in candidate_celltype_cols:
        if col in adata.obs.columns:
            celltype_col_raw = col
            celltype_source = f"provided_obs_metadata ({col})"
            break

    if celltype_col_raw is None and CELLTYPIST_AVAILABLE:
        try:
            from celltypist import models as ct_models
            try:
                ct_models.download_models()
            except Exception as e:
                logger.warning(f"CellTypist model download failed: {e}")
            model_name = "Immune_All_Low.pkl"
            logger.info(
                f"[CELLTYPE-ML] No curated cell-type labels detected; "
                f"applying CellTypist model '{model_name}' ({analysis_name})."
            )
            preds = celltypist.annotate(
                adata,
                model=model_name,
                majority_voting=True,
            )
            plabels = preds.predicted_labels
            if "majority_voting" in plabels.columns:
                adata.obs["celltype_celltypist"] = (
                    plabels["majority_voting"].reindex(adata.obs_names).astype("category")
                )
                celltype_col_raw = "celltype_celltypist"
            elif "predicted_labels" in plabels.columns:
                adata.obs["celltype_celltypist"] = (
                    plabels["predicted_labels"].reindex(adata.obs_names).astype("category")
                )
                celltype_col_raw = "celltype_celltypist"
            if celltype_col_raw is not None:
                celltype_source = "CellTypist (supervised ML cell-type classifier)"
        except Exception as e:
            logger.warning(f"CellTypist failed for {analysis_name}: {e}")

    celltype_col = None
    if celltype_col_raw is not None:
        if not pd.api.types.is_categorical_dtype(adata.obs[celltype_col_raw]):
            adata.obs[celltype_col_raw] = adata.obs[celltype_col_raw].astype("category")
        adata.obs[standard_celltype_col] = adata.obs[celltype_col_raw].astype("category")
        celltype_col = standard_celltype_col

    # ========= Celltype plots & markers =========
    if celltype_col is not None:
        sc.settings.figdir = celltype_anno_dir
        sc.pl.umap(
            adata,
            color=[celltype_col],
            legend_loc="right margin",
            size=15,
            show=False,
            save=f"_{analysis_name}_UMAP_celltypes.png",
        )
        if "X_tsne" in adata.obsm.keys():
            sc.pl.tsne(
                adata,
                color=[celltype_col],
                legend_loc="right margin",
                size=15,
                show=False,
                save=f"_{analysis_name}_TSNE_celltypes.png",
            )
        sc.pl.pca(
            adata,
            color=[celltype_col],
            show=False,
            save=f"_{analysis_name}_PCA_celltypes.png",
        )
        ct_counts = adata.obs[celltype_col].value_counts().sort_values(ascending=False)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.bar(ct_counts.index.astype(str), ct_counts.values)
        ax.set_xticklabels(ct_counts.index.astype(str), rotation=90)
        ax.set_ylabel("cells")
        ax.set_title(f"Cell types ({celltype_col})")
        plt.tight_layout()
        plt.savefig(celltype_anno_dir / f"{analysis_name}_celltype_composition_barplot.png", dpi=300)
        plt.close()

        logger.info(f"[{analysis_name}] Computing celltype-specific marker genes (Scanpy)...")
        compute_celltype_markers(
            adata,
            celltype_col=celltype_col,
            out_dir=celltype_markers_dir,
            analysis_name=analysis_name,
            n_markers_per_type=50,
            reference_dir=reference_dir,
        )

        if "leiden" in adata.obs.columns:
            logger.info(f"[{analysis_name}] Making side-by-side embeddings (leiden vs cell type).")
            sc.settings.figdir = celltype_anno_dir
            sc.pl.pca(
                adata,
                color=["leiden", celltype_col],
                wspace=0.4,
                show=False,
                save=f"_{analysis_name}_PCA_leiden_vs_celltype.png",
            )
            if "X_tsne" in adata.obsm.keys():
                sc.pl.tsne(
                    adata,
                    color=["leiden", celltype_col],
                    wspace=0.4,
                    show=False,
                    save=f"_{analysis_name}_TSNE_leiden_vs_celltype.png",
                )
            sc.pl.umap(
                adata,
                color=["leiden", celltype_col],
                wspace=0.4,
                show=False,
                save=f"_{analysis_name}_UMAP_leiden_vs_celltype.png",
            )

    else:
        logger.info(
            f"[CELLTYPE] No cell type column found for {analysis_name} "
            f"and no ML-based cell-type prediction could be applied."
        )

    # ========= Cluster → celltype mapping =========
    cluster_celltype_map = None
    if celltype_col is not None and "leiden" in adata.obs.columns:
        tmp = adata.obs[["leiden", celltype_col]].dropna()
        if not tmp.empty:
            cluster_celltype_map = {}
            ref_rows = []
            for cl, sub in tmp.groupby("leiden"):
                top_ct = sub[celltype_col].astype(str).value_counts().idxmax()
                safe_ct = str(top_ct).replace(" ", "_").replace("/", "_")
                label = f"{cl}_{safe_ct}"
                cluster_celltype_map[str(cl)] = label
                ref_rows.append(
                    {
                        "leiden_cluster": cl,
                        "major_celltype_label": top_ct,
                        "cluster_celltype_label": label,
                    }
                )
            logger.info(f"[{analysis_name}] Cluster → celltype map: {cluster_celltype_map}")

            reference_dir.mkdir(parents=True, exist_ok=True)
            ref_df = pd.DataFrame(ref_rows)
            ref_map_file = reference_dir / f"{analysis_name}_cluster_to_celltype_map.csv"
            ref_df.to_csv(ref_map_file, index=False)
            logger.info(
                f"[{analysis_name}] Saved cluster→celltype reference table: {ref_map_file}"
            )

    # ========= Cluster markers & enrichment =========
    try:
        logger.info(f"[{analysis_name}] Computing cluster marker genes (leiden, wilcoxon)...")

        intercluster_dir = clustering_dir / "Intercluster_analysis_deg"
        intercluster_dir.mkdir(parents=True, exist_ok=True)
        sc.settings.figdir = intercluster_dir

        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            method="wilcoxon",
            n_genes=50,
        )

        sc.pl.rank_genes_groups(
            adata,
            n_genes=20,
            sharey=False,
            show=False,
            save=f"_{analysis_name}_cluster_markers_rankplot.png",
        )
        heatmap_height_clusters = max(8.0, 0.6 * len(clusters))
        sc.pl.rank_genes_groups_heatmap(
            adata,
            n_genes=10,
            show=False,
            save=f"_{analysis_name}_cluster_markers_heatmap.png",
            figsize=(10, heatmap_height_clusters),
        )
        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=10,
            show=False,
            save=f"_{analysis_name}_cluster_markers_dotplot.png",
        )

        try:
            markers_all = sc.get.rank_genes_groups_df(adata, None)
        except Exception:
            rg = adata.uns["rank_genes_groups"]
            groups = rg["names"].dtype.names
            rows = []

            for g in groups:
                names = rg["names"][g]
                scores = rg["scores"][g] if "scores" in rg else None
                pvals_adj = rg["pvals_adj"][g] if "pvals_adj" in rg else None

                logfc = None
                if "logfoldchanges" in rg:
                    logfc = rg["logfoldchanges"][g]

                for idx, gene in enumerate(names):
                    row = {
                        "group": g,
                        "names": gene,
                        "scores": float(scores[idx]) if scores is not None else np.nan,
                        "pvals_adj": float(pvals_adj[idx]) if pvals_adj is not None else np.nan,
                        "rank": idx + 1,
                    }
                    if logfc is not None:
                        row["logfoldchanges"] = float(logfc[idx])
                    rows.append(row)

            markers_all = pd.DataFrame(rows)


        markers_all["cluster"] = markers_all["group"].astype(str)
        if cluster_celltype_map is not None:
            markers_all["cluster_celltype_label"] = markers_all["cluster"].map(
                cluster_celltype_map
            ).fillna(markers_all["cluster"])

        intercluster_csv = intercluster_dir / "intercluster_cluster_markers.csv"
        markers_all.to_csv(intercluster_csv, index=False)
        logger.info(f"[{analysis_name}] Wrote intercluster markers: {intercluster_csv}")

        for cl, subdf in markers_all.groupby("cluster"):
            label = (
                markers_all.loc[markers_all["cluster"] == cl, "cluster_celltype_label"]
                .iloc[0]
                if "cluster_celltype_label" in markers_all.columns
                else cl
            )
            safe_label = str(label).replace(" ", "_").replace("/", "_")
            out_f = intercluster_dir / f"cluster_{safe_label}_markers.csv"
            subdf.to_csv(out_f, index=False)

        if DO_PATHWAY_CLUSTERING:
            pathway_root_dir.mkdir(parents=True, exist_ok=True)
            cluster_pathway_dir = pathway_root_dir / "cluster_marker_enrichment"
            logger.info(f"[{analysis_name}] Running pathway enrichment for cluster markers...")
            run_cluster_marker_enrichment(
                markers_all,
                out_dir=cluster_pathway_dir,
                analysis_name=analysis_name,
                pval_col="pvals_adj",
                pval_cutoff=0.05,
                top_n=200,
                cluster_celltype_map=cluster_celltype_map,
            )
        else:
            logger.info(
                f"[{analysis_name}] Skipping cluster-marker pathway enrichment "
                f"(DO_PATHWAY_CLUSTERING=False)."
            )

    except Exception as e:
        logger.warning(f"[{analysis_name}] intercluster marker computation failed: {e}")

    # ========= CAR-T SCORING (pre/post, TPEX/TEX/etc) =========
    # ========= CAR-T SCORING (pre/post, TPEX/TEX/etc) =========
    if DO_CART_SCORING:
        try:
            logger.info(f"[CART] Running CAR-T signature scoring for {analysis_name} ...")

            # 1) Compute scores for all CART_SIGNATURES
            compute_cart_scores(adata)   # your function: adds CART_* columns to obs

            # 2) Classify cells into functional CAR-T states
            refine_cart_states(adata)

            # --- NEW: Evidence that signature genes are truly expressed (overall + by CART state) ---
            export_cart_signature_gene_evidence(
                adata,
                cart_dir=cart_dir,
                analysis_name=analysis_name,
                state_col="CART_State_v2",
                use_layer=None,
                expr_threshold=0.0,
            )

            # --- NEW: Patient-level Pre/Post comparison (your main goal) ---
            export_cart_patient_phase_gene_summary(
                adata,
                cart_dir=cart_dir,
                analysis_name=analysis_name,
                patient_col="patient_id",   # comes from metadata.xlsx in MULTI mode
                phase_col="cart_phase",     # comes from metadata.xlsx (recommended)
                fallback_phase_col="group", # if cart_phase missing, uses group
                state_col="CART_State_v2",
                use_layer=None,
                expr_threshold=0.0,
                min_cells=10,
            )

            # --- OPTIONAL: DE markers defining each CART state (Wilcoxon) ---
            export_cart_state_markers(
                adata,
                cart_dir=cart_dir,
                analysis_name=analysis_name,
                state_col="CART_State_v2",
                n_genes=100,
            )

            # 3) Make sure CAR-T folder exists
            cart_dir.mkdir(parents=True, exist_ok=True)

            # 4) Export all CAR-T related obs columns (scores + state)
            cart_cols = [c for c in adata.obs.columns if c.startswith("CART_")]
            if cart_cols:
                cart_scores_file = cart_dir / f"{analysis_name}_CART_scores_obs.csv"
                adata.obs[cart_cols].to_csv(cart_scores_file)
                logger.info(f"[CART] Wrote CAR-T scores table: {cart_scores_file}")

            # 5) Save counts of each CAR-T state
            if "CART_State_v2" in adata.obs.columns:
                cs_counts = adata.obs["CART_State_v2"].value_counts()
                counts_file = cart_dir / f"{analysis_name}_CART_State_v2_counts.csv"
                cs_counts.to_csv(counts_file, header=["n_cells"])
                logger.info(f"[CART] Wrote CAR-T state counts: {counts_file}")

                # Also write a small readable text summary
                lines = ["CART_State_v2 counts:"]
                for st, n in cs_counts.items():
                    lines.append(f"  {st}: {int(n)} cells")
                summary_txt = cart_dir / f"{analysis_name}_CART_summary.txt"
                summary_txt.write_text("\n".join(lines), encoding="utf-8")
                logger.info(f"[CART] Wrote CAR-T summary text: {summary_txt}")

            # 6) UMAP coloured by CAR-T state
            if "CART_State_v2" in adata.obs.columns:
                sc.settings.figdir = cart_dir
                sc.pl.umap(
                    adata,
                    color=["CART_State_v2"],
                    size=15,
                    show=False,
                    save=f"_{analysis_name}_UMAP_CART_State_v2.png",
                )

            # 7) Violin plots of key signatures (TEX, TPEX, Memory, Effector, Proliferation)
            key_scores = [
                "CART_TEX_score",
                "CART_TPEX_score",
                "CART_Effector_TEFF_score",
                "CART_Memory_like_CAR_T_good_score",
                "CART_Proliferation_score",
            ]

            key_scores = [c for c in key_scores if c in adata.obs.columns]
            if key_scores:
                sc.settings.figdir = cart_dir
                # Prefer grouping by group or celltype if available
                groupby_col = None
                if "group" in adata.obs.columns:
                    groupby_col = "group"
                elif "celltype" in adata.obs.columns:
                    groupby_col = "celltype"

                if groupby_col is not None:
                    sc.pl.violin(
                        adata,
                        key_scores,
                        groupby=groupby_col,
                        rotation=90,
                        show=False,
                        save=f"_{analysis_name}_CART_signatures_by_{groupby_col}.png",
                    )
                else:
                    sc.pl.violin(
                        adata,
                        key_scores,
                        rotation=90,
                        show=False,
                        save=f"_{analysis_name}_CART_signatures.png",
                    )

            # 8) Biological legend for TEX/TPEX/etc
            legend_lines = [
                "CAR-T functional state legend (CART_State_v2):",
                "",
                "TEX_terminal: highly exhausted T cells; high inhibitory receptors (PD-1, TIM-3, etc.),",
                "  often reduced proliferative capacity and poor long-term persistence.",
                "",
                "TPEX_like: precursor-exhausted cells; still express TCF7/memory genes plus some",
                "  exhaustion markers; associated with better CAR-T responses and durability.",
                "",
                "TStemCM_like: stem-cell memory–like CAR T cells with strong self-renewal and",
                "  proliferative potential; generally favourable for long-term control.",
                "",
                "Effector_TEFF: cytotoxic effector cells (GZMB, PRF1, IFNG) driving immediate tumour kill.",
                "",
                "Terminal_diff: terminally differentiated KLRG1+/CX3CR1+ effector cells; often short-lived.",
                "",
                "Proliferating_T: actively cycling Ki67+ CAR T cells.",
                "",
                "Naive_or_Tcm: naive/central memory-like cells (CCR7+/CD27+/IL7R+).",
                "",
                "Tem_like: effector memory–like cells.",
                "",
                "Low_signal_T: T cells with low scores across all signatures (unclassified/low activity).",
                "",
                "Non_T_or_other_immune: cells that do not look like T/NK/CAR at all (likely myeloid, B, etc.).",
            ]
            legend_file = cart_dir / "CART_state_legend.txt"
            legend_file.write_text("\n".join(legend_lines), encoding="utf-8")
            logger.info(f"[CART] Wrote CAR-T state legend: {legend_file}")

        except Exception as e:
            logger.warning(f"[CART] Failed to compute/refine CAR-T scores for {analysis_name}: {e}")
    else:
        logger.info(f"[CART] DO_CART_SCORING=False → skipping CAR-T analysis for {analysis_name}.")

    # ========= Save processed AnnData =========
    processed_h5ad = out_dir / f"{analysis_name}_processed_scanpy_output.h5ad"
    adata.write_h5ad(processed_h5ad)

    # ========= Summary report =========
    summary_lines = [
        f"=== {analysis_name} ===",
        f"Output folder: {out_dir}",
        f"Initial cells: {initial_cells}",
        f"Initial genes: {initial_genes}",
        f"Genes after min_cells filter: {genes_after_min_cells}",
        f"Cells after QC filters: {adata.n_obs}",
        f"HVGs used: {hvg_count}",
        f"Final shape (post-HVG selection for embeddings, ALL genes retained for DE): "
        f"{adata.n_obs} cells x {adata.n_vars} genes",
        f"Leiden clusters: {len(clusters)}",
        f"Processed AnnData (Scanpy): {processed_h5ad}",
    ]
    if celltype_col is not None:
        ct_counts2 = adata.obs[celltype_col].value_counts()
        summary_lines.append(f"Celltype column (standard): {celltype_col}")
        if celltype_source is not None:
            summary_lines.append(f"Celltype annotation source: {celltype_source}")
        summary_lines.append("Celltype counts:")
        for ct, n in ct_counts2.items():
            summary_lines.append(f"  {ct}: {int(n)} cells")
    else:
        summary_lines.append("Celltype annotation: not available")

    # CAR-T state summary (supports old and new state columns)
    for col in ["CART_State", "CART_State_v2"]:
        if col in adata.obs.columns:
            cs_counts = adata.obs[col].value_counts()
            summary_lines.append(f"{col} counts:")
            for st, n in cs_counts.items():
                summary_lines.append(f"  {st}: {int(n)} cells")


    summary_file = summary_dir / f"{analysis_name}_analysis_summary.txt"
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))

    logger.info("\n".join(summary_lines))

    # ========= Pseudobulk + DESeq2 + GSEApy =========
    if do_pseudobulk:
        pseudobulk_root_dir.mkdir(parents=True, exist_ok=True)

        pb_mtx_dir = pseudobulk_root_dir / "matrix"
        pb_mtx_dir.mkdir(parents=True, exist_ok=True)

        pb_deseq2_global_dir = pseudobulk_root_dir / "deseq2_global"
        pb_deseq2_global_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"[{analysis_name}] Running pseudobulk aggregation (sample + group) into: {pb_mtx_dir}")
        pb_sample, design, pb_group, can_do_deg = compute_pseudobulk_matrices(
            adata,
            sample_col="sample",
            group_col=group_col,
            out_dir=pb_mtx_dir,
        )

        logger.info(f"[{analysis_name}] Running cell-type pseudobulk aggregation...")
        if celltype_col is not None:
            pb_ct, design_ct = compute_pseudobulk_by_sample_celltype(
                adata,
                sample_col="sample",
                group_col=group_col,
                celltype_col=celltype_col,
                out_dir=pb_mtx_dir,
            )
        else:
            logger.info("[PSEUDOBULK] No celltype column available; skipping cell-type pseudobulk.")
            pb_ct, design_ct = None, None

        if can_do_deg:
            logger.info(f"[{analysis_name}] Calling run_deseq2_and_pathway_analysis on matrices {pb_mtx_dir} → {pb_deseq2_global_dir} ...")
            run_deseq2_and_pathway_analysis(pb_mtx_dir, pb_deseq2_global_dir)

            if DO_PATHWAY_CLUSTERING:
                pathway_root_dir.mkdir(parents=True, exist_ok=True)
                pb_pathway_global_dir = pathway_root_dir / "pseudobulk_global"
                logger.info(f"[{analysis_name}] Running GSEApy enrichment for DESeq2 pseudobulk results...")
                run_gseapy_on_deseq2_results(
                    deseq2_dir=pb_deseq2_global_dir,
                    pathway_out_dir=pb_pathway_global_dir,
                )
        else:
            logger.info(
                f"[{analysis_name}] Skipping DESeq2: pseudobulk requires ≥2 samples per group "
                f"for a valid group design."
            )
    else:
        logger.info(
            f"[{analysis_name}] Pseudobulk + DESeq2 step disabled "
            f"(do_pseudobulk=False)."
        )

    # ========= Group-wise DE & downstream =========
    if do_groupwise_de:
        deg_root_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"[{analysis_name}] Comparing cell type proportions between groups...")
        plot_groupwise_celltype_proportions(
            adata,
            group_col=group_col,
            celltype_col=celltype_col,
            out_dir=celltype_anno_dir,
        )

        ct_deg_pathway_dir = pathway_root_dir / "celltype_DEG_enrichment"

        if celltype_col is not None:
            logger.info(f"[{analysis_name}] Running single-cell DE per cell type (Scanpy)...")
            compute_de_by_celltype(
                adata,
                celltype_col=celltype_col,
                group_col=group_col,
                deg_root_dir=deg_root_dir,
                deg_all_pairs=True,
                pval_cutoff=0.05,
            )
            # for deg diff pairs:
            # compute_de_by_celltype(
            #     adata,
            #     celltype_col=celltype_col,
            #     group_col=group_col,
            #     deg_root_dir=deg_root_dir,
            #     deg_pairs=[("Post", "Pre"), ("D7", "D0")],
            #     pval_cutoff=0.05,
            # )

            try:
                summarize_celltype_degs_markers_pathways(
                    out_dir=out_dir,
                    analysis_name=analysis_name,
                    deg_dir=deg_root_dir / "celltype_specific_deg",
                    celltype_dir=celltype_markers_dir,
                    ct_deg_pathway_dir=ct_deg_pathway_dir,
                )
            except Exception as e:
                logger.warning(f"[SUMMARY-CT-DEG] Failed to summarise celltype DEGs+markers+pathways: {e}")

        logger.info(f"[{analysis_name}] Building group-specific UMAPs...")
        color_col_for_umap = celltype_col if celltype_col is not None else cluster_col
        group_umap_dir = dimred_dir / "groupwise_embeddings"
        group_umap_dir.mkdir(parents=True, exist_ok=True)
        plot_group_specific_umaps(
            adata,
            group_col=group_col,
            color_col=color_col_for_umap,
            out_dir=group_umap_dir,
        )

    logger.info(f"=== Done Scanpy pipeline: {analysis_name} ===")


# =========================
# 11. DRIVER LOGIC (10x ONLY, UPDATED META HANDLING)
# =========================

if SINGLE_MODE:
    adata_single_raw = load_10x_feature_barcode_matrix(SINGLE_10X_DIR)

    adata_single_raw.obs["sample"] = SINGLE_SAMPLE_LABEL
    if SINGLE_GROUP_LABEL is not None:
        adata_single_raw.obs["group"] = SINGLE_GROUP_LABEL

    run_scanpy_pipeline(
        adata_single_raw,
        COMBINED_OUT_DIR,
        analysis_name="single_dataset",
        do_pseudobulk=False,  # single sample → no DESeq2
        batch_key=None,
        integration_method=None,
        do_groupwise_de=False,
        do_dpt=False,
    )

else:
    base_dir = MULTI_BASE_DIR

    # ===== NEW META LOADING (supports pre/post CAR-T, patient_id, etc.) =====
    meta_path = base_dir / MULTI_META_FILENAME
    meta_df = None

    if meta_path.exists():
        try:
            meta_df = pd.read_excel(meta_path)
            if "sample" not in meta_df.columns:
                logger.warning(f"[META] {meta_path.name} has no 'sample' column; ignoring.")
                meta_df = None
            else:
                meta_df["sample"] = meta_df["sample"].astype(str)
                meta_df = meta_df.set_index("sample")
                logger.info(
                    f"[META] Loaded metadata for {meta_df.shape[0]} samples from {meta_path.name}"
                )
        except Exception as e:
            logger.warning(f"[META] Failed to read {meta_path}: {e}")
    else:
        logger.info(f"[META] No meta file found at {meta_path}; using folder names as groups.")

    group_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
    if not group_dirs:
        raise FileNotFoundError(f"No group folders found in: {base_dir}")

    all_adatas_qc = []
    n_samples_total = 0

    for gdir in sorted(group_dirs):
        default_group_name = gdir.name
        logger.info(f"[MULTI-10X] Group folder: {default_group_name}")
        sample_dirs = [d for d in gdir.iterdir() if d.is_dir()]

        if not sample_dirs:
            logger.warning(
                f"[MULTI-10X] no sample subfolders under {default_group_name}, "
                f"treating group folder as a single sample."
            )
            sample_dirs = [gdir]

        for sdir in sorted(sample_dirs):
            sample_label = sdir.name

            # Defaults
            group_name = default_group_name
            extra_meta = {}

            if meta_df is not None and sample_label in meta_df.index:
                row = meta_df.loc[sample_label]
                if "group" in row.index:
                    group_name = str(row["group"])
                for col in ["cart_phase", "cart_compartment", "patient_id", "response"]:
                    if col in row.index:
                        extra_meta[col] = str(row[col])

            logger.info(
                f"[MULTI-10X] Loading sample: sample={sample_label}, group={group_name}, extra_meta={extra_meta}"
            )

            try:
                adata_raw = load_10x_feature_barcode_matrix(sdir)
            except FileNotFoundError as e:
                logger.warning(f"[MULTI-10X] skip sample {sample_label}: {e}")
                continue

            adata_raw.obs["sample"] = sample_label
            adata_raw.obs["group"] = group_name
            for k, v in extra_meta.items():
                adata_raw.obs[k] = v

            adata_qc = simple_qc_filter(adata_raw)
            all_adatas_qc.append(adata_qc)
            n_samples_total += 1

            sample_out_dir = sdir / OUTPUT_FOLDER_NAME
            run_scanpy_pipeline(
                adata_qc.copy(),
                sample_out_dir,
                analysis_name=f"{group_name}_{sample_label}",
                do_pseudobulk=False,
                batch_key=None,
                integration_method=None,
                do_groupwise_de=False,
                do_dpt=False,
            )

    if not all_adatas_qc:
        raise RuntimeError("No valid samples loaded in multi 10x mode.")
    #logger.info(f"[MULTI-10X] Combining {n_samples_total} filtered samples into one AnnData...")

    # 🔴 CRITICAL FIX: make barcodes unique BEFORE concat
    for adata in all_adatas_qc:
        if "sample" not in adata.obs.columns:
            raise RuntimeError("Sample column missing before combine")

        # store original barcode
        adata.obs["barcode"] = adata.obs_names.astype(str)

        # prefix barcode with sample name
        adata.obs_names = (
            adata.obs["sample"].astype(str) + "_" + adata.obs["barcode"].astype(str)
        )

        # safety
        adata.obs_names_make_unique()

    # ✅ NOW safe to concatenate
    adata_all = ad.concat(
        all_adatas_qc,
        join="outer",
        fill_value=0,
    )

    adata_all.var_names_make_unique()

    logger.info(
        f"[MULTI-10X] Combined shape: {adata_all.n_obs} cells x {adata_all.n_vars} genes"
    )

    # logger.info(f"[MULTI-10X] Combining {n_samples_total} filtered samples into one AnnData...")
    # adata_all = ad.concat(
    #     all_adatas_qc,
    #     join="outer",
    #     label=None,
    #     index_unique=None,
    #     fill_value=0,
    # )
    # adata_all.var_names_make_unique()
    # logger.info(f"[MULTI-10X] Combined shape: {adata_all.n_obs} cells x {adata_all.n_vars} genes")

    run_scanpy_pipeline(
        adata_all,
        COMBINED_OUT_DIR,
        analysis_name="combined_all_samples",
        do_pseudobulk=True,
        batch_key="sample",
        integration_method="bbknn",
        do_groupwise_de=True,
        group_col="group",     # e.g. Pre vs Post
        cluster_col="leiden",
        do_dpt=False,
    )

logger.info("DONE — Full single-cell pipeline (10x-only, CAR-T aware, per-sample + combined) finished with structured outputs.")
# =========================
# 11. MULTI/SINGLE ENTRYPOINT
# =========================

def _guess_10x_dirs(base_dir: Path) -> list[Path]:
    """
    Returns a list of folders that look like a 10x folder (contain matrix + barcodes + features/genes).
    This is deliberately permissive.
    """
    base_dir = Path(base_dir)
    candidates = []
    for d in base_dir.rglob("*"):
        if not d.is_dir():
            continue
        m = _find_first_matching(d, ["matrix.mtx", "matrix.mtx.gz", "*.mtx", "*.mtx.gz"])
        b = _find_first_matching(d, ["barcodes.tsv", "barcodes.tsv.gz", "*barcodes.tsv", "*barcodes.tsv.gz"])
        f = _find_first_matching(d, ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"])
        if m and b and f:
            candidates.append(d)
    # dedup + stable order
    candidates = sorted(set(candidates))
    return candidates


def _load_metadata_xlsx(meta_path: Path) -> pd.DataFrame:
    meta = pd.read_excel(meta_path)
    meta.columns = [str(c).strip() for c in meta.columns]
    return meta


def _attach_metadata(adata: ad.AnnData, meta_row: pd.Series):
    """
    Attach known metadata fields to adata.obs. All stored as strings.
    """
    for key in ["sample", "group", "patient_id", "cart_phase"]:
        if key in meta_row.index and pd.notnull(meta_row[key]):
            adata.obs[key] = str(meta_row[key])


def main():
    logger.info("=== PIPELINE MAIN START ===")

    if SINGLE_MODE:
        logger.info("[MODE] SINGLE_MODE")
        adata = load_10x_feature_barcode_matrix(SINGLE_10X_DIR)
        adata.obs["sample"] = SINGLE_SAMPLE_LABEL
        adata.obs["group"] = SINGLE_GROUP_LABEL

        # quick per-sample QC (optional)
        adata = simple_qc_filter(adata)

        run_scanpy_pipeline(
            adata=adata,
            out_dir=COMBINED_OUT_DIR / "single_analysis",
            analysis_name="single_analysis",
            do_pseudobulk=True,
            do_groupwise_de=False,   # single group usually
            group_col="group",
            batch_key=None,
            integration_method=None,
        )

    else:
        logger.info("[MODE] MULTI_MODE")
        meta_path = MULTI_BASE_DIR / MULTI_META_FILENAME
        if not meta_path.exists():
            raise FileNotFoundError(f"metadata.xlsx not found: {meta_path}")

        meta = _load_metadata_xlsx(meta_path)

        # REQUIRED: a column that identifies which 10x folder belongs to this row.
        # Common names you might have: "tenx_dir", "folder", "path", "sample_dir"
        candidate_path_cols = ["tenx_dir", "folder", "path", "sample_dir", "dir"]
        path_col = None
        for c in candidate_path_cols:
            if c in meta.columns:
                path_col = c
                break

        adatas = []

        if path_col is not None:
            logger.info(f"[MULTI] Using metadata path column: {path_col}")

            for _, row in meta.iterrows():
                tenx_dir = Path(row[path_col])
                if not tenx_dir.is_absolute():
                    tenx_dir = (MULTI_BASE_DIR / tenx_dir).resolve()

                if not tenx_dir.exists():
                    logger.warning(f"[MULTI] 10x dir missing, skipping: {tenx_dir}")
                    continue

                a = load_10x_feature_barcode_matrix(tenx_dir)

                # attach metadata columns (sample/group/patient_id/cart_phase if present)
                _attach_metadata(a, row)

                # fallbacks if metadata lacks these
                if "sample" not in a.obs.columns:
                    a.obs["sample"] = tenx_dir.name
                if "group" not in a.obs.columns:
                    a.obs["group"] = "UNKNOWN"

                # pre-QC per sample
                a = simple_qc_filter(a)

                adatas.append(a)

        else:
            logger.info("[MULTI] metadata.xlsx has no path column; auto-discovering 10x dirs...")
            tenx_dirs = _guess_10x_dirs(MULTI_BASE_DIR)
            logger.info(f"[MULTI] Discovered {len(tenx_dirs)} 10x directories.")
            for tenx_dir in tenx_dirs:
                a = load_10x_feature_barcode_matrix(tenx_dir)
                a.obs["sample"] = tenx_dir.name
                a.obs["group"] = "UNKNOWN"
                a = simple_qc_filter(a)
                adatas.append(a)

        if not adatas:
            raise RuntimeError("No samples loaded in MULTI_MODE.")

        logger.info(f"[MULTI] Concatenating {len(adatas)} samples...")
        adata_combined = ad.concat(
            adatas,
            join="outer",
            label="sample",
            keys=[a.obs["sample"].iloc[0] for a in adatas],
            fill_value=0,
            index_unique="-",
        )

        # IMPORTANT: make sure group exists and is string/categorical
        if "group" in adata_combined.obs.columns:
            adata_combined.obs["group"] = adata_combined.obs["group"].astype(str)

        run_scanpy_pipeline(
            adata=adata_combined,
            out_dir=COMBINED_OUT_DIR / "combined_all_samples",
            analysis_name="combined_all_samples",
            do_pseudobulk=True,
            batch_key="sample",
            integration_method="bbknn",     # set None if you don't want BBKNN
            do_groupwise_de=True,           # only works if >=2 groups exist
            group_col="group",
        )

    logger.info("=== PIPELINE MAIN END ===")


if __name__ == "__main__":
    main()
#with car-t
#oncocyrix-multicohort --mode multi --multi-base-dir "C:\Users\shery\Downloads\oncocyrix_multicohort\GSE208653_RAW-v2" --do-cart

#without car-t
#oncocyrix-multicohort --mode multi --multi-base-dir "C:\Users\shery\Downloads\oncocyrix_multicohort\GSE208653_RAW-v2"