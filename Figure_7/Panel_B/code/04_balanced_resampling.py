#!/usr/bin/env python3
"""Balanced, undirected PAGA resampling for the three frozen cohorts."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/emt_paga_numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/emt_paga_mpl")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("emt_invivo_bodyfig_20260901")
TABLE_DIR = ROOT / "tables"
FIGURE_DIR = ROOT / "figures"
SEED = 20260901
N_REPLICATES = int(os.environ.get("EMT_PAGA_N_REPLICATES", "200"))
NEIGHBOR_VALUES = (10, 20)
COHORTS = {"ChinaHCC": "China-HCC", "USAHCC": "USA-HCC", "CCAcombined": "CCA-combined"}
ROLE_ORDER = (
    "normal_cholangiocyte_root", "normal_hepatocyte_root",
    "malignant_epithelial_EMT_root", "qHSC_root",
    "activated_HSC_intermediate", "portal_fibroblast_root",
    "portal_fibroblast_like_root", "fibrotic_MFB_terminal",
)
LABEL = {
    "normal_hepatocyte_root": "Normal hepatocyte",
    "normal_cholangiocyte_root": "Normal cholangiocyte",
    "malignant_epithelial_EMT_root": "Malignant epithelial",
    "qHSC_root": "qHSC", "activated_HSC_intermediate": "Activated HSC",
    "portal_fibroblast_root": "Portal fibroblast",
    "portal_fibroblast_like_root": "Portal fibroblast-like",
    "fibrotic_MFB_terminal": "Fibrotic MFB",
}
COLOR = {
    "normal_hepatocyte_root": "#9E9E9E", "normal_cholangiocyte_root": "#C9B47C",
    "malignant_epithelial_EMT_root": "#E08214", "qHSC_root": "#80C7D8",
    "activated_HSC_intermediate": "#F2A07B", "portal_fibroblast_root": "#43A77B",
    "portal_fibroblast_like_root": "#43A77B",
    "fibrotic_MFB_terminal": "#B65C85",
}


def make_adata(frame: pd.DataFrame, roles: list[str]) -> ad.AnnData:
    harmony = [c for c in frame.columns if c.startswith("Harmony_")]
    if len(harmony) != 30:
        raise ValueError("Expected 30 Harmony dimensions")
    obs = frame[["barcode", "frozen_role"]].copy()
    obs.index = obs.pop("barcode").astype(str)
    obs["frozen_role"] = pd.Categorical(obs["frozen_role"], categories=roles, ordered=True)
    obj = ad.AnnData(X=np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
    obj.obsm["X_harmony"] = frame[harmony].to_numpy(dtype=np.float32)
    return obj


def paga_matrix(frame: pd.DataFrame, roles: list[str], k: int, seed: int) -> np.ndarray:
    obj = make_adata(frame, roles)
    sc.pp.neighbors(obj, n_neighbors=k, use_rep="X_harmony", metric="euclidean",
                    method="umap", random_state=seed)
    sc.tl.paga(obj, groups="frozen_role", model="v1.2")
    return obj.uns["paga"]["connectivities"].toarray()


def summarize(values: pd.Series) -> pd.Series:
    return pd.Series({"mean": values.mean(), "median": values.median(),
                      "q025": values.quantile(0.025), "q975": values.quantile(0.975)})


def plot_distribution(results: pd.DataFrame, cohort: str, k: int, stem: str) -> None:
    selected = results[(results.cohort == cohort) & (results.n_neighbors == k)]
    edges = list(selected.edge.drop_duplicates())
    positions = np.arange(1, len(edges) + 1)
    fig, ax = plt.subplots(figsize=(max(6.6, 1.35 * len(edges)), 5.6), facecolor="white")
    arrays = [selected.loc[selected.edge == edge, "paga_connectivity"].to_numpy() for edge in edges]
    violins = ax.violinplot(arrays, positions=positions, widths=0.72,
                            showmeans=False, showmedians=False, showextrema=False)
    for body, edge in zip(violins["bodies"], edges):
        root = edge.split("--")[0]
        body.set_facecolor(COLOR.get(root, "#777777")); body.set_edgecolor("#444444")
        body.set_linewidth(0.6); body.set_alpha(0.78)
    rng = np.random.default_rng(SEED + k)
    for pos, edge, values in zip(positions, edges, arrays):
        ax.scatter(np.full(len(values), pos) + rng.normal(0, 0.04, len(values)), values,
                   s=7, alpha=0.16, color="#303030", linewidths=0, zorder=2)
        med = float(np.median(values)); lo, hi = np.quantile(values, [0.025, 0.975])
        ax.errorbar(pos, med, yerr=[[med - lo], [hi - med]], fmt="o", markersize=5,
                    color="#202020", capsize=4, linewidth=1.2, zorder=4)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{LABEL.get(edge.split('--')[0], edge.split('--')[0])}\n- MFB" for edge in edges],
                       fontsize=8.5, rotation=0)
    ax.set_ylabel("PAGA connectivity weight", fontsize=11, fontweight="bold")
    ax.set_ylim(-0.03, 1.0)
    ax.set_title(f"{cohort}: balanced PAGA resampling", loc="left", fontsize=13,
                 fontweight="bold", pad=12)
    ax.text(0, 1.01, f"{int(selected.n_per_role.iloc[0])} cells per role x {N_REPLICATES} resamples; n_neighbors={k}",
            transform=ax.transAxes, fontsize=9, color="#4D4D4D")
    ax.text(0, -0.18, "Points are resamples; bars show medians and 95% resampling intervals.\n"
            "Connectivity is undirected and does not establish lineage origin.",
            transform=ax.transAxes, fontsize=8, color="#555555")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.6); ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.86, bottom=0.25)
    fig.savefig(FIGURE_DIR / f"{stem}_k{k}.pdf", bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}_k{k}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def process(tag: str, cohort: str) -> None:
    frame = pd.read_csv(TABLE_DIR / f"S3_{tag}_PAGA_Harmony_input.csv")
    if tag == "USAHCC":
        frame = frame[frame["frozen_role"] != "normal_cholangiocyte_root"].copy()
    roles = [r for r in ROLE_ORDER if r in set(frame.frozen_role)]
    if "fibrotic_MFB_terminal" not in roles or len(roles) < 3:
        raise ValueError(f"Invalid roles for {cohort}: {roles}")
    counts = frame.frozen_role.value_counts().reindex(roles)
    n_per_role = int(counts.min())
    role_indices = {r: np.flatnonzero(frame.frozen_role.to_numpy() == r) for r in roles}
    rows = []
    mfb = roles.index("fibrotic_MFB_terminal")
    source_roles = [r for r in roles if r != "fibrotic_MFB_terminal"]
    for replicate in range(1, N_REPLICATES + 1):
        rng = np.random.default_rng(SEED + replicate + 1000 * list(COHORTS).index(tag))
        sampled = np.concatenate([rng.choice(role_indices[r], n_per_role, replace=False) for r in roles])
        sampled_frame = frame.iloc[sampled].copy()
        for k in NEIGHBOR_VALUES:
            matrix = paga_matrix(sampled_frame, roles, k, SEED + replicate)
            values = {r: float(matrix[roles.index(r), mfb]) for r in source_roles}
            max_value = max(values.values())
            winners = [r for r, v in values.items() if np.isclose(v, max_value)]
            normal_roles = [r for r in source_roles if r.startswith("normal_")]
            normal_max = max([values[r] for r in normal_roles], default=-np.inf)
            normal_is_max = normal_max >= max_value - 1e-12
            portal_is_max = any(r.startswith("portal_") for r in winners)
            malignant_is_max = "malignant_epithelial_EMT_root" in winners
            for root, value in values.items():
                rows.append({"cohort": cohort, "replicate": replicate, "n_per_role": n_per_role,
                             "n_neighbors": k, "edge": f"{root}--MFB",
                             "paga_connectivity": value, "is_max": root in winners,
                             "normalEpi_is_max": normal_is_max,
                             "portal_is_max": portal_is_max,
                             "malignantEMT_is_max": malignant_is_max})
        if replicate % 20 == 0 or replicate == N_REPLICATES:
            print(f"[progress] {cohort}: {replicate}/{N_REPLICATES}", flush=True)
    results = pd.DataFrame(rows)
    results.to_csv(TABLE_DIR / f"S4_{tag}_balanced_resampling_by_iteration.csv", index=False)
    summary = results.groupby(["cohort", "n_neighbors", "edge"], observed=True)["paga_connectivity"].apply(summarize).unstack().reset_index()
    summary.to_csv(TABLE_DIR / f"S4_{tag}_balanced_resampling_edge_summary.csv", index=False)
    flags = results.drop_duplicates(["replicate", "n_neighbors"])
    advantage = flags.groupby(["cohort", "n_neighbors"], observed=True).agg(
        fraction_normalEpi_is_max=("normalEpi_is_max", "mean"),
        fraction_portal_is_max=("portal_is_max", "mean"),
        fraction_malignantEMT_is_max=("malignantEMT_is_max", "mean"),
    ).reset_index()
    # Add the per-root winning frequency to make ties and normal-vs-portal comparisons explicit.
    wins = results[results.is_max].groupby(["cohort", "n_neighbors", "edge"], observed=True).size().rename("n_wins").reset_index()
    wins["fraction_is_max"] = wins["n_wins"] / N_REPLICATES
    wins.to_csv(TABLE_DIR / f"S4_{tag}_MFB_source_advantage.csv", index=False)
    advantage.to_csv(TABLE_DIR / f"S4_{tag}_normalEpi_advantage.csv", index=False)
    for k in NEIGHBOR_VALUES:
        plot_distribution(results, cohort, k, f"S4_{tag}_balanced_resampling")
    print(f"\n[done] {cohort}; n_per_role={n_per_role}")
    print(summary.round(4).to_string(index=False))
    print(advantage.round(4).to_string(index=False))


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True); FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    selected = COHORTS
    cohort_filter = os.environ.get("EMT_COHORT", "")
    if cohort_filter:
        if cohort_filter not in COHORTS:
            raise ValueError(f"Unknown EMT_COHORT: {cohort_filter}")
        selected = {cohort_filter: COHORTS[cohort_filter]}
    for tag, cohort in selected.items():
        process(tag, cohort)


if __name__ == "__main__":
    main()
