#!/usr/bin/env python3
"""Export uncapped observed/expected PAGA edge ratios.

Scanpy PAGA v1.2 caps the displayed ratio at 1.0.  For the expanded
normal-epithelial graph this hides useful differences between strong edges,
so this companion export retains the raw observed-edge count and its random
expectation ratio while using the identical Harmony representation and
30-neighbor graph as the capped analysis.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/emt_paga_numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/emt_paga_mpl")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scanpy.tools._paga import PAGA

ROOT = Path("emt_invivo_bodyfig_20260901")
TABLE_DIR = ROOT / "tables"
SEED = 20260901
COHORTS = {"ChinaHCC": "China-HCC", "USAHCC": "USA-HCC", "CCAcombined": "CCA-combined"}


def process(tag: str, cohort: str) -> None:
    frame = pd.read_csv(TABLE_DIR / f"S3_{tag}_PAGA_Harmony_input.csv")
    roles = list(frame["frozen_role"].drop_duplicates())
    obs = frame[["barcode", "cohort", "frozen_role"]].copy()
    obs.index = obs.pop("barcode").astype(str)
    obs["frozen_role"] = pd.Categorical(obs["frozen_role"], categories=roles, ordered=True)
    harmony = [c for c in frame.columns if c.startswith("Harmony_")]
    if len(harmony) != 30 or obs.index.duplicated().any() or set(frame["cohort"]) != {cohort}:
        raise ValueError(f"Invalid input for {tag}")
    obj = ad.AnnData(X=np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
    obj.obsm["X_harmony"] = frame[harmony].to_numpy(dtype=np.float32)
    sc.pp.neighbors(obj, n_neighbors=30, use_rep="X_harmony", metric="euclidean",
                    method="umap", random_state=SEED)
    paga = PAGA(obj, groups="frozen_role", model="v1.2")
    observed, capped = paga.compute_connectivities()
    observed = observed.toarray()
    expected = paga.expected_n_edges_random.toarray()
    ratio = np.divide(observed, expected, out=np.zeros_like(observed, dtype=float), where=expected > 0)
    wide = pd.DataFrame(ratio, index=roles, columns=roles)
    wide.index.name = "frozen_role"
    wide.to_csv(TABLE_DIR / f"S3_{tag}_PAGA_connectivity_matrix_k30_uncapped.csv")
    long = pd.DataFrame([
        {"cohort": cohort, "n_neighbors": 30, "role_1": roles[i], "role_2": roles[j],
         "observed_edges": observed[i, j], "expected_random_edges": expected[i, j],
         "paga_observed_expected_ratio": ratio[i, j]}
        for i in range(len(roles)) for j in range(i + 1, len(roles))
    ])
    long.to_csv(TABLE_DIR / f"S3_{tag}_PAGA_connectivity_uncapped_long.csv", index=False)
    print(f"[done] {cohort}: uncapped ratio range {ratio[ratio > 0].min():.3f}-{ratio.max():.3f}")
    print(wide.round(3).to_string())


if __name__ == "__main__":
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
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
