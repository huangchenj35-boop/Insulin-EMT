#!/usr/bin/env python3
"""Rebuild uncapped PAGA v1.2 ratios from the frozen balanced samples.

Use exact Euclidean kNN for these small graphs, and the directed-edge-count
null model in Scanpy PAGA v1.2. Re-capping must reproduce EVERY stored result.
No uncapped value is inferred from a capped value.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
TABLES = PROJECT / "emt_invivo_bodyfig_20260901/tables"
OUT = PROJECT / "EMT_Resutls/figures/Balanced_resampling_HIJ/uncapped_data"
SEED = 20260901
COHORTS = {"ChinaHCC": "China-HCC", "USAHCC": "USA-HCC", "CCAcombined": "CCA-combined"}
ROLE_ORDER = ("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root",
              "qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root",
              "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
HASH = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_saved_full_uncapped(tag):
    """Independently validate above-one values against saved full k30 graphs."""
    frame = pd.read_csv(TABLES / f"S3_{tag}_PAGA_Harmony_input.csv")
    reference_path = TABLES / f"S3_{tag}_PAGA_connectivity_matrix_k30_uncapped.csv"
    reference = pd.read_csv(reference_path, index_col=0)
    roles = reference.columns.tolist()
    labels = np.array([roles.index(r) for r in frame.frozen_role])
    coords = frame[[c for c in frame if c.startswith("Harmony_")]].to_numpy(dtype=np.float32).astype(np.float64)
    delta = coords[:, None, :]-coords[None, :, :]
    dist2 = np.einsum("ijk,ijk->ij", delta, delta)
    np.fill_diagonal(dist2, np.inf)
    order = np.argsort(dist2, axis=1, kind="stable")
    neighbors = order[:, :29]
    boundary = np.take_along_axis(dist2, order[:, 28:30], axis=1)
    assert not np.any(boundary[:, 0] == boundary[:, 1]) and not np.any(dist2 == 0)
    directed = np.zeros((len(roles), len(roles)), dtype=np.int64)
    np.add.at(directed, (np.repeat(labels, 29), labels[neighbors].ravel()), 1)
    sizes = np.bincount(labels, minlength=len(roles))
    outgoing = directed.sum(axis=1)
    expected = (outgoing[:,None]*sizes[None,:]+outgoing[None,:]*sizes[:,None])/(len(frame)-1)
    observed = directed+directed.T
    np.fill_diagonal(observed, 0)
    ratio = np.divide(observed, expected, out=np.zeros_like(expected), where=expected>0)
    error = float(np.max(np.abs(ratio-reference.loc[roles,roles].to_numpy())))
    assert error < 1e-12, (tag, error)
    return {"k":30,"max_error":error,"reference_sha256":HASH(reference_path),
        "verified_unique_edges":len(roles)*(len(roles)-1)//2,"above_one_unique_edges":int(np.sum(np.triu(ratio,1)>1))}


def process(tag, cohort, cohort_index):
    path = TABLES / f"S3_{tag}_PAGA_Harmony_input.csv"
    oldpath = TABLES / f"S4_{tag}_balanced_resampling_by_iteration.csv"
    frame = pd.read_csv(path)
    if tag == "USAHCC":
        frame = frame.loc[frame.frozen_role != "normal_cholangiocyte_root"].copy()
    roles = [r for r in ROLE_ORDER if r in set(frame.frozen_role)]
    assert len(roles) == [7, 5, 5][cohort_index]
    assert frame.barcode.is_unique and set(frame.cohort) == {cohort}
    harmony = [c for c in frame if c.startswith("Harmony_")]
    assert len(harmony) == 30
    # Preserve original float32 representation; use float64 accumulation for exact distances.
    coords = frame[harmony].to_numpy(dtype=np.float32).astype(np.float64)
    role_indices = {r: np.flatnonzero(frame.frozen_role.to_numpy() == r) for r in roles}
    n_per = min(len(v) for v in role_indices.values())
    assert n_per == [29, 59, 16][cohort_index]
    labels = np.repeat(np.arange(len(roles)), n_per)
    sizes = np.bincount(labels, minlength=len(roles))
    target = roles.index("fibrotic_MFB_terminal")
    n = len(labels)
    rows, selections = [], []
    ties, duplicates = 0, 0
    start = time.monotonic()
    for replicate in range(1, 201):
        sample_seed = SEED + replicate + 1000 * cohort_index
        rng = np.random.default_rng(sample_seed)
        selected = np.concatenate([rng.choice(role_indices[r], n_per, replace=False) for r in roles])
        barcodes = frame.iloc[selected].barcode.astype(str).tolist()
        selections.append({"cohort": cohort, "replicate": replicate, "sample_seed": sample_seed,
            "selected_barcodes_sha256": hashlib.sha256("\n".join(barcodes).encode()).hexdigest()})
        x = coords[selected]
        diff = x[:, None, :] - x[None, :, :]
        dist2 = np.einsum("ijk,ijk->ij", diff, diff)
        np.fill_diagonal(dist2, np.inf)
        duplicates += int(np.sum(dist2 == 0))
        ordered = np.argsort(dist2, axis=1, kind="stable")
        for k in (10, 20):
            # Scanpy's n_neighbors includes self; its PAGA distance graph has k-1 nonself edges.
            neighbors = ordered[:, :k-1]
            near_d = np.take_along_axis(dist2, neighbors, axis=1)
            boundary = np.take_along_axis(dist2, ordered[:, k-2:k], axis=1)
            ties += int(np.sum(boundary[:, 0] == boundary[:, 1]))
            from_role = np.repeat(labels, k-1)
            to_role = labels[neighbors].ravel()
            positive = near_d.ravel() > 0
            directed = np.zeros((len(roles), len(roles)), dtype=np.int64)
            np.add.at(directed, (from_role[positive], to_role[positive]), 1)
            outgoing = directed.sum(axis=1)  # internal + inter-population directed edges
            observed = directed + directed.T
            np.fill_diagonal(observed, 0)
            expected = (outgoing[:, None]*sizes[None, :] + outgoing[None, :]*sizes[:, None])/(n-1)
            ratio = np.divide(observed, expected, out=np.zeros_like(expected), where=expected > 0)
            for role_idx, role in enumerate(roles):
                if role_idx == target:
                    continue
                value = float(ratio[role_idx, target])
                rows.append({"cohort": cohort, "replicate": replicate, "n_per_role": n_per,
                    "n_neighbors": k, "edge": f"{role}--MFB",
                    "observed_edges": int(observed[role_idx, target]),
                    "expected_random_edges": float(expected[role_idx, target]),
                    "paga_observed_expected_ratio": value, "recapped_connectivity": min(value, 1.0),
                    "sample_seed": sample_seed, "neighbor_seed_original": SEED+replicate})
    result = pd.DataFrame(rows)
    saved = pd.read_csv(oldpath)
    on = ["cohort", "replicate", "n_per_role", "n_neighbors", "edge"]
    check = result.merge(saved[on+["paga_connectivity"]], on=on, validate="one_to_one")
    assert len(check) == len(saved) == len(result)
    check["capped_absolute_error"] = abs(check.recapped_connectivity-check.paga_connectivity)
    err = float(check.capped_absolute_error.max())
    if err > 1e-12 or ties or duplicates:
        OUT.mkdir(parents=True, exist_ok=True)
        check.to_csv(OUT/f"FAILED_{tag}_reconstruction_check.csv", index=False)
        raise RuntimeError(f"Do not use reconstructed output: {tag}: error={err}, ties={ties}, duplicate_coords={duplicates}")
    full_graph_check = verify_saved_full_uncapped(tag)
    audit = {"cohort": cohort, "input_rows": len(frame), "roles_including_MFB": roles,
        "n_per_role": n_per, "repeats_per_k": 200, "k_values": [10,20], "verified_values": len(check),
        "max_recapping_error": err, "knn_boundary_ties": ties, "duplicate_coordinates": duplicates,
        "above_one": int((result.paga_observed_expected_ratio>1).sum()),
        "max_uncapped_ratio": float(result.paga_observed_expected_ratio.max()),
        "harmony_input_sha256": HASH(path), "capped_reference_sha256": HASH(oldpath),
        "full_uncapped_k30_validation":full_graph_check, "elapsed_seconds": time.monotonic()-start}
    print(json.dumps(audit), flush=True)
    return result, pd.DataFrame(selections), check, audit


if __name__ == "__main__":
    batches = [process(tag, cohort, i) for i, (tag, cohort) in enumerate(COHORTS.items())]
    OUT.mkdir(parents=True, exist_ok=True)
    for (tag, _), batch in zip(COHORTS.items(), batches):
        batch[0].to_csv(OUT/f"S4_{tag}_balanced_resampling_uncapped_by_iteration.csv", index=False)
    pd.concat([b[1] for b in batches]).to_csv(OUT/"sample_selection_sha256.csv", index=False)
    pd.concat([b[2] for b in batches]).to_csv(OUT/"all_5600_recapping_checks.csv", index=False)
    provenance = {"method": "Exact Euclidean kNN on original float32 Harmony coordinates; Scanpy PAGA v1.2 directed edge-count null model without the min(value,1) step",
        "local_method_reference": "emt_invivo_bodyfig_20260901/scripts/03_paga_connectivity_uncapped.py",
        "official_formula_reference": "https://github.com/scverse/scanpy/blob/main/src/scanpy/tools/_paga.py",
        "knn_note": "n_neighbors includes self; k-1 nonself neighbors; exact distance order with zero boundary ties",
        "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
        "script_sha256": HASH(__file__), "cohorts": [b[3] for b in batches]}
    (OUT/"uncapped_reconstruction_audit.json").write_text(json.dumps(provenance,indent=2)+"\n")
    print("PASS: all 5600 re-capped values reproduce the original tables within 1e-12.",flush=True)
