#!/usr/bin/env python3
"""Rebuild three EMT final compositions from archived figure data and vector panels.

Historical upstream scripts are provided in each Panel/code folder, but this
entry point rebuilds the published artwork from frozen derived table inputs.
"""
import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TASKS = {
    "Figure_6": {
        "entry": "V20_多库通路组合复核_20261027/code/build_combined_figures.py",
        "result": "V20_多库通路组合复核_20261027/final_combination/Figure_6_combined.pdf",
        "reference": "Figure_6/final/Figure_6_combined.pdf",
        "filename": "Figure_6_combined.pdf"
    },
    "Figure_7": {
        "entry": "V6_最终图_20261024/code/build_figure7.py",
        "result": "JBC_singlecell_20261024_final/Figure_7.pdf",
        "reference": "Figure_7/final/Figure_7.pdf",
        "filename": "Figure_7.pdf"
    },
    "Figure_S1": {
        "entry": "V14_最终图与回复_20261026/code/build_final_figures.py",
        "result": "V14_最终图与回复_20261026/figures/Figure_S1.pdf",
        "reference": "Figure_S1/final/Figure_S1.pdf",
        "filename": "Figure_S1.pdf"
    }
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", choices=["all"] + list(TASKS), default="all")
    ap.add_argument("--output", type=Path, help="New output directory; must not yet exist")
    ap.add_argument("--font-dir", type=Path,
                    default=Path("/System/Library/Fonts/Supplemental"),
                    help="Directory with installed Arial.ttf and Arial Bold.ttf")
    args = ap.parse_args()
    try:
        import fitz
    except ImportError:
        raise SystemExit("Missing PyMuPDF: pip install -r environment/requirements_plotting.txt")
    for nm in ["Arial.ttf", "Arial Bold.ttf"]:
        if not (args.font_dir / nm).is_file():
            raise SystemExit("Font not found: " + str(args.font_dir / nm))
    out = (args.output or ROOT / "outputs" / (
        "rebuild_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") +
        "_" + uuid.uuid4().hex[:6])).resolve()
    if out.exists():
        raise SystemExit("Refusing to overwrite existing directory: " + str(out))
    out.mkdir(parents=True)
    mapping = json.loads((ROOT / "runtime_map.json").read_text())
    selected = list(TASKS) if args.figure == "all" else [args.figure]
    report = []
    for figure in selected:
        task = TASKS[figure]
        with tempfile.TemporaryDirectory(prefix="emt_figure_") as td:
            project = Path(td) / "03_EMT"
            work = project / "figures/EMT_Figure_Final/Weng_Figure"
            for item in mapping:
                name = item["repository_path"]
                if not (name.startswith(figure + "/") or name.startswith("shared/")):
                    continue
                source = ROOT / name
                target = work / item["historical_weng_path"]
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
            builder = work / task["entry"]
            if not builder.is_file():
                raise RuntimeError("Missing assembly code: " + str(builder))
            default_font_path = "/System/Library/Fonts/Supplemental"
            if str(args.font_dir.resolve()) != default_font_path:
                builder.write_text(builder.read_text().replace(
                    default_font_path, str(args.font_dir.resolve())))
            logfile = out / (figure + ".log")
            with logfile.open("w") as stream:
                run = subprocess.run([sys.executable, str(builder)],
                                     cwd=project, stdout=stream,
                                     stderr=subprocess.STDOUT)
            if run.returncode:
                print(logfile.read_text()[-4500:], file=sys.stderr)
                raise RuntimeError("Assembly failed; see " + str(logfile))
            rebuilt = work / task["result"]
            if not rebuilt.is_file():
                raise RuntimeError("Expected PDF was not generated: " + str(rebuilt))
            output = out / task["filename"]
            shutil.copy2(rebuilt, output)
            approved = ROOT / task["reference"]
            doc_a = fitz.open(approved)
            doc_b = fitz.open(output)
            if len(doc_a) != 1 or len(doc_b) != 1:
                raise RuntimeError("Expected a one-page PDF")
            p = doc_a[0].get_pixmap(dpi=144, alpha=False)
            q = doc_b[0].get_pixmap(dpi=144, alpha=False)
            same = (p.width, p.height, p.samples) == (q.width, q.height, q.samples)
            dimensions_equal = doc_a[0].rect == doc_b[0].rect
            doc_b[0].get_pixmap(dpi=100, alpha=False).save(output.with_suffix(".png"))
            report.append({
                "figure": figure,
                "approved_sha256": sha(approved),
                "rebuilt_sha256": sha(output),
                "pixel_identical_144dpi": same,
                "dimensions_equal": dimensions_equal,
                "scope": "Final assembly from frozen derived tables and vector components"
            })
            doc_a.close()
            doc_b.close()
            print(figure + ": " + ("PASS" if same and dimensions_equal else "MISMATCH"),
                  flush=True)
    (out / "rebuild_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    if not all(r["pixel_identical_144dpi"] and r["dimensions_equal"] for r in report):
        raise SystemExit("Rendering comparison failed: " + str(out))
    print("PASS: all selected figures match approved PDFs at 144 dpi")
    print("Output:", out)

if __name__ == "__main__":
    main()
