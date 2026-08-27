#!/usr/bin/env python3
"""Render the combined dataset comparison CSV as include and standalone TeX."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
import pandas as pd


def esc(text: str) -> str:
    return (str(text).replace("\\", r"\textbackslash{}")
            .replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
            .replace("#", r"\#").replace("$", r"\$").replace("{", r"\{")
            .replace("}", r"\}"))


def fmt(value) -> str:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return r"\textemdash{}"
    s = str(value)
    if s.lower() in {"true", "false"}:
        return s.capitalize()
    try:
        x = float(s)
        return f"{int(x):,}" if x.is_integer() else f"{x:,.3g}"
    except ValueError:
        return esc(s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--tables-dir", required=True, type=Path)
    args = ap.parse_args()
    args.tables_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, keep_default_na=False)
    datasets = [c for c in df.columns if c not in {"section", "property"}]
    header = " & ".join(["Property"] + [r"\textbf{" + esc(x) + "}" if x == "Ours" else esc(x) for x in datasets]) + r" \\"
    rows = []
    previous = None
    for r in df.itertuples(index=False):
        if r.section != previous:
            if previous is not None:
                rows.append(r"\addlinespace")
            rows.append(r"\multicolumn{" + str(len(datasets) + 1) + r"}{l}{\textit{" + esc(r.section) + r"}} \\")
            previous = r.section
        values = []
        for ds in datasets:
            value = fmt(getattr(r, ds))
            values.append(r"\textbf{" + value + "}" if ds == "Ours" else value)
        rows.append(esc(r.property) + " & " + " & ".join(values) + r" \\")
    tabular = "\n".join([
        r"\begin{table*}[t]", r"\centering", r"\small",
        r"\begin{threeparttable}",
        r"\caption{Auditable comparison of the current formal database release with external reference datasets.}",
        r"\label{tab:dataset-comparison}",
        r"\begin{adjustbox}{max width=\textwidth}",
        r"\begin{tabular}{l" + "r" * len(datasets) + "}",
        r"\toprule", header, r"\midrule", *rows, r"\bottomrule",
        r"\end{tabular}", r"\end{adjustbox}",
        r"\begin{tablenotes}[flushleft]\footnotesize",
        r"\item Ours is bold and denotes 91,860 unique Filter 4 PASS pairs in the explicit frozen final database-construction release.",
        r"\item \textemdash{} denotes not available from a confirmed source in the current repository; values were not guessed.",
        r"\item Definitions and denominators are recorded in \texttt{qc/metric\_definitions.tsv}.",
        r"\end{tablenotes}", r"\end{threeparttable}", r"\end{table*}", "",
    ])
    (args.tables_dir / "dataset_comparison_table.tex").write_text(tabular, encoding="utf-8", newline="\n")
    standalone = "\n".join([
        r"\documentclass[10pt]{article}", r"\usepackage[margin=0.5in,landscape]{geometry}",
        r"\usepackage{booktabs}", r"\usepackage{threeparttable}", r"\usepackage{adjustbox}",
        r"\begin{document}", tabular, r"\end{document}", "",
    ])
    (args.tables_dir / "dataset_comparison_table_standalone.tex").write_text(standalone, encoding="utf-8", newline="\n")
    qc_dir = args.tables_dir.parent / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    latexmk, pdflatex = shutil.which("latexmk"), shutil.which("pdflatex")
    compiler, status, detail = "", "NOT_RUN", "No LaTeX compiler found on PATH."
    if latexmk:
        compiler = latexmk
        proc = subprocess.run([latexmk, "-pdf", "-interaction=nonstopmode", "dataset_comparison_table_standalone.tex"],
                              cwd=args.tables_dir, text=True, capture_output=True)
        status, detail = ("PASS", "latexmk generated the standalone PDF.") if proc.returncode == 0 else ("FAIL", proc.stdout[-1000:] + proc.stderr[-1000:])
    elif pdflatex:
        compiler = pdflatex
        proc = subprocess.run([pdflatex, "-interaction=nonstopmode", "-halt-on-error", "dataset_comparison_table_standalone.tex"],
                              cwd=args.tables_dir, text=True, capture_output=True)
        status, detail = ("PASS", "pdflatex generated the standalone PDF.") if proc.returncode == 0 else ("FAIL", proc.stdout[-1000:] + proc.stderr[-1000:])
    (qc_dir / "latex_build_status.tsv").write_text(
        "status\tcompiler\tdetail\n" + f"{status}\t{compiler}\t{detail.replace(chr(9), ' ').replace(chr(10), ' ')}\n",
        encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
