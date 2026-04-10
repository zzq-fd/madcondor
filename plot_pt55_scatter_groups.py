#!/usr/bin/env python3
"""Create scatter plots from plot_pt55.py summary table.

This script reads the table produced by plot_pt55.py (pt55_mean_table.txt).
If there are >=3 entries where three of (sintheta, tanbeta, mA, ma) are the same
and the 4th is different, it makes a scatter plot:
  x-axis   = the varying (4th) parameter
  y-axis   = mean_pt

The plot title and filename include the fixed 3 parameters.
All PNGs are saved into the directory containing the input txt file.

Usage:
    python3 plot_pt55_scatter_groups.py <OUTROOT>

Example:
    python3 plot_pt55_scatter_groups.py bbdm_2HDMa_type1_case1_scan

Backward-compatible:
    You can also pass the table file path directly.

Notes:
  - Uses CERN ROOT via the `root` executable (no PyROOT required).
  - Ignores rows where mean_pt is missing (N/A).
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Row:
    sintheta: float
    tanbeta: float
    mA: float
    ma: float
    mchi: float
    cosbma: float
    deltam: float
    mean_pt: float


def _fmt_float_tag(x: float) -> str:
    # filename-friendly tag: 0.9 -> 0p9, 10.0 -> 10
    s = f"{x:g}"
    return s.replace(".", "p")


def _fmt_title_float(x: float) -> str:
    return f"{x:g}"


def read_table(path: Path) -> list[Row]:
    rows: list[Row] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        # Supported formats:
        #   v1 (old): sintheta tanbeta mA ma mean_pt entries dirname
        #   v2 (new): sintheta tanbeta mA ma mchi cosbma deltam mean_pt entries dirname
        if len(parts) < 6:
            continue

        if len(parts) >= 9:
            st_s, tb_s, mA_s, ma_s, mchi_s, cb_s, dm_s, mean_s = (
                parts[0],
                parts[1],
                parts[2],
                parts[3],
                parts[4],
                parts[5],
                parts[6],
                parts[7],
            )
        else:
            st_s, tb_s, mA_s, ma_s, mean_s = parts[0], parts[1], parts[2], parts[3], parts[4]
            # Legacy table had no mchi/cosbma/deltam columns.
            # Assume 0.0 so older outputs remain plottable.
            mchi_s, cb_s, dm_s = "0", "0", "0"

        if (
            mean_s == "N/A"
            or st_s == "N/A"
            or tb_s == "N/A"
            or mA_s == "N/A"
            or ma_s == "N/A"
            or mchi_s == "N/A"
            or cb_s == "N/A"
            or dm_s == "N/A"
        ):
            continue

        try:
            rows.append(
                Row(
                    sintheta=float(st_s),
                    tanbeta=float(tb_s),
                    mA=float(mA_s),
                    ma=float(ma_s),
                    mchi=float(mchi_s),
                    cosbma=float(cb_s),
                    deltam=float(dm_s),
                    mean_pt=float(mean_s),
                )
            )
        except Exception:
            continue

    return rows


def write_root_macro(macro_path: Path) -> None:
    macro_path.write_text(
        r'''#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TROOT.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int make_scatter(const char* data_txt,
                 const char* x_label,
                 const char* plot_title,
                 const char* out_png)
{
  gROOT->SetBatch(kTRUE);

  std::ifstream in(data_txt);
  if (!in) {
    std::cerr << "ERROR: cannot open data file: " << data_txt << std::endl;
    return 2;
  }

  std::vector<double> xs;
  std::vector<double> ys;
  double x, y;
  while (in >> x >> y) {
    xs.push_back(x);
    ys.push_back(y);
  }

  if (xs.size() < 3) {
    std::cerr << "ERROR: need >=3 points, got " << xs.size() << std::endl;
    return 3;
  }

  TGraph g((int)xs.size());
  for (int i = 0; i < (int)xs.size(); ++i) {
    g.SetPoint(i, xs[i], ys[i]);
  }

  g.SetTitle(plot_title);
  g.GetXaxis()->SetTitle(x_label);
    g.GetYaxis()->SetTitle("mean(gen-p_{T}^{miss}=|#vec{p}_{T}(55)+#vec{p}_{T}(36)|) [GeV]");
    g.GetYaxis()->SetTitleOffset(1.25);
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.0);

    TCanvas c("c", "c", 1100, 800);
    c.SetLeftMargin(0.14);
    c.SetBottomMargin(0.12);
  g.Draw("AP");
  c.SaveAs(out_png);
  return 0;
}
'''
    )


def root_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def _in_x_range(x: float, xmin: float | None, xmax: float | None) -> bool:
    if xmin is not None and x < xmin:
        return False
    if xmax is not None and x > xmax:
        return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True, allow_abbrev=False)
    ap.add_argument("input", help="<OUTROOT> or direct path to pt55_mean_table.txt")
    ap.add_argument("--xmin", type=float, default=None, help="Optional x-axis lower bound")
    ap.add_argument("--xmax", type=float, default=None, help="Optional x-axis upper bound")
    args = ap.parse_args()
    if args.xmin is not None and args.xmax is not None and args.xmax < args.xmin:
        raise SystemExit("ERROR: --xmax must be >= --xmin")

    arg = Path(args.input).expanduser()
    if arg.is_dir():
        out_dir = arg
        table_path = out_dir / "pt55_mean_table.txt"
    else:
        table_path = arg
        out_dir = table_path.parent

    if not table_path.exists():
        print(f"ERROR: pt55 mean table not found: {table_path}")
        print("Hint: run plot_pt55.py first to generate pt55_mean_table.txt")
        return 2
    # Output directory layout:
    #   <OUTROOT>/output_scatter_plot/mean_met/
    # Create intermediate directories if missing.
    plots_dir = out_dir / "output_scatter_plot" / "mean_met"
    plots_dir.mkdir(parents=True, exist_ok=True)

    root_exe = shutil.which("root")
    if not root_exe:
        print("ERROR: CERN ROOT executable 'root' not found on PATH.")
        print("Hint: if you use conda, run: conda activate mg5")
        return 2

    rows = read_table(table_path)
    if not rows:
        print("No valid rows found (mean_pt may be missing).")
        return 0

    macro_path = plots_dir / ".pt55_scatter_tmp.C"
    write_root_macro(macro_path)

    params = ["sintheta", "tanbeta", "mA", "ma", "mchi", "cosbma", "deltam"]

    # group[(varying_param, fixed_tuple)] = {x_value: [mean_pt values...]}
    # where fixed_tuple contains all other params in the order of fixed_names.
    groups: dict[tuple[str, tuple[float, ...]], dict[float, list[float]]] = {}

    for r in rows:
        values = {
            "sintheta": r.sintheta,
            "tanbeta": r.tanbeta,
            "mA": r.mA,
            "ma": r.ma,
            "mchi": r.mchi,
            "cosbma": r.cosbma,
            "deltam": r.deltam,
        }

        for vary in params:
            fixed_names = [p for p in params if p != vary]
            fixed_tuple = tuple(values[n] for n in fixed_names)
            key = (vary, fixed_tuple)
            x = values[vary]
            groups.setdefault(key, {}).setdefault(x, []).append(r.mean_pt)

    made = 0
    any_failures = False

    for (vary, fixed_tuple), x_to_ys in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        # Need >=3 distinct x values
        if len(x_to_ys) < 3:
            continue

        # Average duplicates at same x (in case you ran the same point multiple times)
        points = [(x, sum(ys) / len(ys)) for x, ys in x_to_ys.items()]
        points = [(x, y) for x, y in points if _in_x_range(x, args.xmin, args.xmax)]
        if len(points) < 3:
            continue
        points.sort(key=lambda t: t[0])

        fixed_names = [p for p in params if p != vary]
        fixed_map = dict(zip(fixed_names, fixed_tuple))

        # Build labels and filename using ONLY the fixed 3 parameters.
        fixed_label = ", ".join(
            f"{k}={_fmt_title_float(fixed_map[k])}" for k in fixed_names
        )
        title = f"mean(gen-pTmiss=|pT(55)+pT(36)|) vs {vary}  ({fixed_label})"

        # Axis labels
        if vary in ("mA", "ma", "mchi", "deltam"):
            x_label = "m_{#chi} [GeV]" if vary == "mchi" else f"{vary} [GeV]"
        else:
            x_label = vary

        def _tag(k: str) -> str:
            return {
                "sintheta": "st",
                "tanbeta": "tb",
                "mA": "mA",
                "ma": "ma",
                "mchi": "mchi",
                "cosbma": "cbma",
                "deltam": "dm",
            }[k]

        fixed_fname = "_".join(
            f"{_tag(k)}{_fmt_float_tag(fixed_map[k])}" for k in fixed_names
        )
        fname = f"scatter_{vary}__{fixed_fname}.png"
        out_png = plots_dir / fname

        data_txt = plots_dir / f".pt55_scatter_{vary}_{hash((vary, fixed_tuple)) & 0xffffffff:x}.dat"
        with open(data_txt, "w") as f:
            for x, y in points:
                f.write(f"{x} {y}\n")

        stmt = (
            f'gROOT->LoadMacro("{root_quote(macro_path)}+"); '
            f'make_scatter("{root_quote(data_txt)}","{x_label}","{title}","{root_quote(out_png)}");'
        )
        cmd = [root_exe, "-l", "-b", "-q", "-e", stmt]

        log_path = plots_dir / (out_png.stem + ".log")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log_path.write_text(proc.stdout)

        if proc.returncode != 0:
            any_failures = True
            print(f"Warning: failed to make plot {out_png.name} (see {log_path.name})")
            continue

        # Clean up per-plot temp/log files on success to avoid clutter.
        try:
            data_txt.unlink(missing_ok=True)
        except Exception:
            pass
        try:
            log_path.unlink(missing_ok=True)
        except Exception:
            pass

        made += 1

    # Clean up macro + compilation artifacts on success.
    if not any_failures:
        try:
            macro_path.unlink(missing_ok=True)
        except Exception:
            pass

        compiled_base = macro_path.name.replace(".C", "_C")
        for p in plots_dir.glob(f"{compiled_base}*"):
            try:
                if p.is_file():
                    p.unlink(missing_ok=True)
            except Exception:
                pass

    print(f"Wrote {made} scatter plot(s) into: {plots_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
