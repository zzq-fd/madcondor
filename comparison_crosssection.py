#!/usr/bin/env python3
"""Compare cross-sections for type1/type2 scans (auto-discovery like plot_crosssection_scatter.py).

This script scans both folders' `results/Events_*` directories, extracts parameters from directory
names, reads cross sections from `*_banner.txt`, and then *automatically* detects which parameter
is varying (per fixed-parameter slice). It overlays type1/type2 on the same scatter plot.

If you have repeated runs for the same parameter point, you can store them in directories with an
extra suffix like `_job123` (written by runMadscan.sh when using repeated submissions). The script
will aggregate repeats at the same parameter point and plot the mean cross section with a
statistical error bar defined as the standard error on the mean: $\\sigma/\\sqrt{N}$.

Usage:
        python3 comparison_crosssection.py <TYPE1_OUTROOT> <TYPE2_OUTROOT>

Outputs:
    - If exactly one varying group is found: ./comparison/crosssection/<type1>_vs_<type2>_crosssection.png
    - If multiple varying groups exist: ./comparison/crosssection/<type1>_vs_<type2>_crosssection__<vary>__<idx>.png
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

# --- 参数解析 ---
def parse_args():
    parser = argparse.ArgumentParser(description="Compare cross-sections for type1/type2 scans.")
    parser.add_argument("type1", type=str, help="type1 scan folder (e.g. case1_scan2)")
    parser.add_argument("type2", type=str, help="type2 scan folder (e.g. case2_scan2)")
    parser.add_argument("--xmin", type=float, default=None, help="Optional x-axis lower bound (e.g. 0.01)")
    parser.add_argument("--xmax", type=float, default=None, help="Optional x-axis upper bound")
    parser.add_argument(
        "--diagnose",
        action="store_true",
        help="Print repeat-run diagnostics (n>1) and exit.",
    )
    return parser.parse_args()


def _print_repeat_diagnostics(*, label: str, rows: list[Row]) -> None:
    stats = _stats_by_key(rows)
    ns = [v[2] for v in stats.values()]
    nrep = sum(1 for n in ns if n > 1)
    nmax = max(ns) if ns else 0
    print(f"{label}:")
    print(f"  Events_* rows: {len(rows)}")
    print(f"  Unique points: {len(stats)}")
    print(f"  Points with repeats (n>1): {nrep}")
    print(f"  Max repeats n: {nmax}")

    if nrep:
        # Show a few examples keyed by mchi (often the scan parameter).
        examples: list[tuple[float, int, float, float]] = []
        for k, (mean, sem, n) in stats.items():
            if n > 1:
                mchi = k[4]
                examples.append((mchi, n, mean, sem))
        examples.sort(key=lambda t: (t[0], -t[1]))
        for mchi, n, mean, sem in examples[:12]:
            print(f"    mchi={mchi:g}  n={n}  mean={mean:.6g} pb  SEM={sem:.3g} pb")


def _format_repeat_diagnostics(*, label: str, rows: list[Row]) -> str:
    stats = _stats_by_key(rows)
    ns = [v[2] for v in stats.values()]
    nrep = sum(1 for n in ns if n > 1)
    nmax = max(ns) if ns else 0
    lines: list[str] = []
    lines.append(f"{label}:")
    lines.append(f"  Events_* rows: {len(rows)}")
    lines.append(f"  Unique points: {len(stats)}")
    lines.append(f"  Points with repeats (n>1): {nrep}")
    lines.append(f"  Max repeats n: {nmax}")

    if nrep:
        lines.append("  Examples (mchi, n, mean, SEM):")
        examples: list[tuple[float, int, float, float]] = []
        for k, (mean, sem, n) in stats.items():
            if n > 1:
                mchi = k[4]
                examples.append((mchi, n, mean, sem))
        examples.sort(key=lambda t: (t[0], -t[1]))
        for mchi, n, mean, sem in examples[:20]:
            lines.append(f"    mchi={mchi:g}  n={n}  mean={mean:.6g} pb  SEM={sem:.3g} pb")
    return "\n".join(lines) + "\n"

# --- 目录参数匹配 ---
def parse_events_dir_params(dirname):
    name = dirname.replace("Events_", "")
    # Allow repeated runs to be stored under a unique suffix like "_job123".
    # Strip it so it doesn't interfere with parameter matching.
    name = re.sub(r"(?:^|_)job\d+(?:$|_)", "_", name)
    out = {}
    st_match = re.search(r"(?:^|_)st(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if st_match:
        out["sin"] = _tag_to_float(st_match.group(1))
    tb_match = re.search(r"(?:^|_)tb(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if tb_match:
        out["tan"] = _tag_to_float(tb_match.group(1))
    mA_match = re.search(r"mA([0-9]+)", name)
    if mA_match:
        out["mA"] = float(mA_match.group(1))
    ma_match = re.search(r"ma([0-9]+)", name)
    if ma_match:
        out["ma"] = float(ma_match.group(1))
    mchi_match = re.search(r"mchi(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if mchi_match:
        out["mchi"] = _tag_to_float(mchi_match.group(1))
    cbma_match = re.search(r"cbma(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if cbma_match:
        out["cosbma"] = _tag_to_float(cbma_match.group(1))
    dm_match = re.search(r"dm(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if dm_match:
        out["deltam"] = _tag_to_float(dm_match.group(1))
    return out


def _tag_to_float(tag: str) -> float:
    neg = tag.startswith("m")
    if neg:
        tag = tag[1:]
    v = float(tag.replace("p", "."))
    return -v if neg else v


def find_banner_file(events_dir: Path) -> Path | None:
    run_dirs = sorted([p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])
    for run_dir in run_dirs:
        for banner in sorted(run_dir.glob("*_banner.txt")):
            return banner
    return None


def extract_cross_section_pb(banner_path: Path) -> float | None:
    txt = banner_path.read_text(errors="ignore")
    for line in txt.splitlines():
        if "Integrated weight (pb)" in line:
            m = re.search(r":\s*([+\-0-9.eE]+)", line)
            if not m:
                continue
            try:
                return float(m.group(1))
            except Exception:
                return None
    return None


def extract_mchi_from_banner(banner_path: Path) -> float | None:
    """Best-effort extraction of mchi from banner param card.

    We treat MASS(52) as mchi, matching the conventions used elsewhere in this workspace.
    """

    txt = banner_path.read_text(errors="ignore")
    m = re.search(r"^\s*52\s+([+\-0-9.eE]+)\b", txt, flags=re.MULTILINE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


@dataclass(frozen=True)
class Row:
    sin: float
    tan: float
    mA: float
    ma: float
    mchi: float
    cosbma: float
    deltam: float
    xs_pb: float


@dataclass(frozen=True)
class StatRow:
    sin: float
    tan: float
    mA: float
    ma: float
    mchi: float
    cosbma: float
    deltam: float
    xs_pb: float
    xs_err_pb: float
    n: int


def discover_rows(outroot: Path) -> list[Row]:
    results_dir = outroot / "results"
    if not results_dir.is_dir():
        return []
    rows: list[Row] = []
    for ev_dir in sorted([p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith("Events_")]):
        params = parse_events_dir_params(ev_dir.name)
        if "sin" not in params or "tan" not in params or "mA" not in params or "ma" not in params:
            continue
        banner = find_banner_file(ev_dir)
        if not banner:
            continue
        xs = extract_cross_section_pb(banner)
        if xs is None:
            continue

        # Some legacy directory names omit mchi/cosbma/deltam. If mchi is missing,
        # try to infer it from the banner's MASS(52).
        if "mchi" not in params:
            mch = extract_mchi_from_banner(banner)
            if mch is not None:
                params["mchi"] = mch
        rows.append(
            Row(
                sin=float(params["sin"]),
                tan=float(params["tan"]),
                mA=float(params["mA"]),
                ma=float(params["ma"]),
                mchi=float(params.get("mchi", 0.0)),
                cosbma=float(params.get("cosbma", 0.0)),
                deltam=float(params.get("deltam", 0.0)),
                xs_pb=float(xs),
            )
        )
    return rows


def _x_label_for(vary: str) -> str:
        if vary == "deltam":
                return "#Delta m [GeV]"
        if vary == "mchi":
                return "m_{#chi} [GeV]"
        if vary in ("mA", "ma", "mchi"):
                return f"{vary} [GeV]"
        if vary == "sin":
                return "sin#theta"
        if vary == "cosbma":
                return "cos(#beta-#alpha)"
        if vary == "tan":
                return "tan#beta"
        return vary


def _fixed_text_for(vary: str, values: dict[str, float]) -> str:
        fixed_lines: list[str] = []
        if vary != "mA":
                fixed_lines.append(f"m_{{A}} = {values['mA']:g} GeV")
        if vary != "ma":
                fixed_lines.append(f"m_{{a}} = {values['ma']:g} GeV")
        if vary != "mchi":
                fixed_lines.append(f"m_{{#chi}} = {values['mchi']:g} GeV")
        if vary != "sin":
                fixed_lines.append(f"sin#theta = {values['sin']:g}")
        if vary != "tan":
                fixed_lines.append(f"tan#beta = {values['tan']:g}")
        if vary != "cosbma":
                fixed_lines.append(f"cos(#beta-#alpha) = {values['cosbma']:g}")
        if vary != "deltam":
            fixed_lines.append(f"#Delta m = m_{{H^{{+}}}} - m_{{A}} = {values['deltam']:g} GeV")
        return "; ".join(fixed_lines)


def root_quote(path: Path) -> str:
        return str(path).replace("\\", "\\\\").replace('"', '\\"')


def root_quote_str(value: str) -> str:
        return value.replace("\\", "\\\\").replace('"', '\\"')


def write_root_macro(macro_path: Path) -> None:
        macro_path.write_text(
                r'''#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static bool read_xye(const char* data_txt, std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& es, bool& has_err_col) {
    std::ifstream in(data_txt);
    if (!in) return false;
    has_err_col = false;
    double x, y;
    double e;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        std::istringstream ss(line);
        if (!(ss >> x >> y)) continue;
        if (ss >> e) {
            has_err_col = true;
        } else {
            e = 0.0; // legacy 2-column input => no error bars
        }
        xs.push_back(x);
        ys.push_back(y);
        es.push_back(e);
    }
    return true;
}

int make_scatter2(const char* data1_txt,
                                 const char* data2_txt,
                                 const char* x_label,
                                 const char* fixed_text,
                                 const char* out_png)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    // Make error-bar caps visible.
    gStyle->SetEndErrorSize(6);

    std::vector<double> x1, y1, e1, x2, y2, e2;
    bool hasE1 = false;
    bool hasE2 = false;
    if (!read_xye(data1_txt, x1, y1, e1, hasE1)) {
        std::cerr << "ERROR: cannot open data file: " << data1_txt << std::endl;
        return 2;
    }
    if (!read_xye(data2_txt, x2, y2, e2, hasE2)) {
        std::cerr << "ERROR: cannot open data file: " << data2_txt << std::endl;
        return 2;
    }
    if (x1.size() < 3 || x2.size() < 3) {
        std::cerr << "ERROR: need >=3 points for both type1/type2" << std::endl;
        return 3;
    }

    // Decide whether to use log-y based on dynamic range across both.
    double ymin_pos = 0.0;
    double ymax = 0.0;
    auto scan = [&](const std::vector<double>& ys) {
        for (size_t i = 0; i < ys.size(); ++i) {
            const double yv = ys[i];
            if (yv > 0.0) {
                if (ymin_pos == 0.0 || yv < ymin_pos) ymin_pos = yv;
                if (yv > ymax) ymax = yv;
            } else {
                if (yv > ymax) ymax = yv;
            }
        }
    };
    scan(y1);
    scan(y2);
    const double ratio = (ymin_pos > 0.0 && ymax > 0.0) ? (ymax / ymin_pos) : 0.0;
    // Use log-y when dynamic range is large.
    // Threshold lowered to 20x to make type1/type2 overlays easier to read.
    const bool useLogY = (ratio >= 20.0);

    // X axis: keep linear (user-requested).
    const bool useLogX = false;

    std::vector<double> ex1(x1.size(), 0.0);
    std::vector<double> ex2(x2.size(), 0.0);
    TGraphErrors g1e((int)x1.size(), x1.data(), y1.data(), ex1.data(), e1.data());
    TGraphErrors g2e((int)x2.size(), x2.data(), y2.data(), ex2.data(), e2.data());
    TGraph g1((int)x1.size(), x1.data(), y1.data());
    TGraph g2((int)x2.size(), x2.data(), y2.data());

    TGraph* g1p = hasE1 ? (TGraph*)&g1e : &g1;
    TGraph* g2p = hasE2 ? (TGraph*)&g2e : &g2;

    g1p->SetTitle("");
    g1p->GetXaxis()->SetTitle(x_label);
    g1p->GetYaxis()->SetTitle("#sigma [pb]");
    g1p->GetXaxis()->SetTitleSize(0.045);
    g1p->GetYaxis()->SetTitleSize(0.045);
    g1p->GetXaxis()->SetLabelSize(0.040);
    g1p->GetYaxis()->SetLabelSize(0.040);
    g1p->GetYaxis()->SetTitleOffset(1.60);

    g1p->SetMarkerStyle(20);
    g1p->SetMarkerSize(1.1);
    g1p->SetMarkerColor(kBlue+1);
    g1p->SetLineColor(kBlue+1);
    g1p->SetLineWidth(3);

    // Type2: make it visually distinct and hard to miss even if it overlaps.
    g2p->SetMarkerStyle(25);
    g2p->SetMarkerSize(1.35);
    g2p->SetMarkerColor(kRed+1);
    g2p->SetLineColor(kRed+1);
    g2p->SetLineWidth(4);
    g2p->SetLineStyle(2);

    TGaxis::SetExponentOffset(-0.08, 0.0, "y");

    TCanvas c("c", "c", 1400, 1000);
    c.SetLeftMargin(0.16);
    c.SetRightMargin(0.30);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
    if (useLogY) {
        c.SetLogy();
        if (ymin_pos > 0.0 && ymax > 0.0) {
            g1p->GetYaxis()->SetRangeUser(ymin_pos * 0.7, ymax * 1.3);
            if (ratio > 1000.0) {
                // If y spans more than 3 orders of magnitude,
                // keep default decade labels only.
            } else {
                g1p->GetYaxis()->SetMoreLogLabels();
            }
        }
    } else {
        // IMPORTANT: ROOT auto-scales using the first drawn graph only (g1).
        // For overlays, explicitly include both type1/type2 ranges.
        if (ymax > 0.0) {
            g1p->GetYaxis()->SetRangeUser(0.0, ymax * 1.15);
        }
    }

    if (useLogX) {
        c.SetLogx();
    }

    // Draw type1 first to set axes.
    // Use X0 to suppress x-error bars (we only want y-errors).
    if (hasE1) g1p->Draw("APE1X0");
    else g1p->Draw("AP");
    g1p->Draw("L SAME");

    // Draw type2 LAST so it remains visible when overlapping.
    if (hasE2) g2p->Draw("PE1X0 SAME");
    else g2p->Draw("P SAME");
    g2p->Draw("L SAME");

    // Redraw type2 markers/errors on top to make error bars and points pop.
    if (hasE2) g2p->Draw("PE1X0 SAME");
    else g2p->Draw("P SAME");

    // CMS-style header (same placement as plot_crosssection_scatter.py)
    TLatex latex;
    latex.SetNDC(true);
    const double lm = c.GetLeftMargin();
    const double rm = c.GetRightMargin();
    const double tm = c.GetTopMargin();
    const double yCms = (1.0 - tm) + 0.012;
    const double xInsideRight = 1.0 - rm - 0.01;
    const double yInsideTop = 1.0 - tm - 0.004;

    latex.SetTextFont(61);
    latex.SetTextSize(0.055);
    latex.DrawLatex(lm + 0.01, yCms, "CMS Internal");

    latex.SetTextFont(42);
    latex.SetTextSize(0.045);
    latex.SetTextAlign(31);
    latex.DrawLatex(xInsideRight, yInsideTop, "(13.6 TeV)");

    // Manual legend (top-left inside the plot area)
    double lx1 = lm + 0.02;
    double lx2 = lx1 + 0.06;
    double ly = 1.0 - tm - 0.07;
    const double dy = 0.055;
    TLine line;
    line.SetLineWidth(3);
    latex.SetTextFont(42);
    latex.SetTextSize(0.040);
    latex.SetTextAlign(12);
    line.SetLineColor(kBlue+1);
    line.DrawLineNDC(lx1, ly, lx2, ly);
    latex.DrawLatex(lx2 + 0.02, ly - 0.015, "type1");
    ly -= dy;
    line.SetLineColor(kRed+1);
    line.DrawLineNDC(lx1, ly, lx2, ly);
    latex.DrawLatex(lx2 + 0.02, ly - 0.015, "type2");

    // Fixed parameter block (right) (same as plot_crosssection_scatter.py)
    latex.SetTextAlign(13);
    latex.SetTextFont(42);
    latex.SetTextSize(0.038);
    {
        const double x0 = 0.72;
        double y0 = 0.88;
        latex.SetTextSize(0.040);
        latex.DrawLatex(x0, y0, "2HDM+a     b#bar{b} + #chi#bar{#chi}");
        latex.SetTextSize(0.038);

        std::string ft = fixed_text ? fixed_text : "";
        std::vector<std::string> lines;
        size_t start = 0;
        while (start < ft.size()) {
            size_t pos = ft.find(';', start);
            std::string s = (pos == std::string::npos) ? ft.substr(start) : ft.substr(start, pos - start);
            while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(s.begin());
            if (!s.empty()) lines.push_back(s);
            if (pos == std::string::npos) break;
            start = pos + 1;
        }

        double y = y0 - 0.055;
        for (size_t i = 0; i < lines.size(); ++i) {
            latex.DrawLatex(x0, y, lines[i].c_str());
            y -= 0.055;
        }
    }

    c.SaveAs(out_png);
    gSystem->Exit(0);
    return 0;
}
'''
        )


def _run_root_make_scatter(
        *,
        root_exe: str,
        macro: Path,
        data1_txt: Path,
        data2_txt: Path,
        x_label: str,
        fixed_text: str,
        out_png: Path,
) -> tuple[int, str]:
        stmt = (
                f'gROOT->ProcessLine(".L {root_quote(macro)}"); '
                f'make_scatter2("{root_quote(data1_txt)}","{root_quote(data2_txt)}","{root_quote_str(x_label)}","{root_quote_str(fixed_text)}","{root_quote(out_png)}");'
        )
        cmd = [root_exe, "-l", "-b", "-q", "-e", stmt]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return proc.returncode, proc.stdout


def group_rows(rows: list[Row]):
    params = ["sin", "tan", "mA", "ma", "mchi", "cosbma", "deltam"]
    # groups[(vary, fixed_tuple)] = {x: [ys...]}
    groups: dict[tuple[str, tuple[float, ...]], dict[float, list[float]]] = {}
    for r in rows:
        values = {
            "sin": r.sin,
            "tan": r.tan,
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
            groups.setdefault(key, {}).setdefault(x, []).append(r.xs_pb)
    return groups


def group_stat_rows(rows: list[StatRow]):
    params = ["sin", "tan", "mA", "ma", "mchi", "cosbma", "deltam"]
    # groups[(vary, fixed_tuple)] = {x: (mean, err, n)}
    groups: dict[tuple[str, tuple[float, ...]], dict[float, tuple[float, float, int]]] = {}
    for r in rows:
        values = {
            "sin": r.sin,
            "tan": r.tan,
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
            prev = groups.setdefault(key, {}).get(x)
            cur = (r.xs_pb, r.xs_err_pb, r.n)
            if prev is None or cur[2] >= prev[2]:
                groups[key][x] = cur
    return groups


def mean_points(x_to_ys: dict[float, list[float]]) -> list[tuple[float, float]]:
    pts = [(x, sum(ys) / len(ys)) for x, ys in x_to_ys.items()]
    pts.sort(key=lambda t: t[0])
    return pts


def mean_points_err(x_to_stats: dict[float, tuple[float, float, int]]) -> list[tuple[float, float, float, int]]:
    pts = [(x, v[0], v[1], v[2]) for x, v in x_to_stats.items()]
    pts.sort(key=lambda t: t[0])
    return pts


def _in_x_range(x: float, xmin: float | None, xmax: float | None) -> bool:
    if xmin is not None and x < xmin:
        return False
    if xmax is not None and x > xmax:
        return False
    return True


def _key_of_row(r: Row) -> tuple[float, float, float, float, float, float, float]:
    return (r.sin, r.tan, r.mA, r.ma, r.mchi, r.cosbma, r.deltam)


def _avg_by_key(rows: list[Row]) -> dict[tuple[float, float, float, float, float, float, float], float]:
    by: dict[tuple[float, float, float, float, float, float, float], list[float]] = {}
    for r in rows:
        by.setdefault(_key_of_row(r), []).append(r.xs_pb)
    return {k: (sum(v) / len(v)) for k, v in by.items()}


def _mean_std(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    mean = sum(values) / len(values)
    if len(values) < 2:
        return mean, 0.0
    var = sum((x - mean) ** 2 for x in values) / (len(values) - 1)
    if var < 0.0 and var > -1e-18:
        var = 0.0
    return mean, math.sqrt(var)


def _stats_by_key(
    rows: list[Row],
) -> dict[tuple[float, float, float, float, float, float, float], tuple[float, float, int]]:
    """Return mapping key -> (mean, stat_err, n) across repeated runs.

    stat_err is the standard error on the mean: sample_std / sqrt(n).
    """
    by: dict[tuple[float, float, float, float, float, float, float], list[float]] = {}
    for r in rows:
        by.setdefault(_key_of_row(r), []).append(r.xs_pb)
    out: dict[tuple[float, float, float, float, float, float, float], tuple[float, float, int]] = {}
    for k, vals in by.items():
        mean, std = _mean_std(vals)
        n = len(vals)
        stat_err = (std / math.sqrt(n)) if n > 1 else 0.0
        out[k] = (mean, stat_err, n)
    return out


def _overlap_diagnostics(
    k1: set[tuple[float, float, float, float, float, float, float]],
    k2: set[tuple[float, float, float, float, float, float, float]],
) -> str:
    """Return a short hint about which parameter likely prevents overlap."""

    params = ["sin", "tan", "mA", "ma", "mchi", "cosbma", "deltam"]
    if not k1 or not k2:
        return ""

    lines: list[str] = []
    # Show unique mchi sets (often the culprit)
    mchi1 = sorted({k[4] for k in k1})
    mchi2 = sorted({k[4] for k in k2})
    if mchi1 != mchi2:
        lines.append(f"mchi in type1: {mchi1[:6]}{'...' if len(mchi1) > 6 else ''}")
        lines.append(f"mchi in type2: {mchi2[:6]}{'...' if len(mchi2) > 6 else ''}")

    # Check whether ignoring exactly one parameter would create overlap.
    for i, p in enumerate(params):
        proj1 = {k[:i] + k[i + 1 :] for k in k1}
        proj2 = {k[:i] + k[i + 1 :] for k in k2}
        n = len(proj1 & proj2)
        if n > 0:
            lines.append(f"If you ignore '{p}', there would be {n} overlapping points (so '{p}' likely differs).")

    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    if args.xmin is not None and args.xmax is not None and args.xmax < args.xmin:
        print("ERROR: --xmax must be >= --xmin")
        return 2
    out1 = Path(args.type1)
    out2 = Path(args.type2)

    rows1 = discover_rows(out1)
    rows2 = discover_rows(out2)
    if not rows1 and not rows2:
        print("ERROR: no valid Events_* rows found in either folder")
        return 2

    if args.diagnose:
        # Also write to file so diagnostics are accessible even if stdout is suppressed.
        outdir = Path("comparison") / "crosssection"
        outdir.mkdir(parents=True, exist_ok=True)
        diag_path = outdir / f"diagnose__{args.type1}_vs_{args.type2}.txt"
        txt = _format_repeat_diagnostics(label=str(out1), rows=rows1) + "\n" + _format_repeat_diagnostics(
            label=str(out2), rows=rows2
        )
        diag_path.write_text(txt)

        _print_repeat_diagnostics(label=str(out1), rows=rows1)
        _print_repeat_diagnostics(label=str(out2), rows=rows2)
        print(f"Wrote: {diag_path}")
        return 0

    # Compare only common parameter points; otherwise group keys won't match and
    # the plot won't be a true type1-vs-type2 comparison.
    s1 = _stats_by_key(rows1)
    s2 = _stats_by_key(rows2)
    common = sorted(set(s1.keys()) & set(s2.keys()))
    if not common:
        print("ERROR: no overlapping parameter points between type1 and type2.")
        print("Hint: make sure you compare matching scans (e.g. ..._scan2 vs ..._scan2) and that both have completed results.")
        diag = _overlap_diagnostics(set(s1.keys()), set(s2.keys()))
        if diag:
            print(diag)
        return 2

    common_rows1: list[StatRow] = []
    common_rows2: list[StatRow] = []
    for k in common:
        mean1, std1, n1 = s1[k]
        mean2, std2, n2 = s2[k]
        common_rows1.append(
            StatRow(
                sin=k[0],
                tan=k[1],
                mA=k[2],
                ma=k[3],
                mchi=k[4],
                cosbma=k[5],
                deltam=k[6],
                xs_pb=mean1,
                xs_err_pb=std1,
                n=n1,
            )
        )
        common_rows2.append(
            StatRow(
                sin=k[0],
                tan=k[1],
                mA=k[2],
                ma=k[3],
                mchi=k[4],
                cosbma=k[5],
                deltam=k[6],
                xs_pb=mean2,
                xs_err_pb=std2,
                n=n2,
            )
        )

    g1 = group_stat_rows(common_rows1)
    g2 = group_stat_rows(common_rows2)

    keys = sorted(set(g1.keys()) | set(g2.keys()), key=lambda k: (k[0], k[1]))
    # Now both sides are built from the same parameter set; require >=3 points for both.
    good = [k for k in keys if len(g1.get(k, {})) >= 3 and len(g2.get(k, {})) >= 3]
    if not good:
        print("ERROR: no varying-parameter groups with >=3 common points found")
        return 2

    outdir = Path("comparison") / "crosssection"
    outdir.mkdir(parents=True, exist_ok=True)

    root_exe = shutil.which("root")
    if not root_exe:
        print("ERROR: CERN ROOT executable 'root' not found on PATH.")
        return 2

    macro = outdir / ".comparison_crosssection_tmp.C"
    write_root_macro(macro)

    single = len(good) == 1
    made = 0
    try:
        for idx, (vary, fixed_tuple) in enumerate(good, start=1):
            pts1 = mean_points_err(g1.get((vary, fixed_tuple), {}))
            pts2 = mean_points_err(g2.get((vary, fixed_tuple), {}))

            pts1 = [p for p in pts1 if _in_x_range(p[0], args.xmin, args.xmax)]
            pts2 = [p for p in pts2 if _in_x_range(p[0], args.xmin, args.xmax)]
            if len(pts1) < 3 or len(pts2) < 3:
                print(f"Skip group vary={vary}: <3 points after x-range filtering")
                continue

            if single:
                outfile = outdir / f"{args.type1}_vs_{args.type2}_crosssection.png"
            else:
                outfile = outdir / f"{args.type1}_vs_{args.type2}_crosssection__{vary}__{idx}.png"

            data1 = outdir / f".tmp_type1__{vary}__{idx}.txt"
            data2 = outdir / f".tmp_type2__{vary}__{idx}.txt"
            try:
                has_rep_1 = any(n > 1 for _x, _y, _e, n in pts1)
                has_rep_2 = any(n > 1 for _x, _y, _e, n in pts2)

                if has_rep_1:
                    data1.write_text("\n".join([f"{x:g} {y:.12g} {e:.12g}" for x, y, e, _n in pts1]) + "\n")
                else:
                    data1.write_text("\n".join([f"{x:g} {y:.12g}" for x, y, _e, _n in pts1]) + "\n")

                if has_rep_2:
                    data2.write_text("\n".join([f"{x:g} {y:.12g} {e:.12g}" for x, y, e, _n in pts2]) + "\n")
                else:
                    data2.write_text("\n".join([f"{x:g} {y:.12g}" for x, y, _e, _n in pts2]) + "\n")

                fixed_names = [p for p in ["sin", "tan", "mA", "ma", "mchi", "cosbma", "deltam"] if p != vary]
                fixed_values = {fixed_names[i]: fixed_tuple[i] for i in range(len(fixed_names))}
                # Reconstruct full values for the fixed text helper
                values_full = {
                    "sin": fixed_values.get("sin", 0.0),
                    "tan": fixed_values.get("tan", 0.0),
                    "mA": fixed_values.get("mA", 0.0),
                    "ma": fixed_values.get("ma", 0.0),
                    "mchi": fixed_values.get("mchi", 0.0),
                    "cosbma": fixed_values.get("cosbma", 0.0),
                    "deltam": fixed_values.get("deltam", 0.0),
                }
                x_label = _x_label_for(vary)
                fixed_text = _fixed_text_for(vary, values_full)

                rc, out = _run_root_make_scatter(
                    root_exe=root_exe,
                    macro=macro,
                    data1_txt=data1,
                    data2_txt=data2,
                    x_label=x_label,
                    fixed_text=fixed_text,
                    out_png=outfile,
                )
            finally:
                # Always delete intermediate XY text files.
                try:
                    data1.unlink(missing_ok=True)
                    data2.unlink(missing_ok=True)
                except Exception:
                    pass

            if rc != 0:
                log = outdir / (outfile.stem + ".log")
                log.write_text(out)
                print(f"ERROR: ROOT failed (see {log})")
                return 2

            print(f"Saved: {outfile}")
            made += 1

        if made == 0:
            print("ERROR: no groups with >=3 points remain after applying x-range filter")
            return 2
        return 0
    finally:
        # Always delete the temporary ROOT macro.
        try:
            macro.unlink(missing_ok=True)
        except Exception:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
