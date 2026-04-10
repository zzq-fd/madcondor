#!/usr/bin/env python3
"""Make cross-section scatter plots from MadGraph banner files.

Two input modes:

1) Auto-discovery mode (like plot_pt55_scatter_groups.py):
     python3 plot_crosssection_scatter.py <OUTROOT>

     The script scans <OUTROOT>/results/Events_* and parses parameters encoded
     in directory names, then extracts cross sections from *_banner.txt.
     If there are >=3 entries where all-but-one parameters are fixed and the
     remaining one varies, it makes ONE scatter plot per such group.
     If your scan varies multiple parameters across different groups, it will
     generate plots for each varying parameter group.

2) Explicit mode (backward compatible):
     Same as plot_met55_overlay.py:
         --sin/--tan/--mA/--ma/--mchi/--cosbma/--deltam accept comma-separated lists.
     Requirements:
         - exactly one parameter has >=3 values (used as x-axis)
         - y-axis is cross section (pb)

Outputs:
    <OUTROOT>/output_scatter_plot/crosssection/*.png

Notes:
    - mA corresponds to mass(36) input.
    - deltam allows negatives; e.g. --deltam -200,-150,... (no quotes needed).
    - Requires CERN ROOT executable `root` (no PyROOT).
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


def _split_list(value: str) -> list[str]:
    parts = [v.strip() for v in value.split(",")]
    return [p for p in parts if p]


def _to_floats(values: list[str], name: str) -> list[float]:
    out: list[float] = []
    for v in values:
        try:
            out.append(float(v))
        except Exception:
            raise SystemExit(f"ERROR: cannot parse {name} value: {v}")
    return out


def _st_tag(x: float) -> str:
    # runMadscan.sh uses: '-' -> 'm', '.' -> 'p'
    return f"{x:g}".replace("-", "m").replace(".", "p")


def _tag_to_float(tag: str) -> float:
    neg = tag.startswith("m")
    if neg:
        tag = tag[1:]
    v = float(tag.replace("p", "."))
    return -v if neg else v


def parse_events_dir_params(dirname: str) -> dict[str, float]:
    name = dirname.replace("Events_", "")
    out: dict[str, float] = {}

    st_match = re.search(r"(?:^|_)st(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if st_match:
        out["sin"] = _tag_to_float(st_match.group(1))

    tb_match = re.search(r"(?:^|_)tb(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if tb_match:
        # Support tags like: tb3, tb0.9, tb0.01 (written by runMadscan.sh), and also tb0p9/tb0p01.
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


def choose_varying(values_map: dict[str, list[float]]) -> str:
    varying = [k for k, v in values_map.items() if len(v) > 1]
    if len(varying) != 1:
        raise SystemExit(
            "ERROR: 必须且只能有一个参数给多个值（>=3），作为横轴。"
            f" 当前多值参数: {', '.join(varying) if varying else '(none)'}"
        )
    if len(values_map[varying[0]]) < 3:
        raise SystemExit(f"ERROR: 横轴参数 {varying[0]} 的点数必须 >= 3")
    return varying[0]


def find_banner_file(events_dir: Path) -> Path | None:
    run_dirs = sorted([p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])
    for run_dir in run_dirs:
        for banner in sorted(run_dir.glob("*_banner.txt")):
            return banner
    return None


def extract_cross_section_pb(banner_path: Path) -> float | None:
    # Expected line:
    #   #  Integrated weight (pb)  :       0.00017319891580000002
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


def write_root_macro(macro_path: Path) -> None:
    macro_path.write_text(
        r'''#include <TCanvas.h>
    #include <TGraph.h>
    #include <TAxis.h>
    #include <TGaxis.h>
    #include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int make_scatter(const char* data_txt,
                 const char* x_label,
                                 const char* fixed_text,
                 const char* out_png)
{
  gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);

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

    const double logRatioThreshold = 100.0;

    // Decide whether to use log-x based on dynamic range.
    // Require all x > 0 for log scale.
    bool xAllPositive = true;
    double xmin_pos = 0.0;
    double xmax = 0.0;
    for (size_t i = 0; i < xs.size(); ++i) {
        const double xv = xs[i];
        if (xv <= 0.0) {
            xAllPositive = false;
            continue;
        }
        if (xmin_pos == 0.0 || xv < xmin_pos) xmin_pos = xv;
        if (xv > xmax) xmax = xv;
    }
    const double xRatio = (xAllPositive && xmin_pos > 0.0 && xmax > 0.0) ? (xmax / xmin_pos) : 0.0;
    const bool useLogX = (xRatio >= logRatioThreshold);

    // Decide whether to use log-y based on dynamic range.
    // Only consider strictly positive y values.
    double ymin_pos = 0.0;
    double ymax = 0.0;
    for (size_t i = 0; i < ys.size(); ++i) {
        const double yv = ys[i];
        if (yv > 0.0) {
            if (ymin_pos == 0.0 || yv < ymin_pos) ymin_pos = yv;
            if (yv > ymax) ymax = yv;
        } else {
            if (yv > ymax) ymax = yv;
        }
    }
    const double yRatio = (ymin_pos > 0.0 && ymax > 0.0) ? (ymax / ymin_pos) : 0.0;
    const bool useLogY = (yRatio >= logRatioThreshold);

  TGraph g((int)xs.size());
  for (int i = 0; i < (int)xs.size(); ++i) {
    g.SetPoint(i, xs[i], ys[i]);
  }

    g.SetTitle("");
    g.GetXaxis()->SetTitle(x_label);
    g.GetYaxis()->SetTitle("#sigma [pb]");
    g.GetXaxis()->SetTitleSize(0.045);
    g.GetYaxis()->SetTitleSize(0.045);
    g.GetXaxis()->SetLabelSize(0.040);
    g.GetYaxis()->SetLabelSize(0.040);
    // Prevent y-axis title/label overlap
    g.GetYaxis()->SetTitleOffset(1.60);
  g.SetMarkerStyle(20);
    g.SetMarkerSize(1.1);
    g.SetMarkerColor(kBlack);

    // Move scientific-notation exponent (e.g. x10^{-6}) further left.
    TGaxis::SetExponentOffset(-0.08, 0.0, "y");

        TCanvas c("c", "c", 1400, 1000);
    c.SetLeftMargin(0.16);
    c.SetRightMargin(0.30);
    c.SetBottomMargin(0.12);
    c.SetTopMargin(0.08);
        if (useLogX) {
            c.SetLogx();
            // Ensure sensible visible range for log scale.
            if (xmin_pos > 0.0 && xmax > 0.0) {
                g.GetXaxis()->SetLimits(xmin_pos * 0.7, xmax * 1.3);
            }
        }
        if (useLogY) {
            c.SetLogy();
            // Ensure sensible visible range for log scale.
            if (ymin_pos > 0.0 && ymax > 0.0) {
                g.GetYaxis()->SetRangeUser(ymin_pos * 0.7, ymax * 1.3);
                if (yRatio > 1000.0) {
                    // If y spans more than 3 orders of magnitude,
                    // keep default decade labels only.
                } else {
                    g.GetYaxis()->SetMoreLogLabels();
                }
            }
        }

    // X-axis label simplification rules requested:
    // - If all x < 1: show only 3 labels
    // - If xmax >= 10: show only factors 1 and 2 per decade (e.g. 1,2,10,20,...)
    // - Otherwise: if many points, limit to maxLabels labels
    int nPoints = (int)xs.size();
    const int maxLabels = 10;
    double data_xmin = xs.front();
    double data_xmax = xs.back();

    if (data_xmax < 1.0) {
        // small-range decimals: force 3 labels
        g.GetXaxis()->SetNdivisions(3, kTRUE);
        g.GetXaxis()->SetLabelSize(0.035);
    } else if (data_xmax >= 10.0 && xAllPositive) {
        // wide dynamic range on log axis: prefer sparse, human-readable labels
        // Try to show only factor 1 and 2 per decade.
        // Setting Ndivisions to 202 requests 2 main divisions per decade (best-effort).
        g.GetXaxis()->SetNdivisions(202, kTRUE);
        g.GetXaxis()->SetLabelSize(0.035);
        // Avoid ROOT adding many minor log labels
        // (Don't call SetMoreLogLabels here so labels stay sparse.)
    } else if (nPoints > maxLabels) {
        g.GetXaxis()->SetNdivisions(maxLabels, kTRUE);
        g.GetXaxis()->SetLabelSize(0.030);
        g.GetXaxis()->SetLabelOffset(0.01);
    } else {
        g.GetXaxis()->SetNdivisions(510);
    }

  g.Draw("AP");

    // CMS-style header
    TLatex latex;
    latex.SetNDC(true);

    const double lm = c.GetLeftMargin();
    const double rm = c.GetRightMargin();
    const double tm = c.GetTopMargin();
    // Put CMS text just above the top frame line (still in the top margin).
    const double yCms = (1.0 - tm) + 0.012;
    // Put sqrt(s) inside the plot area (not in the right margin).
    const double xInsideRight = 1.0 - rm - 0.01;
    const double yInsideTop = 1.0 - tm - 0.004;

    latex.SetTextFont(61);
    latex.SetTextSize(0.055);
    latex.DrawLatex(lm + 0.01, yCms, "CMS Internal");

    latex.SetTextFont(42);
    latex.SetTextSize(0.045);
    latex.SetTextAlign(31);
    latex.DrawLatex(xInsideRight, yInsideTop, "(13.6 TeV)");

    // Fixed parameter block (right)
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

  // Exit immediately to avoid ROOT shutdown-time segfaults in batch mode.
  gSystem->Exit(0);
  return 0;
}
'''
    )


def root_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def root_quote_str(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


@dataclass(frozen=True)
class FixedParams:
    sin: float
    tan: float
    mA: float
    ma: float
    mchi: float
    cosbma: float
    deltam: float


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


def _x_label_for(param: str) -> str:
    if param == "deltam":
        return "#Delta m [GeV]"
    if param == "mchi":
        return "m_{#chi} [GeV]"
    if param in ("mA", "ma", "mchi"):
        return f"{param} [GeV]"
    if param == "sin":
        return "sin#theta"
    if param == "cosbma":
        return "cos(#beta-#alpha)"
    if param == "tan":
        return "tan#beta"
    return param


def _infer_type_label(outroot: Path) -> str | None:
    name = outroot.name.lower()
    if "type1" in name and "type2" not in name:
        return "type1"
    if "type2" in name and "type1" not in name:
        return "type2"
    return None


def _fixed_text_for(vary: str, values: dict[str, float], type_label: str | None = None) -> str:
    fixed_lines: list[str] = []
    if type_label:
        fixed_lines.append(f"model = {type_label}")
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


def _tag_for_values(values: dict[str, float]) -> str:
    return (
        f"st{_st_tag(values['sin'])}"
        f"_tb{_st_tag(values['tan'])}"
        f"_mA{int(values['mA'])}"
        f"_ma{int(values['ma'])}"
        f"_mchi{_st_tag(values['mchi'])}"
        f"_cbma{_st_tag(values['cosbma'])}"
        f"_dm{_st_tag(values['deltam'])}"
    )


def _discover_rows(results_dir: Path) -> list[Row]:
    event_dirs = sorted([p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith("Events_")])
    rows: list[Row] = []
    for ev_dir in event_dirs:
        params = parse_events_dir_params(ev_dir.name)

        # Required tags to make meaningful grouping. Others default to 0.0 for legacy naming.
        if "sin" not in params or "tan" not in params or "mA" not in params or "ma" not in params:
            # Legacy/unknown dir name; skip.
            continue

        sinv = params["sin"]
        tanv = params["tan"]
        mAv = params["mA"]
        mav = params["ma"]
        mchiv = params.get("mchi", 0.0)
        cbv = params.get("cosbma", 0.0)
        dmv = params.get("deltam", 0.0)

        banner = find_banner_file(ev_dir)
        if not banner:
            continue
        xs_pb = extract_cross_section_pb(banner)
        if xs_pb is None:
            continue

        rows.append(
            Row(
                sin=sinv,
                tan=tanv,
                mA=mAv,
                ma=mav,
                mchi=mchiv,
                cosbma=cbv,
                deltam=dmv,
                xs_pb=xs_pb,
            )
        )
    return rows


def _run_root_make_scatter(*, root_exe: str, macro: Path, data_txt: Path, x_label: str, fixed_text: str, out_png: Path) -> tuple[int, str]:
    stmt = (
        f'gROOT->ProcessLine(".L {root_quote(macro)}"); '
        f'make_scatter("{root_quote(data_txt)}","{root_quote_str(x_label)}","{root_quote_str(fixed_text)}","{root_quote(out_png)}");'
    )
    cmd = [root_exe, "-l", "-b", "-q", "-e", stmt]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc.returncode, proc.stdout


def _in_x_range(x: float, xmin: float | None, xmax: float | None) -> bool:
    if xmin is not None and x < xmin:
        return False
    if xmax is not None and x > xmax:
        return False
    return True


def main_auto(outroot: Path, xmin: float | None = None, xmax: float | None = None) -> int:
    results_dir = outroot / "results"
    if not results_dir.exists():
        print(f"ERROR: results directory not found: {results_dir}")
        return 2

    root_exe = shutil.which("root")
    if not root_exe:
        print("ERROR: CERN ROOT executable 'root' not found on PATH.")
        print("Hint: conda activate mg5")
        return 2

    rows = _discover_rows(results_dir)
    if not rows:
        print("No valid rows found (missing banner or cross section).")
        return 0

    out_dir = outroot / "output_scatter_plot" / "crosssection"
    out_dir.mkdir(parents=True, exist_ok=True)

    macro = out_dir / ".crosssection_tmp.C"
    write_root_macro(macro)
    type_label = _infer_type_label(outroot)

    params = ["sin", "tan", "mA", "ma", "mchi", "cosbma", "deltam"]

    # group[(vary, fixed_tuple)] = {x_value: [xs...]}
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

    made = 0
    any_failures = False

    for (vary, fixed_tuple), x_to_ys in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if len(x_to_ys) < 3:
            continue

        # Average duplicates at same x (multiple runs)
        points = [(x, sum(ys) / len(ys)) for x, ys in x_to_ys.items()]
        points = [(x, y) for x, y in points if _in_x_range(x, xmin, xmax)]
        if len(points) < 3:
            continue
        points.sort(key=lambda t: t[0])

        fixed_names = [p for p in params if p != vary]
        fixed_values = {fixed_names[i]: fixed_tuple[i] for i in range(len(fixed_names))}

        # For tags/labels we still want a complete value dict.
        full_values = dict(fixed_values)
        # Use the first x for the varying param in tags (tag only indicates fixed ones anyway).
        full_values[vary] = points[0][0]

        fixed_text = _fixed_text_for(
            vary,
            {
                "sin": full_values.get("sin", 0.0),
                "tan": full_values.get("tan", 0.0),
                "mA": full_values.get("mA", 0.0),
                "ma": full_values.get("ma", 0.0),
                "mchi": full_values.get("mchi", 0.0),
                "cosbma": full_values.get("cosbma", 0.0),
                "deltam": full_values.get("deltam", 0.0),
            },
            type_label=type_label,
        )

        fixed_tag = _tag_for_values({
            "sin": full_values.get("sin", 0.0),
            "tan": full_values.get("tan", 0.0),
            "mA": full_values.get("mA", 0.0),
            "ma": full_values.get("ma", 0.0),
            "mchi": full_values.get("mchi", 0.0),
            "cosbma": full_values.get("cosbma", 0.0),
            "deltam": full_values.get("deltam", 0.0),
        })

        out_png = out_dir / f"crosssection_{vary}__{fixed_tag}.png"
        x_label = _x_label_for(vary)

        data_txt = out_dir / f".crosssection_{vary}_{hash((vary, fixed_tag)) & 0xffffffff:x}.dat"
        with open(data_txt, "w") as f:
            for x, y in points:
                f.write(f"{x} {y}\n")

        rc, out = _run_root_make_scatter(
            root_exe=root_exe,
            macro=macro,
            data_txt=data_txt,
            x_label=x_label,
            fixed_text=fixed_text,
            out_png=out_png,
        )

        log_path = out_dir / (out_png.stem + ".log")
        log_path.write_text(out)

        if rc != 0:
            any_failures = True
            print(f"ERROR: ROOT failed for {out_png.name} (see {log_path})")
            continue

        # Cleanup temps on success
        try:
            log_path.unlink(missing_ok=True)
        except Exception:
            pass
        try:
            data_txt.unlink(missing_ok=True)
        except Exception:
            pass

        made += 1

    # Clean macro and compiled artifacts
    try:
        macro.unlink(missing_ok=True)
    except Exception:
        pass
    compiled_base = macro.name.replace(".C", "_C")
    compiled_bases = [compiled_base]
    if compiled_base.startswith("."):
        compiled_bases.append("_" + compiled_base[1:])
    for base in compiled_bases:
        for p in out_dir.glob(f"{base}*"):
            try:
                if p.is_file():
                    p.unlink(missing_ok=True)
            except Exception:
                pass

    if made == 0 and not any_failures:
        print("No groups with >=3 distinct x values found.")
        return 0

    print(f"Wrote {made} plot(s) into: {out_dir}")
    return 0 if not any_failures else 1


def main() -> int:
    # Auto mode: first positional arg is OUTROOT and no explicit scan flags are used.
    explicit_flags = {"--sin", "--tan", "--mA", "--ma", "--mchi", "--cosbma", "--deltam"}
    if len(sys.argv) >= 2 and not sys.argv[1].startswith("-"):
        has_explicit_flags = any((tok.split("=")[0] in explicit_flags) for tok in sys.argv[2:])
        if not has_explicit_flags:
            ap_auto = argparse.ArgumentParser(add_help=True, allow_abbrev=False)
            ap_auto.add_argument("outroot", help="Top-level output directory (e.g. bbdm_2HDMa_type1_case1_scan)")
            ap_auto.add_argument("--xmin", type=float, default=None, help="Optional x-axis lower bound")
            ap_auto.add_argument("--xmax", type=float, default=None, help="Optional x-axis upper bound")
            a = ap_auto.parse_args()
            if a.xmin is not None and a.xmax is not None and a.xmax < a.xmin:
                raise SystemExit("ERROR: --xmax must be >= --xmin")
            return main_auto(Path(a.outroot).expanduser(), xmin=a.xmin, xmax=a.xmax)

    ap = argparse.ArgumentParser(add_help=True, allow_abbrev=False)
    ap.add_argument("outroot", help="Top-level output directory (e.g. bbdm_2HDMa_type1_case1_scan)")
    ap.add_argument("--sin", required=True, type=str, help="sin(theta) value(s), comma-separated (e.g. 0.01,0.1)")
    ap.add_argument("--tan", required=True, type=str, help="tan(beta) value(s), comma-separated (e.g. 0.01,0.1)")
    ap.add_argument("--mA", required=True, type=str, help="mA (=mass36 input) value(s) in GeV, comma-separated")
    ap.add_argument("--ma", required=True, type=str, help="ma (m55) value(s) in GeV, comma-separated")
    ap.add_argument("--mchi", required=True, type=str, help="mchi (mass of particle 52) value(s) in GeV, comma-separated")
    ap.add_argument("--cosbma", required=True, type=str, help="cos(beta-alpha) value(s), comma-separated (supports negatives)")
    ap.add_argument("--deltam", required=True, type=str, help="deltam = m35 - m36 value(s) in GeV, comma-separated (supports negatives)")
    ap.add_argument("--xmin", type=float, default=None, help="Optional x-axis lower bound")
    ap.add_argument("--xmax", type=float, default=None, help="Optional x-axis upper bound")

    # 处理负数逗号列表：--deltam -200,-150,... / --cosbma -0.1,-0.05,...
    argv = sys.argv[1:]
    fixed_argv: list[str] = []
    i = 0
    while i < len(argv):
        if argv[i] in ("--deltam", "--cosbma") and i + 1 < len(argv):
            nxt = argv[i + 1]
            if nxt.startswith("-") and ("," in nxt) and not nxt.startswith("--"):
                fixed_argv.append(f"{argv[i]}={nxt}")
                i += 2
                continue
        fixed_argv.append(argv[i])
        i += 1

    args = ap.parse_args(fixed_argv)
    if args.xmin is not None and args.xmax is not None and args.xmax < args.xmin:
        raise SystemExit("ERROR: --xmax must be >= --xmin")

    outroot = Path(args.outroot).expanduser()
    results_dir = outroot / "results"
    if not results_dir.exists():
        print(f"ERROR: results directory not found: {results_dir}")
        return 2

    root_exe = shutil.which("root")
    if not root_exe:
        print("ERROR: CERN ROOT executable 'root' not found on PATH.")
        print("Hint: conda activate mg5")
        return 2

    values_map = {
        "sin": _to_floats(_split_list(args.sin), "sin"),
        "tan": _to_floats(_split_list(args.tan), "tan"),
        "mA": _to_floats(_split_list(args.mA), "mA"),
        "ma": _to_floats(_split_list(args.ma), "ma"),
        "mchi": _to_floats(_split_list(args.mchi), "mchi"),
        "cosbma": _to_floats(_split_list(args.cosbma), "cosbma"),
        "deltam": _to_floats(_split_list(args.deltam), "deltam"),
    }

    vary = choose_varying(values_map)

    def _single(name: str) -> float:
        v = values_map[name]
        if len(v) != 1:
            raise SystemExit(f"ERROR: {name} must be a single value unless it is the varying parameter")
        return v[0]

    # Determine which parameters are encoded in Events_* directory names.
    event_dirs = sorted([p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith("Events_")])
    if not event_dirs:
        print(f"ERROR: no Events_* directories under {results_dir}")
        return 2

    recognized: set[str] = set()
    for d in event_dirs[:50]:
        recognized |= set(parse_events_dir_params(d.name).keys())

    fixed = FixedParams(
        sin=_single("sin") if vary != "sin" else 0.0,
        tan=_single("tan") if vary != "tan" else 0.0,
        mA=_single("mA") if vary != "mA" else 0.0,
        ma=_single("ma") if vary != "ma" else 0.0,
        mchi=_single("mchi") if vary != "mchi" else 0.0,
        cosbma=_single("cosbma") if vary != "cosbma" else 0.0,
        deltam=_single("deltam") if vary != "deltam" else 0.0,
    )

    # If varying parameter isn't encoded in directory names, refuse.
    if vary and vary not in recognized and len(values_map[vary]) > 1:
        raise SystemExit(
            f"ERROR: varying parameter '{vary}' is not encoded in Events_* directory names, cannot match. "
            f"Recognized: {', '.join(sorted(recognized)) or '(none)'}"
        )

    def match_dir_from_inputs(
        *, sinv: float, tanv: float, mAv: float, mav: float, mchiv: float, cosbma_v: float, deltam_v: float
    ) -> Path | None:
        want = {
            "sin": sinv,
            "tan": tanv,
            "mA": mAv,
            "ma": mav,
            "mchi": mchiv,
            "cosbma": cosbma_v,
            "deltam": deltam_v,
        }
        candidates: list[Path] = []
        for d in event_dirs:
            params = parse_events_dir_params(d.name)
            ok = True
            for k in recognized:
                if k not in want:
                    continue
                if k not in params:
                    ok = False
                    break
                if abs(params[k] - want[k]) > 1e-9:
                    ok = False
                    break
            if ok:
                candidates.append(d)
        if not candidates:
            return None
        return sorted(candidates)[0]

    points: list[tuple[float, float]] = []
    varying_values = values_map[vary]

    for vv in varying_values:
        sinv = vv if vary == "sin" else fixed.sin
        tanv = vv if vary == "tan" else fixed.tan
        mAv = vv if vary == "mA" else fixed.mA
        mav = vv if vary == "ma" else fixed.ma
        mchiv = vv if vary == "mchi" else fixed.mchi
        cbv = vv if vary == "cosbma" else fixed.cosbma
        dmv = vv if vary == "deltam" else fixed.deltam

        ev_dir = match_dir_from_inputs(sinv=sinv, tanv=tanv, mAv=mAv, mav=mav, mchiv=mchiv, cosbma_v=cbv, deltam_v=dmv)
        if not ev_dir:
            raise SystemExit(
                f"ERROR: cannot find Events dir for sin={sinv}, tan={tanv}, mA={mAv}, ma={mav}, mchi={mchiv}, cosbma={cbv}, deltam={dmv}"
            )

        banner = find_banner_file(ev_dir)
        if not banner:
            raise SystemExit(f"ERROR: cannot find *_banner.txt under {ev_dir}")

        xs_pb = extract_cross_section_pb(banner)
        if xs_pb is None:
            raise SystemExit(f"ERROR: cannot extract cross section from {banner}")

        points.append((vv, xs_pb))

    points = [(x, y) for x, y in points if _in_x_range(x, args.xmin, args.xmax)]

    if len(points) < 3:
        raise SystemExit(f"ERROR: need >=3 points in selected x range, got {len(points)}")

    # Output directory
    out_dir = outroot / "output_scatter_plot" / "crosssection"
    out_dir.mkdir(parents=True, exist_ok=True)

    fixed_tag = (
        f"st{_st_tag(values_map['sin'][0] if vary != 'sin' else points[0][0])}"
        f"_tb{_st_tag(values_map['tan'][0] if vary != 'tan' else points[0][0])}"
        f"_mA{int(values_map['mA'][0] if vary != 'mA' else points[0][0])}"
        f"_ma{int(values_map['ma'][0] if vary != 'ma' else points[0][0])}"
        f"_mchi{_st_tag(values_map['mchi'][0] if vary != 'mchi' else points[0][0])}"
        f"_cbma{_st_tag(values_map['cosbma'][0] if vary != 'cosbma' else points[0][0])}"
        f"_dm{_st_tag(values_map['deltam'][0] if vary != 'deltam' else points[0][0])}"
    )

    out_png = out_dir / f"crosssection_{vary}__{fixed_tag}.png"

    x_label = _x_label_for(vary)

    fixed_text = _fixed_text_for(
        vary,
        {
            "sin": values_map["sin"][0],
            "tan": values_map["tan"][0],
            "mA": values_map["mA"][0],
            "ma": values_map["ma"][0],
            "mchi": values_map["mchi"][0],
            "cosbma": values_map["cosbma"][0],
            "deltam": values_map["deltam"][0],
        },
        type_label=_infer_type_label(outroot),
    )

    # Write points to txt
    data_txt = out_dir / f".crosssection_{vary}_{hash((vary, fixed_tag)) & 0xffffffff:x}.dat"
    with open(data_txt, "w") as f:
        for x, y in points:
            f.write(f"{x} {y}\n")

    macro = out_dir / ".crosssection_tmp.C"
    write_root_macro(macro)

    rc, out = _run_root_make_scatter(
        root_exe=root_exe,
        macro=macro,
        data_txt=data_txt,
        x_label=x_label,
        fixed_text=fixed_text,
        out_png=out_png,
    )
    log_path = out_dir / (out_png.stem + ".log")
    log_path.write_text(out)

    if rc != 0:
        print(f"ERROR: ROOT failed (see {log_path})")
        return rc

    # Cleanup temps on success
    try:
        log_path.unlink(missing_ok=True)
    except Exception:
        pass
    try:
        data_txt.unlink(missing_ok=True)
    except Exception:
        pass
    try:
        macro.unlink(missing_ok=True)
    except Exception:
        pass

    compiled_base = macro.name.replace(".C", "_C")
    compiled_bases = [compiled_base]
    if compiled_base.startswith("."):
        compiled_bases.append("_" + compiled_base[1:])
    for base in compiled_bases:
        for p in out_dir.glob(f"{base}*"):
            try:
                if p.is_file():
                    p.unlink(missing_ok=True)
            except Exception:
                pass

    print(f"Wrote: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
