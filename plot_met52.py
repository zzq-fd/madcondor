#!/usr/bin/env python3
"""Overlay normalized gen-pT^miss curves from MadGraph ROOT outputs using particle 52 pT.

This script mirrors plot_met55_overlay.py in interface and behavior, but defines
"gen-pTmiss" as the magnitude of the vector sum of transverse momenta of
particle 52 and -52:
    gen-pTmiss = | \vec{p}_{T}(52) + \vec{p}_{T}(-52) |

Scan modes:
  1) Single-parameter scan (default):
     If exactly one of {sin, tan, mA, ma, mchi, cosbma, deltam} is provided with multiple values,
     then that parameter is treated as the scan variable (each value becomes one curve).

  2) Paired two-parameter scan:
     If exactly TWO parameters are provided with multiple values, and the two lists have the
     same length, then values are paired by index (NOT a Cartesian product). Example:
       --ma 200,300 --mA 300,400  => curves: (ma=200,mA=300) and (ma=300,mA=400)

Requirements:
  - CERN ROOT executable `root` available on PATH.
  - Does NOT require PyROOT.

Usage:
    python3 plot_met52.py <OUTROOT> --sin <v[,v2,...]> --tan <v[,v2,...]> --mA <v[,v2,...]> --ma <v[,v2,...]> --mchi <v[,v2,...]> --cosbma <v[,v2,...]> --deltam <v[,v2,...]>

Output:
  - Writes a PNG into <OUTROOT>/output_scatter_plot/met52/.
  - Cleans up temporary ROOT macro artifacts on success.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PointKey:
    sintheta: float
    tanbeta: float
    mA: float
    ma: float
    mchi: float
    cosbma: float
    deltam: float


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
    s = f"{x:g}".replace("-", "m").replace(".", "p")
    return s


def _tag_to_float(tag: str) -> float:
    neg = tag.startswith("m")
    if neg:
        tag = tag[1:]
    v = float(tag.replace("p", "."))
    return -v if neg else v


def parse_events_dir_params(dirname: str) -> dict[str, float]:
    """Parse parameters encoded in Events_* directory name."""

    name = dirname.replace("Events_", "")
    out: dict[str, float] = {}

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


def find_unweighted_root(events_dir: Path) -> Path | None:
    run_dirs = sorted([p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])
    if not run_dirs:
        return None
    preferred = next((p for p in run_dirs if p.name == "run_01"), None)
    chosen = preferred or run_dirs[0]
    f = chosen / "unweighted_events.root"
    return f if f.exists() else None


def choose_varying(values_map: dict[str, list[float]]) -> list[str]:
    varying = [k for k, v in values_map.items() if len(v) > 1]
    if len(varying) > 2:
        raise SystemExit(
            "ERROR: please provide multiple values for at most TWO parameters among sin,tan,mA,ma,mchi,cosbma,deltam. "
            f"Got multiple for: {', '.join(varying)}"
        )
    if len(varying) == 2:
        a, b = varying
        if len(values_map[a]) != len(values_map[b]):
            raise SystemExit(
                "ERROR: when providing multiple values for two parameters, the two lists must have the same length "
                "(paired scan). "
                f"Got {a} has {len(values_map[a])} values, {b} has {len(values_map[b])} values."
            )
    return varying


def write_root_macro(macro_path: Path) -> None:
    macro_path.write_text(
        r'''#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Curve {
  std::string path;
  std::string label;
  int color;
};

static bool read_config(const char* cfg, std::vector<Curve>& curves) {
  std::ifstream in(cfg);
  if (!in) return false;

  std::string path, label;
  int color;
  while (std::getline(in, path)) {
    if (path.empty()) continue;
    if (path[0] == '#') continue;

    // Format: path\tlabel\tcolor
    size_t p1 = path.find('\t');
    if (p1 == std::string::npos) continue;
    size_t p2 = path.find('\t', p1 + 1);
    if (p2 == std::string::npos) continue;

    std::string fpath = path.substr(0, p1);
    std::string flabel = path.substr(p1 + 1, p2 - (p1 + 1));
    std::string scolor = path.substr(p2 + 1);

    try {
      color = std::stoi(scolor);
    } catch (...) {
      color = 1;
    }

    curves.push_back({fpath, flabel, color});
  }

  return !curves.empty();
}

static double compute_max_pt(TTree* t) {
    const char* expr = "sqrt(pow(Sum$(Particle.Px*(abs(Particle.PID)==52)),2) + pow(Sum$(Particle.Py*(abs(Particle.PID)==52)),2))";
    const char* sel_evt = "Sum$(abs(Particle.PID)==52)>0";
    Long64_t n = t->Draw(expr, sel_evt, "goff");
  if (n <= 0) return 0.0;
  const double* v = t->GetV1();
  double vmax = 0.0;
  for (Long64_t i = 0; i < n; ++i) {
    if (v[i] > vmax) vmax = v[i];
  }
  return vmax;
}

static void update_y_range(TH1D* frame, const std::vector<TH1D*>& hs) {
  double ymax = 0.0;
  double ymin_pos = 0.0;
  for (auto* h : hs) {
    if (!h) continue;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double v = h->GetBinContent(b);
      if (v > ymax) ymax = v;
      if (v > 0.0 && (ymin_pos == 0.0 || v < ymin_pos)) ymin_pos = v;
    }
  }

  if (ymax <= 0.0) {
    frame->SetMinimum(1e-4);
    frame->SetMaximum(1.0);
    return;
  }

  double ymin = (ymin_pos > 0.0) ? (ymin_pos / 2.0) : (ymax * 1e-4);
  if (ymin <= 0.0) ymin = ymax * 1e-4;
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax * 5.0);
}

int plot_overlay(const char* cfg_txt,
                 const char* out_png,
                 const char* vary_name,
                 const char* fixed_text,
                 double xmin,
                 double xmax,
                 int nbins)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  std::vector<Curve> curves;
  if (!read_config(cfg_txt, curves)) {
    std::cerr << "ERROR: cannot read config: " << cfg_txt << std::endl;
    return 2;
  }

  // Determine xmax if requested
  double global_max = 0.0;
  for (const auto& c : curves) {
    TFile f(c.path.c_str(), "READ");
    if (f.IsZombie()) continue;
    auto* t = dynamic_cast<TTree*>(f.Get("LHEF"));
    if (!t) continue;
    global_max = std::max(global_max, compute_max_pt(t));
  }
  if (xmax <= 0.0) {
    xmax = (global_max > 0.0) ? (std::ceil(global_max / 50.0) * 50.0) : (xmin + 1000.0);
  }
  if (xmax <= xmin + 1.0) xmax = xmin + 1000.0;

  TCanvas c("c", "c", 1100, 800);
  c.SetLeftMargin(0.14);
  c.SetBottomMargin(0.12);
  c.SetLogy(true);

  c.cd();

  // Avoid TLegend (can crash in some batch setups); draw a simple legend by hand.

  TH1D* hframe = nullptr;
  std::vector<TH1D*> hs;
  std::vector<Curve> drawn;
  hs.reserve(curves.size());
  drawn.reserve(curves.size());

  for (size_t i = 0; i < curves.size(); ++i) {
    const auto& cv = curves[i];

    TFile f(cv.path.c_str(), "READ");
    if (f.IsZombie()) {
      std::cerr << "ERROR: cannot open: " << cv.path << std::endl;
      continue;
    }
    auto* t = dynamic_cast<TTree*>(f.Get("LHEF"));
    if (!t) {
      std::cerr << "ERROR: missing TTree LHEF in: " << cv.path << std::endl;
      continue;
    }

    std::string hname = "h" + std::to_string(i);
    auto* h = new TH1D(hname.c_str(), "", nbins, xmin, xmax);
    h->Sumw2(false);
    h->SetDirectory(nullptr);

        const char* expr = "sqrt(pow(Sum$(Particle.Px*(abs(Particle.PID)==52)),2) + pow(Sum$(Particle.Py*(abs(Particle.PID)==52)),2))";
        const char* sel_evt = "Sum$(abs(Particle.PID)==52)>0";
        Long64_t n = t->Draw(expr, sel_evt, "goff");
    const double* v = t->GetV1();
    for (Long64_t j = 0; j < n; ++j) {
      h->Fill(v[j]);
    }

    // Normalize to unit area over the plotted range.
    double integral = h->Integral();
    if (integral > 0.0) h->Scale(1.0 / integral);

    h->SetLineColor(cv.color);
    h->SetLineWidth(2);
    h->SetFillStyle(0);

    if (!hframe) {
      hframe = (TH1D*)h->Clone("frame");
      hframe->Reset("ICES");
      hframe->SetDirectory(nullptr);
      hframe->SetTitle("");
      hframe->GetXaxis()->SetTitle("gen-p_{T}^{miss} (GeV)");
      hframe->GetYaxis()->SetTitle("arbitrary units");
      hframe->GetYaxis()->SetTitleOffset(1.25);
      hframe->Draw();
    }

    h->Draw("HIST SAME");
    hs.push_back(h);
    drawn.push_back(cv);
  }

  if (hframe) update_y_range(hframe, hs);

  // CMS-style text and manual legend
  TLatex latex;
  latex.SetNDC(true);

  latex.SetTextFont(61);
  latex.SetTextSize(0.055);
    latex.DrawLatex(0.16, 0.92, "CMS Internal");

  latex.SetTextFont(42);
  latex.SetTextSize(0.045);
  latex.SetTextAlign(31);
  latex.DrawLatex(0.90, 0.92, "(13.6 TeV)");

  // Fixed parameter block (left)
  latex.SetTextAlign(13);
  latex.SetTextFont(42);
  latex.SetTextSize(0.038);
  {
    // Header line for the left-side annotation block (placed above mA line).
    latex.SetTextSize(0.040);
    latex.DrawLatex(0.26, 0.88, "2HDM+a     b#bar{b} + #chi#bar{#chi}");
    latex.SetTextSize(0.038);

    std::string ft = fixed_text ? fixed_text : "";
    std::vector<std::string> lines;
    size_t start = 0;
    while (start < ft.size()) {
      size_t pos = ft.find(';', start);
      std::string s = (pos == std::string::npos) ? ft.substr(start) : ft.substr(start, pos - start);
      // trim leading spaces
      while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(s.begin());
      if (!s.empty()) lines.push_back(s);
      if (pos == std::string::npos) break;
      start = pos + 1;
    }
    double y = 0.83;
    for (size_t i = 0; i < lines.size(); ++i) {
      latex.DrawLatex(0.26, y, lines[i].c_str());
      y -= 0.055;
    }
  }

  // Manual legend (top-right)
  double lx1 = 0.64;
  double lx2 = 0.72;
  double ly0 = 0.86;
  double ldy = 0.055;
  latex.SetTextAlign(12);
  latex.SetTextFont(42);
  latex.SetTextSize(0.040);
  for (size_t i = 0; i < drawn.size(); ++i) {
    if (i >= hs.size()) break;
    double y = ly0 - i * ldy;
    auto* line = new TLine();
    line->SetLineColor(drawn[i].color);
    line->SetLineWidth(2);
    line->DrawLineNDC(lx1, y, lx2, y);
    latex.DrawLatex(lx2 + 0.015, y, drawn[i].label.c_str());
  }

  c.Modified();
  c.Update();
  c.SaveAs(out_png);

  gSystem->Exit(0);
  return 0;
}
'''
    )


def root_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True, allow_abbrev=False)
    ap.add_argument("outroot", help="Top-level output directory (e.g. bbdm_2HDMa_type1_case1_scan)")
    ap.add_argument("--sin", required=True, type=str, help="sin(theta) value(s), comma-separated (e.g. 0.01,0.1)")
    ap.add_argument("--tan", required=True, type=str, help="tan(beta) value(s), comma-separated (e.g. 0.01,0.1)")
    ap.add_argument("--mA", required=True, type=str, help="mA value(s) in GeV, comma-separated")
    ap.add_argument("--ma", required=True, type=str, help="ma (m55) value(s) in GeV, comma-separated")
    ap.add_argument("--mchi", required=True, type=str, help="mchi (mass of particle 52) value(s) in GeV, comma-separated")
    ap.add_argument("--cosbma", required=True, type=str, help="cos(beta-alpha) value(s), comma-separated (supports negatives)")
    ap.add_argument("--deltam", required=True, type=str, help="deltam = m35 - m36 value(s) in GeV, comma-separated (supports negatives)")
    ap.add_argument("--xmin", type=float, default=150.0, help="x-axis minimum (default: 150)")
    ap.add_argument("--xmax", type=float, default=0.0, help="x-axis maximum (0 means auto)")
    ap.add_argument("--nbins", type=int, default=20, help="number of bins (default: 20)")

    argv = sys.argv[1:]
    fixed_argv: list[str] = []
    i = 0
    while i < len(argv):
        if argv[i] == "--deltam" and i + 1 < len(argv):
            nxt = argv[i + 1]
            if nxt.startswith("-") and ("," in nxt) and not nxt.startswith("--"):
                fixed_argv.append(f"--deltam={nxt}")
                i += 2
                continue
        fixed_argv.append(argv[i])
        i += 1

    args = ap.parse_args(fixed_argv)

    # Preserve the order of scan variables as provided on CLI (for paired scan labeling).
    flag_to_name = {
        "--sin": "sin",
        "--tan": "tan",
        "--mA": "mA",
        "--ma": "ma",
        "--mchi": "mchi",
        "--cosbma": "cosbma",
        "--deltam": "deltam",
    }
    arg_order: list[str] = []
    for tok in fixed_argv:
        if not tok.startswith("--"):
            continue
        for flag, name in flag_to_name.items():
            if tok == flag or tok.startswith(flag + "="):
                if name not in arg_order:
                    arg_order.append(name)
                break

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

    varying_params = choose_varying(values_map)
    if len(varying_params) == 2:
        varying_params = sorted(varying_params, key=lambda n: (arg_order.index(n) if n in arg_order else 999))

    def _is_varying(name: str) -> bool:
        return name in varying_params

    def _single(name: str) -> float:
        v = values_map[name]
        if len(v) != 1:
            raise SystemExit(f"ERROR: {name} must be a single value unless it is a varying parameter")
        return v[0]

    event_dirs = sorted([p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith("Events_")])
    if not event_dirs:
        print(f"ERROR: no Events_* directories under {results_dir}")
        return 2

    recognized: set[str] = set()
    for d in event_dirs[:20]:
        recognized |= set(parse_events_dir_params(d.name).keys())

    fixed_sin = _single("sin") if not _is_varying("sin") else None
    fixed_tan = _single("tan") if not _is_varying("tan") else None
    fixed_mA = _single("mA") if not _is_varying("mA") else None
    fixed_ma = _single("ma") if not _is_varying("ma") else None
    fixed_mchi = _single("mchi") if not _is_varying("mchi") else None
    fixed_cosbma = _single("cosbma") if not _is_varying("cosbma") else None
    fixed_deltam = _single("deltam") if not _is_varying("deltam") else None

    for vp in varying_params:
        if vp not in recognized and len(values_map[vp]) > 1:
            raise SystemExit(
                f"ERROR: varying parameter '{vp}' is not encoded in Events_* directory names, so it cannot be scanned/overlaid. "
                f"Recognized from directory names: {', '.join(sorted(recognized)) or '(none)'}"
            )

    legend_sym = {
        "sin": "sin#theta",
        "tan": "tan#beta",
        "mA": "m_{A}",
        "ma": "m_{a}",
        "mchi": "m_{#chi}",
        "cosbma": "cos(#beta-#alpha)",
        "deltam": "#Delta m",
    }

    colors = [
        600 + 1,
        800 + 1,
        632 + 1,
        416 + 2,
        880 + 1,
        920 + 1,
        432 + 2,
    ]

    curve_specs: list[tuple[Path, str, int]] = []

    def match_dir_from_inputs(
        sinv: float,
        tanv: float,
        mAv: float,
        mav: float,
        mchiv: float,
        cosbma_v: float,
        deltam_v: float,
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

    def _unit_for(name: str) -> str:
        return " GeV" if name in ("mA", "ma", "mchi", "deltam") else ""

    def _format_pair_label(a_name: str, a_val: float, b_name: str, b_val: float) -> str:
        return (
            f"{legend_sym[a_name]} = {a_val:g}{_unit_for(a_name)}, "
            f"{legend_sym[b_name]} = {b_val:g}{_unit_for(b_name)}"
        )

    def _format_single_label(vary_name: str, vv: float) -> str:
        return f"{legend_sym[vary_name]} = {vv:g}{_unit_for(vary_name)}"

    if len(varying_params) == 0:
        sinv = fixed_sin
        tanv = fixed_tan
        mAv = fixed_mA
        mav = fixed_ma
        mchiv = fixed_mchi
        cosbma_v = fixed_cosbma
        deltam_v = fixed_deltam

        if (
            sinv is None
            or tanv is None
            or mAv is None
            or mav is None
            or mchiv is None
            or cosbma_v is None
            or deltam_v is None
        ):
            raise SystemExit("ERROR: internal: missing fixed values")

        ev_dir = match_dir_from_inputs(
            sinv=sinv,
            tanv=tanv,
            mAv=mAv,
            mav=mav,
            mchiv=mchiv,
            cosbma_v=cosbma_v,
            deltam_v=deltam_v,
        )
        if not ev_dir:
            raise SystemExit(
                f"ERROR: could not find matching Events directory for sin={sinv}, tan={tanv}, mA={mAv}, ma={mav}, mchi={mchiv}, cosbma={cosbma_v}, deltam={deltam_v}"
            )
        root_file = find_unweighted_root(ev_dir)
        if not root_file:
            raise SystemExit(f"ERROR: missing unweighted_events.root under {ev_dir}")
        curve_specs.append((root_file, "point", colors[0]))

    elif len(varying_params) == 1:
        vary = varying_params[0]
        varying_values = values_map[vary]
        for idx, vv in enumerate(varying_values):
            sinv = vv if vary == "sin" else fixed_sin
            tanv = vv if vary == "tan" else fixed_tan
            mAv = vv if vary == "mA" else fixed_mA
            mav = vv if vary == "ma" else fixed_ma
            mchiv = vv if vary == "mchi" else fixed_mchi
            cosbma_v = vv if vary == "cosbma" else fixed_cosbma
            deltam_v = vv if vary == "deltam" else fixed_deltam

            if (
                sinv is None
                or tanv is None
                or mAv is None
                or mav is None
                or mchiv is None
                or cosbma_v is None
                or deltam_v is None
            ):
                raise SystemExit("ERROR: internal: missing fixed values")

            ev_dir = match_dir_from_inputs(
                sinv=sinv,
                tanv=tanv,
                mAv=mAv,
                mav=mav,
                mchiv=mchiv,
                cosbma_v=cosbma_v,
                deltam_v=deltam_v,
            )
            if not ev_dir:
                raise SystemExit(
                    f"ERROR: could not find matching Events directory for sin={sinv}, tan={tanv}, mA={mAv}, ma={mav}, mchi={mchiv}, cosbma={cosbma_v}, deltam={deltam_v}"
                )

            root_file = find_unweighted_root(ev_dir)
            if not root_file:
                raise SystemExit(f"ERROR: missing unweighted_events.root under {ev_dir}")

            label = _format_single_label(vary, vv)
            color = colors[idx % len(colors)]
            curve_specs.append((root_file, label, color))

    else:
        a, b = varying_params
        a_vals = values_map[a]
        b_vals = values_map[b]
        for idx in range(len(a_vals)):
            av = a_vals[idx]
            bv = b_vals[idx]

            sinv = av if a == "sin" else (bv if b == "sin" else fixed_sin)
            tanv = av if a == "tan" else (bv if b == "tan" else fixed_tan)
            mAv = av if a == "mA" else (bv if b == "mA" else fixed_mA)
            mav = av if a == "ma" else (bv if b == "ma" else fixed_ma)
            mchiv = av if a == "mchi" else (bv if b == "mchi" else fixed_mchi)
            cosbma_v = av if a == "cosbma" else (bv if b == "cosbma" else fixed_cosbma)
            deltam_v = av if a == "deltam" else (bv if b == "deltam" else fixed_deltam)

            if (
                sinv is None
                or tanv is None
                or mAv is None
                or mav is None
                or mchiv is None
                or cosbma_v is None
                or deltam_v is None
            ):
                raise SystemExit("ERROR: internal: missing fixed values")

            ev_dir = match_dir_from_inputs(
                sinv=sinv,
                tanv=tanv,
                mAv=mAv,
                mav=mav,
                mchiv=mchiv,
                cosbma_v=cosbma_v,
                deltam_v=deltam_v,
            )
            if not ev_dir:
                raise SystemExit(
                    f"ERROR: could not find matching Events directory for sin={sinv}, tan={tanv}, mA={mAv}, ma={mav}, mchi={mchiv}, cosbma={cosbma_v}, deltam={deltam_v}"
                )

            root_file = find_unweighted_root(ev_dir)
            if not root_file:
                raise SystemExit(f"ERROR: missing unweighted_events.root under {ev_dir}")

            label = _format_pair_label(a, av, b, bv)
            color = colors[idx % len(colors)]
            curve_specs.append((root_file, label, color))

    fixed_lines: list[str] = []
    if "mA" not in varying_params:
        fixed_lines.append(f"m_{{A}} = {fixed_mA:g} GeV")
    if "mchi" not in varying_params:
        fixed_lines.append(f"m_{{#chi}} = {fixed_mchi:g} GeV")
    if "deltam" not in varying_params:
        fixed_lines.append(f"#Delta m = m_{{H^{{+}}}} - m_{{A}} = {fixed_deltam:g} GeV")
    if "cosbma" not in varying_params:
        fixed_lines.append(f"cos(#beta-#alpha) = {fixed_cosbma:g}")
    if "tan" not in varying_params:
        fixed_lines.append(f"tan#beta = {fixed_tan:g}")
    if "sin" not in varying_params:
        fixed_lines.append(f"sin#theta = {fixed_sin:g}")
    if "ma" not in varying_params:
        fixed_lines.append(f"m_{{a}} = {fixed_ma:g} GeV")

    fixed_text = "; ".join(fixed_lines)

    def _tag_or_var(x: float | None, *, kind: str) -> str:
        if x is None:
            return "var"
        if kind == "int":
            return str(int(x))
        return _st_tag(x)

    fixed_tag = (
        f"st{_tag_or_var(fixed_sin, kind='float')}"
        f"_tb{_tag_or_var(fixed_tan, kind='float')}"
        f"_mA{_tag_or_var(fixed_mA, kind='int')}"
        f"_ma{_tag_or_var(fixed_ma, kind='int')}"
        f"_mchi{_tag_or_var(fixed_mchi, kind='float')}"
        f"_cbma{_tag_or_var(fixed_cosbma, kind='float')}"
        f"_dm{_tag_or_var(fixed_deltam, kind='float')}"
    )

    if len(varying_params) == 0:
        vary_tag = "single"
    elif len(varying_params) == 1:
        vary_tag = varying_params[0]
    else:
        vary_tag = f"{varying_params[0]}_{varying_params[1]}_paired"

    out_dir = outroot / "output_scatter_plot" / "met52"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"met52_{vary_tag}__{fixed_tag}.png"

    cfg = out_dir / ".met52_cfg.txt"
    macro = out_dir / ".met52_tmp.C"
    write_root_macro(macro)

    cfg_lines = ["# path\tlabel\tcolor"]
    for fpath, label, color in curve_specs:
        cfg_lines.append(f"{fpath}\t{label}\t{color}")
    cfg.write_text("\n".join(cfg_lines) + "\n")

    fixed_text_escaped = fixed_text.replace("\\", "\\\\").replace('"', '\\"')
    stmt = (
        f'gROOT->ProcessLine(".L {root_quote(macro)}"); '
        f'plot_overlay("{root_quote(cfg)}","{root_quote(out_png)}","{vary_tag}","{fixed_text_escaped}",{args.xmin},{args.xmax},{args.nbins});'
    )
    cmd = [root_exe, "-l", "-b", "-q", "-e", stmt]

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log_path = out_dir / (out_png.stem + ".log")
    log_path.write_text(proc.stdout)

    if proc.returncode != 0:
        print(f"ERROR: ROOT failed (see {log_path})")
        return proc.returncode

    # Cleanup temps on success
    try:
        log_path.unlink(missing_ok=True)
    except Exception:
        pass
    try:
        cfg.unlink(missing_ok=True)
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
