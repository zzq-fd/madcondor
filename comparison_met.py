#!/usr/bin/env python3
"""Overlay normalized gen-pT^miss curves for type1/type2 on one plot.

This follows the definition in plot_met52.py:
        gen-pTmiss = | \vec{p}_T(52) + \vec{p}_T(-52) |

Inputs:
    - Either a single token like `case1_scan2` (the script finds
      `bbdm_2HDMa_type1_<token>` and `bbdm_2HDMa_type2_<token>`),
      or two explicit OUTROOT folders (backward compatible).
    - Parameter selection flags like plot_met52.py, but for comparison you should pass a SINGLE value
        for each parameter (no scan lists), so the plot contains exactly two curves.

Usage example:
    python3 comparison_met.py case1_scan2 \
        --sin 0.9 --tan 3 --mA 200 --ma 300 --mchi 10 --cosbma 0 --deltam 0

Output:
    - ./comparison/met/<case*_scan*>_<param-tags>.png
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

# --- 参数解析 ---
def parse_args():
    parser = argparse.ArgumentParser(description="Overlay gen-pTmiss curves for type1/type2 scans.")
    parser.add_argument(
        "inputs",
        nargs="+",
        help=(
            "Either: <case*_scan*> token (e.g. case1_scan2) OR: <TYPE1_OUTROOT> <TYPE2_OUTROOT>. "
            "If a token is given, folders are resolved as bbdm_2HDMa_type1_<token> and bbdm_2HDMa_type2_<token>."
        ),
    )
    parser.add_argument("--sin", type=str, default=None, help="single sin(theta), e.g. 0.01")
    parser.add_argument("--tan", type=str, default=None, help="single tan(beta), e.g. 0.01")
    parser.add_argument("--mA", type=str, default=None, help="single mA in GeV")
    parser.add_argument("--ma", type=str, default=None, help="single ma in GeV")
    parser.add_argument("--mchi", type=str, default=None, help="single mchi in GeV")
    parser.add_argument("--cosbma", type=str, default=None, help="single cos(beta-alpha)")
    parser.add_argument("--deltam", type=str, default=None, help="single deltam in GeV")
    parser.add_argument("--xmin", type=float, default=150.0, help="x-axis lower bound (default: 150)")
    parser.add_argument("--xmax", type=float, default=0.0, help="x-axis upper bound (0 means auto)")
    return parser.parse_args()

# --- 目录参数匹配 ---
def _tag_to_float(tag: str) -> float:
    neg = tag.startswith("m")
    if neg:
        tag = tag[1:]
    v = float(tag.replace("p", "."))
    return -v if neg else v


def parse_events_dir_params(dirname: str) -> dict[str, float]:
    name = dirname.replace("Events_", "")
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


def find_banner_file(events_dir: Path) -> Path | None:
    run_dirs = sorted([p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])
    for run_dir in run_dirs:
        for banner in sorted(run_dir.glob("*_banner.txt")):
            return banner
    return None


def extract_mchi_from_banner(banner_path: Path) -> float | None:
    txt = banner_path.read_text(errors="ignore")
    m = re.search(r"^\s*52\s+([+\-0-9.eE]+)\b", txt, flags=re.MULTILINE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None

# --- 参数筛选 ---

def _parse_single(value: str | None, name: str, default: float | None = None) -> float:
    if value is None:
        if default is None:
            raise SystemExit(f"ERROR: missing required --{name}")
        return float(default)
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if len(parts) != 1:
        raise SystemExit(f"ERROR: comparison_met expects a single value for --{name} (got {value})")
    return float(parts[0])


def match_params(params: dict[str, float], target: dict[str, float]) -> bool:
    for k in ["sin","tan","mA","ma","mchi","cosbma","deltam"]:
        if k not in params:
            return False
        if float(params[k]) != float(target[k]):
            return False
    return True

# --- 读取 met52 pt 分布 ---
def find_matching_events_dir(outroot: Path, target: dict[str, float]) -> Path | None:
    results_dir = outroot / "results"
    if not results_dir.is_dir():
        return None
    for d in sorted([p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith("Events_")]):
        params = parse_events_dir_params(d.name)
        # Legacy naming may omit some tags; infer mchi from banner when possible.
        if "mchi" not in params:
            banner = find_banner_file(d)
            if banner:
                mch = extract_mchi_from_banner(banner)
                if mch is not None:
                    params["mchi"] = mch

        # If still missing, assume defaults (legacy scans usually used 0 for these).
        params.setdefault("mchi", 0.0)
        params.setdefault("cosbma", 0.0)
        params.setdefault("deltam", 0.0)
        if match_params(params, target):
            return d
    return None


def _st_tag(x: float) -> str:
    return f"{x:g}".replace("-", "m").replace(".", "p")


def _mass_tag(x: float) -> str:
    if float(x).is_integer():
        return str(int(x))
    return _st_tag(float(x))


def _output_param_tag(target: dict[str, float]) -> str:
    return (
        f"st{_st_tag(target['sin'])}"
        f"_tb{_st_tag(target['tan'])}"
        f"_mA{_mass_tag(target['mA'])}"
        f"_ma{_mass_tag(target['ma'])}"
        f"_mchi{_mass_tag(target['mchi'])}"
        f"_cbma{_st_tag(target['cosbma'])}"
        f"_dm{_st_tag(target['deltam'])}"
    )


def find_unweighted_root(events_dir: Path) -> Path | None:
    run_dirs = sorted([p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])
    if not run_dirs:
        return None
    preferred = next((p for p in run_dirs if p.name == "run_01"), None)
    chosen = preferred or run_dirs[0]
    f = chosen / "unweighted_events.root"
    return f if f.exists() else None


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
#include <iostream>
#include <string>

static void update_y_range(TH1D* frame, const TH1D* h1, const TH1D* h2) {
    if (!frame) return;
    double ymax = 0.0;
    double ymin_pos = 0.0;
    const TH1D* hs[2] = {h1, h2};
    for (const auto* h : hs) {
        if (!h) continue;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double v = h->GetBinContent(b);
            if (v > ymax) ymax = v;
            if (v > 0.0 && (ymin_pos == 0.0 || v < ymin_pos)) ymin_pos = v;
        }
    }
    if (ymax <= 0.0) {
        frame->SetMinimum(1e-5);
        frame->SetMaximum(1.0);
        return;
    }
    double ymin = (ymin_pos > 0.0) ? (ymin_pos / 2.0) : (ymax * 1e-5);
    if (ymin <= 0.0) ymin = ymax * 1e-5;
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax * 5.0);
}

static double compute_max_ptmiss(TTree* t) {
  const char* expr = "sqrt(pow(Sum$(Particle.Px*(abs(Particle.PID)==52)),2) + pow(Sum$(Particle.Py*(abs(Particle.PID)==52)),2))";
  const char* sel  = "Sum$(abs(Particle.PID)==52)>0";
  Long64_t n = t->Draw(expr, sel, "goff");
  if (n <= 0) return 0.0;
  const double* v = t->GetV1();
  double vmax = 0.0;
  for (Long64_t i=0;i<n;++i) vmax = std::max(vmax, v[i]);
  return vmax;
}

static TH1D* fill_hist(const char* name, TTree* t, double xmin, double xmax, int nbins) {
  const char* expr = "sqrt(pow(Sum$(Particle.Px*(abs(Particle.PID)==52)),2) + pow(Sum$(Particle.Py*(abs(Particle.PID)==52)),2))";
  const char* sel  = "Sum$(abs(Particle.PID)==52)>0";
  Long64_t n = t->Draw(expr, sel, "goff");
    auto* h = new TH1D(name, "", nbins, xmin, xmax);
    h->Sumw2(false);
    h->SetDirectory(nullptr);
  if (n <= 0) return h;
  const double* v = t->GetV1();
  for (Long64_t i=0;i<n;++i) h->Fill(v[i]);
  const double integ = h->Integral();
  if (integ > 0) h->Scale(1.0/integ);
  return h;
}

int compare_met52(const char* f1, const char* f2, const char* out_png, const char* fixed_text, double xmin_in, double xmax_in) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile a(f1, "READ");
  TFile b(f2, "READ");
  if (a.IsZombie() || b.IsZombie()) {
    std::cerr << "ERROR: cannot open input root files" << std::endl;
    return 2;
  }
  auto* t1 = dynamic_cast<TTree*>(a.Get("LHEF"));
  auto* t2 = dynamic_cast<TTree*>(b.Get("LHEF"));
  if (!t1 || !t2) {
    std::cerr << "ERROR: missing TTree LHEF" << std::endl;
    return 2;
  }

    const double xmax_auto = std::ceil(std::max(compute_max_ptmiss(t1), compute_max_ptmiss(t2)) / 50.0) * 50.0;
    const double xmin = xmin_in;
    double xhi = (xmax_in > 0.0) ? xmax_in : ((xmax_auto > 0.0) ? xmax_auto : (xmin + 1000.0));
    if (xhi <= xmin + 1.0) xhi = xmin + 1000.0;
  const int nbins = 50;

  auto* h1 = fill_hist("h1", t1, xmin, xhi, nbins);
  auto* h2 = fill_hist("h2", t2, xmin, xhi, nbins);
    h1->SetLineColor(kBlue+1);
    h2->SetLineColor(kRed+1);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h1->SetFillStyle(0);
    h2->SetFillStyle(0);

  TCanvas c("c", "c", 1100, 800);
  c.SetLeftMargin(0.14);
  c.SetBottomMargin(0.12);
  c.SetLogy(true);

    // Frame
    auto* frame = (TH1D*)h1->Clone("frame");
    frame->Reset("ICES");
    frame->SetDirectory(nullptr);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("gen-p_{T}^{miss} (GeV)");
    frame->GetYaxis()->SetTitle("arbitrary units");
    frame->GetYaxis()->SetTitleOffset(1.25);
    frame->Draw();

    h1->Draw("HIST SAME");
    h2->Draw("HIST SAME");
    update_y_range(frame, h1, h2);

    // CMS-style text and manual legend (match plot_met55_overlay.py)
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
    latex.SetTextSize(0.040);
    latex.DrawLatex(0.26, 0.88, "2HDM+a     b#bar{b} + #chi#bar{#chi}");
    latex.SetTextSize(0.038);
    {
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
        double y = 0.83;
        for (size_t i = 0; i < lines.size(); ++i) {
            latex.DrawLatex(0.26, y, lines[i].c_str());
            y -= 0.055;
        }
    }

    // Manual legend (top-right)
    double lx1 = 0.64;
    double lx2 = 0.72;
    double ly = 0.86;
    double dy = 0.055;
    TLine line;
    line.SetLineWidth(3);

    latex.SetTextFont(42);
    latex.SetTextSize(0.040);
    latex.SetTextAlign(12);

    line.SetLineColor(h1->GetLineColor());
    line.DrawLineNDC(lx1, ly, lx2, ly);
    latex.DrawLatex(lx2 + 0.02, ly - 0.015, "type1");
    ly -= dy;
    line.SetLineColor(h2->GetLineColor());
    line.DrawLineNDC(lx1, ly, lx2, ly);
    latex.DrawLatex(lx2 + 0.02, ly - 0.015, "type2");

  c.SaveAs(out_png);
  gSystem->Exit(0);
  return 0;
}
'''
    )


def run_root_compare(
    root_exe: str,
    macro: Path,
    f1: Path,
    f2: Path,
    out_png: Path,
    fixed_text: str,
    xmin: float,
    xmax: float,
) -> tuple[int, str]:
    safe_fixed = fixed_text.replace('\\', '\\\\').replace('"', '\\"')
    stmt = (
        f'gROOT->ProcessLine(".L {macro.as_posix()}"); '
        f'compare_met52("{f1.as_posix()}","{f2.as_posix()}","{out_png.as_posix()}","{safe_fixed}",{xmin},{xmax});'
    )
    cmd = [root_exe, "-l", "-b", "-q", "-e", stmt]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc.returncode, proc.stdout

# --- 主流程 ---
def main():
    args = parse_args()
    if args.xmax > 0 and args.xmax < args.xmin:
        raise SystemExit("ERROR: --xmax must be >= --xmin")

    if len(args.inputs) == 1:
        token = args.inputs[0]
        type1_name = f"bbdm_2HDMa_type1_{token}"
        type2_name = f"bbdm_2HDMa_type2_{token}"
        out1 = Path(type1_name)
        out2 = Path(type2_name)
        case_token = token
    elif len(args.inputs) == 2:
        out1 = Path(args.inputs[0])
        out2 = Path(args.inputs[1])
        n1 = out1.name
        n2 = out2.name
        m1 = re.match(r"^bbdm_2HDMa_type1_(.+)$", n1)
        m2 = re.match(r"^bbdm_2HDMa_type2_(.+)$", n2)
        if m1 and m2 and m1.group(1) == m2.group(1):
            case_token = m1.group(1)
        else:
            case_token = f"{n1}_vs_{n2}"
    else:
        raise SystemExit("ERROR: provide either <case*_scan*> OR <type1_outroot> <type2_outroot>")

    if not out1.is_dir():
        raise SystemExit(f"ERROR: folder not found: {out1}")
    if not out2.is_dir():
        raise SystemExit(f"ERROR: folder not found: {out2}")

    target = {
        "sin": _parse_single(args.sin, "sin"),
        "tan": _parse_single(args.tan, "tan"),
        "mA": _parse_single(args.mA, "mA"),
        "ma": _parse_single(args.ma, "ma"),
        "mchi": _parse_single(args.mchi, "mchi", default=0.0),
        "cosbma": _parse_single(args.cosbma, "cosbma", default=0.0),
        "deltam": _parse_single(args.deltam, "deltam", default=0.0),
    }

    ev1 = find_matching_events_dir(out1, target)
    ev2 = find_matching_events_dir(out2, target)
    if not ev1:
        raise SystemExit(f"ERROR: cannot find matching Events_* in {out1}/results for {target}")
    if not ev2:
        raise SystemExit(f"ERROR: cannot find matching Events_* in {out2}/results for {target}")

    f1 = find_unweighted_root(ev1)
    f2 = find_unweighted_root(ev2)
    if not f1:
        raise SystemExit(f"ERROR: unweighted_events.root not found under {ev1}")
    if not f2:
        raise SystemExit(f"ERROR: unweighted_events.root not found under {ev2}")

    root_exe = shutil.which("root")
    if not root_exe:
        raise SystemExit("ERROR: CERN ROOT executable 'root' not found on PATH")

    outdir = Path("comparison") / "met"
    outdir.mkdir(parents=True, exist_ok=True)
    out_png = outdir / f"{case_token}_{_output_param_tag(target)}.png"
    macro = outdir / ".comparison_met52_tmp.C"

    fixed_text = (
        f"m_{{A}}={target['mA']:g} GeV;"
        f"m_{{a}}={target['ma']:g} GeV;"
        f"m_{{#chi}}={target['mchi']:g} GeV;"
        f"sin#theta={target['sin']:g};"
        f"tan#beta={target['tan']:g};"
        f"cos(#beta-#alpha)={target['cosbma']:g};"
        f"#Delta m = m_{{H^{{+}}}} - m_{{A}} = {target['deltam']:g} GeV"
    )
    write_root_macro(macro)
    rc, out = run_root_compare(root_exe, macro, f1, f2, out_png, fixed_text, args.xmin, args.xmax)
    try:
        macro.unlink(missing_ok=True)
    except Exception:
        pass
    if rc != 0:
        log = outdir / (out_png.stem + ".log")
        log.write_text(out)
        raise SystemExit(f"ERROR: ROOT failed (see {log})")

    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
