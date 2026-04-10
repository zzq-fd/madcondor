#!/usr/bin/env python3
"""Plot pT of particle 52 from MadGraph ROOT outputs and summarize mean pT.

Usage:
  python3 plot_pt52.py <OUTROOT>

Example:
  python3 plot_pt52.py bbdm_2HDMa_type1_case1_scan

What it does:
  - Scans <OUTROOT>/results/Events_*/*/unweighted_events.root
  - For each Events_* directory, makes a pT histogram for particle 52 and saves
    a PNG into that Events_* directory.
  - Computes mean pT for particle 52 and writes a summary table into <OUTROOT>/.

Notes:
  - Requires CERN ROOT executable `root` available on PATH.
  - Does not require PyROOT.
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
class PointInfo:
    sintheta: float | None
    tanbeta: float | None
    mA: float | None
    ma: float | None
    mchi: float | None
    cosbma: float | None
    deltam: float | None


def _parse_tag_float(tag: str) -> float | None:
    """Parse a tag-encoded float.

        Supported forms:
            - 0p1 / 0.1 -> 0.1
            - 1         -> 1.0
            - m0p1 / m0.1 -> -0.1
    """

    if not tag:
        return None
    neg = tag.startswith("m")
    if neg:
        tag = tag[1:]
    try:
        val = float(tag.replace("p", "."))
    except Exception:
        return None
    return -val if neg else val


def parse_events_dir_name(dirname: str) -> PointInfo:
    name = dirname.replace("Events_", "")

    st_match = re.search(r"(?:^|_)st(m?[0-9]+(?:[p.][0-9]+)?)", name)
    sintheta = _parse_tag_float(st_match.group(1)) if st_match else None

    tb_match = re.search(r"(?:^|_)tb(m?[0-9]+(?:[p.][0-9]+)?)", name)
    tanbeta = _parse_tag_float(tb_match.group(1)) if tb_match else None

    mA_match = re.search(r"mA([0-9]+)", name)
    mA = float(mA_match.group(1)) if mA_match else None

    ma_match = re.search(r"ma([0-9]+)", name)
    ma = float(ma_match.group(1)) if ma_match else None

    mchi: float | None = None
    mchi_match = re.search(r"mchi(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if mchi_match:
        mchi = _parse_tag_float(mchi_match.group(1))

    cosbma: float | None = None
    cbma_match = re.search(r"cbma(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if cbma_match:
        cosbma = _parse_tag_float(cbma_match.group(1))

    deltam: float | None = None
    dm_match = re.search(r"dm(m?[0-9]+(?:[p.][0-9]+)?)", name)
    if dm_match:
        deltam = _parse_tag_float(dm_match.group(1))

    return PointInfo(sintheta=sintheta, tanbeta=tanbeta, mA=mA, ma=ma, mchi=mchi, cosbma=cosbma, deltam=deltam)


def write_root_macro(macro_path: Path) -> None:
    macro_path.write_text(
        r'''#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

int plot_pt52(const char* input_root,
              const char* output_png,
              const char* output_txt,
              double xmin_in,
              double xmax_in)
{
  gROOT->SetBatch(kTRUE);

  TFile f(input_root, "READ");
  if (f.IsZombie()) {
    std::cerr << "ERROR: cannot open ROOT file: " << input_root << std::endl;
    return 2;
  }

  auto* t = dynamic_cast<TTree*>(f.Get("LHEF"));
  if (!t) {
    std::cerr << "ERROR: TTree 'LHEF' not found in: " << input_root << std::endl;
    return 3;
  }

  // Collect pT values for particle 52.
  // This uses the flat array branches written by MadGraph: Particle.PID and Particle.PT.
  const char* expr = "Particle.PT";
  const char* sel  = "abs(Particle.PID)==52";

  Long64_t n = t->Draw(expr, sel, "goff");
  if (n <= 0) {
    std::ofstream out(output_txt);
    out << "entries 0\n";
    out << "mean_pt nan\n";
    out.close();

    // Still write an empty plot (useful for bookkeeping)
    TCanvas c("c", "c", 1100, 800);
    c.SetLeftMargin(0.14);
    c.SetBottomMargin(0.12);
        double xlow = xmin_in;
        double xhi = (xmax_in > xlow) ? xmax_in : (xlow + 1.0);
        TH1D h("h", ";p_{T}(particle 52) [GeV];Entries", 10, xlow, xhi);
    h.Draw();
    c.SaveAs(output_png);
    return 0;
  }

  const double* v = t->GetV1();
  double sum = 0.0;
  double vmax = 0.0;
  for (Long64_t i = 0; i < n; ++i) {
    const double pt = v[i];
    sum += pt;
    if (pt > vmax) vmax = pt;
  }

    double xlow = xmin_in;
    double xhi = (xmax_in > xlow) ? xmax_in : (vmax * 1.05);
    if (!(xhi > xlow)) xhi = xlow + 1.0;

    TH1D h("h", ";p_{T}(particle 52) [GeV];Entries", 100, xlow, xhi);
    double sum_sel = 0.0;
    Long64_t n_sel = 0;
  for (Long64_t i = 0; i < n; ++i) {
        const double pt = v[i];
        if (pt < xlow || pt > xhi) continue;
        h.Fill(pt);
        sum_sel += pt;
        n_sel += 1;
  }

  TCanvas c("c", "c", 1100, 800);
  c.SetLeftMargin(0.14);
  c.SetBottomMargin(0.12);
  h.GetYaxis()->SetTitleOffset(1.25);
  h.Draw("hist");
  c.SaveAs(output_png);

  std::ofstream out(output_txt);
    if (n_sel > 0) {
        out << "entries " << n_sel << "\n";
        out << "mean_pt " << (sum_sel / double(n_sel)) << "\n";
    } else {
        out << "entries 0\n";
        out << "mean_pt nan\n";
    }
  out.close();

  return 0;
}
'''
    )


def root_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def find_root_file(events_dir: Path) -> Path | None:
    run_dirs = sorted([p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])
    if not run_dirs:
        return None

    preferred = next((p for p in run_dirs if p.name == "run_01"), None)
    chosen = preferred or run_dirs[0]

    root_path = chosen / "unweighted_events.root"
    return root_path if root_path.exists() else None


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True, allow_abbrev=False)
    ap.add_argument("outroot", help="Top-level output directory (e.g. bbdm_2HDMa_type1_case1_scan)")
    ap.add_argument("--xmin", type=float, default=0.0, help="x-axis lower bound (default: 0)")
    ap.add_argument("--xmax", type=float, default=0.0, help="x-axis upper bound (0 means auto)")
    args = ap.parse_args()
    if args.xmax > 0 and args.xmax <= args.xmin:
        raise SystemExit("ERROR: --xmax must be > --xmin")

    outroot = Path(args.outroot).expanduser()
    results_dir = outroot / "results"
    if not results_dir.exists():
        print(f"ERROR: results directory not found: {results_dir}")
        return 2

    root_exe = shutil.which("root")
    if not root_exe:
        print("ERROR: CERN ROOT executable 'root' not found on PATH.")
        print("Hint: if you use conda, run: conda activate mg5")
        return 2

    macro_path = outroot / ".plot_pt52_tmp.C"
    write_root_macro(macro_path)

    event_dirs = sorted([d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith("Events_")])
    print(f"Found {len(event_dirs)} event directories")

    rows: list[tuple[PointInfo, str, float | None, int | None]] = []
    any_failures = False

    for ev_dir in event_dirs:
        root_file = find_root_file(ev_dir)
        if not root_file:
            print(f"Warning: missing run_*/unweighted_events.root under {ev_dir.name}")
            info = parse_events_dir_name(ev_dir.name)
            rows.append((info, ev_dir.name, None, None))
            continue

        out_png = ev_dir / "pt_particle52.png"
        out_txt = ev_dir / "pt_particle52_stats.txt"

        stmt = (
            f'gROOT->LoadMacro("{root_quote(macro_path)}+"); '
            f'plot_pt52("{root_quote(root_file)}","{root_quote(out_png)}","{root_quote(out_txt)}",{args.xmin},{args.xmax});'
        )
        cmd = [root_exe, "-l", "-b", "-q", "-e", stmt]

        log_path = ev_dir / "pt_particle52_root.log"
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log_path.write_text(proc.stdout)

        info = parse_events_dir_name(ev_dir.name)

        mean_pt: float | None = None
        entries: int | None = None
        if out_txt.exists():
            txt = out_txt.read_text().splitlines()
            for line in txt:
                if line.startswith("entries "):
                    try:
                        entries = int(line.split()[1])
                    except Exception:
                        entries = None
                if line.startswith("mean_pt "):
                    try:
                        val = float(line.split()[1])
                        mean_pt = None if (val != val) else val
                    except Exception:
                        mean_pt = None

        if proc.returncode != 0:
            any_failures = True
            print(f"Warning: ROOT failed for {ev_dir.name} (see {log_path})")
        else:
            try:
                log_path.unlink(missing_ok=True)
            except Exception:
                pass

        rows.append((info, ev_dir.name, mean_pt, entries))
        if mean_pt is not None and entries is not None:
            print(f"Processed: {ev_dir.name} -> mean pT={mean_pt:.6g} (entries={entries})")
        else:
            print(f"Processed: {ev_dir.name} -> mean pT=N/A")

    table_path = outroot / "pt52_mean_table.txt"
    with open(table_path, "w") as f:
        f.write("# Mean pT table for particle 52\n")
        f.write("# Columns: sintheta tanbeta mA ma mchi cosbma deltam mean_pt entries\n")
        f.write(
            f"{'sintheta':<10} {'tanbeta':<8} {'mA':<6} {'ma':<6} {'mchi':<10} {'cosbma':<10} {'deltam':<10} {'mean_pt':<14} {'entries':<8} {'dirname'}\n"
        )
        f.write("-" * 130 + "\n")
        for info, dirname, mean_pt, entries in rows:
            st = f"{info.sintheta:.4f}" if info.sintheta is not None else "N/A"
            tb = f"{info.tanbeta:.1f}" if info.tanbeta is not None else "N/A"
            mA = f"{info.mA:.0f}" if info.mA is not None else "N/A"
            ma = f"{info.ma:.0f}" if info.ma is not None else "N/A"
            mchi = f"{info.mchi:g}" if info.mchi is not None else "N/A"
            cosbma = f"{info.cosbma:g}" if info.cosbma is not None else "N/A"
            deltam = f"{info.deltam:g}" if info.deltam is not None else "N/A"
            mp = f"{mean_pt:.8g}" if mean_pt is not None else "N/A"
            en = str(entries) if entries is not None else "N/A"
            f.write(f"{st:<10} {tb:<8} {mA:<6} {ma:<6} {mchi:<10} {cosbma:<10} {deltam:<10} {mp:<14} {en:<8} {dirname}\n")

    print(f"\nTable written to: {table_path}")

    if not any_failures:
        try:
            macro_path.unlink(missing_ok=True)
        except Exception:
            pass

        compiled_base = macro_path.name.replace(".C", "_C")
        for p in outroot.glob(f"{compiled_base}*"):
            try:
                if p.is_file():
                    p.unlink(missing_ok=True)
            except Exception:
                pass

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
