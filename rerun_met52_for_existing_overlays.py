#!/usr/bin/env python3
"""Re-run plot_met52 for existing overlay_met55 plots.

This script scans directories matching:
  bbdm_2HDMa_type1_case*_scan*

For each OUTROOT, it looks for:
  <OUTROOT>/output_scatter_plot/overlay_met/overlay_met55_*.png

It parses the filename to extract the varying tag and fixed-tag parameters.
Then it reconstructs the parameter value lists by inspecting the available
Events_* directories in <OUTROOT>/results, and finally calls:
  python3 plot_met52.py <OUTROOT> ...

Notes:
  - For single-parameter scans, it uses *all* available values in results that
    match the fixed parameters.
  - For paired scans, it uses all available (a,b) pairs matching the fixed
    parameters, sorted by (a,b), and passes paired lists.
  - If the expected output PNG already exists in output_scatter_plot/met52,
    the plot is skipped.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


def _tag_to_float(tag: str) -> float:
    neg = tag.startswith("m")
    if neg:
        tag = tag[1:]
    v = float(tag.replace("p", "."))
    return -v if neg else v


def _st_tag(x: float) -> str:
    return f"{x:g}".replace("-", "m").replace(".", "p")


def parse_events_dir_params(dirname: str) -> dict[str, float]:
    name = dirname.replace("Events_", "")
    out: dict[str, float] = {}

    st_match = re.search(r"st([0-9]+)p([0-9]+)", name)
    if st_match:
        out["sin"] = float(f"{st_match.group(1)}.{st_match.group(2)}")

    tb_match = re.search(r"tb([0-9]+)", name)
    if tb_match:
        out["tan"] = float(tb_match.group(1))

    mA_match = re.search(r"mA([0-9]+)", name)
    if mA_match:
        out["mA"] = float(mA_match.group(1))

    ma_match = re.search(r"ma([0-9]+)", name)
    if ma_match:
        out["ma"] = float(ma_match.group(1))

    mchi_match = re.search(r"mchi(m?[0-9]+(?:p[0-9]+)?)", name)
    if mchi_match:
        out["mchi"] = _tag_to_float(mchi_match.group(1))

    cbma_match = re.search(r"cbma(m?[0-9]+(?:p[0-9]+)?)", name)
    if cbma_match:
        out["cosbma"] = _tag_to_float(cbma_match.group(1))

    dm_match = re.search(r"dm(m?[0-9]+(?:p[0-9]+)?)", name)
    if dm_match:
        out["deltam"] = _tag_to_float(dm_match.group(1))

    return out


def parse_overlay_png_name(png: Path) -> tuple[str, dict[str, float | None]]:
    """Return (vary_tag, fixed_map) parsed from overlay_met55_<vary_tag>__<fixed_tag>.png."""
    m = re.match(r"^overlay_met55_(.+)__([^.]+)\.png$", png.name)
    if not m:
        raise ValueError(f"Unrecognized filename: {png.name}")
    vary_tag = m.group(1)
    fixed_tag = m.group(2)

    # fixed_tag example: st0p9_tb3_mA1500_mavar_mchi10_cbma0_dm0
    parts = fixed_tag.split("_")
    fixed: dict[str, float | None] = {
        "sin": None,
        "tan": None,
        "mA": None,
        "ma": None,
        "mchi": None,
        "cosbma": None,
        "deltam": None,
    }

    for part in parts:
        if part.startswith("st"):
            fixed["sin"] = _tag_to_float(part[2:]) if part[2:] != "var" else None
        elif part.startswith("tb"):
            fixed["tan"] = _tag_to_float(part[2:]) if part[2:] != "var" else None
        elif part.startswith("mA"):
            fixed["mA"] = float(part[2:]) if part[2:] != "var" else None
        elif part.startswith("ma"):
            fixed["ma"] = float(part[2:]) if part[2:] != "var" else None
        elif part.startswith("mchi"):
            fixed["mchi"] = _tag_to_float(part[4:]) if part[4:] != "var" else None
        elif part.startswith("cbma"):
            fixed["cosbma"] = _tag_to_float(part[4:]) if part[4:] != "var" else None
        elif part.startswith("dm"):
            fixed["deltam"] = _tag_to_float(part[2:]) if part[2:] != "var" else None

    return vary_tag, fixed


def _fmt_arg_list(vals: list[float], kind: str) -> str:
    if kind == "int":
        return ",".join(str(int(v)) for v in vals)
    return ",".join(f"{v:g}" for v in vals)


def _expected_met52_png(outroot: Path, vary_tag: str, fixed: dict[str, float | None]) -> Path:
    def _tag_or_var(x: float | None, *, kind: str) -> str:
        if x is None:
            return "var"
        if kind == "int":
            return str(int(x))
        return _st_tag(x)

    fixed_tag = (
        f"st{_tag_or_var(fixed['sin'], kind='float')}"
        f"_tb{_tag_or_var(fixed['tan'], kind='float')}"
        f"_mA{_tag_or_var(fixed['mA'], kind='int')}"
        f"_ma{_tag_or_var(fixed['ma'], kind='int')}"
        f"_mchi{_tag_or_var(fixed['mchi'], kind='float')}"
        f"_cbma{_tag_or_var(fixed['cosbma'], kind='float')}"
        f"_dm{_tag_or_var(fixed['deltam'], kind='float')}"
    )
    out_dir = outroot / "output_scatter_plot" / "met52"
    return out_dir / f"met52_{vary_tag}__{fixed_tag}.png"


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--root", default=".", help="Base directory to search (default: current dir)")
    ap.add_argument("--pattern", default="bbdm_2HDMa_type1_case*_scan*", help="Glob for OUTROOT dirs")
    ap.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    args = ap.parse_args()

    base = Path(args.root).expanduser().resolve()
    outroots = sorted([p for p in base.glob(args.pattern) if p.is_dir()])
    if not outroots:
        print(f"No OUTROOT dirs found under {base} matching {args.pattern}")
        return 2

    script = base / "plot_met52.py"
    if not script.exists():
        print(f"ERROR: plot_met52.py not found at {script}")
        return 2

    total = 0
    ran = 0
    skipped = 0

    for outroot in outroots:
        overlay_dir = outroot / "output_scatter_plot" / "overlay_met"
        if not overlay_dir.exists():
            continue

        pngs = sorted(overlay_dir.glob("overlay_met55_*.png"))
        if not pngs:
            continue

        results_dir = outroot / "results"
        event_dirs = sorted([p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith("Events_")]) if results_dir.exists() else []
        recognized: set[str] = set()
        for d in event_dirs[:50]:
            recognized |= set(parse_events_dir_params(d.name).keys())

        for png in pngs:
            total += 1
            try:
                vary_tag, fixed = parse_overlay_png_name(png)
            except Exception as e:
                print(f"Skip (parse): {png} ({e})")
                continue

            # Determine varying params from vary_tag
            if vary_tag == "single":
                varying = []
            elif vary_tag.endswith("_paired"):
                core = vary_tag[: -len("_paired")]
                # e.g. mA_ma_paired
                varying = core.split("_")
            else:
                varying = [vary_tag]

            # Legacy overlay_met55 filenames from older script versions may not use 'var' in fixed_tag.
            # If a parameter is varying, ignore any fixed value parsed from the filename.
            for vp in varying:
                if vp in fixed:
                    fixed[vp] = None

            # Reconstruct value lists from results
            # Filter Events_* by fixed values (only for recognized params)
            def _match_fixed(params: dict[str, float]) -> bool:
                for key, want in fixed.items():
                    if want is None:
                        continue
                    if key not in recognized:
                        continue
                    if key not in params:
                        return False
                    if abs(params[key] - want) > 1e-9:
                        return False
                return True

            matched_params: list[dict[str, float]] = []
            for ev in event_dirs:
                params = parse_events_dir_params(ev.name)
                if _match_fixed(params):
                    matched_params.append(params)

            if not matched_params:
                print(f"Skip (no matching Events_*): {png}")
                continue

            # Build argument lists
            value_lists: dict[str, list[float]] = {}
            all_params = ("sin", "tan", "mA", "ma", "mchi", "cosbma", "deltam")
            for name in all_params:
                if name in varying:
                    continue
                if fixed.get(name) is not None:
                    value_lists[name] = [fixed[name]]
                    continue

                # Older overlay_met55 filenames might not encode some parameters.
                # Infer from the matching Events_* directories; require a unique value.
                uniq = sorted({p[name] for p in matched_params if name in p})
                if len(uniq) == 1:
                    value_lists[name] = [uniq[0]]
                    continue
                if len(uniq) == 0:
                    # If a parameter is not encoded in Events_* names, it also won't be inferable here.
                    # For legacy scans, cosbma and deltam were typically fixed to 0.
                    if name in ("cosbma", "deltam"):
                        value_lists[name] = [0.0]
                        continue
                    raise SystemExit(f"ERROR: cannot infer required fixed parameter '{name}' from results for {png}")
                raise SystemExit(
                    f"ERROR: parameter '{name}' is not encoded in filename and is not unique in results (values={uniq}) for {png}"
                )

            if len(varying) == 0:
                pass
            elif len(varying) == 1:
                v = varying[0]
                vals = sorted({p[v] for p in matched_params if v in p})
                value_lists[v] = vals
            else:
                a, b = varying[0], varying[1]
                pairs = sorted({(p.get(a), p.get(b)) for p in matched_params if a in p and b in p})
                a_vals = [pa for pa, _ in pairs]
                b_vals = [pb for _, pb in pairs]
                value_lists[a] = a_vals
                value_lists[b] = b_vals

            # Sanity: ensure required keys exist
            for req in all_params:
                if req not in value_lists:
                    raise SystemExit(f"ERROR: internal: missing argument list for '{req}' for {png}")

            # Expected output
            out_png = _expected_met52_png(outroot, vary_tag, fixed)
            if out_png.exists():
                skipped += 1
                continue

            cmd = [
                sys.executable,
                str(script),
                str(outroot),
                "--sin",
                _fmt_arg_list(value_lists["sin"], "float"),
                "--tan",
                _fmt_arg_list(value_lists["tan"], "float"),
                "--mA",
                _fmt_arg_list(value_lists["mA"], "int"),
                "--ma",
                _fmt_arg_list(value_lists["ma"], "int"),
                "--mchi",
                _fmt_arg_list(value_lists["mchi"], "float"),
                "--cosbma",
                _fmt_arg_list(value_lists["cosbma"], "float"),
                "--deltam",
                _fmt_arg_list(value_lists["deltam"], "float"),
            ]

            print("RUN:", " ".join(cmd))
            if not args.dry_run:
                out_png.parent.mkdir(parents=True, exist_ok=True)
                proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                if proc.returncode != 0:
                    print(proc.stdout)
                    print(f"FAILED (rc={proc.returncode}): {png}")
                    continue
            ran += 1

    print(f"Done. found={total}, ran={ran}, skipped_existing={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
