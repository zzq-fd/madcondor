"""Generate HTCondor job lists for MadGraph scans.

This script supports three modes:

1) User top-config mode (default, no --preset/--from-results)
   Edit USER_SCAN_INPUT near the top of this file.
   Syntax examples:
     - Fixed value: 3
     - Independent scan: "scan:[1,2,3]"
     - Synchronized scan group: "scan1:[1,2,3]"
       (parameters with the same scan label change together by index)

2) Preset mode (--preset)
   Use a built-in grid such as type1case1scan1, tan0p01to50, legacy_default.

3) Sync-from-results mode (--from-results)
   Parse existing Events_* directory names and regenerate a matching joblist.

Output format per line (space-separated):
    SINTHETA TANBETA M36 M55 MCHI COSBMA DELTAM
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path
from typing import Iterable, Optional


# ---------------------------------------------------------------------------
# User-editable quick config (default mode)
# ---------------------------------------------------------------------------
#
# Edit only this dictionary for fast scan setup.
#
# Rules:
#   - Non-scan parameter: write a single value (int/float), e.g. 0.9
#   - Independent scan: "scan:[1,2,3]"
#   - Synchronized scan: "scan1:[...]", "scan2:[...]", ...
#       Parameters with the same scan label are zipped by index.
#
# Example:
#   "M36": "scan1:[600,700,800]"
#   "DELTAM": "scan1:[0,50,100]"
# produces aligned pairs:
#   (600,0), (700,50), (800,100)
#
USER_SCAN_INPUT: dict[str, object] = {
    "SINTHETA": 0.9,
    "TANBETA": 3,
    "M36": 1500,
    "M55": 300,
    "MCHI": "scan:[140,140.5,141,141.5,142,142.5,143,143.5,144,144.5,145,145.5,146,146.5,147,147.5,148,148.5,149,149.5,150]",
    "COSBMA": 0,
    "DELTAM": 0,
}


PARAM_ORDER = ["SINTHETA", "TANBETA", "M36", "M55", "MCHI", "COSBMA", "DELTAM"]
_SCAN_SPEC_RE = re.compile(r"^scan(?P<label>[A-Za-z0-9_]*)\s*:\s*(?P<values>\[.*\])\s*$")


def _decimal_range(start: float, stop: float, step: float) -> list[float]:
    """Generate an inclusive float range using decimal arithmetic."""
    from decimal import Decimal

    if step <= 0:
        raise ValueError("step must be > 0")

    start_d = Decimal(str(start))
    stop_d = Decimal(str(stop))
    step_d = Decimal(str(step))

    out: list[float] = []
    x = start_d
    while x <= stop_d:
        out.append(float(x))
        x += step_d
    return out


def _to_float(value: object, *, context: str) -> float:
    try:
        return float(value)
    except Exception as exc:
        raise ValueError(f"{context}: expected numeric value, got {valuerm -rf /home/zzq/madcondor/bbdm_2HDMa_type1_case1_scan_tan.tar.gz}") from exc


def _parse_scan_values(text: str, *, context: str) -> list[float]:
    try:
        parsed = ast.literal_eval(text)
    except Exception as exc:
        raise ValueError(f"{context}: cannot parse scan list: {text}") from exc
    if not isinstance(parsed, (list, tuple)):
        raise ValueError(f"{context}: scan part must be a list, got {type(parsed).__name__}")
    if len(parsed) == 0:
        raise ValueError(f"{context}: scan list cannot be empty")
    return [_to_float(v, context=context) for v in parsed]


def build_jobs_from_user_scan(user_cfg: dict[str, object]) -> list[tuple[float, float, float, float, float, float, float]]:
    """Build rows from USER_SCAN_INPUT.

    Behavior:
      - fixed value => one value
      - scan:[...] => independent scan dimension
      - scanN:[...] => synchronized scan group by label N
    """

    missing = [p for p in PARAM_ORDER if p not in user_cfg]
    if missing:
        raise ValueError(f"USER_SCAN_INPUT missing required keys: {', '.join(missing)}")

    extras = [k for k in user_cfg.keys() if k not in PARAM_ORDER]
    if extras:
        raise ValueError(f"USER_SCAN_INPUT has unknown keys: {', '.join(extras)}")

    fixed: dict[str, float] = {}
    independent_scans: list[tuple[str, list[float]]] = []
    sync_groups: dict[str, list[tuple[str, list[float]]]] = {}

    for param in PARAM_ORDER:
        raw = user_cfg[param]

        if isinstance(raw, (int, float)):
            fixed[param] = float(raw)
            continue

        if isinstance(raw, (list, tuple)):
            values = [_to_float(v, context=f"{param}") for v in raw]
            if not values:
                raise ValueError(f"{param}: list scan cannot be empty")
            independent_scans.append((param, values))
            continue

        if isinstance(raw, str):
            s = raw.strip()
            # Allow plain numeric strings as fixed values.
            try:
                fixed[param] = float(s)
                continue
            except Exception:
                pass

            m = _SCAN_SPEC_RE.fullmatch(s)
            if not m:
                raise ValueError(
                    f"{param}: invalid spec {raw!r}. Use number, scan:[...], or scanN:[...]"
                )

            label = m.group("label")
            values_text = m.group("values")
            values = _parse_scan_values(values_text, context=param)

            if label == "":
                # Plain 'scan:[...]' is independent per parameter.
                independent_scans.append((param, values))
            else:
                sync_groups.setdefault(label, []).append((param, values))
            continue

        raise ValueError(f"{param}: unsupported value type {type(raw).__name__}")

    combos: list[dict[str, float]] = [dict(fixed)]

    # Independent scan dimensions (Cartesian product)
    for param, values in independent_scans:
        new_combos: list[dict[str, float]] = []
        for combo in combos:
            for v in values:
                c = dict(combo)
                c[param] = float(v)
                new_combos.append(c)
        combos = new_combos

    # Synchronized scan groups (zip by index)
    for label, entries in sync_groups.items():
        lengths = {len(values) for _, values in entries}
        if len(lengths) != 1:
            detail = ", ".join(f"{p}:{len(vs)}" for p, vs in entries)
            raise ValueError(f"scan{label}: all synchronized lists must have equal length, got {detail}")

        n = next(iter(lengths))
        new_combos = []
        for combo in combos:
            for i in range(n):
                c = dict(combo)
                for param, values in entries:
                    c[param] = float(values[i])
                new_combos.append(c)
        combos = new_combos

    if not combos:
        raise ValueError("No scan points generated from USER_SCAN_INPUT")

    rows = [
        (
            float(c["SINTHETA"]),
            float(c["TANBETA"]),
            float(c["M36"]),
            float(c["M55"]),
            float(c["MCHI"]),
            float(c["COSBMA"]),
            float(c["DELTAM"]),
        )
        for c in combos
    ]

    # De-duplicate and keep stable ordering.
    rows_unique = sorted(set(rows), key=lambda r: (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
    return rows_unique


def _tag_to_float(tag: str) -> float:
    """Convert directory tag encoding to float.

    Conventions used in this repo:
      - '.' can be encoded as 'p' (e.g. 0p9 -> 0.9)
      - negative numbers can encode '-' as leading 'm' (e.g. m0p1 -> -0.1)
    """

    sign = 1.0
    s = tag
    if s.startswith("m"):
        sign = -1.0
        s = s[1:]
    s = s.replace("p", ".")
    return sign * float(s)


def _maybe_get(pattern: str, text: str) -> Optional[str]:
    m = re.search(pattern, text)
    if not m:
        return None
    return m.group(1)


def parse_events_dirname(dirname: str) -> dict:
    """Parse Events_* directory name into parameters."""

    name = dirname
    if name.startswith("Events_"):
        name = name[len("Events_") :]

    # Allow both 'p' and '.' decimal formats.
    st_tag = _maybe_get(r"(?:^|_)st(m?[0-9]+(?:[p.][0-9]+)?)", name)
    tb_tag = _maybe_get(r"(?:^|_)tb(m?[0-9]+(?:[p.][0-9]+)?)", name)
    mA_tag = _maybe_get(r"(?:^|_)mA(m?[0-9]+(?:[p.][0-9]+)?)", name)
    ma_tag = _maybe_get(r"(?:^|_)ma(m?[0-9]+(?:[p.][0-9]+)?)", name)

    # Optional fields
    mchi_tag = _maybe_get(r"(?:^|_)mchi(m?[0-9]+(?:[p.][0-9]+)?)", name)
    cbma_tag = _maybe_get(r"(?:^|_)cbma(m?[0-9]+(?:[p.][0-9]+)?)", name)
    dm_tag = _maybe_get(r"(?:^|_)dm(m?[0-9]+(?:[p.][0-9]+)?)", name)

    if not (st_tag and tb_tag and mA_tag and ma_tag):
        raise ValueError(f"Unrecognized Events directory name: {dirname}")

    out = {
        "sintheta": _tag_to_float(st_tag),
        "tanbeta": _tag_to_float(tb_tag),
        "mA": _tag_to_float(mA_tag),
        "ma": _tag_to_float(ma_tag),
        "mchi": _tag_to_float(mchi_tag) if mchi_tag else None,
        "cosbma": _tag_to_float(cbma_tag) if cbma_tag else None,
        "deltam": _tag_to_float(dm_tag) if dm_tag else None,
        "dirname": dirname,
    }
    return out


def collect_jobs_from_results(
    base: Path,
    default_mchi: float,
    default_cosbma: float,
    default_deltam: float,
) -> list[tuple[float, float, float, float, float, float, float]]:
    """Scan an existing output folder and return job rows."""

    base = base.expanduser().resolve()

    results_dir = base
    # Accept either <OUTROOT> or <OUTROOT>/results
    if (base / "results").is_dir():
        results_dir = base / "results"

    if not results_dir.is_dir():
        raise FileNotFoundError(f"results directory not found: {results_dir}")

    jobs: set[tuple[float, float, float, float, float, float, float]] = set()
    for child in results_dir.iterdir():
        if not child.is_dir():
            continue
        if not child.name.startswith("Events_"):
            continue
        params = parse_events_dirname(child.name)
        st = float(params["sintheta"])
        tb = float(params["tanbeta"])
        mA = float(params["mA"])
        ma = float(params["ma"])
        mchi = float(params["mchi"]) if params["mchi"] is not None else float(default_mchi)
        cosbma = float(params["cosbma"]) if params["cosbma"] is not None else float(default_cosbma)
        deltam = float(params["deltam"]) if params["deltam"] is not None else float(default_deltam)
        jobs.add((st, tb, mA, ma, mchi, cosbma, deltam))

    return sorted(jobs, key=lambda r: (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))


def build_preset_jobs(preset: str) -> list[tuple[float, float, float, float, float, float, float]]:
    """Return unique rows for a built-in preset."""

    if preset == "type1case1scan1":
        mchi_values = _decimal_range(140.0, 150.0, 0.5)
        return [(0.9, 3.0, 1500.0, 300.0, mchi, 0.0, 0.0) for mchi in mchi_values]

    if preset == "tan0p01to50":
        tanbeta_values = [
            0.01,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0,
            2.0,
            5.0,
            10.0,
            20.0,
            30.0,
            40.0,
            50.0,
        ]
        return [(0.9, float(tb), 600.0, 300.0, 1.0, 0.0, 0.0) for tb in tanbeta_values]

    if preset == "legacy_default":
        mA_fixed = 1500.0
        mchi_fixed = 10.0
        cosbma_fixed = 0.0
        deltam_fixed = 0.0

        sin_scan = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        tan_fixed = 3.0
        ma_fixed = 300.0

        sin_fixed = 0.9
        ma_scan = [
            50,
            100,
            150,
            200,
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
        ]
        tan_scan = [3, 5, 10, 20, 30, 40, 50]

        jobs_set: set[tuple[float, float, float, float, float, float, float]] = set()
        for st in sin_scan:
            jobs_set.add(
                (
                    float(st),
                    float(tan_fixed),
                    float(mA_fixed),
                    float(ma_fixed),
                    float(mchi_fixed),
                    float(cosbma_fixed),
                    float(deltam_fixed),
                )
            )
        for ma in ma_scan:
            jobs_set.add(
                (
                    float(sin_fixed),
                    float(tan_fixed),
                    float(mA_fixed),
                    float(ma),
                    float(mchi_fixed),
                    float(cosbma_fixed),
                    float(deltam_fixed),
                )
            )
        for tb in tan_scan:
            jobs_set.add(
                (
                    float(sin_fixed),
                    float(tb),
                    float(mA_fixed),
                    float(ma_fixed),
                    float(mchi_fixed),
                    float(cosbma_fixed),
                    float(deltam_fixed),
                )
            )
        return sorted(jobs_set, key=lambda r: (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))

    raise ValueError(f"Unknown preset: {preset}")


def write_joblist(path: Path, jobs: Iterable[tuple[float, float, float, float, float, float, float]]) -> int:
    path = path.expanduser()
    count = 0
    with open(path, "w") as f:
        for row in jobs:
            f.write("{} {} {} {} {} {} {}\n".format(*[f"{x:g}" for x in row]))
            count += 1
    return count


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Generate HTCondor joblist for MadGraph scans")
    parser.add_argument(
        "--preset",
        type=str,
        default=None,
        choices=["type1case1scan1", "tan0p01to50", "legacy_default"],
        help=(
            "Use a built-in grid. "
            "type1case1scan1: mchi=140..150 step 0.5; "
            "tan0p01to50: tan list 0.01..50 at fixed masses; "
            "legacy_default: previous hard-coded union grid."
        ),
    )
    parser.add_argument(
        "--from-results",
        type=str,
        default=None,
        help=(
            "Parse existing Events_* directories to regenerate joblist. "
            "Accepts either <OUTROOT> or <OUTROOT>/results."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="joblist.txt",
        help="Output joblist file (default: joblist.txt)",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Repeat each parameter point N times (default: 1).",
    )
    parser.add_argument("--default-mchi", type=float, default=10.0, help="Default mchi if missing in dirname")
    parser.add_argument("--default-cosbma", type=float, default=0.0, help="Default cosbma if missing in dirname")
    parser.add_argument("--default-deltam", type=float, default=0.0, help="Default deltam if missing in dirname")

    args = parser.parse_args(argv)

    if args.repeat < 1:
        raise SystemExit("ERROR: --repeat must be >= 1")

    def _apply_repeat(rows: list[tuple[float, float, float, float, float, float, float]]) -> list[
        tuple[float, float, float, float, float, float, float]
    ]:
        if args.repeat == 1:
            return rows
        repeated: list[tuple[float, float, float, float, float, float, float]] = []
        for _ in range(args.repeat):
            repeated.extend(rows)
        return repeated

    if args.from_results:
        base = Path(args.from_results)
        jobs_unique = collect_jobs_from_results(
            base=base,
            default_mchi=args.default_mchi,
            default_cosbma=args.default_cosbma,
            default_deltam=args.default_deltam,
        )
        unique_n = len(jobs_unique)
        jobs = _apply_repeat(jobs_unique)
        out_path = Path(args.output)
        n = write_joblist(out_path, jobs)
        if args.repeat == 1:
            print(f"Found {unique_n} points from {base}, wrote {out_path}")
        else:
            print(f"Found {unique_n} points from {base}, wrote {out_path} (repeat={args.repeat}, total_lines={n})")
        return 0

    if args.preset:
        jobs_unique = build_preset_jobs(args.preset)
        unique_n = len(jobs_unique)
        jobs = _apply_repeat(jobs_unique)
        out_path = Path(args.output)
        n = write_joblist(out_path, jobs)
        if args.repeat == 1:
            print(f"Wrote {unique_n} lines to {out_path} (preset={args.preset})")
        else:
            print(f"Wrote {unique_n}x{args.repeat}={n} lines to {out_path} (preset={args.preset})")
        return 0

    # Default mode: use the top editable USER_SCAN_INPUT block.
    try:
        jobs_unique = build_jobs_from_user_scan(USER_SCAN_INPUT)
    except Exception as exc:
        raise SystemExit(f"ERROR in USER_SCAN_INPUT: {exc}") from exc

    unique_n = len(jobs_unique)
    jobs = _apply_repeat(jobs_unique)
    out_path = Path(args.output)
    n = write_joblist(out_path, jobs)

    if args.repeat == 1:
        print(f"Wrote {unique_n} points to {out_path} (mode=user_top_config)")
    else:
        print(f"Wrote {unique_n} points to {out_path} (repeat={args.repeat}, total_lines={n}, mode=user_top_config)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
