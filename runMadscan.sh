#!/usr/bin/env bash

# Exit immediately on error and echo every command
set -e
set -x

ulimit -s unlimited

# HTCondor execute nodes sometimes provide a very minimal environment.
# MadGraph expects PATH (and often HOME) to exist.
export PATH="${PATH:-/usr/local/bin:/usr/bin:/bin}"
# Ensure /bin is present even if PATH was set but incomplete
case ":${PATH}:" in
  *:/bin:*) ;;
  *) export PATH="${PATH}:/bin";;
esac
export HOME="${HOME:-$PWD}"

# -----------------------------
# Process package autodetect
# -----------------------------
#############################
# Process package autodetect
#############################
# Respect PROCESS_TAR/OUTROOT passed by HTCondor via environment.
# This prevents mismatches where Condor expects $(OUTROOT)/results but the script
# writes to a different folder derived from a different tarball name.

# Default tarball name if none specified.
PROCESS_TAR=${PROCESS_TAR:-bbdm_2HDMa_type1_case1_scan.tar.gz}

# Output root folder name. Prefer env OUTROOT (set by submit file); otherwise derive.
OUTROOT=${OUTROOT:-${PROCESS_TAR%.tar.gz}}

# Always create the output dir early so output transfer doesn't fail even if we exit early.
mkdir -p "${OUTROOT}/results"

# If the requested tarball isn't present, try falling back to any matching tarball.
if [[ ! -f "${PROCESS_TAR}" ]]; then
  # Also allow tarballs like ..._scan1.tar.gz
  alt=$(ls -1 bbdm_2HDMa_type1_case1_scan*.tar.gz 2>/dev/null | head -n 1 || true)
  if [[ -n "${alt}" ]]; then
    PROCESS_TAR="${alt}"
  fi
fi

if [[ -z "${PROCESS_TAR}" || ! -f "${PROCESS_TAR}" ]]; then
  echo "ERROR: Could not find process tarball. Expected '${PROCESS_TAR}' or a file matching 'bbdm_2HDMa_type1_case1_scan*.tar.gz' in $PWD" >&2
  exit 2
fi

# Detect the top-level directory name inside the tarball
# (redirect stderr to avoid tar SIGPIPE warnings when head exits early)
PROCESS_DIR=$(tar -tf "${PROCESS_TAR}" 2>/dev/null | head -n 1 | cut -d/ -f1)
if [[ -z "${PROCESS_DIR}" ]]; then
  echo "ERROR: Could not detect top-level directory inside ${PROCESS_TAR}" >&2
  exit 2
fi

# -----------------------------
# Input arguments
# -----------------------------
JOBID=$1
SINTHETA=$2
TANBETA=$3
M36_IN=$4
M55=$5
MCHI=${6:-1}
COSBMA=${7:-0}
DELTAM=${8:-0}

# -----------------------------
# Random seed (per job)
# -----------------------------
# For repeated submissions, ensure each job uses a distinct random seed so
# event generation and integration are statistically independent.
# You can override via env var ISEED.
ISEED=${ISEED:-$((1000 + JOBID))}

# Use a python helper to handle floating point math robustly.
PYTHON_BIN=$(command -v python3 || command -v python || true)
if [[ -z "${PYTHON_BIN}" ]]; then
  echo "ERROR: python3/python not found on execute node; required for float computations (sinbma, m36)." >&2
  exit 3
fi

read -r M35 M36 M37 SINBMA < <(M36_IN="${M36_IN}" COSBMA="${COSBMA}" DELTAM="${DELTAM}" "${PYTHON_BIN}" - <<'PY'
import math
import os

m36_in = float(os.environ["M36_IN"])
cosbma = float(os.environ["COSBMA"])
deltam = float(os.environ["DELTAM"])

# DELTAM is defined as m35 - m36
m36 = m36_in
m35 = m36 + deltam
m37 = m35

# User provides cos(beta-alpha); param card expects sin(beta-alpha).
# Take the positive branch and clamp for numerical safety.
arg = 1.0 - cosbma * cosbma
if arg < 0 and arg > -1e-12:
    arg = 0.0
if arg < 0:
    raise SystemExit(f"cosbma={cosbma} is outside [-1,1]; cannot compute sinbma")
sinbma = math.sqrt(arg)

print(f"{m35:g} {m36:g} {m37:g} {sinbma:.12g}")
PY
)
SINTHETA_TAG=${SINTHETA//./p}
MCHI_TAG=${MCHI//./p}
COSBMA_TAG=${COSBMA//-/m}
COSBMA_TAG=${COSBMA_TAG//./p}
DELTAM_TAG=${DELTAM//-/m}
DELTAM_TAG=${DELTAM_TAG//./p}

# -----------------------------
# Job info
# -----------------------------
echo "===================================="
echo "JOBID        = ${JOBID}"
echo "iseed        = ${ISEED}"
echo "sinθ         = ${SINTHETA}"
echo "tanβ         = ${TANBETA}"
echo "cos(b-a)     = ${COSBMA}"
echo "Δm=m35-m36   = ${DELTAM}"
echo "sin(b-a)     = ${SINBMA}"
echo "mA (=m36)    = ${M36_IN}"
echo "m35          = ${M35}"
echo "m36          = ${M36}"
echo "m37          = ${M37}"
echo "m55          = ${M55}"
echo "mchi         = ${MCHI}"
echo "===================================="

# -----------------------------
# Untar production directory
# -----------------------------
echo "# Untar the Production Dir"
tar -xaf "${PROCESS_TAR}"
echo "Done"

# -----------------------------
# Prepare launch card
# -----------------------------
LAUNCH_CARD="${PROCESS_DIR}/launch_card_st${SINTHETA_TAG}_tb${TANBETA}_mA${M36_IN}_ma${M55}_mchi${MCHI_TAG}_cbma${COSBMA_TAG}_dm${DELTAM_TAG}.dat"

cp launch_card_template.dat "${LAUNCH_CARD}"

sed -i \
  -e "s/__TANBETA__/${TANBETA}/g" \
  -e "s/__SINTHETA__/${SINTHETA}/g" \
  -e "s/__SINBMA__/${SINBMA}/g" \
  -e "s/__M35__/${M35}/g" \
  -e "s/__M36__/${M36}/g" \
  -e "s/__M37__/${M37}/g" \
  -e "s/__M55__/${M55}/g" \
  -e "s/__MCHI__/${MCHI}/g" \
  "${LAUNCH_CARD}"

# Enforce per-job random seed (important for repeat jobs).
{
  echo ""
  echo "###################################"
  echo "## Random seed (per job)"
  echo "###################################"
  echo "set run_card iseed ${ISEED}"
} >> "${LAUNCH_CARD}"

# -----------------------------
# Ensure auto decays are set
# -----------------------------
# Some template variants may omit these; enforce them here to be safe.
need_decay_fix=0
for pid in 25 35 36 37 55; do
  if ! grep -qE "^set param_card[[:space:]]+DECAY[[:space:]]+${pid}[[:space:]]+auto\b" "${LAUNCH_CARD}"; then
    need_decay_fix=1
    break
  fi
done

if [[ ${need_decay_fix} -eq 1 ]]; then
  {
    echo ""
    echo "###################################"
    echo "## Decays (enforced)"
    echo "###################################"
    echo "set param_card DECAY 25 auto"
    echo "set param_card DECAY 35 auto"
    echo "set param_card DECAY 36 auto"
    echo "set param_card DECAY 37 auto"
    echo "set param_card DECAY 55 auto"
  } >> "${LAUNCH_CARD}"
fi

# -----------------------------
# Run MadEvent
# -----------------------------
cd "${PROCESS_DIR}"

# Use an explicit Python interpreter for MadGraph if available.
# This avoids relying on the execute node's default python (often too old).
MG5_PYTHON_DEFAULT="${HOME}/miniconda3/envs/mg5/bin/python3"
if [[ ! -x "${MG5_PYTHON_DEFAULT}" ]]; then
  MG5_PYTHON_DEFAULT=""
fi

# Allow explicit override via environment variable MG5_PYTHON.
# Otherwise try a common conda env location under $HOME, and finally fall back
# to whatever python was detected earlier for float computations.
MG5_PYTHON=${MG5_PYTHON:-${MG5_PYTHON_DEFAULT}}
if [[ -z "${MG5_PYTHON}" ]]; then
  MG5_PYTHON="${PYTHON_BIN}"
fi

# If we're using a conda env python, prepend its bin/ to PATH so tools like
# gfortran installed in that env can be found.
if [[ -x "${MG5_PYTHON}" ]]; then
  MG5_ENV_BIN=$(dirname "${MG5_PYTHON}")
  export PATH="${MG5_ENV_BIN}:${PATH}"
fi

# MadEvent needs a Fortran compiler to build the process code.
echo "[env] PATH=${PATH}"
echo "[env] HOME=${HOME}"
echo "[env] MG5_PYTHON=${MG5_PYTHON}"
command -v gfortran >/dev/null 2>&1 && echo "[env] gfortran=$(command -v gfortran)" || echo "[env] gfortran=NOT_FOUND"
command -v gfortran >/dev/null 2>&1 && gfortran --version | head -n 1 || true
if ! command -v gfortran >/dev/null 2>&1; then
  echo "ERROR: No Fortran compiler (gfortran) found in PATH on the execute node." >&2
  echo "Install gfortran on the execute nodes, or install it into your conda env and ensure that env's bin/ is on PATH." >&2
  echo "Current PATH=${PATH}" >&2
  exit 4
fi

if [[ -x "${MG5_PYTHON}" ]]; then
  "${MG5_PYTHON}" ./bin/madevent "launch_card_st${SINTHETA_TAG}_tb${TANBETA}_mA${M36_IN}_ma${M55}_mchi${MCHI_TAG}_cbma${COSBMA_TAG}_dm${DELTAM_TAG}.dat"
else
  "${PYTHON_BIN}" ./bin/madevent "launch_card_st${SINTHETA_TAG}_tb${TANBETA}_mA${M36_IN}_ma${M55}_mchi${MCHI_TAG}_cbma${COSBMA_TAG}_dm${DELTAM_TAG}.dat"
fi
wait
# -----------------------------
# Copy events
# -----------------------------
RESULT_DIR="../${OUTROOT}/results/Events_st${SINTHETA_TAG}_tb${TANBETA}_mA${M36_IN}_ma${M55}_mchi${MCHI_TAG}_cbma${COSBMA_TAG}_dm${DELTAM_TAG}_job${JOBID}"
mkdir -p "${RESULT_DIR}"

# Sanity check: MadGraph should have produced at least one banner file.
if [[ ! -d "Events" ]]; then
  echo "ERROR: No Events/ directory produced by MadGraph; job failed." >&2
  exit 5
fi

if ! find "Events" -type f -name "*_banner.txt" -print -quit | grep -q .; then
  echo "ERROR: No *_banner.txt found under Events/. MadGraph likely failed before event generation." >&2
  exit 5
fi

cp -r Events/. "${RESULT_DIR}/"

echo "===================================="
echo "Job finished successfully"
echo "===================================="
