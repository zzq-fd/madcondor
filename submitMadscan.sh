#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")"

usage() {
  cat <<'EOF'
Usage:
  ./submitMadscan.sh <PROCESS_TAR>

Example:
  ./submitMadscan.sh bbdm_2HDMa_type1_case1_sint0p7_scan.tar.gz

Notes:
  - The output folder will be the tarball name without the .tar.gz suffix.
  - stdout/stderr/condor.log and results will be placed under that folder.
EOF
}

# Pick the process tarball (required).
if [[ $# -lt 1 ]]; then
  echo "ERROR: Missing required argument <PROCESS_TAR>" >&2
  usage >&2
  exit 2
fi

PROCESS_TAR=$1

if [[ "${PROCESS_TAR}" == "-h" || "${PROCESS_TAR}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ ! -f "${PROCESS_TAR}" ]]; then
  echo "ERROR: Tarball not found: ${PROCESS_TAR}" >&2
  echo "In directory: $PWD" >&2
  echo "Available tarballs:" >&2
  ls -1 *.tar.gz 2>/dev/null >&2 || true
  echo >&2
  usage >&2
  exit 2
fi

if [[ -z "${PROCESS_TAR}" || ! -f "${PROCESS_TAR}" ]]; then
  echo "ERROR: Could not find process tarball in $PWD" >&2
  echo "Expected '${PROCESS_TAR}' or a file matching 'bbdm_2HDMa_type1_case1*_scan.tar.gz'" >&2
  exit 2
fi

OUTROOT="${PROCESS_TAR%.tar.gz}"

# Create output folder structure on the submit machine.
# HTCondor will write stdout/stderr/log here; directories must exist.
mkdir -p \
  "${OUTROOT}/logs/output" \
  "${OUTROOT}/logs/error" \
  "${OUTROOT}/logs/log" \
  "${OUTROOT}/results"

# Submit with explicit macro overrides so the submit file doesn't need manual edits.
exec condor_submit \
  -append "PROCESS_TAR=${PROCESS_TAR}" \
  -append "OUTROOT=${OUTROOT}" \
  subMadscan.sub
