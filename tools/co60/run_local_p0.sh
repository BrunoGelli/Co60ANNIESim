#!/usr/bin/env bash
# Geometry/navigation and synthetic ROOT I/O only. No transport, no RAT baseline.
set -euo pipefail
if [[ $# != 2 ]]; then
  echo 'usage: bash run_local_p0.sh PRIVATE_UPLOADED_APP NEW_PRIVATE_RESULT_DIR' >&2
  exit 2
fi
co60_repo=$(cd -- "$(dirname -- "$0")/../.." && pwd)
co60_app=$(realpath -- "$1")
co60_out=$(realpath -m -- "$2")
case "$co60_out/" in "$co60_repo/"*) echo 'Results must be outside public repo' >&2; exit 2;; esac
if [[ -e "$co60_out" ]]; then echo 'Choose a new output directory' >&2; exit 2; fi
mkdir -p -- "$co60_out"
exec > >(tee "$co60_out/commands.log") 2>&1
set -x
cmake -S "$co60_repo/tests/co60" -B "$co60_out/build"
cmake --build "$co60_out/build" -j2
python3 "$co60_repo/tools/co60/prepare_geometry.py" --app "$co60_app" --out "$co60_out/geometry" \
  --disk-thickness-mm 2 --radial-clearance-mm .25 --axial-clearance-mm .1 \
  --active-radius-mm 1 --active-thickness-mm .01 --active-z-mm 0 \
  --bgo-bore-radius-mm 1 --plastic-material polystyrene --source-position-mm 0 750 0
python3 "$co60_repo/tools/co60/probe_inputs.py" "$co60_out/geometry"
"$co60_out/build/geometry_probe" "$co60_out/geometry/assembly.tsv" \
  "$co60_out/geometry/vessel.txt" "$co60_out/geometry/probes.tsv"
python3 "$co60_repo/tests/co60/geometry_scan.py" --app "$co60_app" \
  --probe "$co60_out/build/geometry_probe" --out "$co60_out/scans"
python3 "$co60_repo/tests/co60/root_smoke.py" "$co60_out/SYNTHETIC_schema_smoke.root"
c++ -std=c++17 "$co60_app/utilities/pmtinfo_to_ratdb.cpp" -o "$co60_out/pmtinfo_to_ratdb"
"$co60_out/pmtinfo_to_ratdb" "$co60_app/ratdb/ANNIE/PMTPositions_Scan_Glass.txt" \
  "$co60_out/PMTINFO_regenerated.ratdb"
python3 "$co60_repo/tools/co60/check_pmt_tables.py" \
  "$co60_app/ratdb/ANNIE/PMTINFO_run2_config14_scan.ratdb" \
  "$co60_out/PMTINFO_regenerated.ratdb" "$co60_app/ratdb/ANNIE/ANNIE_run2_config14_scan.geo"
