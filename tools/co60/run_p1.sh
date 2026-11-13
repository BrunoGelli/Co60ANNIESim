#!/usr/bin/env bash
# Local standalone transport only. Explicitly authorized bounded samples, no optics.
set -euo pipefail
if [[ $# != 5 ]]; then
  echo 'usage: bash run_p1.sh PRIVATE_UPLOADED_APP PASSED_P0 VISUALIZATION NEW_PRIVATE_OUTPUT FOAM_DENSITY_G_CM3' >&2
  exit 2
fi
co60_repo=$(cd -- "$(dirname -- "$0")/../.." && pwd)
co60_app=$(realpath -- "$1")
co60_p0=$(realpath -- "$2")
co60_vis=$(realpath -- "$3")
co60_out=$(realpath -m -- "$4")
co60_foam=$5
case "$co60_out/" in "$co60_repo/"*) echo 'Output must be outside public repository' >&2; exit 2;; esac
if [[ -e "$co60_out" ]]; then echo 'Choose a new output directory' >&2; exit 2; fi
mkdir -p "$co60_out/logs"
exec > >(tee "$co60_out/commands.log") 2>&1
set -x
python3 "$co60_repo/tools/co60/prepare_transport.py" --app "$co60_app" --p0 "$co60_p0" \
  --visualization "$co60_vis" --out "$co60_out/inputs" --foam-density-g-cm3 "$co60_foam"
cmake -S "$co60_repo/transport/co60" -B "$co60_out/build" -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "$co60_out/build" -j2
for co60_mode in full pair; do
  if [[ "$co60_mode" == full ]]; then co60_seed=601001; else co60_seed=602001; fi
  "$co60_out/build/co60_transport" "$co60_out/inputs" "$co60_mode" 10000 "$co60_seed" .01 \
    "$co60_out/${co60_mode}_10000.root" 1 > "$co60_out/logs/${co60_mode}_10000.log" 2>&1
  python3 "$co60_repo/tools/co60/analyze_transport.py" validate "$co60_out/${co60_mode}_10000.root" \
    --strict --out "$co60_out/${co60_mode}_10000_validation.json"
done
# Larger files are created only if BOTH 10k strict validations pass.
for co60_mode in full pair; do
  if [[ "$co60_mode" == full ]]; then co60_seed=601002; else co60_seed=602002; fi
  "$co60_out/build/co60_transport" "$co60_out/inputs" "$co60_mode" 100000 "$co60_seed" .01 \
    "$co60_out/${co60_mode}_100000.root" 1 > "$co60_out/logs/${co60_mode}_100000.log" 2>&1
  python3 "$co60_repo/tools/co60/analyze_transport.py" validate "$co60_out/${co60_mode}_100000.root" \
    --strict --out "$co60_out/${co60_mode}_100000_validation.json"
done
"$co60_out/build/co60_transport" "$co60_out/inputs" full 10000 601001 .005 \
  "$co60_out/full_cut005_10000.root" 1 > "$co60_out/logs/full_cut005_10000.log" 2>&1
python3 "$co60_repo/tools/co60/analyze_transport.py" validate "$co60_out/full_cut005_10000.root" \
  --strict --out "$co60_out/full_cut005_10000_validation.json"
"$co60_out/build/co60_transport" "$co60_out/inputs" natural_clock_test 10 600103 .01 \
  "$co60_out/natural_clock_10.root" 1 > "$co60_out/logs/natural_clock_10.log" 2>&1
"$co60_out/build/co60_transport" "$co60_out/inputs" ni60_timing_test 5000 600104 .01 \
  "$co60_out/ni60_timing_5000.root" 1 > "$co60_out/logs/ni60_timing_5000.log" 2>&1
python3 "$co60_repo/tools/co60/analyze_transport.py" plot --full "$co60_out/full_100000.root" \
  --pair "$co60_out/pair_100000.root" --out "$co60_out/analysis"
python3 "$co60_repo/tools/co60/finish_p1.py" --results "$co60_out" --app "$co60_app"
