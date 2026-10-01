#!/usr/bin/env bash
# Full reproduction on Linux/macOS.
set -e; cd "$(dirname "$0")"
python3 -m pip install -r requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 TIDE_REPS=3
python3 scripts/machine_info.py
python3 scripts/run_experiments.py
python3 scripts/experiments/pipeline_study.py
python3 scripts/experiments/radial_shortcut.py
python3 scripts/experiments/timestep_check.py
python3 scripts/experiments/newton_steps.py 1000
for f in fig_schematic fig_convergence fig_maps fig_cost fig_pipeline; do python3 scripts/plotting/$f.py; done
python3 scripts/paper_numbers.py
