# TIDE – Trochoid Inverse-Distance Envelope

Reference implementation and replication code for

> *TIDE: A Discretisation-Free Point-Wise Evaluation of the Kinematic Surface in Face Milling with Round Inserts*  
> Duong Duc Tri, Nguyen Quoc Chi – Faculty of Mechanical Engineering, HCMUT.

TIDE computes the height of the kinematic surface of face milling with round inserts **at individual points**, without a grid, time stepping or edge discretisation: per insert and point, at most six analytically seeded contact times are refined by Newton's method (see [`METHOD.md`](METHOD.md)). The repository also contains a C++ implementation of the forward Z-map method (FSM) of Bakhshan et al. (*An efficient open-source framework for high-fidelity 3D surface topography and roughness prediction in milling*, arXiv:2603.28160), used as the comparison method, and every experiment, figure script and result file of the study.

## 1. Quick start (about 5 minutes)

Requirements: CMake ≥ 3.15, a C++17 compiler (GCC/Clang, or MSVC from Visual Studio 2022), Python ≥ 3.10.

```bash
pip install -r requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release      # Windows/MSVC: cmake -S . -B build -A x64
cmake --build build --config Release -j
OMP_NUM_THREADS=1 ./build/tide_tools tide --config benchmark/test_cases/paper_A1.json --grid-base 100
```

Expected (up to the timing `t_s`): `"Sa":1.58088191701734,"Sq":1.843714875967,"Sz":9.93395141174691,"n_points":20301`.
On Windows with MinGW use `-G "MinGW Makefiles"`; with MSVC the executables are in `build/Release/`.

## 2. Reproduce everything

| Platform | Command |
|---|---|
| Windows (MSVC) | `run_all.bat` |
| Linux / macOS | `bash run_all.sh` |
| Any (Python driver with logs and resume, build the code first) | `python run_local.py --fresh` |

All scripts build with `-O3`/`/O2` and **no fast-math**, run every experiment **single-threaded** (`OMP_NUM_THREADS=1`) and report the **median of three runs** (`TIDE_REPS=3`). A full run takes about 1–2 h on a modern 8-core desktop/laptop CPU (most of it FSM at 1600–3200 cells and the hybrid-model study). A single step: `python scripts/run_experiments.py <step>` with step ∈ `verify sampling ladder cost resolution surfaces scaling`; the hybrid-model study is `python scripts/experiments/pipeline_study.py [cloud ref cores learn]` (`learn` reads the shipped `data/raw/pipeline_reference.csv` and `pipeline_cores.csv` and runs in seconds; `ref` is the long step, about 1 h on 2 cores).

**What to expect.** Accuracy results (heights, Sa, Sz, errors, convergence orders, verification against the brute-force search, the hybrid-model errors) are deterministic and reproduce to round-off (the largest relative difference we observed between two platforms and compilers is 5×10⁻¹², in Sa). Wall-clock times depend on the machine, compiler and C library; the ratios between TIDE and FSM change with them (see `results/machine_info.json` for the machine of the shipped results).

## 3. Where each result comes from

| Item | Step | Output (shipped in this repository) |
|---|---|---|
| Schematic of the method | – | `scripts/plotting/fig_schematic.py` |
| Cases A1–A3 | – | `benchmark/test_cases/paper_A*.json` |
| Verification against a brute-force crossing search | `verify` | `data/processed/edge_verification.txt` |
| Newton steps needed (residual of the crossing after 0–4 steps) | `python scripts/experiments/newton_steps.py 1000` | `data/processed/newton_steps.txt` |
| Radial approximation vs edge-ray contact (supplementary) | `experiments/radial_shortcut.py` | `data/processed/radial_shortcut.json` |
| TIDE runtime vs cutting speed and feed | `scaling` | `data/raw/scaling.csv` |
| TIDE sampling study, runtime per point | `sampling` | `data/raw/tide_sampling.csv` |
| Sampling condition r < r_max | `resolution` | `data/raw/resolution_condition.csv` |
| Refinement of FSM and TIDE (released and tied time step) | `ladder` | `data/raw/ladder_dt_regimes.csv` (`plotting/fig_convergence.py`) |
| Sixteen-fold smaller time step | `experiments/timestep_check.py` | `data/processed/timestep_check.json` |
| Height maps and profile of A1 | `surfaces` | `data/raw/surfaces/A1_fsm_vs_tide.npz` (`plotting/fig_maps.py`) |
| FSM over 64 numerical settings | `experiments/pipeline_study.py cloud` | `data/raw/fsm_settings_cloud.csv` (`plotting/fig_cost.py`) |
| Runtime vs error on doubled grids | `cost` | `data/raw/cost_accuracy_dense.csv` (`plotting/fig_cost.py`) |
| Hybrid model: reference, cores, learning | `experiments/pipeline_study.py ref cores learn` | `data/raw/pipeline_reference.csv`, `pipeline_cores.csv`, `data/processed/pipeline_study.json` (`plotting/fig_pipeline.py`) |
| All numbers quoted in the text | `python scripts/paper_numbers.py` | `results/paper_numbers.txt` |

The figures are redrawn by `scripts/plotting/*.py` into `figures/` (PDF, SVG, 600-dpi PNG; not shipped, created on demand).

## 4. Layout

```
include/tide_core.hpp    TIDE: seeds -> realignment -> 4 Newton steps on f(t) -> lower envelope (header only)
include/fsm_core.hpp     forward Z-map (FSM), C++ implementation of the algorithm of Bakhshan et al.
include/tide_json.hpp    case-file reader
src/tide/tide_tools.cpp  driver: modes fsm | tide | points
src/benchmarks/          TIDE runtime vs grid, cutting speed and feed
src/reference_python/    Python TIDE (verification, schematic)
src/reference_efsm/      Python FSM reference
benchmark/test_cases/    cases A1–A3
scripts/                 run_experiments.py, experiments/, plotting/, paper_numbers.py, machine_info.py
data/, results/          results of the shipped run
compare_bakhshan/        comparison with the original released FSM code (section 5)
run_all.bat, run_all.sh, run_local.py   one-shot reproduction
```

Timings are single-thread wall-clock times (`std::chrono::steady_clock`) measured inside `tide_tools`; process start-up and file I/O are not included. TIDE has no reduction across points, so its result is bitwise identical with and without OpenMP; `OMP_NUM_THREADS=n` runs it on `n` threads. FSM uses a time step tied to the grid (r = v_c Δt/Δy = 0.5) unless marked `released` (Δt = 8·10⁻⁷ s).

## 5. Comparison with the original FSM code (optional)

`compare_bakhshan/` checks that `include/fsm_core.hpp` reproduces the original code node by node, and times both on one and on many threads. The original code is **not** redistributed here. To run it:

```bash
git clone https://github.com/HadiBakhshan/surf-topo.git && git -C surf-topo checkout 1f5d3e0
# compile the unmodified original source with a thin driver (needs Eigen 3.4 headers)
g++ -std=c++17 -O3 -march=native -fopenmp -DEIGEN_NO_DEBUG -I surf-topo/src -I /path/to/eigen -I include \
    compare_bakhshan/bakhshan_driver.cpp surf-topo/src/simulation.cpp -o compare_bakhshan/bakhshan_gcc_native
g++ -std=c++17 -O3 -fopenmp -DEIGEN_NO_DEBUG -I surf-topo/src -I /path/to/eigen -I include \
    compare_bakhshan/bakhshan_driver.cpp surf-topo/src/simulation.cpp -o compare_bakhshan/bakhshan_gcc_generic
# optional: build its Python module (see its README) and point to it
export SURFTOPO_BIN=$PWD/surf-topo/bin PYD_PYTHON=python
python compare_bakhshan/run_compare.py            # steps: acc time threads mt3 (about 35 min)
```

The results of the shipped run are in `compare_bakhshan/results/*.csv`.

## License and acknowledgements

MIT License (see `LICENSE`). The forward method follows the kinematic model and test cases released by Bakhshan et al. (MIT licensed, <https://github.com/HadiBakhshan/surf-topo>).
