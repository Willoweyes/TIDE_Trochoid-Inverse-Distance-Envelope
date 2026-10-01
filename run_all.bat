@echo off
REM Full reproduction on Windows. Needs: CMake, a C++17 compiler (VS 2022 or MinGW), Python 3.10+.
cd /d %~dp0
python -m pip install -r requirements.txt || goto :err
cmake -S . -B build -A x64 || goto :err
cmake --build build --config Release || goto :err
set OMP_NUM_THREADS=1
set OPENBLAS_NUM_THREADS=1
set MKL_NUM_THREADS=1
set TIDE_REPS=3
python scripts\machine_info.py || goto :err
python scripts\run_experiments.py || goto :err
python scripts\experiments\pipeline_study.py || goto :err
python scripts\experiments\radial_shortcut.py || goto :err
python scripts\experiments\timestep_check.py || goto :err
python scripts\experiments\newton_steps.py 1000 || goto :err
for %%f in (fig_schematic fig_convergence fig_maps fig_cost fig_pipeline) do python scripts\plotting\%%f.py || goto :err
python scripts\paper_numbers.py || goto :err
echo DONE. See results\paper_numbers.txt and figures\
exit /b 0
:err
echo FAILED & exit /b 1
