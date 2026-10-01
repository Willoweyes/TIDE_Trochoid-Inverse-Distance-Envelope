"""Record the hardware/software used for all timings -> results/machine_info.json/.txt"""
import json, platform, subprocess, os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
def sh(cmd):
    try: return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception: return ""
info = dict(os=platform.platform(), machine=platform.machine(), python=sys.version.split()[0],
            logical_cores=os.cpu_count())
if os.name == "nt":
    ps = lambda q: sh(f'powershell -NoProfile -Command "{q}"')
    info["cpu"] = ps("(Get-CimInstance Win32_Processor).Name")
    info["cpu_max_clock_MHz"] = ps("(Get-CimInstance Win32_Processor).MaxClockSpeed")
    info["physical_cores"] = ps("(Get-CimInstance Win32_Processor).NumberOfCores")
    info["ram_GB"] = ps("[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)")
    info["power_plan"] = ps("powercfg /getactivescheme")
else:
    lc = sh("lscpu")
    get = lambda k: next((l.split(":", 1)[1].strip() for l in lc.splitlines() if l.startswith(k)), "")
    info["cpu"] = get("Model name"); info["cpu_max_clock_MHz"] = get("CPU max MHz")
    info["physical_cores"] = get("Core(s) per socket"); info["ram_GB"] = sh("free -g | awk '/Mem:/{print $2}'")
cf = ROOT/"build"/"compiler.txt"
info["compiler"] = cf.read_text().strip() if cf.exists() else ""
for m in ("numpy", "scipy", "sklearn", "matplotlib"):
    try: info[m] = __import__(m).__version__
    except Exception: info[m] = "missing"
(ROOT/"results").mkdir(exist_ok=True)
(ROOT/"results/machine_info.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info, indent=1))
