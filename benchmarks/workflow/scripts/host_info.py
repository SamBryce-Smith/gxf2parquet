#!/usr/bin/env python3
"""Emit a small JSON blob describing the host: CPU model, arch, cores,
memory, kernel, OS. Written once per pipeline run and joined onto every
row of the bench TSVs at aggregation time.

Adapted (AWS probe removed) from fulcrumgenomics/riker benchmark-pipeline (MIT License,
Copyright (c) 2026 Fulcrum Genomics LLC)."""

import json
import os
import platform
import re
import subprocess
import sys


def cpu_model() -> str:
    if sys.platform == "darwin":
        try:
            return subprocess.check_output(
                ["sysctl", "-n", "machdep.cpu.brand_string"], text=True
            ).strip()
        except Exception:
            return "unknown"
    # Linux
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith(("model name", "Model", "cpu model")):
                    return line.split(":", 1)[1].strip()
            # aarch64 /proc/cpuinfo often only has implementer/part — probe.
            fh.seek(0)
            for line in fh:
                if line.startswith("CPU part"):
                    return f"aarch64 CPU part {line.split(':')[1].strip()}"
    except FileNotFoundError:
        pass
    return "unknown"


def total_mem_bytes() -> int:
    if sys.platform == "darwin":
        try:
            return int(subprocess.check_output(["sysctl", "-n", "hw.memsize"]).strip())
        except Exception:
            return 0
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(re.search(r"(\d+)", line).group(1)) * 1024
    except FileNotFoundError:
        pass
    return 0


def main():
    info = {
        "hostname": platform.node(),
        "os": platform.system(),
        "os_release": platform.release(),
        "arch": platform.machine(),
        "cpu_model": cpu_model(),
        "cpu_count_logical": os.cpu_count(),
        "total_mem_bytes": total_mem_bytes(),
        "python": platform.python_version(),
    }
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
