"""Tests for the portable package launcher shell script."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def test_conda_runner_prefers_environment_cpp_runtime(tmp_path: Path) -> None:
    """Conda's C++ libraries precede inherited system libraries."""

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_prefix = tmp_path / "conda/envs/orthofinder_results"
    (fake_prefix / "lib").mkdir(parents=True)
    capture = tmp_path / "arguments.txt"
    conda = fake_bin / "conda"
    conda.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "if [[ $# -eq 4 && $1 == run && $2 == --name && $4 == env ]]; then\n"
        "  printf 'CONDA_PREFIX=%s\\n' \"$FAKE_CONDA_PREFIX\"\n"
        "  exit 0\n"
        "fi\n"
        "printf '%s\\n' \"$@\" > \"$CAPTURE_FILE\"\n",
        encoding="utf-8",
    )
    conda.chmod(0o755)
    repository = Path(__file__).resolve().parents[1]
    environment = {
        **os.environ,
        "PATH": f"{fake_bin}:{os.environ['PATH']}",
        "FAKE_CONDA_PREFIX": str(fake_prefix),
        "CAPTURE_FILE": str(capture),
        "LD_LIBRARY_PATH": "/system/lib",
    }
    completed = subprocess.run(
        [
            "bash",
            str(repository / "run_orthofinder_results.sh"),
            "--conda-env",
            "orthofinder_results",
            "--python-executable",
            "python3",
            "--",
            "--action",
            "inspect",
        ],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    arguments = capture.read_text(encoding="utf-8").splitlines()
    assert arguments[:5] == [
        "run",
        "--no-capture-output",
        "--name",
        "orthofinder_results",
        "env",
    ]
    assert arguments[5] == (
        f"LD_LIBRARY_PATH={fake_prefix}/lib:/system/lib"
    )
    assert arguments[6:] == [
        "python3",
        "-m",
        "orthofinder_results",
        "--action",
        "inspect",
    ]
