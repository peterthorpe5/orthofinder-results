"""End-to-end tests for the matched-dispersion Slurm scratch wrapper."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _write_fake_package(*, root: Path) -> Path:
    """Create a fake runner and both bundled biological authorities."""

    data = root / "src/orthofinder_interrogation_app/data"
    data.mkdir(parents=True)
    (data / "e3_seed_catalogue.tsv").write_text(
        "seed_id\nE3_A\n",
        encoding="utf-8",
    )
    (data / "arabidopsis_dispersion_benchmarks.tsv").write_text(
        "marker_id\nAT1G00001\n",
        encoding="utf-8",
    )
    runner = root / "run_orthofinder_results.sh"
    table_names = (
        "e3_cluster_results",
        "benchmark_marker_audit",
        "benchmark_group_profiles",
        "benchmark_matched_controls",
        "benchmark_cluster_results",
        "benchmark_background_statistics",
        "benchmark_contrasts",
        "benchmark_individual_comparisons",
        "benchmark_cluster_classifications",
        "pairwise_distances",
    )
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "printf '%s\\n' \"$@\" > \"$CAPTURE_FILE\"\n"
        "output=''\n"
        "while (($#)); do\n"
        "  if [[ \"$1\" == '--output-dir' ]]; then output=$2; shift 2; else shift; fi\n"
        "done\n"
        "[[ -n \"$output\" ]]\n"
        "mkdir -p \"$output/tables\" \"$output/duckdb\" \"$output/qc\"\n"
        + "".join(
            f"printf 'header\\nrow\\n' | gzip -c > "
            f'"$output/tables/{name}.tsv.gz"\n'
            for name in table_names
        )
        + "printf 'duckdb\\n' > \"$output/duckdb/orthofinder_results.duckdb\"\n"
        "printf 'check_name\\tstatus\\nall\\tPASS\\n' > "
        '"$output/qc/validation_checks.tsv"\n'
        "printf '{\"status\": \"complete\"}\\n' > \"$output/run_manifest.json\"\n",
        encoding="utf-8",
    )
    runner.chmod(0o755)
    return root


def _arguments(*, tmp_path: Path) -> tuple[list[str], Path, Path]:
    """Return valid wrapper arguments, formal output and capture path."""

    package = _write_fake_package(root=tmp_path / "package")
    results = tmp_path / "results"
    results.mkdir()
    (results / "Log.txt").write_text(
        "OrthoFinder version 2.5.5\n",
        encoding="utf-8",
    )
    output = tmp_path / "persistent/dispersion_benchmark"
    capture = tmp_path / "captured_arguments.txt"
    return (
        [
            str(package),
            "orthofinder_results",
            str(results),
            "-",
            "-",
            "benchmark-run",
            str(output),
            "--distance-max-members",
            "50",
            "--benchmark-controls-per-group",
            "3",
        ],
        output,
        capture,
    )


def test_benchmark_wrapper_stages_inputs_and_keeps_formal_output_persistent(
    tmp_path: Path,
) -> None:
    """Heavy inputs and work use TMPDIR while final publication uses GPFS-like storage."""

    arguments, output, capture = _arguments(tmp_path=tmp_path)
    node_tmp = tmp_path / "node_tmp"
    node_tmp.mkdir()
    environment = {
        **os.environ,
        "TMPDIR": str(node_tmp),
        "USER": "test_user",
        "SLURM_JOB_ID": "13579",
        "SLURM_CPUS_PER_TASK": "9",
        "CAPTURE_FILE": str(capture),
    }
    script = Path(__file__).resolve().parents[1] / "slurm/dispersion_benchmark.sbatch"
    completed = subprocess.run(
        ["bash", str(script), *arguments],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert (output / "tables/benchmark_contrasts.tsv.gz").is_file()
    assert "Profile contrasts" in completed.stdout
    captured = capture.read_text(encoding="utf-8").splitlines()
    for option in (
        "--results-dir",
        "--focus-proteins",
        "--benchmark-proteins",
        "--work-dir",
    ):
        index = captured.index(option)
        assert captured[index + 1].startswith(str(node_tmp))
    assert captured[captured.index("--action") + 1] == "dispersion-benchmark"
    assert captured[captured.index("--output-dir") + 1] == str(output)
    assert not tuple(node_tmp.iterdir())
    assert not tuple(output.parent.glob(".*.publish.lock"))


def test_benchmark_wrapper_is_mac_testable_and_protects_managed_options() -> None:
    """Portable date/move syntax and fixed selection controls prevent unsafe overrides."""

    script = Path(__file__).resolve().parents[1] / "slurm/dispersion_benchmark.sbatch"
    source = script.read_text(encoding="utf-8")
    assert "date --iso-8601" not in source
    assert "mv -T" not in source
    assert "${TMPDIR:-/tmp}" not in source
    assert "date -u '+%Y-%m-%dT%H:%M:%SZ'" in source
    assert '--output-dir "$PERSISTENT_OUTPUT"' in source
    assert "--benchmark-proteins" in source
    assert "completed_resource" not in source


def test_benchmark_wrapper_requires_scheduler_tmpdir(tmp_path: Path) -> None:
    """Missing scheduler scratch fails before any source staging or publication."""

    arguments, output, capture = _arguments(tmp_path=tmp_path)
    environment = {
        **os.environ,
        "USER": "test_user",
        "CAPTURE_FILE": str(capture),
    }
    environment.pop("TMPDIR", None)
    script = Path(__file__).resolve().parents[1] / "slurm/dispersion_benchmark.sbatch"
    completed = subprocess.run(
        ["bash", str(script), *arguments],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 2
    assert "scheduler-managed TMPDIR" in completed.stderr
    assert not output.exists()
    assert not capture.exists()
