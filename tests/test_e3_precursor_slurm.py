"""End-to-end tests for the E3 precursor scheduler scratch wrapper."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _write_fake_package(*, root: Path) -> Path:
    """Create a fake package runner and bundled E3 authority.

    Args:
        root: New fake checkout directory.

    Returns:
        Fake package root.
    """

    data = root / "src/orthofinder_interrogation_app/data"
    data.mkdir(parents=True)
    (data / "e3_seed_catalogue.tsv").write_text(
        "seed_id\nE3_A\n",
        encoding="utf-8",
    )
    runner = root / "run_orthofinder_results.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'printf \'%s\\n\' "$@" > "$CAPTURE_FILE"\n'
        "output=''\n"
        "while (($#)); do\n"
        "  if [[ \"$1\" == '--output-dir' ]]; then output=$2; shift 2; else shift; fi\n"
        "done\n"
        '[[ -n "$output" ]]\n'
        'mkdir -p "$output/tables" "$output/duckdb" "$output/qc"\n'
        "for name in e3_cluster_results e3_seed_catalogue_audit "
        "e3_seed_matches pairwise_distances; do\n"
        "  printf 'header\\nrow\\n' | gzip -c > \"$output/tables/${name}.tsv.gz\"\n"
        "done\n"
        "printf 'duckdb\\n' > \"$output/duckdb/orthofinder_results.duckdb\"\n"
        "printf 'check_name\\tstatus\\nall\\tPASS\\n' > "
        '"$output/qc/validation_checks.tsv"\n'
        "printf '{\"status\": \"complete\"}\\n' > \"$output/run_manifest.json\"\n",
        encoding="utf-8",
    )
    runner.chmod(0o755)
    return root


def _wrapper_arguments(*, tmp_path: Path) -> tuple[list[str], Path, Path]:
    """Return valid wrapper arguments, output and capture paths.

    Args:
        tmp_path: Pytest temporary directory.

    Returns:
        Argument list, final output and captured runner arguments.
    """

    package = _write_fake_package(root=tmp_path / "package")
    results = tmp_path / "results"
    results.mkdir()
    (results / "Log.txt").write_text("OrthoFinder version 2.5.5\n", encoding="utf-8")
    output = tmp_path / "persistent/e3_precursor"
    capture = tmp_path / "captured_arguments.txt"
    return (
        [
            str(package),
            "orthofinder_results",
            str(results),
            "-",
            "e3-run",
            str(output),
            "--distance-max-members",
            "50",
        ],
        output,
        capture,
    )


def test_e3_wrapper_stages_verifies_and_atomically_publishes(
    tmp_path: Path,
) -> None:
    """Source reads and initial analytical outputs remain under scheduler TMPDIR."""

    arguments, output, capture = _wrapper_arguments(tmp_path=tmp_path)
    node_tmp = tmp_path / "node_tmp"
    node_tmp.mkdir()
    environment = {
        **os.environ,
        "TMPDIR": str(node_tmp),
        "USER": "test_user",
        "SLURM_JOB_ID": "24680",
        "SLURM_CPUS_PER_TASK": "9",
        "CAPTURE_FILE": str(capture),
    }
    script = Path(__file__).resolve().parents[1] / "slurm/e3_precursor.sbatch"
    completed = subprocess.run(
        ["bash", str(script), *arguments],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert (output / "tables/e3_cluster_results.tsv.gz").is_file()
    assert "Primary cluster table" in completed.stdout
    captured = capture.read_text(encoding="utf-8").splitlines()
    for option in ("--results-dir", "--focus-proteins", "--output-dir", "--work-dir"):
        index = captured.index(option)
        assert captured[index + 1].startswith(str(node_tmp))
    assert captured[captured.index("--action") + 1] == "e3-precursor"
    assert captured[captured.index("--run-id") + 1] == "e3-run"
    assert not tuple((tmp_path / "persistent").glob(".*.incoming.*"))
    assert not tuple((tmp_path / "persistent").glob(".*.publish.lock"))
    assert not tuple(node_tmp.iterdir())


def test_e3_wrapper_uses_portable_date_and_move_commands() -> None:
    """The wrapper remains testable on macOS despite executing on Linux Slurm."""

    script = Path(__file__).resolve().parents[1] / "slurm/e3_precursor.sbatch"
    source = script.read_text(encoding="utf-8")
    assert "date --iso-8601" not in source
    assert "mv -T" not in source
    assert "${TMPDIR:-/tmp}" not in source
    assert "date -u '+%Y-%m-%dT%H:%M:%SZ'" in source


def test_e3_wrapper_preserves_existing_publication_lock(tmp_path: Path) -> None:
    """A concurrent publisher's lock is reported and never removed."""

    arguments, output, capture = _wrapper_arguments(tmp_path=tmp_path)
    node_tmp = tmp_path / "node_tmp"
    node_tmp.mkdir()
    output.parent.mkdir()
    lock = output.parent / f".{output.name}.publish.lock"
    lock.mkdir()
    environment = {
        **os.environ,
        "TMPDIR": str(node_tmp),
        "USER": "test_user",
        "SLURM_JOB_ID": "24680",
        "CAPTURE_FILE": str(capture),
    }
    script = Path(__file__).resolve().parents[1] / "slurm/e3_precursor.sbatch"
    completed = subprocess.run(
        ["bash", str(script), *arguments],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 2
    assert "another publication holds the output lock" in completed.stderr
    assert lock.is_dir()
    assert not output.exists()
    assert not tuple(output.parent.glob(".*.incoming.*"))
    assert not tuple(node_tmp.iterdir())


def test_e3_wrapper_requires_scheduler_managed_tmpdir(tmp_path: Path) -> None:
    """The cluster wrapper fails before reading sources when TMPDIR is absent."""

    arguments, output, capture = _wrapper_arguments(tmp_path=tmp_path)
    environment = {**os.environ, "USER": "test_user", "CAPTURE_FILE": str(capture)}
    environment.pop("TMPDIR", None)
    script = Path(__file__).resolve().parents[1] / "slurm/e3_precursor.sbatch"
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
