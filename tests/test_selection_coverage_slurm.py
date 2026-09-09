"""End-to-end tests for the selection-coverage cluster scratch wrapper."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _write_fake_package(*, root: Path) -> Path:
    """Create a package runner that publishes one checksummed fixture artifact."""

    root.mkdir()
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
        'mkdir "$output"\n'
        "printf 'fixture\\n' > \"$output/result.tsv\"\n"
        '(cd "$output" && sha256sum result.tsv > checksums.sha256)\n',
        encoding="utf-8",
    )
    runner.chmod(0o755)
    return root


def _wrapper_arguments(*, tmp_path: Path) -> tuple[list[str], Path, Path]:
    """Return valid positional inputs, final output and captured runner arguments."""

    package = _write_fake_package(root=tmp_path / "package")
    resource = tmp_path / "resource"
    resource.mkdir()
    (resource / "run_manifest.json").write_text("{}\n", encoding="utf-8")
    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text("mapping\n", encoding="utf-8")
    expected = tmp_path / "expected.tsv"
    expected.write_text("expected\n", encoding="utf-8")
    focus = tmp_path / "focus.tsv"
    focus.write_text("protein_identifier\nQ9SA03\n", encoding="utf-8")
    output = tmp_path / "persistent" / "selection_result"
    capture = tmp_path / "captured_arguments.txt"
    arguments = [
        str(package),
        "orthofinder_results",
        str(resource),
        str(taxonomy),
        str(expected),
        str(focus),
        str(output),
        "--include-clade-tax-id",
        "33090",
    ]
    return arguments, output, capture


def test_selection_wrapper_stages_verifies_and_atomically_publishes(
    tmp_path: Path,
) -> None:
    """All heavy reads and initial output occur below scheduler-managed TMPDIR."""

    arguments, output, capture = _wrapper_arguments(tmp_path=tmp_path)
    node_tmp = tmp_path / "node_tmp"
    node_tmp.mkdir()
    environment = {
        **os.environ,
        "TMPDIR": str(node_tmp),
        "USER": "test_user",
        "SLURM_JOB_ID": "12345",
        "CAPTURE_FILE": str(capture),
    }
    script = Path(__file__).resolve().parents[1] / "slurm/selection_coverage_tree.sbatch"
    completed = subprocess.run(
        ["bash", str(script), *arguments],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert (output / "result.tsv").read_text(encoding="utf-8") == "fixture\n"
    assert "result.tsv: OK" in completed.stdout
    captured = capture.read_text(encoding="utf-8").splitlines()
    resource_index = captured.index("--resource-dir")
    output_index = captured.index("--output-dir")
    assert captured[resource_index + 1].startswith(str(node_tmp))
    assert captured[output_index + 1].startswith(str(node_tmp))
    assert not tuple((tmp_path / "persistent").glob(".*.incoming.*"))
    assert not tuple((tmp_path / "persistent").glob(".*.publish.lock"))
    assert not tuple(node_tmp.iterdir())


def test_selection_wrapper_uses_portable_date_and_move_commands() -> None:
    """The wrapper avoids GNU-only options so its end-to-end test runs on macOS."""

    script = Path(__file__).resolve().parents[1] / "slurm/selection_coverage_tree.sbatch"
    source = script.read_text(encoding="utf-8")
    assert "date --iso-8601" not in source
    assert "mv -T" not in source
    assert "date -u '+%Y-%m-%dT%H:%M:%SZ'" in source


def test_selection_wrapper_preserves_an_existing_publication_lock(
    tmp_path: Path,
) -> None:
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
        "SLURM_JOB_ID": "12345",
        "CAPTURE_FILE": str(capture),
    }
    script = Path(__file__).resolve().parents[1] / "slurm/selection_coverage_tree.sbatch"
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


def test_selection_wrapper_requires_scheduler_managed_tmpdir(tmp_path: Path) -> None:
    """The cluster wrapper fails before work when TMPDIR is absent."""

    arguments, output, capture = _wrapper_arguments(tmp_path=tmp_path)
    environment = {**os.environ, "USER": "test_user", "CAPTURE_FILE": str(capture)}
    environment.pop("TMPDIR", None)
    script = Path(__file__).resolve().parents[1] / "slurm/selection_coverage_tree.sbatch"
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
