"""Portable shell-wrapper tests and optional unmodified upstream extraction checks."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

SCRIPT = Path(__file__).parents[1] / "scripts" / "extract_deepclust_selected.sh"


@pytest.fixture
def extraction_arguments(tmp_path: Path) -> list[str]:
    """Create a small private database layout and an argv-recording upstream stand-in."""
    repository = tmp_path / "retrieval"
    repository.mkdir()
    (repository / "cluster_index.py").write_text(
        "import sys\nprint('\\n'.join(sys.argv[1:]))\n", encoding="utf-8"
    )
    centroids = tmp_path / "centroids.txt"
    centroids.write_text("NR_UGC79169.1\n", encoding="utf-8")
    database = tmp_path / "database.parquet"
    database.touch()
    index = tmp_path / "index.parquet"
    index.touch()
    return [
        "bash",
        str(SCRIPT),
        "--retrieval-repo",
        str(repository),
        "--centroid-file",
        str(centroids),
        "--database",
        str(database),
        "--index",
        str(index),
        "--output",
        str(tmp_path / "output"),
        "--python",
        sys.executable,
    ]


def test_wrapper_preserves_upstream_named_arguments(extraction_arguments):
    """Single-dash source options and per-cluster output reach the exact upstream interface."""
    result = subprocess.run(args=extraction_arguments, check=True, capture_output=True, text=True)
    arguments = result.stdout.splitlines()
    assert "-path_to_DCD" in arguments and "--path_to_DCD" not in arguments
    assert arguments[arguments.index("--per-clust-output") + 1] == "1"
    assert arguments[arguments.index("--threads") + 1] == "1"
    assert arguments[arguments.index("--centroid_file") + 1].endswith("centroids.txt")
    assert subprocess.run(args=["bash", "-n", str(SCRIPT)], check=False).returncode == 0
    assert subprocess.run(args=["bash", str(SCRIPT), "--help"], capture_output=True).returncode == 0


@pytest.mark.parametrize(
    "extra", [["--threads", "0"], ["--threads", "257"], ["--wrong"], ["--index"]]
)
def test_wrapper_rejects_invalid_controls(extraction_arguments, extra):
    """Malformed options never start the upstream process."""
    result = subprocess.run(args=extraction_arguments + extra, capture_output=True, text=True)
    assert result.returncode == 2


@pytest.mark.parametrize("centroid", ["", "NR_X\tY", "../escape", "NR_X' OR 1=1", ".."])
def test_wrapper_rejects_unsafe_or_nonexact_centroids(extraction_arguments, centroid):
    """Upstream SQL interpolation and raw filename construction receive only safe exact IDs."""
    path = Path(extraction_arguments[extraction_arguments.index("--centroid-file") + 1])
    path.write_text(centroid + "\n", encoding="utf-8")
    result = subprocess.run(args=extraction_arguments, capture_output=True, text=True)
    assert result.returncode == 2


def test_wrapper_refuses_nonempty_output(extraction_arguments):
    """Existing scientific exports remain intact when an output path is reused."""
    output = Path(extraction_arguments[extraction_arguments.index("--output") + 1])
    output.mkdir()
    (output / "keep.txt").write_text("unchanged", encoding="utf-8")
    result = subprocess.run(args=extraction_arguments, capture_output=True, text=True)
    assert result.returncode == 2 and (output / "keep.txt").read_text() == "unchanged"


@pytest.mark.skipif(
    not os.environ.get("DEEPCLUST_RETRIEVAL_REPO"),
    reason="Set DEEPCLUST_RETRIEVAL_REPO to a separately cloned unmodified upstream checkout",
)
def test_unmodified_upstream_extracts_spanning_cluster(tmp_path):
    """A real upstream script extracts only the requested cluster across physical row groups."""
    database = tmp_path / "sequences.parquet"
    pq.write_table(
        table=pa.table(
            {
                "f0": ["NR_TARGET", "NR_MEMBER1", "NR_MEMBER2", "NR_OTHER"],
                "f1": ["ACDE", "ACDF", "ACDG", "MMMM"],
                "f2": ["NR_TARGET", "NR_TARGET", "NR_TARGET", "NR_OTHER"],
                "f3": [0, 1, 2, 3],
            }
        ),
        where=database,
        row_group_size=2,
    )
    index = tmp_path / "index.parquet"
    pq.write_table(
        table=pa.table(
            {
                "RowStart": [0, 3],
                "NMEMBER": [3, 1],
                "CLUSTER": ["NR_TARGET", "NR_OTHER"],
                "RowGroup": ["0,1", "1"],
            }
        ),
        where=index,
    )
    centroids = tmp_path / "centroids.txt"
    centroids.write_text("NR_TARGET\n", encoding="utf-8")
    output = tmp_path / "export"
    result = subprocess.run(
        args=[
            "bash",
            str(SCRIPT),
            "--retrieval-repo",
            os.environ["DEEPCLUST_RETRIEVAL_REPO"],
            "--centroid-file",
            str(centroids),
            "--database",
            str(database),
            "--index",
            str(index),
            "--output",
            str(output),
            "--python",
            sys.executable,
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    exported = (output / "NR_TARGET.fa").read_text()
    assert exported.count(">") == 3 and ">NR_TARGET" in exported and "NR_OTHER" not in exported
