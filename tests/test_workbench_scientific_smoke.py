"""CLI validation and optional real external-tool workflow tests."""

from __future__ import annotations

import os
import runpy
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parents[1] / "scripts" / "smoke_protein_workbench.py"


def test_scientific_smoke_reports_missing_tools_and_preserves_output(tmp_path, monkeypatch):
    """Missing executables and non-empty destinations return controlled errors."""
    namespace = runpy.run_path(path_name=str(SCRIPT))
    main = namespace["main"]
    monkeypatch.setattr(namespace["shutil"], "which", lambda **kwargs: None)
    assert main(argv=["--output", str(tmp_path / "new")]) == 2
    (tmp_path / "keep.txt").write_text("keep", encoding="utf-8")
    assert main(argv=["--output", str(tmp_path)]) == 2
    assert (tmp_path / "keep.txt").read_text() == "keep"


@pytest.mark.skipif(
    not os.environ.get("WORKBENCH_SCIENTIFIC_CONFIG"),
    reason="Set WORKBENCH_SCIENTIFIC_CONFIG to a profile with DIAMOND >=2.2.1, MAFFT and FastTree",
)
def test_real_scientific_smoke(tmp_path):
    """The installed tools produce checksummed clustering, alignment and Auspice artifacts."""
    main = runpy.run_path(path_name=str(SCRIPT))["main"]
    assert (
        main(
            argv=[
                "--config",
                os.environ["WORKBENCH_SCIENTIFIC_CONFIG"],
                "--output",
                str(tmp_path / "real"),
            ]
        )
        == 0
    )
