"""Test independent queues sharing the installed server database and tools."""

import json
from dataclasses import fields

import pytest

from orthofinder_interrogation_app.research.cli import main
from orthofinder_interrogation_app.research.common import ResearchConfig, ResearchError, load_config
from orthofinder_interrogation_app.research.profile_copy import copy_profile


def test_profile_copy_shares_data_without_copying_database_or_jobs(tmp_path, capsys):
    """Only workspace changes; arrays, escapes and every optional route survive."""
    source = tmp_path / "e3.toml"
    source_workspace = tmp_path / "e3-jobs"
    destination = tmp_path / "orthofinder" / "research.toml"
    workspace = tmp_path / "orthofinder-jobs"
    source.write_text(
        f"workspace={json.dumps(str(source_workspace))}\n"
        'index="/mounted/deepclust/index.parquet"\n'
        'sequence_map="/mounted/deepclust/seq2cluster.parquet"\n'
        'sequences="/mounted/deepclust/sequence data/*.parquet"\n'
        'identifier_prefixes=["NR_", "UNIPROT_"]\n'
        'diamond="/existing/e3/bin/diamond"\nthreads=8\n'
        'email="operator@example.org"\nesm_endpoint="https://endpoint.invalid"\n',
        encoding="utf-8",
    )
    before = source.read_bytes()
    assert (
        main(
            argv=[
                "copy-profile",
                "--source",
                str(source),
                "--output",
                str(destination),
                "--workspace",
                str(workspace),
            ]
        )
        == 0
    )
    assert str(destination) in capsys.readouterr().out
    original = load_config(path=source)
    copied = load_config(path=destination)
    for field in fields(ResearchConfig):
        expected = workspace if field.name == "workspace" else getattr(original, field.name)
        assert getattr(copied, field.name) == expected
    assert source.read_bytes() == before
    assert not workspace.exists() and not source_workspace.exists()
    assert not (tmp_path / "mounted").exists()


def test_profile_copy_never_replaces_files_or_shares_queue(tmp_path):
    """Source and destination profiles remain unchanged on invalid requests."""
    source = tmp_path / "source.toml"
    workspace = tmp_path / "existing-jobs"
    source.write_text(f"workspace={json.dumps(str(workspace))}\n", encoding="utf-8")
    output = tmp_path / "output.toml"
    with pytest.raises(ResearchError, match="separate workspace"):
        copy_profile(source=source, output=output, workspace=workspace)
    assert not output.exists()
    with pytest.raises(ResearchError, match="already exists"):
        copy_profile(source=source, output=source, workspace=tmp_path / "different")
    output.write_text("retained", encoding="utf-8")
    with pytest.raises(ResearchError, match="already exists"):
        copy_profile(source=source, output=output, workspace=tmp_path / "different")
    assert output.read_text() == "retained"
    assert (
        main(
            argv=[
                "copy-profile",
                "--source",
                str(tmp_path / "absent"),
                "--output",
                str(output),
                "--workspace",
                str(tmp_path / "different"),
            ]
        )
        == 2
    )
