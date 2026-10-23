"""Tests for exact metadata inheritance and topology-preserving Auspice stream labels."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from orthofinder_interrogation_app.research.analysis import auspice_dataset
from orthofinder_interrogation_app.research.common import (
    ResearchError,
    atomic_json,
    file_sha256,
)
from orthofinder_interrogation_app.research.jobs import (
    claim_job,
    finish_job,
    job_directory,
    submit_job,
)
from orthofinder_interrogation_app.research.tree_metadata import inherited_tip_metadata


def complete_metadata_job(*, workspace: Path, parameters: dict, tables: dict[str, str]) -> str:
    """Publish small checksummed metadata tables through the real queue contracts.

    Args:
        workspace: Test queue root.
        parameters: Explicit parent/source controls.
        tables: Named table content.

    Returns:
        A completed job ID with verified request and table artifacts.
    """
    job = submit_job(workspace=workspace, operation="subcluster", parameters=parameters)
    claim_job(workspace=workspace, job_id=job)
    directory = job_directory(workspace=workspace, job_id=job)
    artifacts = {}
    for name, content in tables.items():
        (directory / name).write_text(content, encoding="utf-8")
    for name in ["request.json", *tables]:
        path = directory / name
        artifacts[name] = {
            "path": name,
            "sha256": file_sha256(path=path),
            "bytes": path.stat().st_size,
        }
    atomic_json(path=directory / "manifest.json", value={"job_id": job, "artifacts": artifacts})
    finish_job(workspace=workspace, job_id=job)
    return job


def test_metadata_inherits_species_and_most_recent_subclusters(tmp_path):
    """Local species survive reclustering while newer subcluster assignments take priority."""
    parent = complete_metadata_job(
        workspace=tmp_path,
        parameters={},
        tables={
            "source_context.tsv": "input_id\tspecies_label\na\tPlant A\nb\tPlant B\n",
            "membership.tsv": "sequence_id\tsubcluster_id\na\tOLD\nb\tOLD\n",
        },
    )
    child = complete_metadata_job(
        workspace=tmp_path,
        parameters={"source_job": parent},
        tables={
            "membership.tsv": "sequence_id\tsubcluster_id\na\tNEW\nb\tNEW\n",
        },
    )
    assert inherited_tip_metadata(
        workspace=tmp_path, parameters={"source_job": child}, identifiers=["a"]
    ) == {"a": {"species": "Plant A", "subcluster_id": "NEW"}}
    assert inherited_tip_metadata(workspace=tmp_path, parameters={}, identifiers=["a"]) == {}
    assert (
        inherited_tip_metadata(
            workspace=tmp_path, parameters={"source_job": child}, identifiers=["z"]
        )
        == {}
    )


@pytest.mark.parametrize(
    "content,message",
    [
        ("wrong\tcolumns\na\tX\n", "required columns"),
        ("sequence_id\tsubcluster_id\na\tX\na\tY\n", "Conflicting"),
    ],
)
def test_inheritance_refuses_ambiguous_or_malformed_tables(tmp_path, content, message):
    """Ambiguous metadata do not become asserted phylogenetic traits."""
    job = complete_metadata_job(
        workspace=tmp_path, parameters={}, tables={"membership.tsv": content}
    )
    with pytest.raises(ResearchError, match=message):
        inherited_tip_metadata(
            workspace=tmp_path, parameters={"source_job": job}, identifiers=["a"]
        )


def test_inheritance_rejects_cycles_and_corrupt_requests(tmp_path):
    """Checksummed source traversal detects invalid request structure and dependency cycles."""
    job = complete_metadata_job(workspace=tmp_path, parameters={}, tables={})
    directory = job_directory(workspace=tmp_path, job_id=job)
    request = json.loads((directory / "request.json").read_text())
    for parameters, message in [({"source_job": job}, "source chain"), ([], "Corrupt")]:
        request["parameters"] = parameters
        atomic_json(path=directory / "request.json", value=request)
        manifest = json.loads((directory / "manifest.json").read_text())
        manifest["artifacts"]["request.json"]["sha256"] = file_sha256(
            path=directory / "request.json"
        )
        atomic_json(path=directory / "manifest.json", value=manifest)
        with pytest.raises(ResearchError, match=message):
            inherited_tip_metadata(
                workspace=tmp_path, parameters={"source_job": job}, identifiers=["a"]
            )


def test_auspice_labels_only_supported_multi_tip_subcluster_branches():
    """Mixed branches and isolated members remain ordinary branches in the original topology."""
    dataset = auspice_dataset(
        newick="((a:0.1,b:0.2):0.1,(c:0.1,d:0.2):0.2,e:0.5);",
        metadata={
            "a": {"subcluster_id": "S1"},
            "b": {"subcluster_id": "S1"},
            "c": {"subcluster_id": "S2"},
            "d": {"subcluster_id": "S2"},
        },
        title="Synthetic protein family",
    )
    assert dataset["meta"]["stream_labels"] == ["subcluster"]
    root = dataset["tree"]
    assert "branch_attrs" not in root
    labels = [
        child.get("branch_attrs", {}).get("labels", {}).get("subcluster")
        for child in root["children"]
    ]
    assert labels == ["S1__part_1", "S2__part_1", None]
    assert root["children"][0]["children"][0]["node_attrs"]["div"] == pytest.approx(0.2)


def test_auspice_keeps_nonmonophyletic_parts_distinct_and_unknowns_unlabelled():
    """One sequence-similarity subcluster can legitimately occupy several inferred clades."""
    newick = "((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1,u:0.1);"
    metadata = {name: {"subcluster_id": "S"} for name in "abcd"}
    dataset = auspice_dataset(newick=newick, metadata=metadata, title="Family")
    labels = [
        child.get("branch_attrs", {}).get("labels", {}).get("subcluster")
        for child in dataset["tree"]["children"]
    ]
    assert labels == ["S__part_1", "S__part_2", None]
    unknown = auspice_dataset(newick=newick, metadata={}, title="Unknown traits")
    assert "stream_labels" not in unknown["meta"]
