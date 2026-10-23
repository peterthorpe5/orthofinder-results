"""Trace exact input traits and label supported protein subclades for Auspice."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

from .common import ResearchError
from .jobs import artifact_path


def inherited_tip_metadata(
    *, workspace: Path, parameters: Mapping[str, object], identifiers: Sequence[str]
) -> dict[str, dict[str, str]]:
    """Recover audited species and latest subcluster assignments from checksummed ancestors.

    Args:
        workspace: Shared queue workspace.
        parameters: Current job controls with an optional source_job reference.
        identifiers: Exact alignment tip IDs.

    Returns:
        Available source traits for the selected tips; absent traits remain unknown.

    Raises:
        ResearchError: If ancestor data are corrupt, cyclic or exceed 32 processing stages.
    """
    selected = set(identifiers)
    result: dict[str, dict[str, str]] = {}
    source = parameters.get("source_job")
    seen: set[str] = set()
    while source:
        if not isinstance(source, str) or source in seen or len(seen) >= 32:
            raise ResearchError("Invalid or excessively deep research source chain")
        seen.add(source)
        request_path = artifact_path(workspace=workspace, job_id=source, artifact="request.json")
        manifest_path = request_path.parent / "manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            request = json.loads(request_path.read_text(encoding="utf-8"))
            artifacts = manifest["artifacts"]
            controls = request["parameters"]
            if not isinstance(artifacts, dict) or not isinstance(controls, dict):
                raise ValueError("Invalid ancestor contract")
        except (OSError, KeyError, ValueError, TypeError) as error:
            raise ResearchError("Corrupt source metadata in the research chain") from error
        for artifact, key, trait, column in (
            ("membership.tsv", "sequence_id", "subcluster_id", "subcluster_id"),
            ("source_context.tsv", "input_id", "species", "species_label"),
        ):
            if artifact not in artifacts:
                continue
            path = artifact_path(workspace=workspace, job_id=source, artifact=artifact)
            values: dict[str, set[str]] = {}
            with path.open(encoding="utf-8") as handle:
                reader = csv.DictReader(f=handle, delimiter="\t")
                if not {key, column}.issubset(reader.fieldnames or ()):
                    raise ResearchError("Source metadata table lacks its required columns")
                for row in reader:
                    if row[key] in selected and row[column]:
                        values.setdefault(row[key], set()).add(row[column])
            for identifier, observed in values.items():
                if len(observed) != 1:
                    raise ResearchError("Conflicting source traits for one exact protein ID")
                result.setdefault(identifier, {}).setdefault(trait, next(iter(observed)))
        source = controls.get("source_job")
    return result


def label_streamtree_partitions(*, dataset: dict[str, object]) -> None:
    """Label maximal multi-tip clades whose tips share one audited subcluster assignment.

    Args:
        dataset: Valid Auspice v2 dataset, modified in place.

    Notes:
        Non-monophyletic subclusters may occupy several distinctly labelled parts. This
        labels actual inferred branches without modifying topology or inventing ancestry.
        Mixed and unassigned clades remain ordinary tree branches.
    """
    root = dataset["tree"]
    stack = [(root, False)]
    summaries: dict[int, tuple[str | None, int]] = {}
    while stack:
        node, visited = stack.pop()
        children = node.get("children", [])
        if children and not visited:
            stack.append((node, True))
            stack.extend((child, False) for child in children)
            continue
        if not children:
            value = node.get("node_attrs", {}).get("subcluster_id", {}).get("value")
            summaries[id(node)] = (str(value) if value else None, 1)
        else:
            labels = {summaries[id(child)][0] for child in children}
            label = next(iter(labels)) if len(labels) == 1 and None not in labels else None
            summaries[id(node)] = (label, sum(summaries[id(child)][1] for child in children))
    pending = [root]
    parts: dict[str, int] = {}
    while pending:
        node = pending.pop()
        label, count = summaries[id(node)]
        if label and count >= 2:
            parts[label] = parts.get(label, 0) + 1
            node.setdefault("branch_attrs", {}).setdefault("labels", {})["subcluster"] = (
                f"{label}__part_{parts[label]}"
            )
        else:
            pending.extend(reversed(node.get("children", [])))
    if parts:
        dataset["meta"]["stream_labels"] = ["subcluster"]
