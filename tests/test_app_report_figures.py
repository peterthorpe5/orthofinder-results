"""Tests for bounded report-data parsing and evolutionary Plotly figures."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from orthofinder_interrogation_app import figures, report_data
from orthofinder_interrogation_app.figures import (
    distance_matrix_array,
    distance_matrix_figure,
    linked_member_ids,
    nearest_neighbour_figure,
    pcoa_figure,
    phylogram_figure,
    shepard_figure,
)
from orthofinder_interrogation_app.report_data import load_visualisation_catalog
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_results.errors import InputValidationError


def _catalog(application_resource: Path):
    """Load the real synthetic offline report through resource validation."""

    resource = open_resource(path=application_resource)
    assert resource.report_path is not None
    return load_visualisation_catalog(
        report_path=resource.report_path, expected_run_id=resource.run_id
    )


def _entry(application_resource: Path) -> dict:
    """Return the fixture's only bounded visual record."""

    return _catalog(application_resource).get(label="HOG|N0|N0.HOG1")


def _write_payload(*, path: Path, payload: object) -> Path:
    """Write one compact report payload for defensive-parser tests."""

    path.write_text(
        '<script id="orthofinder-results-data" type="application/json">'
        + json.dumps(payload, separators=(",", ":"))
        + "</script>",
        encoding="utf-8",
    )
    return path


def test_catalog_identity_labels_bounds_and_missing_group(
    application_resource: Path,
) -> None:
    """A valid report exposes deterministic group labels and requested graph k."""

    catalog = _catalog(application_resource)
    assert catalog.run_id == "test_run"
    assert catalog.labels() == ("HOG|N0|N0.HOG1",)
    assert catalog.nearest_neighbours == 3
    assert len(catalog.get(label=catalog.labels()[0])["members"]) == 3
    with pytest.raises(InputValidationError, match="not available"):
        catalog.get(label="missing")


def test_linked_selection_and_all_evolutionary_figures(
    application_resource: Path,
) -> None:
    """Member/species selections propagate across every requested visual type."""

    entry = _entry(application_resource)
    assert linked_member_ids(
        entry=entry,
        selected_members=("beta_1", "missing"),
        selected_species=("Species_A",),
    ) == frozenset({"alpha_1", "alpha_2", "beta_1"})
    selected = linked_member_ids(entry=entry, selected_members=("alpha_1",))
    pcoa = pcoa_figure(entry=entry, selected_members=selected)
    assert pcoa.layout.title.text.startswith("Complete-distance PCoA")
    assert len(pcoa.data) == 2
    shepard = shepard_figure(entry=entry)
    assert len(shepard.data) == 2
    assert shepard.layout.xaxis.title.text.startswith("Exact input")
    phylogram = phylogram_figure(entry=entry, selected_members=selected, show_all_labels=True)
    assert len(phylogram.data) == 3
    assert list(phylogram.data[-1].text) == ["alpha_1", "alpha_2", "beta_1"]
    matrix = distance_matrix_figure(entry=entry, selected_members=("alpha_1", "beta_1"))
    assert np.asarray(matrix.data[0].z).shape == (2, 2)
    topology = nearest_neighbour_figure(
        entry=entry, selected_members=selected, include_connectors=False
    )
    topology_with_connectors = nearest_neighbour_figure(entry=entry, include_connectors=True)
    assert len(topology.data) == 2
    assert len(topology_with_connectors.data) == 3


def test_exact_distance_matrix_reconstruction_and_validation(
    application_resource: Path,
) -> None:
    """The compact upper triangle reconstructs exact symmetric displayed values."""

    entry = _entry(application_resource)
    labels, matrix = distance_matrix_array(entry=entry)
    assert labels == ("alpha_1", "alpha_2", "beta_1")
    assert matrix.tolist() == [[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]]
    with pytest.raises(InputValidationError, match="absent"):
        distance_matrix_figure(entry=entry, selected_members=("missing",))
    malformed = copy.deepcopy(entry)
    malformed["distanceMatrix"]["upperTriangle"] = [0.1]
    with pytest.raises(InputValidationError, match="requires 3"):
        distance_matrix_array(entry=malformed)
    malformed["distanceMatrix"]["status"] = "UNAVAILABLE"
    with pytest.raises(InputValidationError, match="no exact"):
        distance_matrix_array(entry=malformed)


@pytest.mark.parametrize(
    ("mutation", "call", "message"),
    [
        (
            lambda entry: entry.update(nodes="wrong"),
            lambda entry: pcoa_figure(entry=entry),
            "node array",
        ),
        (
            lambda entry: entry["nodes"][0].pop("projectionX"),
            lambda entry: pcoa_figure(entry=entry),
            "PCoA coordinates",
        ),
        (
            lambda entry: entry["distanceProjection"].update(shepard_points=[]),
            lambda entry: shepard_figure(entry=entry),
            "no Shepard",
        ),
        (
            lambda entry: entry["distanceProjection"].update(shepard_points=[[0.1]]),
            lambda entry: shepard_figure(entry=entry),
            "malformed point",
        ),
        (
            lambda entry: entry["distanceProjection"].update(shepard_points=[[0.1, float("nan")]]),
            lambda entry: shepard_figure(entry=entry),
            "non-finite",
        ),
        (
            lambda entry: entry["phylogram"].update(status="UNAVAILABLE"),
            lambda entry: phylogram_figure(entry=entry),
            "no usable",
        ),
        (
            lambda entry: entry["phylogram"].update(nodes="wrong"),
            lambda entry: phylogram_figure(entry=entry),
            "lacks node",
        ),
        (
            lambda entry: entry.update(edges="wrong"),
            lambda entry: nearest_neighbour_figure(entry=entry),
            "edge array",
        ),
        (
            lambda entry: entry.update(nodes=[]),
            lambda entry: nearest_neighbour_figure(entry=entry),
            "no nodes",
        ),
        (
            lambda entry: entry.update(members="wrong"),
            lambda entry: linked_member_ids(entry=entry),
            "member array",
        ),
        (
            lambda entry: entry.update(members=["wrong"]),
            lambda entry: linked_member_ids(entry=entry),
            "malformed member",
        ),
        (
            lambda entry: entry.update(
                members=[
                    {"member_id": "same", "species_label": "a"},
                    {"member_id": "same", "species_label": "b"},
                ]
            ),
            lambda entry: linked_member_ids(entry=entry),
            "duplicate/empty",
        ),
    ],
)
def test_figure_payload_errors_are_controlled(
    application_resource: Path, mutation, call, message: str
) -> None:
    """Malformed visual subsets produce user-facing errors instead of blank charts."""

    entry = copy.deepcopy(_entry(application_resource))
    mutation(entry)
    with pytest.raises(InputValidationError, match=message):
        call(entry)


def test_report_parser_rejects_identity_json_structure_and_entry_errors(
    tmp_path: Path, application_resource: Path
) -> None:
    """Only a matching, bounded and structurally valid report payload is accepted."""

    entry = _entry(application_resource)
    with pytest.raises(InputValidationError, match="unavailable"):
        load_visualisation_catalog(
            report_path=tmp_path / "missing.html", expected_run_id="test_run"
        )
    no_marker = tmp_path / "no_marker.html"
    no_marker.write_text("not a report")
    with pytest.raises(InputValidationError, match="lacks"):
        load_visualisation_catalog(report_path=no_marker, expected_run_id="test_run")
    invalid_json = tmp_path / "invalid.html"
    invalid_json.write_text(
        '<script id="orthofinder-results-data" type="application/json">{bad}</script>'
    )
    with pytest.raises(InputValidationError, match="invalid JSON"):
        load_visualisation_catalog(report_path=invalid_json, expected_run_id="test_run")
    mismatch = _write_payload(
        path=tmp_path / "mismatch.html",
        payload={"run": {"run_id": "other"}, "networks": {}, "limits": {}},
    )
    with pytest.raises(InputValidationError, match="disagree"):
        load_visualisation_catalog(report_path=mismatch, expected_run_id="test_run")
    bad_entry = _write_payload(
        path=tmp_path / "bad_entry.html",
        payload={
            "run": {"run_id": "test_run"},
            "networks": {"group": {**entry, "members": "bad"}},
            "limits": {},
        },
    )
    with pytest.raises(InputValidationError, match="lacks member"):
        load_visualisation_catalog(report_path=bad_entry, expected_run_id="test_run")
    unsafe_limit = _write_payload(
        path=tmp_path / "limit.html",
        payload={
            "run": {"run_id": "test_run"},
            "networks": {},
            "limits": {"nearestNeighbours": 101},
        },
    )
    with pytest.raises(InputValidationError, match="unsafe"):
        load_visualisation_catalog(report_path=unsafe_limit, expected_run_id="test_run")


def test_additional_matrix_topology_and_colour_defences(
    application_resource: Path,
) -> None:
    """Matrix metadata, irrelevant edges and report colours are treated defensively."""

    entry = _entry(application_resource)
    malformed = copy.deepcopy(entry)
    malformed["distanceMatrix"]["memberOrder"] = "wrong"
    with pytest.raises(InputValidationError, match="malformed"):
        distance_matrix_array(entry=malformed)
    duplicate = copy.deepcopy(entry)
    duplicate["distanceMatrix"]["memberOrder"] = ["same", "same", "third"]
    with pytest.raises(InputValidationError, match="order is invalid"):
        distance_matrix_array(entry=duplicate)
    negative = copy.deepcopy(entry)
    negative["distanceMatrix"]["upperTriangle"][0] = -1
    with pytest.raises(InputValidationError, match="invalid distance"):
        distance_matrix_array(entry=negative)
    no_projection = copy.deepcopy(entry)
    no_projection["distanceProjection"] = None
    with pytest.raises(InputValidationError, match="distance projection"):
        pcoa_figure(entry=no_projection)
    odd_edges = copy.deepcopy(entry)
    odd_edges["edges"] = [
        "malformed",
        {"edgeType": "OTHER", "from": "alpha_1", "to": "alpha_2"},
        {"edgeType": "NEAREST_NEIGHBOUR", "from": "missing", "to": "alpha_2"},
    ]
    assert len(nearest_neighbour_figure(entry=odd_edges).data) == 1
    assert figures._safe_colour("red; background:url(x)") == "#5b6475"


def test_report_payload_defensive_collection_bounds(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Visual group/member/edge limits and duplicate identities are enforced."""

    entry = _entry(application_resource)
    with pytest.raises(InputValidationError, match="JSON object"):
        report_data._validate_payload(payload=[], expected_run_id="test_run")
    with pytest.raises(InputValidationError, match="lacks run"):
        report_data._validate_payload(payload={}, expected_run_id="test_run")
    with pytest.raises(InputValidationError, match="label and object"):
        report_data._validate_payload(
            payload={
                "run": {"run_id": "test_run"},
                "networks": {"": entry},
                "limits": {},
            },
            expected_run_id="test_run",
        )
    monkeypatch.setattr(report_data, "MAX_VISUAL_GROUPS", 0)
    with pytest.raises(InputValidationError, match="visual groups"):
        report_data._validate_payload(
            payload={
                "run": {"run_id": "test_run"},
                "networks": {"group": entry},
                "limits": {},
            },
            expected_run_id="test_run",
        )
    monkeypatch.setattr(report_data, "MAX_VISUAL_GROUPS", 100)
    duplicate = copy.deepcopy(entry)
    duplicate["members"][1]["member_id"] = duplicate["members"][0]["member_id"]
    with pytest.raises(InputValidationError, match="duplicate or malformed"):
        report_data._validate_entry(label="group", entry=duplicate)
    monkeypatch.setattr(report_data, "MAX_VISUAL_MEMBERS", 2)
    with pytest.raises(InputValidationError, match="member safety"):
        report_data._validate_entry(label="group", entry=entry)
    monkeypatch.setattr(report_data, "MAX_VISUAL_MEMBERS", 1_000)
    monkeypatch.setattr(report_data, "MAX_VISUAL_EDGES", 1)
    with pytest.raises(InputValidationError, match="edge safety"):
        report_data._validate_entry(label="group", entry=entry)


def test_report_unterminated_payload_and_nonnumeric_neighbour_limit(
    tmp_path: Path,
) -> None:
    """Truncated scripts fail, while malformed optional k safely becomes unknown."""

    truncated = tmp_path / "truncated.html"
    truncated.write_text('<script id="orthofinder-results-data" type="application/json">{}')
    with pytest.raises(InputValidationError, match="not terminated"):
        load_visualisation_catalog(report_path=truncated, expected_run_id="test_run")
    catalog = report_data._validate_payload(
        payload={
            "run": {"run_id": "test_run"},
            "networks": {},
            "limits": {"nearestNeighbours": "unknown"},
        },
        expected_run_id="test_run",
    )
    assert catalog.nearest_neighbours == 0
