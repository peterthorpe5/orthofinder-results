"""Exact OrthoFinder group selection and optional sequence transfer to the workbench."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import duckdb
import pandas as pd
import pyarrow as pa

from .models import GroupKey, GroupSearchFilters
from .queries import MAX_GROUP_MEMBER_ROWS, OrthoFinderQueryService
from .research.common import ResearchError, bounded_number, normalise_protein
from .resource import connect_read_only
from .terminal_motif import validate_sequence_sidecar

LOGGER = logging.getLogger(__name__)
INPUT_COLUMNS = (
    "input_id",
    "species_label",
    "source_run_id",
    "source_group_type",
    "source_hierarchy_node",
    "source_group_id",
    "selection_rank",
    "source_member_count",
)


def select_group_inputs(
    *,
    service: OrthoFinderQueryService,
    group_type: str,
    hierarchy_node: str,
    group_ids: Sequence[str] = (),
    maximum_groups: int = 200,
    included_species: Sequence[str] = (),
    require_all_species: bool = False,
    member_id_contains: str = "",
    maximum_inputs: int = 10_000,
) -> pd.DataFrame:
    """Select complete group memberships by exact IDs or decreasing local group size.

    Args:
        service: Read-only service for a validated run.
        group_type: HOG or LEGACY_ORTHOGROUP authority.
        hierarchy_node: Exact species-tree level; empty string selects the flat ROOT set.
        group_ids: Exact group IDs overriding top-N selection.
        maximum_groups: Maximum groups in the bounded selection, up to 500.
        included_species: Exact stored species labels.
        require_all_species: Require every supplied species rather than any one.
        member_id_contains: Optional literal protein-ID substring used to find groups.
        maximum_inputs: Maximum distinct input IDs; larger requests fail in full.

    Returns:
        Complete member contexts with exact group keys and local size-selection rank.

    Raises:
        ResearchError: If IDs are missing, input bounds are exceeded or membership is incomplete.
    """
    maximum = bounded_number(
        value=maximum_groups, label="Selected groups", minimum=1, maximum=500, integer=True
    )
    bound = bounded_number(
        value=maximum_inputs, label="Selected proteins", minimum=1, maximum=10000, integer=True
    )
    if isinstance(group_ids, (str, bytes)) or any(
        not isinstance(value, str) or not value.strip() for value in group_ids
    ):
        raise ResearchError("Supply exact group IDs as a list")
    identifiers = tuple(dict.fromkeys(group_ids))
    if identifiers:
        if len(identifiers) > maximum:
            raise ResearchError("Explicit group selection exceeds the selected group bound")
        groups = [
            service.get_group(
                key=GroupKey(
                    run_id=service.resource.run_id,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node,
                    group_id=identifier,
                )
            )
            for identifier in identifiers
        ]
        groups.sort(key=lambda row: (-int(row["member_count"]), str(row["group_id"])))
    else:
        groups = service.search_groups(
            filters=GroupSearchFilters(
                group_type=group_type,
                hierarchy_node=hierarchy_node,
                included_species=tuple(included_species),
                include_mode="ALL" if require_all_species else "ANY",
                member_id_contains=member_id_contains,
                page_size=maximum,
                sort_mode="MEMBER_COUNT_DESC",
            )
        ).rows
    records = []
    inputs: set[str] = set()
    for rank, group in enumerate(groups, start=1):
        count = int(group["member_count"])
        if count > MAX_GROUP_MEMBER_ROWS or count > bound:
            raise ResearchError(
                "A selected group exceeds the complete input bound; use a smaller selection "
                "or an explicit FASTA export. No members have been silently omitted."
            )
        key = GroupKey(
            **{
                name: str(group[name])
                for name in ("run_id", "group_type", "hierarchy_node", "group_id")
            }
        )
        members = service.get_group_members(key=key, maximum=MAX_GROUP_MEMBER_ROWS)
        if len(members) != count:
            raise ResearchError("Stored group membership count disagrees with group statistics")
        for member in members:
            inputs.add(str(member["member_id"]))
            if len(inputs) > bound:
                raise ResearchError(
                    "Selection exceeds 10000 distinct proteins; select fewer groups"
                )
            records.append(
                {
                    "input_id": str(member["member_id"]),
                    "species_label": str(member["species_label"]),
                    "source_run_id": key.run_id,
                    "source_group_type": key.group_type,
                    "source_hierarchy_node": key.hierarchy_node,
                    "source_group_id": key.group_id,
                    "selection_rank": rank,
                    "source_member_count": count,
                }
            )
    LOGGER.info("Selected %d complete local groups and %d protein IDs", len(groups), len(inputs))
    return pd.DataFrame(data=records, columns=INPUT_COLUMNS)


def selected_group_fasta(
    *, service: OrthoFinderQueryService, inputs: pd.DataFrame, sidecar: Path | None = None
) -> str:
    """Export every selected protein using exact member and species joins.

    Args:
        service: Validated immutable OrthoFinder resource.
        inputs: Complete member contexts from select_group_inputs.
        sidecar: Optional installed protein-sequence Parquet from the existing motif workflow.

    Returns:
        Validated complete FASTA with unchanged canonical identifiers.

    Raises:
        ResearchError: If sequences are unavailable, ambiguous, missing or too large for the UI.
    """
    if inputs.empty or not {"input_id", "species_label"}.issubset(inputs.columns):
        raise ResearchError("Prepare a non-empty group membership selection first")
    requested = inputs[["input_id", "species_label"]].drop_duplicates()
    if requested["input_id"].duplicated().any():
        raise ResearchError("Protein IDs repeat across species; supply a uniquely named FASTA")
    if sidecar is None and not service.has_relation(relation="protein_sequences"):
        raise ResearchError(
            "Protein sequences are unavailable. Configure --terminal-motif-parquet or upload FASTA."
        )
    source = validate_sequence_sidecar(path=sidecar) if sidecar else None
    try:
        with connect_read_only(database_path=service.resource.database_path) as connection:
            connection.register(
                view_name="research_requested",
                python_object=pa.Table.from_pandas(df=requested, preserve_index=False),
            )
            relation = "read_parquet(?)" if source else "protein_sequences"
            rows = connection.execute(
                "SELECT r.input_id,r.species_label,s.sequence FROM research_requested r "
                f"LEFT JOIN {relation} s ON r.input_id=s.member_id "
                "AND r.species_label=s.species_label ORDER BY r.input_id",
                [str(source)] if source else [],
            ).fetchall()
    except duckdb.Error as error:
        raise ResearchError(
            "Could not join the selected proteins to their sequence authority"
        ) from error
    if len(rows) != len(requested) or any(row[2] is None for row in rows):
        raise ResearchError(
            "Missing or duplicated sequence evidence; complete FASTA export refused"
        )
    records = []
    byte_count = 0
    for identifier, _species, sequence in rows:
        record = f">{identifier}\n{normalise_protein(sequence=str(sequence))}\n"
        byte_count += len(record.encode("utf-8"))
        if byte_count > 20_000_000:
            raise ResearchError("Selected FASTA exceeds the web limit; use the large-file CLI")
        records.append(record)
    return "".join(records)
