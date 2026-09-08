"""Immutable models for application queries and resource identity."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from orthofinder_results.errors import InputValidationError

INCLUDE_MODES = frozenset({"ANY", "ALL", "EXACT_SET"})
DISTANCE_AVAILABILITY_MODES = frozenset({"ANY", "CALCULATED", "NOT_CALCULATED"})
SORT_MODES = frozenset(
    {
        "MEMBER_COUNT_DESC",
        "SPECIES_COUNT_DESC",
        "GROUP_ID_ASC",
        "MEAN_DISTANCE_ASC",
        "MEAN_DISTANCE_DESC",
        "DISTANCE_SD_ASC",
    }
)
TAXONOMY_SEARCH_MODES = frozenset({"CONTAINS", "ENRICHED", "SAMPLED_EXCLUSIVE", "NEAR_EXCLUSIVE"})


@dataclass(frozen=True)
class ResourceIdentity:
    """Identity and capabilities of one immutable resource.

    Attributes:
        resource_path: Completed resource directory or direct DuckDB path.
        database_path: Physical DuckDB opened by the application.
        report_path: Optional bundled offline report.
        run_id: Immutable run identifier.
        schema_version: Resource schema version.
        resource_package_version: Package version that created the resource.
        orthofinder_version: Source OrthoFinder version.
        adapter_name: Version-specific parser adapter.
        primary_group_authority: Primary HOG or legacy group authority.
        counts: Manifest row and capability counts.
        relations: Exact DuckDB relations available to capability-driven pages.
    """

    resource_path: Path
    database_path: Path
    report_path: Path | None
    run_id: str
    schema_version: int
    resource_package_version: str
    orthofinder_version: str
    adapter_name: str
    primary_group_authority: str
    counts: dict[str, int] = field(default_factory=dict)
    relations: frozenset[str] = field(default_factory=frozenset)


@dataclass(frozen=True)
class GroupKey:
    """Composite identifier for one run-scoped OrthoFinder group."""

    run_id: str
    group_type: str
    hierarchy_node: str
    group_id: str

    def display_label(self) -> str:
        """Return an unambiguous user-facing label.

        Returns:
            Composite group label suitable for a selector.
        """

        node = self.hierarchy_node or "ROOT"
        return f"{self.group_type} | {node} | {self.group_id}"


@dataclass(frozen=True)
class GroupSearchFilters:
    """Validated, bounded filters for lazy group searching."""

    group_type: str = ""
    hierarchy_node: str | None = None
    group_id_contains: str = ""
    member_id_contains: str = ""
    included_species: tuple[str, ...] = ()
    include_mode: str = "ALL"
    excluded_species: tuple[str, ...] = ()
    minimum_member_count: int | None = None
    maximum_member_count: int | None = None
    minimum_species_count: int | None = None
    maximum_species_count: int | None = None
    distance_availability: str = "ANY"
    maximum_mean_distance: float | None = None
    maximum_distance_sd: float | None = None
    sort_mode: str = "MEMBER_COUNT_DESC"
    page_size: int = 100
    page_number: int = 1

    def __post_init__(self) -> None:
        """Reject contradictory or unsafe search controls."""

        if self.include_mode not in INCLUDE_MODES:
            raise InputValidationError(f"Unsupported include mode: {self.include_mode}")
        if self.distance_availability not in DISTANCE_AVAILABILITY_MODES:
            raise InputValidationError(
                f"Unsupported distance availability: {self.distance_availability}"
            )
        if self.sort_mode not in SORT_MODES:
            raise InputValidationError(f"Unsupported sort mode: {self.sort_mode}")
        if not 1 <= self.page_size <= 500:
            raise InputValidationError("page_size must be between 1 and 500.")
        if self.page_number < 1:
            raise InputValidationError("page_number must be at least one.")
        _validate_optional_range(
            minimum=self.minimum_member_count,
            maximum=self.maximum_member_count,
            label="member count",
        )
        _validate_optional_range(
            minimum=self.minimum_species_count,
            maximum=self.maximum_species_count,
            label="species count",
        )
        for label, value in (
            ("maximum_mean_distance", self.maximum_mean_distance),
            ("maximum_distance_sd", self.maximum_distance_sd),
        ):
            if value is not None and value < 0:
                raise InputValidationError(f"{label} must not be negative.")
        included = _normalise_species(values=self.included_species, label="included")
        excluded = _normalise_species(values=self.excluded_species, label="excluded")
        overlap = sorted(set(included) & set(excluded))
        if overlap:
            raise InputValidationError(
                "Species cannot be both required and rejected: " + "; ".join(overlap)
            )
        object.__setattr__(self, "included_species", included)
        object.__setattr__(self, "excluded_species", excluded)

    @property
    def offset(self) -> int:
        """Return the zero-based SQL offset for the requested page."""

        return (self.page_number - 1) * self.page_size


@dataclass(frozen=True)
class SearchPage:
    """One bounded page of group-search results."""

    rows: tuple[dict[str, Any], ...]
    total_rows: int
    page_number: int
    page_size: int


@dataclass(frozen=True)
class TaxonomySearchFilters:
    """Validated controls for one descendant-aware taxonomic search."""

    group_type: str
    hierarchy_node: str
    target_taxon_id: int
    mode: str = "CONTAINS"
    minimum_target_species_count: int = 1
    minimum_target_coverage: float = 0.0
    minimum_mapped_purity: float = 0.0
    maximum_outside_species_count: int = 1
    maximum_unresolved_species_count: int = 0
    maximum_q_value: float = 0.05
    minimum_odds_ratio: float = 1.0
    page_size: int = 100
    page_number: int = 1

    def __post_init__(self) -> None:
        """Reject ambiguous authorities and unsafe statistical controls."""

        if not self.group_type.strip():
            raise InputValidationError("Taxonomy search requires one exact group type.")
        if self.mode not in TAXONOMY_SEARCH_MODES:
            raise InputValidationError(f"Unsupported taxonomy search mode: {self.mode}")
        if self.target_taxon_id <= 0:
            raise InputValidationError("Target NCBI taxon ID must be positive.")
        for label, value in (
            ("minimum_target_species_count", self.minimum_target_species_count),
            ("maximum_outside_species_count", self.maximum_outside_species_count),
            ("maximum_unresolved_species_count", self.maximum_unresolved_species_count),
        ):
            if value < 0:
                raise InputValidationError(f"{label} must not be negative.")
        for label, value in (
            ("minimum_target_coverage", self.minimum_target_coverage),
            ("minimum_mapped_purity", self.minimum_mapped_purity),
            ("maximum_q_value", self.maximum_q_value),
        ):
            if not 0 <= value <= 1:
                raise InputValidationError(f"{label} must be between zero and one.")
        if self.minimum_odds_ratio < 0:
            raise InputValidationError("minimum_odds_ratio must not be negative.")
        if not 1 <= self.page_size <= 500:
            raise InputValidationError("page_size must be between 1 and 500.")
        if self.page_number < 1:
            raise InputValidationError("page_number must be at least one.")

    @property
    def offset(self) -> int:
        """Return the zero-based result offset."""

        return (self.page_number - 1) * self.page_size


@dataclass(frozen=True)
class TaxonomySearchPage:
    """One bounded taxonomic-search page and its statistical universe."""

    rows: tuple[dict[str, Any], ...]
    total_rows: int
    tested_group_count: int
    target_species_count: int
    outside_species_count: int
    unresolved_species_count: int
    page_number: int
    page_size: int


def _normalise_species(*, values: tuple[str, ...], label: str) -> tuple[str, ...]:
    """Strip and deterministically deduplicate selected species labels.

    Args:
        values: Raw selected labels.
        label: Selection role for validation errors.

    Returns:
        Sorted unique labels.

    Raises:
        InputValidationError: If a label is empty.
    """

    stripped = tuple(value.strip() for value in values)
    if any(not value for value in stripped):
        raise InputValidationError(f"The {label} species selection contains an empty label.")
    return tuple(sorted(set(stripped)))


def _validate_optional_range(*, minimum: int | None, maximum: int | None, label: str) -> None:
    """Validate an optional non-negative inclusive integer range.

    Args:
        minimum: Optional lower bound.
        maximum: Optional upper bound.
        label: User-facing range name.

    Raises:
        InputValidationError: If either bound is negative or the range is reversed.
    """

    if minimum is not None and minimum < 0:
        raise InputValidationError(f"Minimum {label} must not be negative.")
    if maximum is not None and maximum < 0:
        raise InputValidationError(f"Maximum {label} must not be negative.")
    if minimum is not None and maximum is not None and minimum > maximum:
        raise InputValidationError(f"Minimum {label} must not exceed maximum {label}.")
