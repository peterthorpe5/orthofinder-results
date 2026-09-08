"""Exact distance geometry and within-group dispersion summaries."""

from __future__ import annotations

import math
import statistics
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from orthofinder_results.errors import InputValidationError


@dataclass(frozen=True)
class PcoaGeometry:
    """Classical PCoA coordinates with explicit preservation diagnostics.

    Attributes:
        member_ids: Matrix and coordinate row order.
        species_labels: Species aligned to ``member_ids``.
        coordinates: Up to three positive-axis coordinates per member.
        axis_fractions: Positive-inertia fractions for the displayed axes.
        cumulative_fraction: Positive inertia retained by all displayed axes.
        negative_inertia_fraction: Absolute negative share of total absolute inertia.
        stress_2d: Normalised raw stress in the first two axes.
        stress_3d: Normalised raw stress in the first three axes.
        correlation_2d: Pearson correlation of exact and 2D projected distances.
        correlation_3d: Pearson correlation of exact and 3D projected distances.
    """

    member_ids: tuple[str, ...]
    species_labels: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    axis_fractions: tuple[float, float, float]
    cumulative_fraction: float
    negative_inertia_fraction: float
    stress_2d: float
    stress_3d: float
    correlation_2d: float | None
    correlation_3d: float | None


def exact_distance_matrix(
    *,
    rows: Sequence[Mapping[str, Any]],
    members: Sequence[Mapping[str, Any]],
) -> tuple[tuple[str, ...], tuple[str, ...], np.ndarray]:
    """Construct and validate one complete symmetric distance matrix.

    Args:
        rows: Long-form pairwise distance rows.
        members: Displayed member-to-species records.

    Returns:
        Lexical member order, aligned species labels and a symmetric matrix.

    Raises:
        InputValidationError: If identifiers, values or pair coverage are invalid.
    """

    species_by_member: dict[str, str] = {}
    for row in members:
        member_id = str(row.get("member_id", ""))
        species = str(row.get("species_label", ""))
        if not member_id or not species:
            raise InputValidationError("Every displayed member requires an identifier and species.")
        previous = species_by_member.get(member_id)
        if previous is not None and previous != species:
            raise InputValidationError(f"Member {member_id!r} maps to multiple species.")
        species_by_member[member_id] = species
    member_ids = tuple(sorted(species_by_member))
    if len(member_ids) < 2:
        raise InputValidationError("At least two unique members are required for dispersion.")
    index = {member_id: position for position, member_id in enumerate(member_ids)}
    matrix = np.zeros((len(member_ids), len(member_ids)), dtype=float)
    observed: set[tuple[str, str]] = set()
    for row in rows:
        left, right = str(row.get("member_a", "")), str(row.get("member_b", ""))
        if left == right or left not in index or right not in index:
            raise InputValidationError("A distance row has an invalid displayed endpoint.")
        try:
            distance = float(row.get("distance"))
        except (TypeError, ValueError) as error:
            raise InputValidationError("A displayed distance is not numeric.") from error
        if not math.isfinite(distance) or distance < 0:
            raise InputValidationError("A displayed distance is non-finite or negative.")
        pair = tuple(sorted((left, right)))
        if pair in observed:
            raise InputValidationError(f"Duplicate displayed distance pair: {pair!r}.")
        observed.add(pair)
        left_index, right_index = index[left], index[right]
        matrix[left_index, right_index] = distance
        matrix[right_index, left_index] = distance
    expected = len(member_ids) * (len(member_ids) - 1) // 2
    if len(observed) != expected:
        raise InputValidationError(
            f"Complete dispersion requires {expected:,} pairs; observed {len(observed):,}."
        )
    species = tuple(species_by_member[member_id] for member_id in member_ids)
    return member_ids, species, matrix


def classical_pcoa(
    *,
    rows: Sequence[Mapping[str, Any]],
    members: Sequence[Mapping[str, Any]],
) -> PcoaGeometry:
    """Calculate deterministic three-axis classical PCoA and fit diagnostics.

    Args:
        rows: Complete exact pairwise distance rows.
        members: Displayed member-to-species records.

    Returns:
        Three-axis coordinates and separate 2D/3D preservation measures.

    Raises:
        InputValidationError: If the matrix has no positive PCoA axis.
    """

    member_ids, species, matrix = exact_distance_matrix(rows=rows, members=members)
    squared = matrix**2
    row_means = squared.mean(axis=1)
    centred = -0.5 * (
        squared - row_means[:, None] - row_means[None, :] + squared.mean()
    )
    try:
        eigenvalues, eigenvectors = np.linalg.eigh(centred)
    except np.linalg.LinAlgError as error:
        raise InputValidationError(f"PCoA eigendecomposition failed: {error}") from error
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    tolerance = max(float(np.max(np.abs(eigenvalues))), 1.0) * 1e-12
    positive_indices = np.flatnonzero(eigenvalues > tolerance)
    if not positive_indices.size:
        raise InputValidationError("The exact matrix produced no positive PCoA axes.")
    coordinates = np.zeros((len(member_ids), 3), dtype=float)
    positive_sum = float(np.sum(eigenvalues[eigenvalues > tolerance]))
    fractions = [0.0, 0.0, 0.0]
    for axis, eigen_index in enumerate(positive_indices[:3]):
        coordinate = eigenvectors[:, eigen_index] * math.sqrt(float(eigenvalues[eigen_index]))
        anchor = int(np.argmax(np.abs(coordinate)))
        if coordinate[anchor] < 0:
            coordinate *= -1
        coordinates[:, axis] = coordinate
        fractions[axis] = float(eigenvalues[eigen_index] / positive_sum)
    upper = np.triu_indices(len(member_ids), k=1)
    exact = matrix[upper]
    projected_2d = _projected_distances(coordinates=coordinates[:, :2], upper=upper)
    projected_3d = _projected_distances(coordinates=coordinates, upper=upper)
    absolute_sum = float(np.sum(np.abs(eigenvalues)))
    negative_fraction = (
        float(np.sum(np.abs(eigenvalues[eigenvalues < -tolerance])) / absolute_sum)
        if absolute_sum > 0
        else 0.0
    )
    return PcoaGeometry(
        member_ids=member_ids,
        species_labels=species,
        coordinates=tuple(tuple(float(value) for value in row) for row in coordinates),
        axis_fractions=tuple(fractions),
        cumulative_fraction=sum(fractions),
        negative_inertia_fraction=negative_fraction,
        stress_2d=_normalised_stress(exact=exact, projected=projected_2d),
        stress_3d=_normalised_stress(exact=exact, projected=projected_3d),
        correlation_2d=_distance_correlation(exact=exact, projected=projected_2d),
        correlation_3d=_distance_correlation(exact=exact, projected=projected_3d),
    )


def member_dispersion_rows(
    *,
    rows: Sequence[Mapping[str, Any]],
    members: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    """Summarise every displayed member's distance to all other members.

    Args:
        rows: Complete pairwise distance rows.
        members: Displayed member-to-species records.

    Returns:
        Rows sorted from the sample medoid towards the most peripheral member.
    """

    member_ids, species, matrix = exact_distance_matrix(rows=rows, members=members)
    output = []
    for index, member_id in enumerate(member_ids):
        values = np.delete(matrix[index], index)
        nearest_index = min(
            (other for other in range(len(member_ids)) if other != index),
            key=lambda other: (matrix[index, other], member_ids[other]),
        )
        output.append(
            {
                "member_id": member_id,
                "species_label": species[index],
                "mean_distance": float(np.mean(values)),
                "median_distance": float(np.median(values)),
                "population_stddev_distance": float(np.std(values)),
                "minimum_distance": float(np.min(values)),
                "maximum_distance": float(np.max(values)),
                "nearest_member_id": member_ids[nearest_index],
                "nearest_species_label": species[nearest_index],
                "nearest_distance": float(matrix[index, nearest_index]),
            }
        )
    medoid_id = min(
        output,
        key=lambda row: (float(row["mean_distance"]), str(row["member_id"])),
    )["member_id"]
    for row in output:
        row["is_sample_medoid"] = row["member_id"] == medoid_id
    output.sort(key=lambda row: (float(row["mean_distance"]), str(row["member_id"])))
    return tuple(output)


def distance_class_rows(
    *, rows: Sequence[Mapping[str, Any]]
) -> tuple[dict[str, Any], ...]:
    """Label every pair as within-species or between-species for plotting.

    Args:
        rows: Species-decorated distance records.

    Returns:
        Concise distance, class and endpoint records.
    """

    output = []
    for row in rows:
        species_a, species_b = str(row.get("species_a", "")), str(row.get("species_b", ""))
        if not species_a or not species_b:
            raise InputValidationError("Distance-class analysis requires endpoint species.")
        output.append(
            {
                "member_a": str(row["member_a"]),
                "member_b": str(row["member_b"]),
                "species_a": species_a,
                "species_b": species_b,
                "distance": float(row["distance"]),
                "pair_class": (
                    "Within species" if species_a == species_b else "Between species"
                ),
            }
        )
    return tuple(output)


def species_dispersion_rows(
    *, rows: Sequence[Mapping[str, Any]]
) -> tuple[dict[str, Any], ...]:
    """Summarise pair-distance distributions by endpoint species pair.

    Args:
        rows: Species-decorated exact distance rows.

    Returns:
        Deterministic species-pair count, mean, median, SD and range records.
    """

    values_by_pair: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        species_a, species_b = str(row.get("species_a", "")), str(row.get("species_b", ""))
        if not species_a or not species_b:
            raise InputValidationError("Species dispersion requires endpoint species.")
        pair = tuple(sorted((species_a, species_b)))
        values_by_pair[pair].append(float(row["distance"]))
    output = []
    for (species_a, species_b), values in sorted(values_by_pair.items()):
        output.append(
            {
                "species_a": species_a,
                "species_b": species_b,
                "pair_class": (
                    "Within species" if species_a == species_b else "Between species"
                ),
                "pair_count": len(values),
                "mean_distance": statistics.fmean(values),
                "median_distance": statistics.median(values),
                "population_stddev_distance": statistics.pstdev(values),
                "minimum_distance": min(values),
                "maximum_distance": max(values),
            }
        )
    return tuple(output)


def _projected_distances(
    *, coordinates: np.ndarray, upper: tuple[np.ndarray, np.ndarray]
) -> np.ndarray:
    """Return upper-triangle Euclidean distances for projected coordinates."""

    differences = coordinates[:, None, :] - coordinates[None, :, :]
    return np.sqrt(np.sum(differences**2, axis=2))[upper]


def _normalised_stress(*, exact: np.ndarray, projected: np.ndarray) -> float:
    """Return normalised raw stress for one projected distance vector."""

    denominator = float(np.sum(exact**2))
    return (
        math.sqrt(float(np.sum((exact - projected) ** 2)) / denominator)
        if denominator > 0
        else 0.0
    )


def _distance_correlation(*, exact: np.ndarray, projected: np.ndarray) -> float | None:
    """Return Pearson preservation correlation when it is mathematically defined."""

    if exact.size <= 1 or np.std(exact) == 0 or np.std(projected) == 0:
        return None
    return float(np.corrcoef(exact, projected)[0, 1])
