"""Plain-language labels for machine-readable dispersion benchmark codes."""

_EXACT_PROFILE_LABELS = {
    "E3_ALL": "All E3 clusters",
    "HOUSEKEEPING_ALL": "All housekeeping-reference clusters",
    "R_NLR_ALL": "All R/NLR-reference clusters",
    "MATCHED_NON_FOCUS": "Matched non-focus controls",
}
_PROFILE_PREFIX_LABELS = {
    "E3_CATEGORY::": "E3 category — ",
    "HOUSEKEEPING_SUBCLASS::": "Housekeeping panel — ",
    "R_NLR_SUBCLASS::": "R/NLR subclass — ",
}
_CLASS_LABELS = {
    "E3": "E3",
    "HOUSEKEEPING": "Housekeeping reference",
    "R_NLR": "R/NLR reference",
    "NON_FOCUS_CONTROL": "Matched non-focus control",
}


def benchmark_profile_label(*, profile_id: str) -> str:
    """Return a readable label while preserving unfamiliar profile information.

    Args:
        profile_id: Exact machine-readable benchmark profile identifier.

    Returns:
        Stable user-facing profile label.
    """

    exact = str(profile_id).strip()
    if exact in _EXACT_PROFILE_LABELS:
        return _EXACT_PROFILE_LABELS[exact]
    for prefix, label in _PROFILE_PREFIX_LABELS.items():
        if exact.startswith(prefix):
            return label + _readable_code(value=exact[len(prefix):])
    return _readable_code(value=exact)


def benchmark_class_label(*, profile_classes: str) -> str:
    """Return a readable label for one or several benchmark classes.

    Args:
        profile_classes: Semicolon-delimited exact benchmark classes.

    Returns:
        A single label or an explicit combined-profile label.
    """

    classes = tuple(
        dict.fromkeys(value.strip() for value in str(profile_classes).split(";") if value.strip())
    )
    if not classes:
        return "Unlabelled profile"
    labels = tuple(_CLASS_LABELS.get(value, _readable_code(value=value)) for value in classes)
    return labels[0] if len(labels) == 1 else "Mixed: " + " + ".join(labels)


def _readable_code(*, value: str) -> str:
    """Make an unfamiliar controlled code readable without discarding its text."""

    text = str(value).replace("::", " — ").replace("_", " ").strip()
    return text or "Unlabelled profile"
