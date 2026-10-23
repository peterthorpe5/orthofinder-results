"""Reuse server database settings with an independent OrthoFinder job workspace."""

from __future__ import annotations

import json
import logging
from dataclasses import fields
from pathlib import Path

from .common import ResearchConfig, ResearchError, load_config

LOGGER = logging.getLogger(__name__)


def copy_profile(*, source: Path, output: Path, workspace: Path) -> Path:
    """Create a profile sharing data and tool settings but using a different queue.

    Args:
        source: Existing validated E3 or OrthoFinder research TOML profile.
        output: New TOML file; existing files are never overwritten.
        workspace: Independent writable OrthoFinder workspace.

    Returns:
        Absolute path to the created profile.

    Raises:
        ResearchError: If the workspace is shared or the destination already exists.
        OSError: If the new file cannot be written.
    """
    original = source.expanduser().resolve()
    destination = output.expanduser().resolve()
    target_workspace = workspace.expanduser().resolve()
    config = load_config(path=original)
    if target_workspace == config.workspace:
        raise ResearchError("Choose a separate workspace from the source app's job queue")
    if destination == original or destination.exists():
        raise ResearchError(f"Profile already exists; it was not changed: {destination}")
    lines = [
        "# Database and tool settings reused from: " + str(original).replace("\n", " "),
        "# Only workspace differs. No databases or source jobs were copied.",
    ]
    for field in fields(ResearchConfig):
        value = target_workspace if field.name == "workspace" else getattr(config, field.name)
        serialisable = str(value) if isinstance(value, Path) else value
        lines.append(f"{field.name} = {json.dumps(serialisable, ensure_ascii=True)}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation makes concurrent setup attempts safe without replacing a profile.
    with destination.open(mode="x", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    LOGGER.info("Created independent research profile %s from %s", destination, original)
    return destination
