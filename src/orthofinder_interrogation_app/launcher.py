"""Named launcher for the standalone Streamlit application."""

from __future__ import annotations

import argparse
import logging
import os
import socket
import subprocess
import sys
from importlib.util import find_spec
from pathlib import Path
from typing import Sequence

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError
from orthofinder_results.io_utils import configure_logging

from .distance_data import default_cache_directory
from .resource import open_resource

_LOGGER = logging.getLogger("orthofinder_interrogation_app.launcher")
RESOURCE_ENVIRONMENT_VARIABLE = "ORTHOFINDER_RESULTS_RESOURCE"
LOG_ENVIRONMENT_VARIABLE = "ORTHOFINDER_RESULTS_APP_LOG"
TAXONOMY_ENVIRONMENT_VARIABLE = "ORTHOFINDER_RESULTS_TAXONOMY"
CACHE_ENVIRONMENT_VARIABLE = "ORTHOFINDER_RESULTS_CACHE"
DEFAULT_PORT = 8501
MAX_AUTOMATIC_PORT_ATTEMPTS = 100


def build_parser() -> argparse.ArgumentParser:
    """Build the all-named-option application launcher parser."""

    parser = argparse.ArgumentParser(
        prog="orthofinder-interrogation-app",
        description="Open a completed OrthoFinder resource in a read-only local application.",
    )
    parser.add_argument(
        "--resource-dir",
        required=True,
        type=Path,
        help="Completed resource directory or orthofinder_results.duckdb file.",
    )
    parser.add_argument("--server-address", default="127.0.0.1")
    parser.add_argument(
        "--server-port",
        type=int,
        help=(
            "Exact Streamlit port. When omitted, select the first available port from "
            f"{DEFAULT_PORT} upwards."
        ),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        help="Optional persistent application log outside the immutable resource.",
    )
    parser.add_argument(
        "--taxonomy-map",
        type=Path,
        help="Optional reviewed taxonomy TSV sidecar; the resource is never modified.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help=(
            "Persistent sidecar cache for on-demand analyses. Defaults to the macOS "
            "Library cache or the Linux/XDG user cache; /tmp is never assumed."
        ),
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="Do not ask Streamlit to open a browser automatically.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Validate the resource and launch Streamlit.

    Args:
        argv: Optional named command-line arguments.

    Returns:
        Streamlit subprocess exit status, or two for a controlled error.
    """

    parser = build_parser()
    arguments = parser.parse_args(argv)
    if arguments.server_port is not None and not 1 <= arguments.server_port <= 65535:
        parser.error("--server-port must be between 1 and 65535.")
    configure_logging(log_path=arguments.log_file, verbose=arguments.verbose)
    if find_spec("streamlit") is None:
        _LOGGER.error(
            "Streamlit is not installed. Install the application dependencies with "
            "python -m pip install --editable '.[app]'."
        )
        return 2
    try:
        resource = open_resource(path=arguments.resource_dir)
        selected_port = select_server_port(
            address=arguments.server_address,
            requested_port=arguments.server_port,
        )
    except OrthoFinderResultsError as error:
        _LOGGER.error("Cannot launch application: %s", error)
        return 2
    environment = os.environ.copy()
    environment[RESOURCE_ENVIRONMENT_VARIABLE] = str(resource.resource_path)
    cache_dir = (
        arguments.cache_dir.expanduser().resolve()
        if arguments.cache_dir is not None
        else default_cache_directory().resolve()
    )
    environment[CACHE_ENVIRONMENT_VARIABLE] = str(cache_dir)
    if arguments.log_file is not None:
        environment[LOG_ENVIRONMENT_VARIABLE] = str(arguments.log_file.expanduser().resolve())
    if arguments.taxonomy_map is not None:
        taxonomy_path = arguments.taxonomy_map.expanduser().resolve()
        if not taxonomy_path.is_file():
            _LOGGER.error("Taxonomy mapping is not a file: %s", taxonomy_path)
            return 2
        environment[TAXONOMY_ENVIRONMENT_VARIABLE] = str(taxonomy_path)
    app_path = Path(__file__).with_name("app.py").resolve()
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_path),
        "--server.address",
        arguments.server_address,
        "--server.port",
        str(selected_port),
        "--server.headless",
        str(arguments.headless).lower(),
        "--browser.gatherUsageStats",
        "false",
    ]
    _LOGGER.info(
        "Launching read-only application: run=%s, address=%s, port=%s, cache=%s",
        resource.run_id,
        arguments.server_address,
        selected_port,
        cache_dir,
    )
    try:
        completed = subprocess.run(command, check=False, env=environment)
    except OSError as error:
        _LOGGER.error("Could not start Streamlit: %s", error)
        return 2
    return completed.returncode


def select_server_port(*, address: str, requested_port: int | None) -> int:
    """Return an exact free port or find a bounded automatic alternative.

    Args:
        address: Interface address passed to Streamlit.
        requested_port: User-requested exact port, or ``None`` for automatic selection.

    Returns:
        Available TCP port.

    Raises:
        InputValidationError: If the address is invalid or no acceptable port is free.
    """

    if requested_port is not None:
        if _port_available(address=address, port=requested_port):
            return requested_port
        raise InputValidationError(
            f"Requested port {requested_port} is already in use on {address}."
        )
    for port in range(DEFAULT_PORT, DEFAULT_PORT + MAX_AUTOMATIC_PORT_ATTEMPTS):
        if _port_available(address=address, port=port):
            if port != DEFAULT_PORT:
                _LOGGER.warning(
                    "Default port %s is unavailable; selected port %s.",
                    DEFAULT_PORT,
                    port,
                )
            return port
    raise InputValidationError(
        f"No free port was found on {address} from {DEFAULT_PORT} to "
        f"{DEFAULT_PORT + MAX_AUTOMATIC_PORT_ATTEMPTS - 1}."
    )


def _port_available(*, address: str, port: int) -> bool:
    """Return whether at least one resolved TCP address can be bound."""

    try:
        candidates = socket.getaddrinfo(
            address,
            port,
            type=socket.SOCK_STREAM,
            flags=socket.AI_PASSIVE,
        )
    except socket.gaierror as error:
        raise InputValidationError(f"Server address could not be resolved: {address}") from error
    for family, socket_type, protocol, _, socket_address in candidates:
        try:
            with socket.socket(family, socket_type, protocol) as probe:
                probe.bind(socket_address)
        except OSError:
            continue
        return True
    return False


if __name__ == "__main__":  # pragma: no cover - exercised through the console script
    raise SystemExit(main())
