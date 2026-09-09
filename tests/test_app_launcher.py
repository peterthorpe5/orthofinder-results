"""Tests for the named Streamlit application launcher."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock, Mock

import pytest

from orthofinder_interrogation_app import launcher
from orthofinder_results.errors import InputValidationError


def test_launcher_validates_resource_and_builds_streamlit_command(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The launcher passes only validated named settings to Streamlit."""

    completed = Mock(returncode=7)
    run = Mock(return_value=completed)
    monkeypatch.setattr(launcher.subprocess, "run", run)
    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: True)
    log_file = tmp_path / "logs" / "app.log"
    expected_file = tmp_path / "expected.tsv"
    expected_file.write_text("taxon_id\treason\tsource\tincluded\n", encoding="utf-8")
    focus_file = tmp_path / "focus.tsv"
    focus_file.write_text("protein_identifier\nQ9SA03\n", encoding="utf-8")
    status = launcher.main(
        [
            "--resource-dir",
            str(application_resource),
            "--server-address",
            "0.0.0.0",
            "--server-port",
            "8765",
            "--log-file",
            str(log_file),
            "--taxonomy-map",
            str(taxonomy_mapping_file),
            "--expected-taxa",
            str(expected_file),
            "--focus-proteins",
            str(focus_file),
            "--headless",
            "--verbose",
        ]
    )
    assert status == 7
    command = run.call_args.args[0]
    assert command[:4] == [launcher.sys.executable, "-m", "streamlit", "run"]
    assert "--server.address" in command and "0.0.0.0" in command
    assert "--server.port" in command and "8765" in command
    environment = run.call_args.kwargs["env"]
    assert environment[launcher.RESOURCE_ENVIRONMENT_VARIABLE] == str(
        application_resource.resolve()
    )
    assert environment[launcher.LOG_ENVIRONMENT_VARIABLE] == str(log_file.resolve())
    assert environment[launcher.TAXONOMY_ENVIRONMENT_VARIABLE] == str(
        taxonomy_mapping_file.resolve()
    )
    assert environment[launcher.EXPECTED_TAXA_ENVIRONMENT_VARIABLE] == str(
        expected_file.resolve()
    )
    assert environment[launcher.FOCUS_ENVIRONMENT_VARIABLE] == str(
        focus_file.resolve()
    )
    assert log_file.is_file()


def test_launcher_returns_controlled_errors(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Missing dependencies, resources and subprocess failures return status two."""

    monkeypatch.setattr(launcher, "find_spec", lambda name: None)
    assert launcher.main(["--resource-dir", str(application_resource)]) == 2
    monkeypatch.setattr(launcher, "find_spec", lambda name: object())
    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: True)
    assert launcher.main(["--resource-dir", str(application_resource / "missing")]) == 2
    monkeypatch.setattr(
        launcher.subprocess,
        "run",
        Mock(side_effect=OSError("cannot execute")),
    )
    assert launcher.main(["--resource-dir", str(application_resource)]) == 2


@pytest.mark.parametrize("port", (0, 65536))
def test_launcher_rejects_invalid_ports(application_resource: Path, port: int) -> None:
    """Invalid network ports fail in argument parsing before launch."""

    with pytest.raises(SystemExit) as error:
        launcher.main(["--resource-dir", str(application_resource), "--server-port", str(port)])
    assert error.value.code == 2


def test_launcher_without_log_file_preserves_environment(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An optional app log is not invented inside the resource."""

    monkeypatch.delenv(launcher.LOG_ENVIRONMENT_VARIABLE, raising=False)
    run = Mock(return_value=Mock(returncode=0))
    monkeypatch.setattr(launcher.subprocess, "run", run)
    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: True)
    assert launcher.main(["--resource-dir", str(application_resource)]) == 0
    environment = run.call_args.kwargs["env"]
    assert launcher.LOG_ENVIRONMENT_VARIABLE not in environment
    assert environment[launcher.RESOURCE_ENVIRONMENT_VARIABLE]
    assert os.path.isabs(environment[launcher.RESOURCE_ENVIRONMENT_VARIABLE])


def test_automatic_port_selection_skips_busy_ports(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An omitted port selects the first available local Streamlit port."""

    def available(*, address: str, port: int) -> bool:
        del address
        return port == 8503

    monkeypatch.setattr(launcher, "_port_available", available)
    run = Mock(return_value=Mock(returncode=0))
    monkeypatch.setattr(launcher.subprocess, "run", run)
    assert launcher.main(["--resource-dir", str(application_resource)]) == 0
    command = run.call_args.args[0]
    assert command[command.index("--server.port") + 1] == "8503"


def test_requested_busy_port_and_missing_taxonomy_are_controlled(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exact occupied ports and absent optional sidecars fail before Streamlit starts."""

    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: False)
    assert (
        launcher.main(
            [
                "--resource-dir",
                str(application_resource),
                "--server-port",
                "8501",
            ]
        )
        == 2
    )
    assert (
        launcher.main(
            [
                "--resource-dir",
                str(application_resource),
                "--expected-taxa",
                str(application_resource / "missing_expected.tsv"),
            ]
        )
        == 2
    )
    assert (
        launcher.main(
            [
                "--resource-dir",
                str(application_resource),
                "--focus-proteins",
                str(application_resource / "missing_focus.tsv"),
            ]
        )
        == 2
    )
    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: True)
    assert (
        launcher.main(
            [
                "--resource-dir",
                str(application_resource),
                "--taxonomy-map",
                str(application_resource / "missing.tsv"),
            ]
        )
        == 2
    )


def test_port_selection_exhaustion_and_invalid_address(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bounded port search and address resolution errors remain actionable."""

    real_port_available = launcher._port_available
    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: False)
    with pytest.raises(InputValidationError, match="No free port"):
        launcher.select_server_port(address="127.0.0.1", requested_port=None)
    monkeypatch.setattr(launcher, "_port_available", real_port_available)
    with pytest.raises(InputValidationError, match="could not be resolved"):
        launcher._port_available(address="invalid host name !", port=8501)


def test_port_probe_tries_resolved_addresses_and_closes_sockets(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A failed address family does not hide a later bindable address."""

    candidates = [
        (2, 1, 6, "", ("127.0.0.1", 8501)),
        (2, 1, 6, "", ("127.0.0.2", 8501)),
    ]
    monkeypatch.setattr(launcher.socket, "getaddrinfo", lambda *args, **kwargs: candidates)
    failed = MagicMock()
    failed.__enter__.return_value.bind.side_effect = OSError("busy")
    available = MagicMock()
    probes = iter((failed, available))
    monkeypatch.setattr(launcher.socket, "socket", lambda *args: next(probes))
    assert launcher._port_available(address="localhost", port=8501)

    only_failed = MagicMock()
    only_failed.__enter__.return_value.bind.side_effect = OSError("busy")
    monkeypatch.setattr(launcher.socket, "socket", lambda *args: only_failed)
    assert not launcher._port_available(address="localhost", port=8501)
