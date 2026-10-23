"""Exercise SSH transport and remote service ownership without contacting a server."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

LAUNCHER = Path(__file__).resolve().parents[1] / "scripts" / "launch_orthofinder_denbi.sh"

STUB = r"""
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

name = Path(sys.argv[0]).name
args = sys.argv[1:]
root = Path(os.environ["OF_TEST_ROOT"])
with (root / "calls.jsonl").open("a") as handle:
    handle.write(json.dumps([name, args]) + "\n")
state_path = root / "state.json"
state = (
    json.loads(state_path.read_text())
    if state_path.exists() else {"windows": {}, "master": False}
)

def save():
    state_path.write_text(json.dumps(state))

if name == "ssh":
    if "-T" in args:
        result = subprocess.run(["bash", "-c", args[-1]], input=sys.stdin.read(), text=True)
        sys.exit(result.returncode)
    if "-M" in args:
        state["master"] = True
    elif "-O" in args:
        action = args[args.index("-O") + 1]
        if action == "check":
            sys.exit(0 if state["master"] else 1)
        if action == "exit":
            state["master"] = False
        if action == "forward" and os.environ.get("OF_TEST_FORWARD_FAIL"):
            sys.exit(1)
    save()
elif name == "curl":
    sys.exit(1 if os.environ.get("OF_TEST_HEALTH_FAIL") else 0)
elif name == "orthofinder-interrogation-app":
    sys.exit(2 if os.environ.get("OF_TEST_VALIDATION_FAIL") else 0)
elif name == "orthofinder-protein-workbench":
    print("[]")
elif name == "tmux":
    command = args[0]
    target = args[args.index("-t") + 1] if "-t" in args else ""
    windows = state["windows"]
    if command == "list-panes":
        if target not in windows:
            sys.exit(1)
        key = "command" if args[-1] == "#{pane_start_command}" else "dead"
        print(windows[target][key])
    elif command == "show-options":
        print(windows[target]["options"].get(args[-1], ""))
    elif command == "set-option":
        windows[target]["options"][args[-2]] = args[-1]
    elif command in {"new-session", "new-window", "respawn-pane"}:
        if command == "respawn-pane":
            argv = args[args.index("-t") + 2:]
        else:
            position = args.index("-n")
            target = "orthofinder-research:=" + args[position + 1]
            argv = args[position + 2:]
        windows[target] = {"dead": 0, "command": shlex.join(argv), "options": {}}
    elif command == "kill-pane":
        del windows[target]
    elif command in {"list-windows", "has-session"}:
        if not windows:
            sys.exit(1)
        if command == "list-windows":
            print("app worker")
    elif command == "capture-pane":
        print("Fixture console diagnostics")
    else:
        raise AssertionError(args)
    save()
"""


@pytest.fixture
def denbi(tmp_path):
    """Create isolated transport, tmux and service executables with a persistent call log."""
    binaries = tmp_path / "environment" / "bin"
    binaries.mkdir(parents=True)
    for name in (
        "ssh",
        "curl",
        "tmux",
        "orthofinder-interrogation-app",
        "orthofinder-protein-workbench",
    ):
        path = binaries / name
        path.write_text(f"#!{sys.executable}\n" + STUB, encoding="utf-8")
        path.chmod(0o755)
    profile = tmp_path / "research.toml"
    profile.write_text('workspace="/server/orthofinder/jobs"\n', encoding="utf-8")
    resource = tmp_path / "resource with spaces and 'quotes'"
    resource.mkdir()
    environment = dict(
        os.environ, PATH=f"{binaries}:{os.environ['PATH']}", OF_TEST_ROOT=str(tmp_path)
    )
    arguments = [
        "--no-browser",
        "--env-prefix",
        str(binaries.parent),
        "--research-config",
        str(profile),
        "--resource-dir",
        str(resource),
        "--state-dir",
        str(tmp_path / "research state"),
    ]
    return SimpleNamespace(
        root=tmp_path,
        environment=environment,
        arguments=arguments,
        profile=profile,
        resource=resource,
    )


def launch(*, denbi, arguments=()):
    """Run the unmodified shell with controlled executables and a bounded timeout."""
    return subprocess.run(
        ["bash", str(LAUNCHER), *denbi.arguments, *arguments],
        env=denbi.environment,
        capture_output=True,
        text=True,
        timeout=15,
    )


def calls(*, denbi):
    """Read observed executable arguments from the stub boundary."""
    return [json.loads(line) for line in (denbi.root / "calls.jsonl").read_text().splitlines()]


def test_launch_reconnect_status_and_disconnect_preserve_remote_services(denbi):
    """Remote analysis stays running across repeated connections and tunnel shutdown."""
    initial = launch(denbi=denbi)
    assert initial.returncode == 0, initial.stderr
    assert "http://127.0.0.1:8510" in initial.stdout
    state = json.loads((denbi.root / "state.json").read_text())
    assert set(state["windows"]) == {"orthofinder-research:=app", "orthofinder-research:=worker"}
    app_call = next(
        args for name, args in calls(denbi=denbi) if name == "orthofinder-interrogation-app"
    )
    assert app_call[app_call.index("--resource-dir") + 1] == str(denbi.resource)
    assert "--validate-only" in app_call and "8506" in app_call
    assert launch(denbi=denbi).returncode == 0
    assert launch(denbi=denbi, arguments=["--status"]).returncode == 0
    after = calls(denbi=denbi)
    assert sum(args[0] == "new-session" for name, args in after if name == "tmux") == 1
    assert sum(args[0] == "new-window" for name, args in after if name == "tmux") == 1
    assert launch(denbi=denbi, arguments=["--stop-tunnel"]).returncode == 0
    final = json.loads((denbi.root / "state.json").read_text())
    assert final["windows"] == state["windows"] and not final["master"]
    forwards = [args for name, args in after if name == "ssh" and "forward" in args]
    assert all("ClearAllForwardings=yes" not in args and "/dev/null" in args for args in forwards)


def test_remote_preflight_failure_does_not_start_services_or_tunnel(denbi):
    """An invalid resource is rejected before a worker or SSH master is created."""
    denbi.environment["OF_TEST_VALIDATION_FAIL"] = "1"
    result = launch(denbi=denbi)
    assert result.returncode == 1 and "Resource validation failed" in result.stderr
    assert not any(name == "ssh" and "-M" in args for name, args in calls(denbi=denbi))
    assert not any(
        name == "tmux" and args[0].startswith("new-") for name, args in calls(denbi=denbi)
    )


def test_app_restart_preserves_worker_and_changed_profile_is_not_silently_reused(denbi):
    """Only an explicitly requested GUI restart replaces its owned pane."""
    assert launch(denbi=denbi).returncode == 0
    changed_port = launch(denbi=denbi, arguments=["--remote-port", "8512"])
    assert changed_port.returncode == 1 and "App settings changed" in changed_port.stderr
    restarted = launch(denbi=denbi, arguments=["--remote-port", "8512", "--restart-app"])
    assert restarted.returncode == 0, restarted.stderr
    killed = [
        args for name, args in calls(denbi=denbi) if name == "tmux" and args[0] == "kill-pane"
    ]
    assert killed == [["kill-pane", "-t", "orthofinder-research:=app"]]
    denbi.profile.write_text('workspace="/server/different/jobs"\n', encoding="utf-8")
    changed_profile = launch(denbi=denbi)
    assert (
        changed_profile.returncode == 1 and "live worker uses different" in changed_profile.stderr
    )


def test_launcher_refuses_unrecognised_windows_and_cleans_failed_tunnel(denbi):
    """An unrelated pane cannot be replaced; failed forwards preserve server processes."""
    assert launch(denbi=denbi).returncode == 0
    path = denbi.root / "state.json"
    state = json.loads(path.read_text())
    state["windows"]["orthofinder-research:=app"]["options"] = {}
    path.write_text(json.dumps(state))
    refused = launch(denbi=denbi, arguments=["--restart-app"])
    assert refused.returncode == 1 and "unrecognised window" in refused.stderr
    state["windows"]["orthofinder-research:=app"]["options"] = {
        "@orthofinder_owner": "orthofinder-launcher-1"
    }
    state["windows"] = {}
    state["master"] = False
    path.write_text(json.dumps(state))
    denbi.environment["OF_TEST_FORWARD_FAIL"] = "1"
    failed = launch(denbi=denbi)
    assert failed.returncode == 1 and "Tunnel failed" in failed.stderr
    final = json.loads(path.read_text())
    assert len(final["windows"]) == 2 and not final["master"]


@pytest.mark.parametrize(
    "arguments",
    [
        ["--local-port", "22"],
        ["--remote-port", "65536"],
        ["--remote-port", "bad"],
        ["--remote-port"],
        ["--resource-dir", "relative/path"],
        ["--state-dir", ""],
        ["--ssh-host", "-oProxyCommand=bad"],
        ["--status", "--stop-tunnel"],
        ["--restart-app", "--status"],
        ["--unknown"],
        ["--resource-dir", "/path\nline"],
    ],
)
def test_invalid_options_fail_before_network(denbi, arguments):
    """Invalid ports, remote paths and conflicting modes cause no SSH connection."""
    result = launch(denbi=denbi, arguments=arguments)
    assert result.returncode == 2
    assert not (denbi.root / "calls.jsonl").exists()
