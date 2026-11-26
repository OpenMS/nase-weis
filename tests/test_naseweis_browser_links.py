"""Linked panels of Explore results in a real browser.

Starts the app on a fresh workspace seeded from the demo and runs
tests/browser/linked_panels.js, which clicks through the Abundances, Spectra
and Identification quality views and checks that each panel follows a pick.
Skipped where Node or Playwright is not installed.
"""
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tests" / "browser" / "linked_panels.js"


def _node_env():
    env = dict(os.environ)
    if shutil.which("npm"):
        root = subprocess.run(["npm", "root", "-g"], capture_output=True, text=True).stdout.strip()
        env["NODE_PATH"] = os.pathsep.join(p for p in (env.get("NODE_PATH"), root) if p)
    return env


def _playwright_available(env) -> bool:
    if not shutil.which("node"):
        return False
    probe = subprocess.run(["node", "-e", "require('playwright')"], env=env, capture_output=True)
    return probe.returncode == 0


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_linked_panels_follow_picks_in_the_browser(tmp_path):
    env = _node_env()
    if not _playwright_available(env):
        pytest.skip("Node with Playwright is not installed")
    port = _free_port()
    app_env = dict(os.environ, LOCAL_WORKSPACES_DIR=str(tmp_path / "workspaces"))
    server = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "app.py", "--server.port", str(port),
         "--server.headless", "true"],
        cwd=ROOT, env=app_env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        url = f"http://localhost:{port}"
        for _ in range(120):
            try:
                urllib.request.urlopen(url + "/_stcore/health", timeout=2)
                break
            except OSError:
                time.sleep(1)
        else:
            pytest.fail("the app did not start")
        result = subprocess.run(["node", str(SCRIPT), url], env=env, capture_output=True,
                                text=True, timeout=900)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "linked panels: ok" in result.stdout
    finally:
        server.terminate()
        server.wait(timeout=30)
