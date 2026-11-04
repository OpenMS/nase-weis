"""Run the regression gate locally or in a freshly built, pinned test container."""

import argparse
import os
import shutil
import subprocess
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], *, env=None) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True, env=env)


def smoke_test_image(engine: str) -> None:
    """Exercise the production image's actual entrypoint without host ports/data."""
    name = f"naseweis-regression-smoke-{uuid.uuid4().hex[:10]}"
    try:
        run(
            [
                engine,
                "run",
                "--rm",
                "-d",
                "--network",
                "none",
                "--name",
                name,
                "-e",
                "LOCAL_WORKSPACES_DIR=/tmp/regression-workspaces",
                "naseweis-regression-app:local",
            ]
        )
        deadline = time.monotonic() + 45
        while time.monotonic() < deadline:
            health = subprocess.run(
                [
                    engine,
                    "exec",
                    name,
                    "/root/miniforge3/envs/streamlit-env/bin/python",
                    "-c",
                    (
                        "import urllib.request; r=urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=2); "
                        "assert r.status == 200 and r.read().strip() == b'ok'"
                    ),
                ],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if health.returncode == 0:
                print(
                    "Production container entrypoint and health endpoint passed.",
                    flush=True,
                )
                return
            time.sleep(0.25)
        subprocess.run([engine, "logs", name], check=False)
        raise SystemExit(
            "The production app image did not become healthy within 45 seconds."
        )
    finally:
        subprocess.run(
            [engine, "rm", "-f", name],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--container",
        action="store_true",
        help="Full gate, including scientific references, in the pinned runtime.",
    )
    parser.add_argument("--engine", default="docker", choices=["docker", "podman"])
    parser.add_argument(
        "--baseline",
        action="store_true",
        help="Require the reference runtime (used inside the container).",
    )
    parser.add_argument("--report-dir", type=Path, default=ROOT / ".local/regressions")
    args = parser.parse_args()
    reports = args.report_dir.resolve()
    reports.mkdir(parents=True, exist_ok=True)
    if args.container:
        if not shutil.which(args.engine):
            parser.error(f"{args.engine} is not installed")
        run(
            [
                args.engine,
                "build",
                "-f",
                "docker/Dockerfile.app",
                "-t",
                "naseweis-regression-app:local",
                ".",
            ]
        )
        smoke_test_image(args.engine)
        run(
            [
                args.engine,
                "build",
                "-f",
                "docker/Dockerfile.test",
                "-t",
                "naseweis-regression-test:local",
                ".",
            ]
        )
        run(
            [
                args.engine,
                "run",
                "--rm",
                "--network",
                "none",
                "-v",
                f"{reports}:/reports",
                "--entrypoint",
                "/root/miniforge3/envs/streamlit-env/bin/python",
                "naseweis-regression-test:local",
                "scripts/check_regressions.py",
                "--baseline",
                "--report-dir",
                "/reports",
            ]
        )
        return
    env = os.environ | {"COVERAGE_FILE": str(reports / ".coverage")}
    run([sys.executable, "-m", "utils.app_spec"])
    pytest = [
        sys.executable,
        "-m",
        "coverage",
        "run",
        "-m",
        "pytest",
        "test.py",
        "test_gui.py",
        "tests/",
        "--strict-regressions",
        "--strict-markers",
        "-o",
        "xfail_strict=True",
        "-q",
        "-ra",
        f"--junitxml={reports / 'junit.xml'}",
    ]
    if args.baseline:
        pytest.append("--naseweis-baseline")
    run(pytest, env=env)
    run([sys.executable, "-m", "coverage", "report"], env=env)
    run(
        [sys.executable, "-m", "coverage", "xml", "-o", str(reports / "coverage.xml")],
        env=env,
    )
    run([sys.executable, "-m", "pip", "check"])
    print(f"Regression gate passed. Reports: {reports}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        raise SystemExit(error.returncode) from error
