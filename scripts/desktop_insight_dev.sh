#!/usr/bin/env bash
# TEMPORARY: until openms-insight 0.7.10 is on PyPI, the desktop build installs it
# from its development branch (t0mdavid-m/OpenMS-Insight#151 and #155), frontend included.
# Called by build-desktop-app.yaml as: bash <script> <payload python> <payload dir>.
set -euo pipefail
py="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
commit=e98352380f7b28083bd406662b96e7b8fa221a86
work="${RUNNER_TEMP:-${TMPDIR:-/tmp}}/openms-insight-dev"
rm -rf "$work"
git clone --quiet https://github.com/t0mdavid-m/OpenMS-Insight "$work"
git -C "$work" checkout --quiet "$commit"
(cd "$work/js-component" && npm ci --no-audit --no-fund && npm run build)
mkdir -p "$work/openms_insight/js-component"
cp -R "$work/js-component/dist" "$work/openms_insight/js-component/"
"$py" -m pip install --no-deps --force-reinstall "$work"
"$py" -c "import openms_insight, pathlib; print('openms-insight', openms_insight.__version__); from openms_insight import SequenceCoverage; assert SequenceCoverage._supports_settings; assert (pathlib.Path(openms_insight.__file__).parent / 'js-component/dist').is_dir()"
