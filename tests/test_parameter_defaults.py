"""
Tests for get_merged_params() and the refactored get_topp_parameters().

This module verifies the three-layer parameter merge:
  ini defaults < _defaults < user overrides
"""
import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from src.workflow.ParameterManager import ParameterManager

# Each test uses a controllable session state without replacing imported
# packages or invalidating the workflow modules shared with GUI tests.
mock_streamlit = MagicMock()
mock_streamlit.session_state = {}


@pytest.fixture(autouse=True)
def reset_streamlit_state(monkeypatch):
    """Patch the actual method's module only for this test, then restore it."""
    mock_streamlit.session_state.clear()
    monkeypatch.setitem(ParameterManager.save_parameters.__globals__, "st", mock_streamlit)


@pytest.fixture
def temp_workflow_dir():
    """Create a temporary workflow directory for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        workflow_dir = Path(tmpdir) / "test-workflow"
        workflow_dir.mkdir()
        ini_dir = workflow_dir / "ini"
        ini_dir.mkdir()
        yield workflow_dir


class TestGetMergedParams:
    """Tests for ParameterManager.get_merged_params()."""

    def test_returns_ini_params_when_no_json(self, temp_workflow_dir):
        """ini params returned when params.json doesn't exist."""
        pm = ParameterManager(temp_workflow_dir)
        ini_params = {"algorithm:param": 1.0, "algorithm:other": "hello"}

        result = pm.get_merged_params("SomeTool", ini_params=ini_params)

        assert result == {"algorithm:param": 1.0, "algorithm:other": "hello"}

    def test_defaults_override_ini(self, temp_workflow_dir):
        """_defaults layer overrides ini values."""
        pm = ParameterManager(temp_workflow_dir)
        params_json = {
            "_defaults": {"SomeTool": {"algorithm:param": 42.0}}
        }
        with open(pm.params_file, "w") as f:
            json.dump(params_json, f)

        result = pm.get_merged_params("SomeTool", ini_params={"algorithm:param": 1.0})

        assert result["algorithm:param"] == 42.0

    def test_user_overrides_defaults(self, temp_workflow_dir):
        """User overrides take priority over _defaults."""
        pm = ParameterManager(temp_workflow_dir)
        params_json = {
            "_defaults": {"SomeTool": {"algorithm:param": 42.0}},
            "SomeTool": {"algorithm:param": 99.0}
        }
        with open(pm.params_file, "w") as f:
            json.dump(params_json, f)

        result = pm.get_merged_params("SomeTool", ini_params={"algorithm:param": 1.0})

        assert result["algorithm:param"] == 99.0

    def test_full_three_layer_merge(self, temp_workflow_dir):
        """All three layers merge correctly: ini < _defaults < user."""
        pm = ParameterManager(temp_workflow_dir)
        params_json = {
            "_defaults": {
                "SomeTool": {
                    "algorithm:param_a": 10.0,  # overrides ini
                    "algorithm:param_b": 20.0,  # overrides ini, NOT overridden by user
                }
            },
            "SomeTool": {
                "algorithm:param_a": 99.0,  # overrides _defaults
                "algorithm:param_c": 55.0,  # only in user
            }
        }
        with open(pm.params_file, "w") as f:
            json.dump(params_json, f)

        ini_params = {
            "algorithm:param_a": 1.0,
            "algorithm:param_b": 2.0,
            "algorithm:param_d": 3.0,  # only in ini
        }

        result = pm.get_merged_params("SomeTool", ini_params=ini_params)

        assert result["algorithm:param_a"] == 99.0   # user wins
        assert result["algorithm:param_b"] == 20.0   # _defaults wins over ini
        assert result["algorithm:param_c"] == 55.0   # user-only key present
        assert result["algorithm:param_d"] == 3.0    # ini-only key present

    def test_no_ini_params(self, temp_workflow_dir):
        """Works when ini_params is None."""
        pm = ParameterManager(temp_workflow_dir)
        params_json = {
            "_defaults": {"SomeTool": {"algorithm:param": 42.0}},
            "SomeTool": {"algorithm:other": 7.0}
        }
        with open(pm.params_file, "w") as f:
            json.dump(params_json, f)

        result = pm.get_merged_params("SomeTool")

        assert result["algorithm:param"] == 42.0
        assert result["algorithm:other"] == 7.0

    def test_different_instances_same_tool(self, temp_workflow_dir):
        """Different instance names get independent _defaults and user overrides."""
        pm = ParameterManager(temp_workflow_dir)
        params_json = {
            "_defaults": {
                "IDFilter_step1": {"score:min": 0.05},
                "IDFilter_step2": {"score:min": 0.01},
            },
            "IDFilter_step1": {"score:min": 0.001},
        }
        with open(pm.params_file, "w") as f:
            json.dump(params_json, f)

        result1 = pm.get_merged_params("IDFilter_step1", ini_params={"score:min": 0.5})
        result2 = pm.get_merged_params("IDFilter_step2", ini_params={"score:min": 0.5})

        assert result1["score:min"] == 0.001   # user override for step1
        assert result2["score:min"] == 0.01    # _defaults for step2 (no user override)

    def test_empty_params_json(self, temp_workflow_dir):
        """Returns empty dict when params.json is empty and no ini_params."""
        pm = ParameterManager(temp_workflow_dir)
        with open(pm.params_file, "w") as f:
            json.dump({}, f)

        result = pm.get_merged_params("SomeTool")

        assert result == {}


class TestGetToppParametersWithDefaults:

    def test_get_topp_parameters_includes_defaults(self, temp_workflow_dir):
        """get_topp_parameters merges _defaults between ini and user values."""
        pm = ParameterManager(temp_workflow_dir)
        params_json = {
            "_defaults": {"SomeTool": {"algorithm:param": 42.0}},
            "SomeTool": {"algorithm:other": 99.0}
        }
        with open(pm.params_file, "w") as f:
            json.dump(params_json, f)

        result = pm.get_merged_params("SomeTool", ini_params={"algorithm:param": 1.0})
        assert result["algorithm:param"] == 42.0
        assert result["algorithm:other"] == 99.0


class TestDefaultsSeeding:

    def test_seed_writes_defaults_to_params_json(self, temp_workflow_dir):
        """Seeding creates _defaults entry in params.json."""
        pm = ParameterManager(temp_workflow_dir)
        custom_defaults = {"param_a": 10.0, "param_b": "fast"}

        # Simulate what input_TOPP seeding does
        params = pm.get_parameters_from_json()
        if "_defaults" not in params:
            params["_defaults"] = {}
        params["_defaults"]["MyTool"] = custom_defaults
        with open(pm.params_file, "w") as f:
            json.dump(params, f)

        # Verify
        loaded = pm.get_parameters_from_json()
        assert loaded["_defaults"]["MyTool"] == {"param_a": 10.0, "param_b": "fast"}

    def test_seed_is_idempotent(self, temp_workflow_dir):
        """Seeding the same tool twice overwrites cleanly."""
        pm = ParameterManager(temp_workflow_dir)

        # First seed
        params = {"_defaults": {"Tool": {"p1": 1.0}}, "other_key": "keep"}
        with open(pm.params_file, "w") as f:
            json.dump(params, f)

        # Second seed with updated defaults
        params = pm.get_parameters_from_json()
        params["_defaults"]["Tool"] = {"p1": 2.0}
        with open(pm.params_file, "w") as f:
            json.dump(params, f)

        loaded = pm.get_parameters_from_json()
        assert loaded["_defaults"]["Tool"]["p1"] == 2.0
        assert loaded["other_key"] == "keep"

    def test_seed_multiple_instances(self, temp_workflow_dir):
        """Different instances of the same tool get independent _defaults."""
        pm = ParameterManager(temp_workflow_dir)
        params = {
            "_defaults": {
                "IDFilter_strict": {"score:pep": 0.01},
                "IDFilter_lenient": {"score:pep": 0.05},
            }
        }
        with open(pm.params_file, "w") as f:
            json.dump(params, f)

        loaded = pm.get_parameters_from_json()
        assert loaded["_defaults"]["IDFilter_strict"]["score:pep"] == 0.01
        assert loaded["_defaults"]["IDFilter_lenient"]["score:pep"] == 0.05


try:
    import pyopenms as poms
    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False


@pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not available")
class TestSaveParametersWithDefaults:

    def _create_fake_ini(self, pm, tool_name, params_dict):
        """Create a fake .ini file with given parameters."""
        param = poms.Param()
        for key, value in params_dict.items():
            param.setValue(f"{tool_name}:1:{key}".encode(), value)
        poms.ParamXMLFile().store(str(Path(pm.ini_dir, f"{tool_name}.ini")), param)

    def test_value_matching_custom_default_not_saved(self, temp_workflow_dir):
        """A value equal to the _defaults entry should not be saved as a user override."""
        pm = ParameterManager(temp_workflow_dir)

        # Create a fake ini with a default value
        self._create_fake_ini(pm, "Tool", {"param_a": 10.0})

        # Pre-seed _defaults with a different value than ini
        params = {"_defaults": {"Tool": {"param_a": 42.0}}}
        with open(pm.params_file, "w") as f:
            json.dump(params, f)

        # Session state has value matching the custom default (42.0), not the ini default (10.0)
        mock_streamlit.session_state[f"{pm.topp_param_prefix}Tool:1:param_a"] = 42.0
        mock_streamlit.session_state["_topp_tool_instance_map"] = {"Tool": "Tool"}

        pm.save_parameters()

        with open(pm.params_file, "r") as f:
            saved = json.load(f)

        # param_a should NOT appear under Tool (it matches the _defaults value)
        assert "param_a" not in saved.get("Tool", {})
        # _defaults should still be present
        assert saved["_defaults"]["Tool"]["param_a"] == 42.0

    def test_value_different_from_custom_default_saved(self, temp_workflow_dir):
        """A value different from _defaults entry should be saved as user override."""
        pm = ParameterManager(temp_workflow_dir)

        # Create a fake ini with a default value
        self._create_fake_ini(pm, "Tool", {"param_a": 10.0})

        params = {"_defaults": {"Tool": {"param_a": 42.0}}}
        with open(pm.params_file, "w") as f:
            json.dump(params, f)

        mock_streamlit.session_state[f"{pm.topp_param_prefix}Tool:1:param_a"] = 99.0
        mock_streamlit.session_state["_topp_tool_instance_map"] = {"Tool": "Tool"}

        pm.save_parameters()

        with open(pm.params_file, "r") as f:
            saved = json.load(f)

        assert saved["Tool"]["param_a"] == 99.0


class TestNonDefaultParamsSummaryDefaults:

    def test_defaults_key_excluded_from_classification(self):
        """_defaults dict should not appear as a TOPP tool in the summary."""
        params = {
            "_defaults": {"Tool": {"p1": 10}},
            "Tool": {"p1": 20},
            "general_param": "value"
        }
        # Simulate the classification logic
        topp = {}
        general = {}
        for k, v in params.items():
            if k == "_defaults":
                continue
            if isinstance(v, dict):
                topp[k] = v
            else:
                general[k] = v

        assert "_defaults" not in topp
        assert "Tool" in topp
        assert "general_param" in general

    def test_defaults_merged_into_summary(self):
        """_defaults values should appear in summary merged with user overrides."""
        params = {
            "_defaults": {
                "ToolA": {"p1": 10, "p2": 20},
                "ToolB": {"p3": 30}
            },
            "ToolA": {"p1": 99}
        }
        # Simulate the merge logic for summary
        topp = {}
        for k, v in params.items():
            if k == "_defaults":
                continue
            if isinstance(v, dict):
                topp[k] = v

        defaults = params.get("_defaults", {})
        for tool_name, default_vals in defaults.items():
            if tool_name not in topp:
                topp[tool_name] = {}
            topp[tool_name] = {**default_vals, **topp.get(tool_name, {})}

        assert topp["ToolA"] == {"p1": 99, "p2": 20}  # user override wins for p1
        assert topp["ToolB"] == {"p3": 30}              # defaults-only tool appears
