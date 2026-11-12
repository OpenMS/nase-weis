"""Reuse native workflow validation before creating a run or writing its logs."""

from copy import copy, deepcopy


class _StartupLog:
    def log(self, message: str, level: int = 0) -> None:
        """Turn a worker's validation diagnostic into an actionable UI error."""
        raise ValueError(message.removeprefix("ERROR:").strip())


def validate_startup(workflow, params, execute) -> None:
    """Validate a private settings copy without outputs, logs or a worker."""
    candidate = copy(workflow)
    candidate.params = deepcopy(params)
    candidate.logger = _StartupLog()
    if not execute(candidate, validate_only=True):
        raise ValueError("Review the input files and settings on Configure before starting.")
