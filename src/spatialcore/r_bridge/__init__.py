"""R/Python integration via subprocess and h5ad interchange.

This module provides utilities for running R scripts from Python,
with proper configuration for reticulate to use the correct Python environment.

Key features:
- Subprocess-based R execution (no rpy2 dependency required)
- Automatic RETICULATE_PYTHON configuration
- JSON output parsing from R scripts
- Timeout handling for long-running scripts

Example usage:
    >>> from spatialcore.r_bridge import check_r_available, run_r_script
    >>> if check_r_available():
    ...     result = run_r_script("my_script.R", args=["input.h5ad"])
    ...     print(result)

Notes:
    - R must be installed and Rscript accessible in PATH
    - Install R via: mamba install -c conda-forge r-base
    - For reticulate-based scripts, RETICULATE_PYTHON is set automatically
    - Windows support is best-effort. If sf-based workflows crash (e.g. access
      violation 1073741819), use WSL or Linux for domain detection.
"""

from spatialcore.r_bridge.subprocess_runner import (
    RBridgeError,
    RExecutionError,
    RNotFoundError,
    RTimeoutError,
    check_r_available,
    get_r_version,
    run_r_code,
    run_r_script,
)

__all__ = [
    # Functions
    "check_r_available",
    "get_r_version",
    "run_r_script",
    "run_r_code",
    # Exceptions
    "RBridgeError",
    "RNotFoundError",
    "RExecutionError",
    "RTimeoutError",
]
