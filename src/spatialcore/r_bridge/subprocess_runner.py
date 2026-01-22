"""R subprocess runner with reticulate configuration.

This module provides utilities for running R scripts from Python via subprocess,
with proper configuration for reticulate to use the correct Python environment.

Key features:
- Automatic RETICULATE_PYTHON configuration
- Timeout handling for long-running R scripts
- JSON output parsing from R script stdout
- Graceful error handling with descriptive exceptions

References:
    - reticulate: https://rstudio.github.io/reticulate/
    - HIERATYPE_PIPELINE_SPECIFICATION.md for implementation details
"""

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Exception Classes
# ============================================================================


class RBridgeError(Exception):
    """Base exception for R bridge errors."""

    pass


class RNotFoundError(RBridgeError):
    """R is not installed or not in PATH."""

    pass


class RExecutionError(RBridgeError):
    """R script execution failed."""

    pass


class RTimeoutError(RBridgeError):
    """R script execution timed out."""

    pass


# ============================================================================
# R Availability Check
# ============================================================================


def _find_mamba_or_conda() -> Optional[str]:
    """
    Find mamba or conda executable for running R in proper environment.

    On Windows, mamba/conda handle DLL paths that Rscript alone cannot.

    Returns
    -------
    str or None
        Path to mamba.bat or conda.bat, or None if not found.
    """
    # Check for CONDA_EXE environment variable first
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe and Path(conda_exe).exists():
        # Prefer mamba if available in same location
        mamba_path = Path(conda_exe).parent / "mamba.exe"
        if mamba_path.exists():
            return str(mamba_path)
        return conda_exe

    # Try to find mamba/conda from sys.prefix (current environment)
    env_prefix = Path(sys.prefix)

    # Go up to find base conda/mamba installation
    # e.g., C:/SpatialCore/miniforge3/envs/spatialcore-dev -> C:/SpatialCore/miniforge3
    if "envs" in env_prefix.parts:
        idx = env_prefix.parts.index("envs")
        base_prefix = Path(*env_prefix.parts[:idx])
    else:
        base_prefix = env_prefix

    # Look for mamba/conda in condabin or Scripts
    search_paths = [
        base_prefix / "condabin" / "mamba.bat",
        base_prefix / "condabin" / "conda.bat",
        base_prefix / "Scripts" / "mamba.exe",
        base_prefix / "Scripts" / "conda.exe",
        base_prefix / "bin" / "mamba",
        base_prefix / "bin" / "conda",
    ]

    for path in search_paths:
        if path.exists():
            return str(path)

    return None


def _get_current_env_name() -> Optional[str]:
    """Get the name of the current conda/mamba environment."""
    # Check CONDA_DEFAULT_ENV first
    env_name = os.environ.get("CONDA_DEFAULT_ENV")
    if env_name:
        return env_name

    # Try to extract from sys.prefix
    # e.g., C:/SpatialCore/miniforge3/envs/spatialcore-dev -> spatialcore-dev
    env_prefix = Path(sys.prefix)
    if "envs" in env_prefix.parts:
        idx = env_prefix.parts.index("envs")
        if idx + 1 < len(env_prefix.parts):
            return env_prefix.parts[idx + 1]

    return None


def _find_rscript() -> Optional[str]:
    """
    Find Rscript executable in conda/mamba environment or system PATH.

    On Windows with conda/mamba, Rscript is typically at:
    - {env}/Lib/R/bin/Rscript.exe
    - {env}/Scripts/Rscript.exe

    Returns
    -------
    str or None
        Path to Rscript executable, or None if not found.
    """
    import shutil

    # First check system PATH
    rscript_path = shutil.which("Rscript")
    if rscript_path:
        return rscript_path

    # On Windows, check conda/mamba environment locations
    env_prefix = Path(sys.prefix)

    # Typical Windows conda R locations
    windows_paths = [
        env_prefix / "Lib" / "R" / "bin" / "Rscript.exe",
        env_prefix / "Lib" / "R" / "bin" / "x64" / "Rscript.exe",
        env_prefix / "Scripts" / "Rscript.exe",
        env_prefix / "Library" / "mingw-w64" / "bin" / "Rscript.exe",
    ]

    # Unix/Linux conda paths
    unix_paths = [
        env_prefix / "lib" / "R" / "bin" / "Rscript",
        env_prefix / "bin" / "Rscript",
    ]

    # Try Windows paths first, then Unix
    for path in windows_paths + unix_paths:
        if path.exists():
            return str(path)

    return None


def _build_rscript_command(rscript_args: List[str]) -> Optional[List[str]]:
    """
    Build the appropriate command for running Rscript.

    If in a conda/mamba environment, uses 'mamba/conda run -n env_name Rscript ...'
    to ensure libraries are properly loaded. Otherwise, calls Rscript directly.

    Parameters
    ----------
    rscript_args : list of str
        Arguments to pass to Rscript (e.g., ["--version"] or ["script.R", "arg1"])

    Returns
    -------
    list of str or None
        Complete command list ready for subprocess.run(), or None if Rscript not found.
    """
    rscript = _find_rscript()
    if not rscript:
        return None

    # Check if we're in a conda/mamba environment
    mamba_conda = _find_mamba_or_conda()
    env_name = _get_current_env_name()

    if mamba_conda and env_name:
        # Use mamba/conda run to properly set up library paths
        return [mamba_conda, "run", "-n", env_name, "Rscript"] + rscript_args
    else:
        # Direct Rscript execution (system R or activated conda env)
        return [rscript] + rscript_args


def check_r_available() -> bool:
    """
    Check if R is available in conda environment or system PATH.

    Returns
    -------
    bool
        True if Rscript is found and executable, False otherwise.

    Examples
    --------
    >>> from spatialcore.r_bridge import check_r_available
    >>> if check_r_available():
    ...     print("R is available")
    ... else:
    ...     print("R not found - install R first")
    """
    cmd = _build_rscript_command(["--version"])
    if not cmd:
        return False

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False
    except subprocess.TimeoutExpired:
        return False
    except Exception:
        return False


def get_r_version() -> Optional[str]:
    """
    Get the installed R version.

    Returns
    -------
    str or None
        R version string if available, None otherwise.

    Examples
    --------
    >>> from spatialcore.r_bridge import get_r_version
    >>> version = get_r_version()
    >>> if version:
    ...     print(f"R version: {version}")
    """
    cmd = _build_rscript_command(["--version"])
    if not cmd:
        return None

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0:
            # R --version outputs to stderr
            output = result.stderr or result.stdout
            # First line usually contains version
            for line in output.split("\n"):
                if "version" in line.lower():
                    return line.strip()
            return output.split("\n")[0].strip()
        return None
    except Exception:
        return None


# ============================================================================
# R Script Execution
# ============================================================================


def run_r_script(
    script_path: Union[str, Path],
    args: Optional[List[str]] = None,
    env_vars: Optional[Dict[str, str]] = None,
    timeout: int = 3600,
    cwd: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Run an R script via subprocess with proper reticulate configuration.

    This function sets RETICULATE_PYTHON to the current Python interpreter,
    ensuring that R's reticulate package uses the same Python environment
    (with all installed packages like anndata, scipy, etc.).

    Parameters
    ----------
    script_path : str or Path
        Path to the R script to execute.
    args : list of str, optional
        Command-line arguments to pass to the script.
    env_vars : dict, optional
        Additional environment variables to set.
        RETICULATE_PYTHON is set automatically.
    timeout : int, default 3600
        Maximum execution time in seconds (default 1 hour).
    cwd : str or Path, optional
        Working directory for script execution.
        Defaults to the script's parent directory.

    Returns
    -------
    Dict[str, Any]
        Parsed JSON output from the R script.
        The R script should output valid JSON on its last stdout line.
        If no valid JSON is found, returns {"stdout": ..., "stderr": ...}.

    Raises
    ------
    FileNotFoundError
        If the R script file does not exist.
    RNotFoundError
        If Rscript is not found in PATH.
    RExecutionError
        If the R script returns a non-zero exit code.
    RTimeoutError
        If the script exceeds the timeout.

    Examples
    --------
    >>> from spatialcore.r_bridge import run_r_script
    >>> result = run_r_script(
    ...     "my_script.R",
    ...     args=["input.h5ad", "output.h5ad"],
    ...     timeout=600,
    ... )
    >>> print(result)
    {'status': 'success', 'n_cells': 10000}

    Notes
    -----
    The R script should output its result as JSON on the last line of stdout.
    All other stdout lines are logged at DEBUG level.

    RETICULATE_PYTHON is set to sys.executable to ensure R's reticulate
    uses the same Python environment that called this function.
    """
    script_path = Path(script_path)
    if not script_path.exists():
        raise FileNotFoundError(f"R script not found: {script_path}")

    rscript_args = [str(script_path)]
    if args:
        rscript_args.extend([str(a) for a in args])

    cmd = _build_rscript_command(rscript_args)
    if not cmd:
        raise RNotFoundError(
            "Rscript not found in conda environment or PATH. "
            "Install R via: mamba install -c conda-forge r-base"
        )

    # Build environment with RETICULATE_PYTHON
    env = os.environ.copy()
    env["RETICULATE_PYTHON"] = sys.executable
    if env_vars:
        env.update(env_vars)

    uses_conda_run = len(cmd) >= 4 and cmd[1] == "run"
    if uses_conda_run:
        runner_name = Path(cmd[0]).name.lower()
        if "mamba" in runner_name and not env.get("MAMBA_ROOT_PREFIX"):
            env_prefix = Path(sys.prefix)
            if "envs" in env_prefix.parts:
                idx = env_prefix.parts.index("envs")
                base_prefix = Path(*env_prefix.parts[:idx])
            else:
                base_prefix = env_prefix
            env["MAMBA_ROOT_PREFIX"] = str(base_prefix)
            logger.info(f"Setting MAMBA_ROOT_PREFIX={base_prefix}")
    if uses_conda_run:
        logger.info("Using conda/mamba run for Rscript execution")
    else:
        logger.info("Using direct Rscript execution")

    # Determine working directory
    working_dir = str(cwd) if cwd else str(script_path.parent)

    logger.info(f"Running R script: {script_path.name}")
    logger.debug(f"  Command: {' '.join(cmd)}")
    logger.debug(f"  Working directory: {working_dir}")
    logger.debug(f"  RETICULATE_PYTHON: {env['RETICULATE_PYTHON']}")
    logger.debug(f"  Timeout: {timeout}s")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=working_dir,
            env=env,
        )
    except subprocess.TimeoutExpired as e:
        raise RTimeoutError(
            f"R script timed out after {timeout} seconds. "
            "Consider increasing timeout for large datasets."
        ) from e

    # Log stdout (excluding JSON line)
    if result.stdout:
        lines = result.stdout.strip().split("\n")
        for line in lines:
            # Skip empty lines and the JSON output line
            if line.strip() and not line.strip().startswith("{"):
                logger.debug(f"R: {line}")

    # Check for errors
    if result.returncode != 0:
        error_msg = result.stderr.strip() if result.stderr else "Unknown error"
        # Also include stdout as it may contain error messages
        if result.stdout:
            error_msg = f"{error_msg}\n\nStdout:\n{result.stdout}"
        if "there is no package called" in error_msg:
            error_msg = (
                f"{error_msg}\n\n"
                "Hint: Missing R package(s). Install the required R dependencies "
                "for this function."
            )
        if sys.platform == "win32" and (
            "1073741819" in error_msg or "0xC0000005" in error_msg
        ):
            error_msg = (
                f"{error_msg}\n\n"
                "Hint: R crashed on Windows (access violation). "
                "For sf-based spatial workflows, use WSL or Linux."
            )
        raise RExecutionError(
            f"R script failed (exit code {result.returncode}):\n{error_msg}"
        )

    # Parse JSON from last line of stdout
    output_lines = result.stdout.strip().split("\n") if result.stdout else []
    if not output_lines:
        logger.warning("R script produced no output")
        return {}

    # Find the last line that looks like JSON
    json_line = None
    for line in reversed(output_lines):
        line = line.strip()
        if line.startswith("{") and line.endswith("}"):
            json_line = line
            break

    if json_line:
        try:
            return json.loads(json_line)
        except json.JSONDecodeError as e:
            logger.warning(f"Failed to parse JSON from R output: {e}")
            return {"stdout": result.stdout, "stderr": result.stderr, "parse_error": str(e)}
    else:
        logger.debug("R script did not return JSON on last line")
        return {"stdout": result.stdout, "stderr": result.stderr}


def run_r_code(
    code: str,
    env_vars: Optional[Dict[str, str]] = None,
    timeout: int = 60,
) -> Dict[str, Any]:
    """
    Run R code directly via Rscript -e.

    This is a convenience function for running short R code snippets
    without creating a script file.

    Parameters
    ----------
    code : str
        R code to execute.
    env_vars : dict, optional
        Additional environment variables to set.
    timeout : int, default 60
        Maximum execution time in seconds.

    Returns
    -------
    Dict[str, Any]
        Parsed JSON output or stdout/stderr.

    Raises
    ------
    RNotFoundError
        If Rscript is not found.
    RExecutionError
        If the R code fails.
    RTimeoutError
        If execution times out.

    Examples
    --------
    >>> from spatialcore.r_bridge import run_r_code
    >>> result = run_r_code('cat(jsonlite::toJSON(list(x = 1, y = 2)))')
    >>> print(result)
    {'x': 1, 'y': 2}
    """
    cmd = _build_rscript_command(["-e", code])
    if not cmd:
        raise RNotFoundError(
            "Rscript not found in conda environment or PATH. "
            "Install R via: mamba install -c conda-forge r-base"
        )

    # Build environment
    env = os.environ.copy()
    env["RETICULATE_PYTHON"] = sys.executable
    if env_vars:
        env.update(env_vars)

    uses_conda_run = len(cmd) >= 4 and cmd[1] == "run"
    if uses_conda_run:
        runner_name = Path(cmd[0]).name.lower()
        if "mamba" in runner_name and not env.get("MAMBA_ROOT_PREFIX"):
            env_prefix = Path(sys.prefix)
            if "envs" in env_prefix.parts:
                idx = env_prefix.parts.index("envs")
                base_prefix = Path(*env_prefix.parts[:idx])
            else:
                base_prefix = env_prefix
            env["MAMBA_ROOT_PREFIX"] = str(base_prefix)
            logger.info(f"Setting MAMBA_ROOT_PREFIX={base_prefix}")
    if uses_conda_run:
        logger.info("Using conda/mamba run for R execution")
    else:
        logger.info("Using direct Rscript execution")
    logger.debug(f"Running R code: {code[:100]}...")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
    except subprocess.TimeoutExpired as e:
        raise RTimeoutError(f"R code timed out after {timeout} seconds.") from e

    if result.returncode != 0:
        error_msg = result.stderr.strip() if result.stderr else "Unknown error"
        if "there is no package called" in error_msg:
            error_msg = (
                f"{error_msg}\n\n"
                "Hint: Missing R package(s). Install the required R dependencies "
                "for this function."
            )
        if sys.platform == "win32" and (
            "1073741819" in error_msg or "0xC0000005" in error_msg
        ):
            error_msg = (
                f"{error_msg}\n\n"
                "Hint: R crashed on Windows (access violation). "
                "For sf-based spatial workflows, use WSL or Linux."
            )
        raise RExecutionError(f"R code failed:\n{error_msg}")

    # Try to parse JSON from output
    output = result.stdout.strip()
    if output.startswith("{") or output.startswith("["):
        try:
            return json.loads(output)
        except json.JSONDecodeError:
            pass

    return {"stdout": result.stdout, "stderr": result.stderr}
