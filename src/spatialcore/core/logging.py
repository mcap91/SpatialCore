"""Logging utilities for SpatialCore."""

import logging
import sys
from pathlib import Path
from typing import Optional, Union

_LOGGER_NAME = "spatialcore"
_DEFAULT_FORMAT = "[%(levelname)s] %(name)s: %(message)s"
_INITIALIZED = False


def _ensure_initialized() -> None:
    """
    Ensure the root spatialcore logger has at least one handler.

    Called automatically by get_logger() to prevent silent log drops.
    """
    global _INITIALIZED
    if _INITIALIZED:
        return

    root_logger = logging.getLogger(_LOGGER_NAME)

    # Only add default handler if none exist
    if not root_logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
        handler.setLevel(logging.INFO)
        root_logger.addHandler(handler)
        root_logger.setLevel(logging.INFO)
        root_logger.propagate = False

    _INITIALIZED = True


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """
    Get a logger instance for SpatialCore.

    Automatically initializes a default stdout handler if none exists,
    preventing silent log message drops.

    Parameters
    ----------
    name
        Optional submodule name. If provided, returns logger for
        'spatialcore.{name}', otherwise returns the root spatialcore logger.

    Returns
    -------
    logging.Logger
        Configured logger instance with at least one handler.
    """
    _ensure_initialized()

    if name:
        return logging.getLogger(f"{_LOGGER_NAME}.{name}")
    return logging.getLogger(_LOGGER_NAME)


def setup_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
) -> None:
    """
    Configure logging for SpatialCore.

    Clears existing handlers and configures fresh logging. Safe to call
    multiple times (idempotent after first call with same parameters).

    Parameters
    ----------
    level
        Logging level (e.g., logging.INFO, logging.DEBUG).
    format_string
        Custom format string. Defaults to '[%(levelname)s] %(name)s: %(message)s'.
    """
    global _INITIALIZED

    if format_string is None:
        format_string = _DEFAULT_FORMAT

    logger = logging.getLogger(_LOGGER_NAME)

    # Clear existing handlers to prevent duplicates
    logger.handlers.clear()

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(format_string))
    handler.setLevel(level)

    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False

    _INITIALIZED = True


def setup_file_logging(
    log_path: Union[str, Path],
    level: int = logging.INFO,
    format_string: Optional[str] = None,
) -> logging.FileHandler:
    """
    Add file handler to SpatialCore logger.

    Parameters
    ----------
    log_path
        Path to log file.
    level
        Logging level.
    format_string
        Custom format string.

    Returns
    -------
    logging.FileHandler
        The file handler (can be removed later with logger.removeHandler).
    """
    if format_string is None:
        format_string = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    handler = logging.FileHandler(log_path, mode="w")
    handler.setFormatter(logging.Formatter(format_string))
    handler.setLevel(level)

    logger = logging.getLogger(_LOGGER_NAME)
    logger.addHandler(handler)

    return handler
