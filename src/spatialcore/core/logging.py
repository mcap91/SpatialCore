"""Logging utilities for SpatialCore."""

import logging
import sys
from typing import Optional

_LOGGER_NAME = "spatialcore"


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """
    Get a logger instance for SpatialCore.

    Parameters
    ----------
    name
        Optional submodule name. If provided, returns logger for
        'spatialcore.{name}', otherwise returns the root spatialcore logger.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    if name:
        return logging.getLogger(f"{_LOGGER_NAME}.{name}")
    return logging.getLogger(_LOGGER_NAME)


def setup_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
) -> None:
    """
    Configure logging for SpatialCore.

    Parameters
    ----------
    level
        Logging level (e.g., logging.INFO, logging.DEBUG).
    format_string
        Custom format string. Defaults to '[%(levelname)s] %(name)s: %(message)s'.
    """
    if format_string is None:
        format_string = "[%(levelname)s] %(name)s: %(message)s"

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(format_string))

    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
