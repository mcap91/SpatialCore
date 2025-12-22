"""Caching utilities for intermediate results."""

import hashlib
from functools import wraps
from pathlib import Path
from typing import Callable, Optional, Union

import anndata as ad

from spatialcore.core.logging import get_logger

logger = get_logger("cache")

_CACHE_DIR = Path(".cache")


def get_cache_path(name: str, cache_dir: Optional[Union[str, Path]] = None) -> Path:
    """
    Get path for a cached file.

    Parameters
    ----------
    name
        Name of the cached file (without extension).
    cache_dir
        Optional cache directory. Defaults to '.cache/'.

    Returns
    -------
    Path
        Full path to the cache file.
    """
    cache_dir = Path(cache_dir) if cache_dir else _CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{name}.h5ad"


def cache_result(
    name: Optional[str] = None,
    cache_dir: Optional[Union[str, Path]] = None,
) -> Callable:
    """
    Decorator to cache AnnData results.

    Parameters
    ----------
    name
        Optional cache name. If not provided, uses function name + hash of args.
    cache_dir
        Optional cache directory.

    Returns
    -------
    Callable
        Decorated function.

    Examples
    --------
    >>> @cache_result(name="my_analysis")
    ... def expensive_computation(adata):
    ...     # ... computation ...
    ...     return adata
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            cache_name = name or _generate_cache_name(func.__name__, args, kwargs)
            cache_path = get_cache_path(cache_name, cache_dir)

            if cache_path.exists():
                logger.info(f"Loading cached result from {cache_path}")
                return ad.read_h5ad(cache_path)

            result = func(*args, **kwargs)

            if isinstance(result, ad.AnnData):
                logger.info(f"Caching result to {cache_path}")
                result.write_h5ad(cache_path)

            return result
        return wrapper
    return decorator


def clear_cache(cache_dir: Optional[Union[str, Path]] = None) -> int:
    """
    Clear all cached files.

    Parameters
    ----------
    cache_dir
        Optional cache directory. Defaults to '.cache/'.

    Returns
    -------
    int
        Number of files removed.
    """
    cache_dir = Path(cache_dir) if cache_dir else _CACHE_DIR
    if not cache_dir.exists():
        return 0

    count = 0
    for f in cache_dir.glob("*.h5ad"):
        f.unlink()
        count += 1

    logger.info(f"Cleared {count} cached files from {cache_dir}")
    return count


def _generate_cache_name(func_name: str, args: tuple, kwargs: dict) -> str:
    """Generate a unique cache name based on function and arguments."""
    hasher = hashlib.md5()
    hasher.update(func_name.encode())
    hasher.update(str(args).encode())
    hasher.update(str(sorted(kwargs.items())).encode())
    return f"{func_name}_{hasher.hexdigest()[:8]}"
