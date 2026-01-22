# R Bridge

The R Bridge executes R code and scripts from Python via subprocess while
ensuring the correct Python environment is available to R (reticulate).

## Modules

- `spatialcore.r_bridge`: Public helpers and exceptions.
- `spatialcore.r_bridge.subprocess_runner`: Low-level subprocess execution,
  environment detection, and error handling.

## Environment Handling

### Windows (dev)

When running inside a conda or mamba environment on Windows, R is executed via
`conda/mamba run -n <env> Rscript ...` to ensure DLLs and libraries are found.
If conda/mamba is not available, the bridge uses direct `Rscript`.

Windows support is best-effort. If you see access violations (for example
`1073741819`), treat Windows as unsupported for sf-based spatial workflows and
run those steps in WSL or Linux instead.

### Linux (production)

On Linux, `Rscript` is used directly when available. If running inside a conda
or mamba environment and the executable is found, `conda/mamba run` may be used
to ensure environment consistency.
If `mamba run` is used and `MAMBA_ROOT_PREFIX` is not set, the bridge sets it
automatically based on the active environment prefix.

## Errors and Diagnostics

- Missing R: raises `RNotFoundError`.
- Execution failure: raises `RExecutionError` with stderr/stdout details.
- Timeout: raises `RTimeoutError`.
- Missing R packages: detected from R error output and surfaced with hints.

Use `check_r_available()` to verify R is available before calling functions
that rely on R.

## Manual Check

To confirm that the bridge can execute R and that sf settings are applied,
run the manual test:

```
SPATIALCORE_MANUAL_R=1 python -m pytest tests/spatialcore_celltypist_test_suite/test_r_bridge.py -k sf_use_s2_windows_manual -v
```
