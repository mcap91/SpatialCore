"""
YAML-based test runner for cross-language validation.

This module reads test definitions from tests.yaml and executes them
in both Python and R to ensure consistent results across languages.
"""

import importlib
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import yaml


def load_test_definitions(yaml_path: str = "tests.yaml") -> List[Dict[str, Any]]:
    """Load test definitions from YAML file."""
    path = Path(__file__).parent / yaml_path
    with open(path) as f:
        data = yaml.safe_load(f)
    return data.get("tests", [])


def run_python_test(test_def: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run a single Python test from definition.

    Parameters
    ----------
    test_def
        Test definition dictionary containing module, function, inputs, etc.

    Returns
    -------
    Dict containing 'passed', 'result', and optionally 'error'.
    """
    module_name = f"spatialcore.{test_def['module']}"
    func_name = test_def["function"]

    try:
        module = importlib.import_module(module_name)
        func = getattr(module, func_name)
    except (ImportError, AttributeError) as e:
        return {"passed": False, "error": f"Could not import {module_name}.{func_name}: {e}"}

    results = []
    tolerance = test_def.get("tolerance", 1e-6)

    for inp, expected in zip(test_def["inputs"], test_def["expected_outputs"]):
        try:
            # Handle AnnData loading
            if "adata" in inp and isinstance(inp["adata"], str):
                import anndata as ad
                adata_path = Path(__file__).parent / inp["adata"]
                inp = {**inp, "adata": ad.read_h5ad(adata_path)}

            result = func(**inp)
            passed = compare_outputs(result, expected, tolerance)
            results.append({"passed": passed, "result": result, "expected": expected})
        except Exception as e:
            results.append({"passed": False, "error": str(e)})

    all_passed = all(r.get("passed", False) for r in results)
    return {"passed": all_passed, "results": results}


def compare_outputs(
    result: Any,
    expected: Dict[str, Any],
    tolerance: float = 1e-6,
) -> bool:
    """
    Compare function output to expected values.

    Parameters
    ----------
    result
        Actual output from function.
    expected
        Dictionary of expected values.
    tolerance
        Tolerance for numeric comparisons.

    Returns
    -------
    True if outputs match within tolerance.
    """
    if isinstance(result, dict):
        for key, exp_val in expected.items():
            if key not in result:
                return False
            if isinstance(exp_val, (int, float)):
                if not np.isclose(result[key], exp_val, rtol=tolerance):
                    return False
            elif result[key] != exp_val:
                return False
        return True
    return False


def run_all_python_tests() -> Dict[str, Any]:
    """Run all Python tests from YAML definitions."""
    tests = load_test_definitions()
    python_tests = [t for t in tests if t["language"] == "python"]

    results = {}
    for test in python_tests:
        test_name = f"{test['module']}.{test['function']}"
        results[test_name] = run_python_test(test)

    return results


if __name__ == "__main__":
    results = run_all_python_tests()
    for name, result in results.items():
        status = "PASS" if result["passed"] else "FAIL"
        print(f"[{status}] {name}")
        if not result["passed"] and "error" in result:
            print(f"       Error: {result['error']}")
