import importlib
from pathlib import Path


def test_package_imports_successfully():
    import nmr_calculator

    assert nmr_calculator.__version__ == "0.12.0"


def test_placeholder_modules_import_successfully():
    module_names = [
        "nmr_calculator.database",
        "nmr_calculator.molecule",
        "nmr_calculator.plotting",
        "nmr_calculator.predictor",
        "nmr_calculator.report",
        "nmr_calculator.spectrum",
        "nmr_calculator.visualization",
    ]

    for module_name in module_names:
        assert importlib.import_module(module_name)


def test_pyproject_exists():
    assert Path("pyproject.toml").exists()
