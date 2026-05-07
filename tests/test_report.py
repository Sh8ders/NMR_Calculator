import pytest

from nmr_calculator.report import generate_1h_report


EXPECTED_REPORT_FILES = {
    "predictions_csv": "predictions.csv",
    "peaks_csv": "peaks.csv",
    "molecule_png": "molecule.png",
    "spectrum_png": "spectrum.png",
    "report_txt": "report.txt",
}


def test_ethanol_report_creates_expected_files(tmp_path):
    paths = generate_1h_report("CCO", str(tmp_path))

    assert paths["output_dir"] == str(tmp_path)
    for key, filename in EXPECTED_REPORT_FILES.items():
        path = tmp_path / filename
        assert paths[key] == str(path)
        assert path.exists()
        assert path.stat().st_size > 0

    report_text = (tmp_path / "report.txt").read_text(encoding="utf-8")
    assert "Input SMILES: CCO" in report_text
    assert "Predicted 1H NMR Summary" in report_text
    assert "approximate rule-based estimates" in report_text


def test_report_accepts_rules_method(tmp_path):
    paths = generate_1h_report("CCO", str(tmp_path), method="rules")

    assert (tmp_path / "report.txt").exists()
    assert paths["report_txt"] == str(tmp_path / "report.txt")


def test_benzene_report_creates_expected_files(tmp_path):
    paths = generate_1h_report("c1ccccc1", str(tmp_path))

    for key, filename in EXPECTED_REPORT_FILES.items():
        path = tmp_path / filename
        assert paths[key] == str(path)
        assert path.exists()
        assert path.stat().st_size > 0


def test_generate_1h_report_raises_value_error_for_invalid_smiles(tmp_path):
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        generate_1h_report("not-a-smiles", str(tmp_path))
