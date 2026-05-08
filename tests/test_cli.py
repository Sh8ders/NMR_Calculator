from typer.testing import CliRunner

from nmr_calculator.cli import app
from nmr_calculator.database import load_shift_database_csv


def test_predict_command_reports_rule_based_shifts():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Predicted 1H NMR chemical shifts:" in result.output
    assert "ppm" in result.output
    assert "3H" in result.output
    assert "2H" in result.output
    assert "1H" in result.output
    old_carbon_nmr_phrase = "13" + "C NMR"
    assert old_carbon_nmr_phrase not in result.output


def test_predict_command_accepts_rules_method():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO", "--method", "rules"])

    assert result.exit_code == 0
    assert "Predicted 1H NMR chemical shifts:" in result.output


def test_inspect_command_reports_proton_environments():
    runner = CliRunner()

    result = runner.invoke(app, ["inspect", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Hydrogen-bearing atoms:" in result.output
    assert "alkyl CH3" in result.output
    assert "heteroatom-adjacent alkyl proton" in result.output
    assert "alcohol exchangeable proton" in result.output


def test_groups_command_reports_grouped_proton_environments():
    runner = CliRunner()

    result = runner.invoke(app, ["groups", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Proton environment groups:" in result.output
    assert "3H" in result.output
    assert "2H" in result.output
    assert "1H" in result.output


def test_peaks_command_reports_predicted_peak_list():
    runner = CliRunner()

    result = runner.invoke(app, ["peaks", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Predicted 1H NMR peak list:" in result.output
    assert "ppm" in result.output
    assert "3H" in result.output
    assert "2H" in result.output
    assert "1H" in result.output


def test_peaks_command_accepts_rules_method():
    runner = CliRunner()

    result = runner.invoke(app, ["peaks", "CCO", "--method", "rules"])

    assert result.exit_code == 0
    assert "Predicted 1H NMR peak list:" in result.output


def test_plot_command_saves_predicted_spectrum(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "ethanol_1h_nmr.png"

    result = runner.invoke(
        app, ["plot", "CCO", "--output", str(output_path)]
    )

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Saved predicted 1H NMR spectrum" in result.output
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_command_accepts_rules_method(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "ethanol_1h_nmr.png"

    result = runner.invoke(
        app, ["plot", "CCO", "--output", str(output_path), "--method", "rules"]
    )

    assert result.exit_code == 0
    assert "Saved predicted 1H NMR spectrum" in result.output
    assert output_path.exists()


def test_molecule_image_command_saves_png(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "ethanol_structure.png"

    result = runner.invoke(
        app, ["molecule-image", "CCO", "--output", str(output_path)]
    )

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Saved molecule image" in result.output
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_report_command_generates_report_files(tmp_path):
    runner = CliRunner()
    output_dir = tmp_path / "ethanol_report"

    result = runner.invoke(
        app, ["report", "CCO", "--output-dir", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Generated 1H NMR report" in result.output
    for filename in [
        "predictions.csv",
        "peaks.csv",
        "molecule.png",
        "spectrum.png",
        "report.txt",
    ]:
        assert filename in result.output
        assert (output_dir / filename).exists()


def test_report_command_accepts_rules_method(tmp_path):
    runner = CliRunner()
    output_dir = tmp_path / "ethanol_report"

    result = runner.invoke(
        app, ["report", "CCO", "--output-dir", str(output_dir), "--method", "rules"]
    )

    assert result.exit_code == 0
    assert "Generated 1H NMR report" in result.output
    assert (output_dir / "report.txt").exists()


def test_predict_command_shows_database_placeholder_message():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO", "--method", "database"])

    assert result.exit_code != 0
    assert "Database path is required" in result.output


def test_predict_command_shows_hybrid_placeholder_message():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO", "--method", "hybrid"])

    assert result.exit_code != 0
    assert "Hybrid prediction is not implemented yet" in result.output


def test_predict_command_accepts_database_method_with_database_path():
    runner = CliRunner()

    result = runner.invoke(
        app,
        [
            "predict",
            "CCO",
            "--method",
            "database",
            "--database",
            "tests/fixtures/test_1h_shift_database.csv",
        ],
    )

    assert result.exit_code == 0
    assert "database" in result.output
    assert "ppm" in result.output


def test_peaks_command_accepts_database_method_with_database_path():
    runner = CliRunner()

    result = runner.invoke(
        app,
        [
            "peaks",
            "CCO",
            "--method",
            "database",
            "--database",
            "tests/fixtures/test_1h_shift_database.csv",
        ],
    )

    assert result.exit_code == 0
    assert "Predicted 1H NMR peak list:" in result.output


def test_plot_command_accepts_database_method_with_database_path(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "ethanol_database.png"

    result = runner.invoke(
        app,
        [
            "plot",
            "CCO",
            "--output",
            str(output_path),
            "--method",
            "database",
            "--database",
            "tests/fixtures/test_1h_shift_database.csv",
        ],
    )

    assert result.exit_code == 0
    assert output_path.exists()


def test_report_command_accepts_database_method_with_database_path(tmp_path):
    runner = CliRunner()
    output_dir = tmp_path / "ethanol_database_report"

    result = runner.invoke(
        app,
        [
            "report",
            "CCO",
            "--output-dir",
            str(output_dir),
            "--method",
            "database",
            "--database",
            "tests/fixtures/test_1h_shift_database.csv",
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / "report.txt").exists()


def test_database_check_command_validates_fixture():
    runner = CliRunner()

    result = runner.invoke(
        app, ["database-check", "tests/fixtures/test_1h_shift_database.csv"]
    )

    assert result.exit_code == 0
    assert "Records: 6" in result.output
    assert "Unique molecules: 3" in result.output
    assert "Database validation passed." in result.output


def test_database_cache_command_writes_cache(tmp_path):
    runner = CliRunner()
    cache_path = tmp_path / "test_cache.csv"

    result = runner.invoke(
        app,
        [
            "database-cache",
            "tests/fixtures/test_1h_shift_database.csv",
            "--output",
            str(cache_path),
        ],
    )

    assert result.exit_code == 0
    assert "Saved normalized database cache" in result.output
    assert cache_path.exists()
    assert cache_path.stat().st_size > 0


def test_nmrshiftdb2_import_command_writes_compatible_csv(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "nmrshiftdb2_1h.csv"

    result = runner.invoke(
        app,
        [
            "nmrshiftdb2-import",
            "tests/fixtures/tiny_nmrshiftdb2_1h.sdf",
            "--output",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert "Input path: tests/fixtures/tiny_nmrshiftdb2_1h.sdf" in result.output
    assert "Molecules scanned: 2" in result.output
    assert "Molecules with 1H-like properties: 2" in result.output
    assert "Records extracted: 5" in result.output
    assert "Skipped records: 0" in result.output
    assert f"Output path: {output_path}" in result.output
    assert "Extracted 1H shift records: 5" in result.output
    assert "Unique molecules: 2" in result.output
    assert output_path.exists()

    imported_database = load_shift_database_csv(str(output_path))
    assert len(imported_database) == 5


def test_nmrshiftdb2_inspect_command_reports_property_previews():
    runner = CliRunner()

    result = runner.invoke(
        app,
        [
            "nmrshiftdb2-inspect",
            "tests/fixtures/tiny_nmrshiftdb2_1h.sdf",
            "--max-molecules",
            "1",
        ],
    )

    assert result.exit_code == 0
    assert "Input path: tests/fixtures/tiny_nmrshiftdb2_1h.sdf" in result.output
    assert "Molecules inspected: 1" in result.output
    assert "Molecule 0" in result.output
    assert "Name: ethanol" in result.output
    assert "Canonical SMILES: CCO" in result.output
    assert "NMREDATA_1D_1H: 0,1.23,CH3" in result.output
