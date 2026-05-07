from typer.testing import CliRunner

from nmr_calculator.cli import app


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
    assert "Database prediction is not implemented yet" in result.output


def test_predict_command_shows_hybrid_placeholder_message():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO", "--method", "hybrid"])

    assert result.exit_code != 0
    assert "Hybrid prediction is not implemented yet" in result.output


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
