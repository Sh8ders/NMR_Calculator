from typer.testing import CliRunner

from nmr_calculator.cli import app


def test_predict_command_reports_phase_1_status():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Phase 1 project setup complete." in result.output
    assert "13C NMR prediction will be implemented in later phases." in result.output

