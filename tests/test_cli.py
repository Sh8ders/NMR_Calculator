from typer.testing import CliRunner

from nmr_calculator.cli import app


def test_predict_command_reports_phase_3_status():
    runner = CliRunner()

    result = runner.invoke(app, ["predict", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Phase 3 proton environment grouping ready." in result.output
    assert "1H NMR prediction will be implemented in later phases." in result.output
    old_carbon_nmr_phrase = "13" + "C NMR"
    assert old_carbon_nmr_phrase not in result.output


def test_inspect_command_reports_proton_environments():
    runner = CliRunner()

    result = runner.invoke(app, ["inspect", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Hydrogen-bearing atoms:" in result.output
    assert "alkyl CH3" in result.output
    assert "heteroatom-adjacent alkyl proton" in result.output
    assert "alcohol/amine/thiol exchangeable proton" in result.output


def test_groups_command_reports_grouped_proton_environments():
    runner = CliRunner()

    result = runner.invoke(app, ["groups", "CCO"])

    assert result.exit_code == 0
    assert "Input SMILES: CCO" in result.output
    assert "Proton environment groups:" in result.output
    assert "3H" in result.output
    assert "2H" in result.output
    assert "1H" in result.output
