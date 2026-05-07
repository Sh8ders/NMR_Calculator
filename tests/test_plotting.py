import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nmr_calculator.plotting import plot_1h_spectrum_from_smiles


def test_ethanol_plot_creates_png(tmp_path):
    output_path = tmp_path / "ethanol_1h_nmr.png"

    fig, ax = plot_1h_spectrum_from_smiles("CCO", output_path=str(output_path))
    plt.close(fig)

    assert output_path.exists()
    assert output_path.stat().st_size > 0
    assert ax.get_xlabel() == "Chemical shift (ppm)"
    assert ax.get_ylabel() == "Relative intensity"
    assert ax.get_xlim()[0] > ax.get_xlim()[1]


def test_benzene_plot_creates_png(tmp_path):
    output_path = tmp_path / "benzene_1h_nmr.png"

    fig, _ = plot_1h_spectrum_from_smiles(
        "c1ccccc1", output_path=str(output_path)
    )
    plt.close(fig)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_raises_value_error_for_invalid_smiles(tmp_path):
    output_path = tmp_path / "invalid.png"

    with pytest.raises(ValueError, match="Invalid SMILES string"):
        plot_1h_spectrum_from_smiles("not-a-smiles", output_path=str(output_path))
