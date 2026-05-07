"""Plotting helpers for predicted 1H NMR spectra."""

import os
import tempfile
from pathlib import Path

matplotlib_config_dir = Path(tempfile.gettempdir()) / "nmr_calculator_matplotlib"
matplotlib_config_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_config_dir))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from nmr_calculator.spectrum import generate_1h_peak_list  # noqa: E402


def plot_1h_spectrum(
    peaks: pd.DataFrame,
    smiles: str | None = None,
    output_path: str | None = None,
):
    """Plot a simple simulated 1H NMR spectrum from a predicted peak list."""
    fig, ax = plt.subplots(figsize=(10, 4.8))

    max_integration = max(peaks["integration"].max(), 1) if not peaks.empty else 1
    for peak in peaks.to_dict("records"):
        ppm = float(peak["ppm"])
        integration = int(peak["integration"])
        height = integration / max_integration

        ax.vlines(ppm, 0, height, color="#1f4e79", linewidth=2.2)
        ax.text(
            ppm,
            height + 0.04,
            f"{ppm:.1f} ppm",
            ha="center",
            va="bottom",
            fontsize=9,
        )
        ax.text(
            ppm,
            max(height * 0.5, 0.08),
            f"{integration}H",
            ha="center",
            va="center",
            fontsize=9,
            rotation=90,
        )

    title = "Predicted 1H NMR Spectrum"
    if smiles:
        title = f"{title}\nSMILES: {smiles}"

    ax.set_title(title)
    ax.set_xlabel("Chemical shift (ppm)")
    ax.set_ylabel("Relative intensity")
    ax.set_xlim(12, 0)
    ax.set_ylim(0, 1.25)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()

    if output_path:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_file, dpi=200, format="png")

    return fig, ax


def plot_1h_spectrum_from_smiles(
    smiles: str,
    output_path: str | None = None,
    method: str = "rules",
    database_path: str | None = None,
):
    """Generate and plot a simple predicted 1H NMR spectrum from SMILES."""
    peak_list = generate_1h_peak_list(
        smiles, method=method, database_path=database_path
    )
    return plot_1h_spectrum(peak_list, smiles=smiles, output_path=output_path)
