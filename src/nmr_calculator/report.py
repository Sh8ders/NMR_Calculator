"""Complete report generation for predicted 1H NMR outputs."""

from pathlib import Path

import matplotlib.pyplot as plt

from nmr_calculator.plotting import plot_1h_spectrum
from nmr_calculator.predictor import get_predictor
from nmr_calculator.spectrum import format_peak_label, generate_peak_list
from nmr_calculator.visualization import save_molecule_image


def generate_1h_report(
    smiles: str,
    output_dir: str,
    method: str = "rules",
    database_path: str | None = None,
) -> dict[str, str]:
    """Generate a complete 1H NMR report folder for a SMILES string."""
    report_dir = Path(output_dir)
    report_dir.mkdir(parents=True, exist_ok=True)

    predictions = get_predictor(method, database_path=database_path).predict(smiles)
    peak_list = generate_peak_list(predictions)

    predictions_csv = report_dir / "predictions.csv"
    peaks_csv = report_dir / "peaks.csv"
    molecule_png = report_dir / "molecule.png"
    spectrum_png = report_dir / "spectrum.png"
    report_txt = report_dir / "report.txt"

    predictions.to_csv(predictions_csv, index=False)
    peak_list.to_csv(peaks_csv, index=False)
    save_molecule_image(smiles, str(molecule_png))

    fig, _ = plot_1h_spectrum(peak_list, smiles=smiles, output_path=str(spectrum_png))
    plt.close(fig)

    report_txt.write_text(_build_report_text(smiles, peak_list), encoding="utf-8")

    return {
        "output_dir": str(report_dir),
        "predictions_csv": str(predictions_csv),
        "peaks_csv": str(peaks_csv),
        "molecule_png": str(molecule_png),
        "spectrum_png": str(spectrum_png),
        "report_txt": str(report_txt),
    }


def _build_report_text(smiles: str, peak_list) -> str:
    peak_summary = "\n".join(
        f"- {format_peak_label(peak)}" for peak in peak_list.to_dict("records")
    )
    total_integration = int(peak_list["integration"].sum()) if not peak_list.empty else 0

    return (
        "NMR Calculator Report\n"
        f"Input SMILES: {smiles}\n\n"
        "Predicted 1H NMR Summary:\n"
        f"- Proton environments: {len(peak_list)}\n"
        f"- Total integration: {total_integration}H\n\n"
        "Peak List:\n"
        f"{peak_summary}\n\n"
        "Notes:\n"
        "These predictions are approximate rule-based estimates. "
        "They are not database-validated or ML-generated yet. "
        "Multiplicity is currently placeholder-level.\n"
    )
