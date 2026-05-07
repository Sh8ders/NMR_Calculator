"""Command-line interface for NMR Calculator."""

import typer

from nmr_calculator.molecule import (
    get_hydrogen_bearing_atoms,
    get_proton_environment_groups,
    prepare_molecule,
)
from nmr_calculator.plotting import plot_1h_spectrum_from_smiles
from nmr_calculator.predictor import predict_1h_shifts
from nmr_calculator.spectrum import format_peak_label, generate_1h_peak_list
from nmr_calculator.visualization import save_molecule_image

app = typer.Typer(help="1H NMR prediction command-line tools.")


@app.callback()
def main() -> None:
    """Run NMR Calculator commands."""


@app.command()
def predict(smiles: str) -> None:
    """Predict approximate rule-based 1H NMR chemical shifts."""
    predictions = predict_1h_shifts(smiles)

    typer.echo(f"Input SMILES: {smiles}")
    typer.echo("Predicted 1H NMR chemical shifts:")
    for prediction in predictions.to_dict("records"):
        typer.echo(
            "  "
            f"Group {prediction['group_id']}: "
            f"atoms {prediction['atom_indices']}, "
            f"{prediction['proton_count']}H, "
            f"{prediction['environment_label']}, "
            f"~{prediction['predicted_shift_ppm']:.1f} ppm "
            f"({prediction['shift_min_ppm']:.1f}-"
            f"{prediction['shift_max_ppm']:.1f} ppm)"
        )


@app.command()
def inspect(smiles: str) -> None:
    """Inspect hydrogen-bearing atoms and simple proton environments."""
    mol = prepare_molecule(smiles)
    environments = get_hydrogen_bearing_atoms(mol)

    typer.echo(f"Input SMILES: {smiles}")
    typer.echo("Hydrogen-bearing atoms:")
    for environment in environments:
        typer.echo(
            "  "
            f"Atom {environment['atom_index']} "
            f"({environment['atom_symbol']}H{environment['total_hydrogens']}): "
            f"{environment['environment_label']}"
        )


@app.command()
def groups(smiles: str) -> None:
    """Print grouped equivalent proton environments."""
    mol = prepare_molecule(smiles)
    proton_groups = get_proton_environment_groups(mol)

    typer.echo(f"Input SMILES: {smiles}")
    typer.echo("Proton environment groups:")
    for group in proton_groups:
        typer.echo(
            "  "
            f"Group {group['group_id']}: "
            f"atoms {group['atom_indices']}, "
            f"{group['proton_count']}H, "
            f"{group['environment_label']}"
        )


@app.command()
def peaks(smiles: str) -> None:
    """Print a sorted predicted 1H NMR peak list."""
    peak_list = generate_1h_peak_list(smiles)

    typer.echo(f"Input SMILES: {smiles}")
    typer.echo("Predicted 1H NMR peak list:")
    for peak in peak_list.to_dict("records"):
        typer.echo(f"  Peak {peak['peak_id']}: {format_peak_label(peak)}")


@app.command()
def plot(
    smiles: str,
    output: str = typer.Option(..., "--output", "-o", help="PNG output path."),
) -> None:
    """Save a simple predicted 1H NMR spectrum plot."""
    fig, _ = plot_1h_spectrum_from_smiles(smiles, output_path=output)
    # Close the figure after saving so repeated CLI/test calls do not leak figures.
    import matplotlib.pyplot as plt

    plt.close(fig)

    typer.echo(f"Input SMILES: {smiles}")
    typer.echo(f"Saved predicted 1H NMR spectrum to {output}")


@app.command()
def molecule_image(
    smiles: str,
    output: str = typer.Option(..., "--output", "-o", help="PNG output path."),
    show_atom_indices: bool = typer.Option(
        True,
        "--atom-indices/--no-atom-indices",
        help="Show RDKit atom index labels.",
    ),
) -> None:
    """Save a molecule structure image."""
    saved_path = save_molecule_image(
        smiles, output, show_atom_indices=show_atom_indices
    )

    typer.echo(f"Input SMILES: {smiles}")
    typer.echo(f"Saved molecule image to {saved_path}")


if __name__ == "__main__":
    app()
