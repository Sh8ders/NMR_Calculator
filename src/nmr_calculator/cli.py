"""Command-line interface for NMR Calculator."""

import typer

from nmr_calculator.molecule import get_hydrogen_bearing_atoms, prepare_molecule

app = typer.Typer(help="1H NMR prediction command-line tools.")


@app.callback()
def main() -> None:
    """Run NMR Calculator commands."""


@app.command()
def predict(smiles: str) -> None:
    """Accept a SMILES string and report Phase 2 scaffold status."""
    typer.echo(f"Input SMILES: {smiles}")
    typer.echo("Phase 2 molecule parsing ready.")
    typer.echo("1H NMR prediction will be implemented in later phases.")


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


if __name__ == "__main__":
    app()
