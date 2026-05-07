"""Command-line interface for NMR Calculator."""

import typer

app = typer.Typer(help="13C NMR prediction command-line tools.")


@app.callback()
def main() -> None:
    """Run NMR Calculator commands."""


@app.command()
def predict(smiles: str) -> None:
    """Accept a SMILES string and report Phase 1 scaffold status."""
    typer.echo(f"Input SMILES: {smiles}")
    typer.echo("Phase 1 project setup complete.")
    typer.echo("13C NMR prediction will be implemented in later phases.")


if __name__ == "__main__":
    app()
