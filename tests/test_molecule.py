import pytest

from nmr_calculator.molecule import (
    get_hydrogen_bearing_atoms,
    label_proton_environment,
    parse_molecule,
    prepare_molecule,
)


def test_valid_smiles_parses_successfully():
    mol = parse_molecule("CCO")

    assert mol.GetNumAtoms() == 3


def test_invalid_smiles_raises_value_error():
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        parse_molecule("not-a-smiles")


def test_ethanol_has_expected_proton_bearing_environments():
    mol = prepare_molecule("CCO")
    environments = get_hydrogen_bearing_atoms(mol)
    labels = {environment["environment_label"] for environment in environments}

    assert len(environments) == 3
    assert "alkyl CH3" in labels
    assert "heteroatom-adjacent alkyl proton" in labels
    assert "alcohol/amine/thiol exchangeable proton" in labels


def test_benzene_has_aromatic_proton_environments():
    mol = prepare_molecule("c1ccccc1")
    environments = get_hydrogen_bearing_atoms(mol)

    assert len(environments) == 6
    assert {
        environment["environment_label"] for environment in environments
    } == {"aromatic proton"}


def test_acetone_has_methyl_proton_environments():
    mol = prepare_molecule("CC(=O)C")
    environments = get_hydrogen_bearing_atoms(mol)
    methyl_environments = [
        environment
        for environment in environments
        if environment["environment_label"] == "alkyl CH3"
    ]

    assert len(methyl_environments) == 2
    assert all(environment["total_hydrogens"] == 3 for environment in methyl_environments)


def test_label_proton_environment_handles_atom_without_hydrogens():
    mol = prepare_molecule("CC(=O)C")

    assert label_proton_environment(mol, 1) == "unknown proton environment"
