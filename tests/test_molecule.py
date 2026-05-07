import pytest

from nmr_calculator.molecule import (
    canonicalize_smiles,
    get_hydrogen_bearing_atoms,
    get_proton_environment_groups,
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


def test_canonicalize_smiles_returns_stable_non_empty_smiles():
    assert canonicalize_smiles("CCO")
    assert canonicalize_smiles("OCC") == canonicalize_smiles("CCO")


def test_canonicalize_smiles_raises_value_error_for_invalid_smiles():
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        canonicalize_smiles("not-a-smiles")


def test_ethanol_has_expected_proton_bearing_environments():
    mol = prepare_molecule("CCO")
    environments = get_hydrogen_bearing_atoms(mol)
    labels = {environment["environment_label"] for environment in environments}

    assert len(environments) == 3
    assert "simple alkyl CH3" in labels
    assert "heteroatom-adjacent alkyl proton" in labels
    assert "alcohol exchangeable proton" in labels


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
        if environment["environment_label"] == "alkyl alpha to carbonyl"
    ]

    assert len(methyl_environments) == 2
    assert all(environment["total_hydrogens"] == 3 for environment in methyl_environments)


def test_label_proton_environment_handles_atom_without_hydrogens():
    mol = prepare_molecule("CC(=O)C")

    assert label_proton_environment(mol, 1) == "unknown proton environment"


def test_ethanol_has_three_proton_environment_groups():
    mol = prepare_molecule("CCO")
    groups = get_proton_environment_groups(mol)

    assert len(groups) == 3
    assert sum(group["proton_count"] for group in groups) == 6
    assert _find_group(groups, "simple alkyl CH3")["proton_count"] == 3
    assert _find_group(groups, "heteroatom-adjacent alkyl proton")[
        "proton_count"
    ] == 2

    exchangeable_group = _find_group(groups, "alcohol exchangeable proton")
    assert exchangeable_group["proton_count"] == 1
    assert exchangeable_group["is_exchangeable"] is True


def test_benzene_has_one_aromatic_proton_environment_group():
    mol = prepare_molecule("c1ccccc1")
    groups = get_proton_environment_groups(mol)

    assert len(groups) == 1
    assert groups[0]["proton_count"] == 6
    assert groups[0]["is_aromatic"] is True
    assert groups[0]["environment_label"] == "aromatic proton"


def test_acetone_has_one_methyl_proton_environment_group():
    mol = prepare_molecule("CC(=O)C")
    groups = get_proton_environment_groups(mol)

    assert len(groups) == 1
    assert groups[0]["proton_count"] == 6
    assert groups[0]["environment_label"] == "alkyl alpha to carbonyl"


def test_toluene_groups_include_methyl_and_aromatic_protons():
    mol = prepare_molecule("Cc1ccccc1")
    groups = get_proton_environment_groups(mol)
    labels = {group["environment_label"] for group in groups}

    assert "benzylic proton" in labels
    assert "aromatic proton" in labels
    assert _find_group(groups, "benzylic proton")["proton_count"] == 3


def test_phase_9_environment_labels_cover_common_functional_groups():
    examples = {
        "CC=O": {"alkyl alpha to carbonyl", "aldehyde proton"},
        "CC(=O)O": {
            "alkyl alpha to carbonyl",
            "carboxylic acid exchangeable proton",
        },
        "C=C": {"vinylic/alkene proton"},
        "CC#C": {"simple alkyl CH3", "terminal alkyne proton"},
    }

    for smiles, expected_labels in examples.items():
        mol = prepare_molecule(smiles)
        labels = {
            group["environment_label"]
            for group in get_proton_environment_groups(mol)
        }
        assert expected_labels <= labels


def _find_group(groups, environment_label):
    return next(
        group for group in groups if group["environment_label"] == environment_label
    )
