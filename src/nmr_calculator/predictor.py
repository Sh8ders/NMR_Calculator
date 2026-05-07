"""Rule-based 1H NMR proton chemical shift prediction."""

from typing import NamedTuple

import pandas as pd
from rdkit.Chem.rdchem import Atom, BondType, Mol

from nmr_calculator.database import filter_database_by_smiles, load_shift_database_csv
from nmr_calculator.molecule import (
    get_attached_heteroatom_type,
    get_proton_environment_groups,
    is_aldehyde_proton,
    is_alpha_to_carbonyl,
    is_allylic,
    is_benzylic,
    is_carboxylic_acid_proton,
    is_terminal_alkyne_proton,
    is_vinylic,
    prepare_molecule,
)


class ShiftRule(NamedTuple):
    predicted_shift_ppm: float
    shift_range_ppm: tuple[float, float]
    confidence: str
    notes: str


class BaseProtonPredictor:
    """Base interface for 1H NMR proton shift predictors."""

    name = "base"

    def predict(self, smiles: str) -> pd.DataFrame:
        """Predict 1H NMR shifts for a SMILES string."""
        raise NotImplementedError


class RuleBasedProtonPredictor(BaseProtonPredictor):
    """Rule-based 1H NMR proton shift predictor."""

    name = "rules"

    def predict(self, smiles: str) -> pd.DataFrame:
        """Predict approximate 1H NMR shifts using local rules."""
        return _predict_1h_shifts_with_rules(smiles)


class DatabaseProtonPredictor(BaseProtonPredictor):
    """Exact-match local CSV database predictor."""

    name = "database"

    def __init__(self, database_path: str):
        self.database_path = database_path

    def predict(self, smiles: str) -> pd.DataFrame:
        """Predict shifts from exact local database molecule matches."""
        database = load_shift_database_csv(self.database_path)
        matching_records = filter_database_by_smiles(database, smiles)
        if matching_records.empty:
            raise ValueError(
                "No database records found for this molecule. "
                "Try --method hybrid to fall back to rules."
            )

        mol = prepare_molecule(smiles)
        groups = get_proton_environment_groups(mol)
        rows = []
        for group in groups:
            group_records = matching_records[
                matching_records["atom_index"].isin(group["atom_indices"])
            ]
            if group_records.empty:
                raise ValueError(
                    "Database records do not cover all proton environments. "
                    "Try --method hybrid."
                )

            shift_values = group_records["shift_ppm"]
            rows.append(
                {
                    "group_id": group["group_id"],
                    "atom_indices": group["atom_indices"],
                    "proton_count": group["proton_count"],
                    "environment_label": group["environment_label"],
                    "predicted_shift_ppm": float(shift_values.mean()),
                    "shift_min_ppm": float(shift_values.min()),
                    "shift_max_ppm": float(shift_values.max()),
                    "confidence": "high",
                    "prediction_method": "database",
                    "notes": "Exact molecule database match.",
                    "database_match_count": int(len(group_records)),
                    "database_source": ", ".join(
                        sorted(set(group_records["source"].astype(str)))
                    ),
                }
            )

        return pd.DataFrame(
            rows,
            columns=[
                "group_id",
                "atom_indices",
                "proton_count",
                "environment_label",
                "predicted_shift_ppm",
                "shift_min_ppm",
                "shift_max_ppm",
                "confidence",
                "prediction_method",
                "notes",
                "database_match_count",
                "database_source",
            ],
        )


class HybridProtonPredictor(BaseProtonPredictor):
    """Placeholder for future database-first hybrid prediction."""

    name = "hybrid"

    def predict(self, smiles: str) -> pd.DataFrame:
        """Raise until hybrid prediction is implemented."""
        raise NotImplementedError(
            "Hybrid prediction is not implemented yet. "
            "It will be added after database prediction."
        )


SUPPORTED_PREDICTION_METHODS = {
    "rules": RuleBasedProtonPredictor,
    "rule": RuleBasedProtonPredictor,
    "rule-based": RuleBasedProtonPredictor,
    "database": DatabaseProtonPredictor,
    "db": DatabaseProtonPredictor,
    "hybrid": HybridProtonPredictor,
}


def get_predictor(method: str, database_path: str | None = None) -> BaseProtonPredictor:
    """Return a proton predictor for a supported method string."""
    normalized_method = method.strip().lower()
    predictor_class = SUPPORTED_PREDICTION_METHODS.get(normalized_method)
    if predictor_class is None:
        supported_methods = ", ".join(sorted(SUPPORTED_PREDICTION_METHODS))
        raise ValueError(
            f"Unknown prediction method: {method!r}. "
            f"Supported methods: {supported_methods}."
        )

    if predictor_class is DatabaseProtonPredictor:
        if database_path is None:
            raise ValueError("Database path is required for database prediction.")
        return predictor_class(database_path)

    return predictor_class()


def predict_1h_shifts(
    smiles: str, method: str = "rules", database_path: str | None = None
) -> pd.DataFrame:
    """Predict approximate 1H NMR shifts for grouped proton environments."""
    return get_predictor(method, database_path=database_path).predict(smiles)


def _predict_1h_shifts_with_rules(smiles: str) -> pd.DataFrame:
    """Predict approximate 1H NMR shifts with the rule-based implementation."""
    mol = prepare_molecule(smiles)
    predictions = [
        predict_proton_shift_for_group(mol, group)
        for group in get_proton_environment_groups(mol)
    ]

    rows = []
    for prediction in predictions:
        shift_min_ppm, shift_max_ppm = prediction["shift_range_ppm"]
        rows.append(
            {
                "group_id": prediction["group_id"],
                "atom_indices": prediction["atom_indices"],
                "proton_count": prediction["proton_count"],
                "environment_label": prediction["environment_label"],
                "predicted_shift_ppm": prediction["predicted_shift_ppm"],
                "shift_min_ppm": shift_min_ppm,
                "shift_max_ppm": shift_max_ppm,
                "confidence": prediction["confidence"],
                "prediction_method": prediction["prediction_method"],
                "notes": prediction["notes"],
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "group_id",
            "atom_indices",
            "proton_count",
            "environment_label",
            "predicted_shift_ppm",
            "shift_min_ppm",
            "shift_max_ppm",
            "confidence",
            "prediction_method",
            "notes",
        ],
    )


def predict_proton_shift_for_group(mol: Mol, group: dict) -> dict:
    """Predict an approximate 1H NMR shift for a proton environment group."""
    rule = get_shift_rule_for_group(mol, group)
    return {
        "group_id": group["group_id"],
        "atom_indices": group["atom_indices"],
        "proton_count": group["proton_count"],
        "environment_label": group["environment_label"],
        "predicted_shift_ppm": rule.predicted_shift_ppm,
        "shift_range_ppm": rule.shift_range_ppm,
        "confidence": rule.confidence,
        "prediction_method": "phase_9_rule_based",
        "notes": rule.notes,
    }


def get_shift_rule_for_group(mol: Mol, group: dict) -> ShiftRule:
    """Select the first applicable approximate shift rule for a proton group."""
    environment_label = str(group["environment_label"])
    atom_indices = [int(atom_index) for atom_index in group["atom_indices"]]

    if "aldehyde proton" in environment_label or any(
        is_aldehyde_proton(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            9.7, (9.0, 10.5), "medium", "Matched aldehyde proton rule."
        )

    if "aromatic proton" in environment_label:
        return ShiftRule(
            7.2, (6.0, 8.5), "medium", "Matched aromatic proton rule."
        )

    if "vinylic/alkene proton" in environment_label or any(
        is_vinylic(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            5.5, (4.5, 6.5), "medium", "Matched vinylic alkene proton rule."
        )

    if "carboxylic acid exchangeable proton" in environment_label or any(
        is_carboxylic_acid_proton(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            11.5,
            (10.0, 13.0),
            "low",
            "Matched exchangeable carboxylic acid proton rule.",
        )

    if "alcohol exchangeable proton" in environment_label:
        return ShiftRule(
            2.5,
            (0.5, 5.5),
            "low",
            "Matched exchangeable alcohol proton rule.",
        )

    if "amine exchangeable proton" in environment_label:
        return ShiftRule(
            2.0,
            (0.5, 5.0),
            "low",
            "Matched exchangeable amine proton rule.",
        )

    if "thiol exchangeable proton" in environment_label:
        return ShiftRule(
            2.0,
            (1.0, 4.0),
            "low",
            "Matched exchangeable thiol proton rule.",
        )

    if "alkyl alpha to carbonyl" in environment_label or any(
        is_alpha_to_carbonyl(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            2.2,
            (2.0, 3.0),
            "medium",
            "Matched alpha-to-carbonyl alkyl proton rule.",
        )

    if "heteroatom-adjacent alkyl proton" in environment_label:
        return ShiftRule(
            3.5,
            (3.0, 4.5),
            "medium",
            "Matched heteroatom-adjacent alkyl proton rule.",
        )

    if "benzylic proton" in environment_label or any(
        is_benzylic(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            2.5, (2.2, 3.0), "medium", "Matched benzylic proton rule."
        )

    if "allylic proton" in environment_label or any(
        is_allylic(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            2.0, (1.6, 2.4), "medium", "Matched allylic proton rule."
        )

    if "terminal alkyne proton" in environment_label or any(
        is_terminal_alkyne_proton(mol, atom_index) for atom_index in atom_indices
    ):
        return ShiftRule(
            2.5, (2.0, 3.0), "low", "Matched terminal alkyne proton rule."
        )

    if environment_label.startswith("simple alkyl CH"):
        return ShiftRule(
            1.2, (0.8, 1.8), "medium", "Matched simple alkyl proton rule."
        )

    if any(get_attached_heteroatom_type(mol, atom_index) for atom_index in atom_indices):
        return ShiftRule(
            3.5,
            (3.0, 4.5),
            "low",
            "Matched heteroatom-adjacent fallback rule.",
        )

    return ShiftRule(
        1.5, (0.5, 8.5), "low", "Matched fallback unknown environment rule."
    )


def is_adjacent_to_carbonyl(mol: Mol, atom_index: int) -> bool:
    """Return whether an atom is bonded to a carbonyl carbon."""
    return is_alpha_to_carbonyl(mol, atom_index)


def is_adjacent_to_aromatic_or_alkene(mol: Mol, atom_index: int) -> bool:
    """Return whether an atom is next to an aromatic atom or C=C pi system."""
    atom = mol.GetAtomWithIdx(atom_index)
    if atom.GetAtomicNum() == 1 or atom.GetIsAromatic():
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:
            continue
        bond = mol.GetBondBetweenAtoms(atom_index, neighbor.GetIdx())
        if neighbor.GetIsAromatic() or bond.GetBondType() == BondType.DOUBLE:
            return True
    return False


def _is_carbonyl_carbon(atom: Atom) -> bool:
    if atom.GetSymbol() != "C":
        return False

    mol = atom.GetOwningMol()
    return any(
        neighbor.GetSymbol() == "O"
        and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()
        == BondType.DOUBLE
        for neighbor in atom.GetNeighbors()
    )


def _is_carboxylic_acid_oxygen(atom: Atom) -> bool:
    if atom.GetSymbol() != "O":
        return False

    return any(_is_carbonyl_carbon(neighbor) for neighbor in atom.GetNeighbors())


def _is_terminal_alkyne_carbon(atom: Atom) -> bool:
    if atom.GetSymbol() != "C":
        return False

    mol = atom.GetOwningMol()
    return any(
        neighbor.GetSymbol() == "C"
        and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()
        == BondType.TRIPLE
        for neighbor in atom.GetNeighbors()
    )
