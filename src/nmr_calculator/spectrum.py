"""1H NMR peak list data structure utilities."""

import pandas as pd

from nmr_calculator.predictor import get_predictor


def generate_peak_list(predictions: pd.DataFrame) -> pd.DataFrame:
    """Convert predicted 1H NMR shift groups into a sorted peak list."""
    peaks = []
    for prediction in predictions.to_dict("records"):
        peaks.append(
            {
                "group_id": prediction["group_id"],
                "ppm": prediction["predicted_shift_ppm"],
                "shift_min_ppm": prediction["shift_min_ppm"],
                "shift_max_ppm": prediction["shift_max_ppm"],
                "integration": prediction["proton_count"],
                "multiplicity": infer_basic_multiplicity(prediction),
                "environment_label": prediction["environment_label"],
                "atom_indices": prediction["atom_indices"],
                "confidence": prediction["confidence"],
                "prediction_method": prediction["prediction_method"],
                "notes": prediction["notes"],
            }
        )

    peak_list = pd.DataFrame(
        peaks,
        columns=[
            "group_id",
            "ppm",
            "shift_min_ppm",
            "shift_max_ppm",
            "integration",
            "multiplicity",
            "environment_label",
            "atom_indices",
            "confidence",
            "prediction_method",
            "notes",
        ],
    )

    if peak_list.empty:
        peak_list.insert(0, "peak_id", [])
        return peak_list

    peak_list = peak_list.sort_values("ppm", ascending=False).reset_index(drop=True)
    peak_list.insert(0, "peak_id", range(1, len(peak_list) + 1))
    return peak_list


def generate_1h_peak_list(smiles: str, method: str = "rules") -> pd.DataFrame:
    """Generate a predicted 1H NMR peak list from a SMILES string."""
    predictions = get_predictor(method).predict(smiles)
    return generate_peak_list(predictions)


def infer_basic_multiplicity(row: dict) -> str:
    """Infer placeholder multiplicity for an early peak-list scaffold."""
    environment_label = str(row["environment_label"])
    if "exchangeable proton" in environment_label:
        return "br s"
    if "aromatic proton" in environment_label:
        return "m"
    return "unknown"


def format_peak_label(row: dict) -> str:
    """Format a compact human-readable peak label."""
    return (
        f"{row['ppm']:.1f} ppm, "
        f"{row['integration']}H, "
        f"{row['multiplicity']}, "
        f"{row['environment_label']}"
    )
