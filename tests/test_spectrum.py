import pytest

from nmr_calculator.predictor import predict_1h_shifts
from nmr_calculator.spectrum import generate_1h_peak_list, generate_peak_list


def test_generate_peak_list_from_predictions():
    predictions = predict_1h_shifts("CCO")
    peak_list = generate_peak_list(predictions)

    assert list(peak_list.columns) == [
        "peak_id",
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
    ]
    assert len(peak_list) == 3


def test_ethanol_peak_list_has_expected_peaks():
    peak_list = generate_1h_peak_list("CCO")

    assert len(peak_list) == 3
    assert peak_list["integration"].sum() == 6
    assert peak_list["ppm"].tolist() == sorted(peak_list["ppm"], reverse=True)
    assert 3 in set(peak_list["integration"])
    assert 2 in set(peak_list["integration"])

    exchangeable_peak = peak_list[
        peak_list["environment_label"] == "alcohol/amine/thiol exchangeable proton"
    ].iloc[0]
    assert exchangeable_peak.integration == 1
    assert exchangeable_peak.multiplicity == "br s"


def test_benzene_peak_list_has_one_aromatic_multiplet():
    peak_list = generate_1h_peak_list("c1ccccc1")

    assert len(peak_list) == 1
    assert peak_list.iloc[0].integration == 6
    assert peak_list.iloc[0].multiplicity == "m"
    assert 6.0 <= peak_list.iloc[0].ppm <= 8.5


def test_acetone_peak_list_has_one_carbonyl_adjacent_peak():
    peak_list = generate_1h_peak_list("CC(=O)C")

    assert len(peak_list) == 1
    assert peak_list.iloc[0].integration == 6
    assert 2.0 <= peak_list.iloc[0].ppm <= 3.0


def test_generate_1h_peak_list_raises_value_error_for_invalid_smiles():
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        generate_1h_peak_list("not-a-smiles")
