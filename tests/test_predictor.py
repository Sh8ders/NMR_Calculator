import pytest

from nmr_calculator.predictor import predict_1h_shifts


def test_ethanol_predictions_include_expected_shift_ranges():
    predictions = predict_1h_shifts("CCO")

    assert len(predictions) == 3
    assert predictions["proton_count"].sum() == 6

    ch3 = _row_for_label(predictions, "simple alkyl CH3")
    assert ch3.shift_min_ppm == 0.8
    assert ch3.shift_max_ppm == 1.8
    assert 0.8 <= ch3.predicted_shift_ppm <= 1.8

    ch2_o = _row_for_label(predictions, "heteroatom-adjacent alkyl proton")
    assert ch2_o.shift_min_ppm == 3.0
    assert ch2_o.shift_max_ppm == 4.5
    assert 3.0 <= ch2_o.predicted_shift_ppm <= 4.5

    exchangeable = _row_for_label(predictions, "alcohol exchangeable proton")
    assert exchangeable.shift_min_ppm == 0.5
    assert exchangeable.shift_max_ppm == 5.5
    assert 0.5 <= exchangeable.predicted_shift_ppm <= 5.5


def test_benzene_prediction_is_aromatic():
    predictions = predict_1h_shifts("c1ccccc1")

    assert len(predictions) == 1
    assert predictions.iloc[0].proton_count == 6
    assert 6.0 <= predictions.iloc[0].predicted_shift_ppm <= 8.5
    assert predictions.iloc[0].environment_label == "aromatic proton"


def test_acetone_prediction_is_carbonyl_adjacent_methyl():
    predictions = predict_1h_shifts("CC(=O)C")

    assert len(predictions) == 1
    assert predictions.iloc[0].proton_count == 6
    assert 2.0 <= predictions.iloc[0].predicted_shift_ppm <= 3.0
    assert "carbonyl" in predictions.iloc[0].notes


def test_predict_1h_shifts_raises_value_error_for_invalid_smiles():
    with pytest.raises(ValueError, match="Invalid SMILES string"):
        predict_1h_shifts("not-a-smiles")


def test_toluene_predictions_include_benzylic_and_aromatic_groups():
    predictions = predict_1h_shifts("Cc1ccccc1")

    assert predictions["proton_count"].sum() == 8
    benzylic = _row_for_label(predictions, "benzylic proton")
    assert benzylic.proton_count == 3
    assert 2.2 <= benzylic.predicted_shift_ppm <= 3.0

    aromatic = predictions[predictions["environment_label"] == "aromatic proton"]
    assert aromatic["proton_count"].sum() == 5
    assert aromatic["predicted_shift_ppm"].between(6.0, 8.5).all()


def test_acetaldehyde_predictions_include_aldehyde_and_alpha_carbonyl():
    predictions = predict_1h_shifts("CC=O")

    assert predictions["proton_count"].sum() == 4
    aldehyde = _row_for_label(predictions, "aldehyde proton")
    assert 9.0 <= aldehyde.predicted_shift_ppm <= 10.5

    methyl = _row_for_label(predictions, "alkyl alpha to carbonyl")
    assert methyl.proton_count == 3
    assert 2.0 <= methyl.predicted_shift_ppm <= 3.0


def test_acetic_acid_predictions_include_carboxylic_acid_proton():
    predictions = predict_1h_shifts("CC(=O)O")

    assert predictions["proton_count"].sum() == 4
    acid = _row_for_label(predictions, "carboxylic acid exchangeable proton")
    assert 10.0 <= acid.predicted_shift_ppm <= 13.0

    methyl = _row_for_label(predictions, "alkyl alpha to carbonyl")
    assert 2.0 <= methyl.predicted_shift_ppm <= 3.0


def test_ethene_predictions_are_vinylic():
    predictions = predict_1h_shifts("C=C")

    assert len(predictions) == 1
    assert predictions.iloc[0].proton_count == 4
    assert predictions.iloc[0].environment_label == "vinylic/alkene proton"
    assert 4.5 <= predictions.iloc[0].predicted_shift_ppm <= 6.5


def test_propyne_predictions_include_terminal_alkyne_and_methyl():
    predictions = predict_1h_shifts("CC#C")

    assert predictions["proton_count"].sum() == 4
    alkyne = _row_for_label(predictions, "terminal alkyne proton")
    assert 2.0 <= alkyne.predicted_shift_ppm <= 3.0

    methyl = _row_for_label(predictions, "simple alkyl CH3")
    assert methyl.proton_count == 3


def _row_for_label(predictions, label):
    matching_rows = predictions[predictions["environment_label"] == label]
    assert len(matching_rows) == 1
    return matching_rows.iloc[0]
