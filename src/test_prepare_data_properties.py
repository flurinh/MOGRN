"""Checks that data preparation consumes the current ST5 property workbook."""

from pathlib import Path

from prepare_data import (
    PROPERTY_FILE,
    create_structure_mapping,
    load_property_data,
)


def test_st5_is_the_workflow_property_authority() -> None:
    properties = load_property_data()

    assert PROPERTY_FILE.name == "mo_exp_ST5_HEK1.xlsx"
    assert len(properties) == 128
    assert properties["short_name"].is_unique
    assert "8wew" not in set(properties["exp_structure_id"])

    expected = {
        "HwMR": "9jws",
        "MacR": "8rso",
        "KnChR": "9j7w",
        "SrRhPDE": "7cj3",
    }
    observed = properties.set_index("short_name")["exp_structure_id"].to_dict()
    for short_name, structure_id in expected.items():
        assert observed[short_name] == structure_id


def test_new_st5_experimental_structure_is_in_held_out_set() -> None:
    properties = load_property_data()
    knchr = properties.set_index("short_name").loc["KnChR"]

    assert knchr["experimentally_determined"] == 1
    assert knchr["dataset_split"] == "B"


def test_real_knchr_structure_maps_to_existing_prediction() -> None:
    properties = load_property_data()
    available = {
        "9j7w": Path("9j7w.cif"),
        "KnChR_J444_refine5": Path("KnChR_J444_refine5.cif"),
        "KnChR_J444_refine5_model_0": Path(
            "KnChR_J444_refine5_model_0.cif"
        ),
    }

    mapping = create_structure_mapping(properties, available)

    assert mapping["9j7w"] == "KnChR_J444_refine5_model_0"
    assert "KnChR_J444_refine5" not in mapping
