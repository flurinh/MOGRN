"""Tests for pre-analysis tandem structure splitting."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from src.tandem_structure_preprocessing import (
    apply_structure_replacements,
    dataset_signature,
    expand_structure_mapping,
    preprocess_registered_tandem_structures,
    update_tandem_helix_definitions,
)


class FakeProcessor:
    def __init__(self, entities):
        self.entities = entities
        self.saved = {}

    def load_entity(self, structure_id):
        return self.entities.get(structure_id)

    def save_entity(self, structure_id, frame, metadata):
        self.saved[structure_id] = {"frame": frame.copy(), "metadata": metadata}


def _parent_frame(length=500):
    rows = []
    for residue in range(1, length + 1):
        rows.append(
            {
                "structure_id": "parent",
                "group": "ATOM",
                "atom_name": "CA",
                "res_atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": residue + 50,
                "label_seq_id": residue,
                "insertion": "",
                "res_name3l": "ALA",
                "x": float(residue),
                "y": 0.0,
                "z": 0.0,
            }
        )
    for atom in range(3):
        rows.append(
            {
                "structure_id": "parent",
                "group": "HETATM",
                "atom_name": f"C{atom}",
                "res_atom_name": f"C{atom}",
                "auth_chain_id": "A",
                "auth_seq_id": 900,
                "label_seq_id": pd.NA,
                "insertion": "",
                "res_name3l": "RET",
                "x": 375.0 + atom * 0.1,
                "y": 1.0,
                "z": 0.0,
            }
        )
    return pd.DataFrame(rows)


def _write_config(path, minimum=500):
    helices = {str(number): [number, number + 1] for number in range(1, 8)}
    config = {
        "schema": "mogrn.tandem-structure-domains/v1",
        "minimum_parent_residues": minimum,
        "parents": {
            "parent": {
                "chain": "A",
                "domains": [
                    {
                        "id": "parent_A",
                        "label": "A",
                        "ordinal_start": 1,
                        "ordinal_end": 250,
                        "helices": helices,
                    },
                    {
                        "id": "parent_B",
                        "label": "B",
                        "ordinal_start": 251,
                        "ordinal_end": 500,
                        "helices": helices,
                    },
                ],
            }
        },
    }
    path.write_text(json.dumps(config))


def test_preprocessor_registers_two_ordinary_locally_numbered_entities(tmp_path):
    config = tmp_path / "domains.json"
    _write_config(config)
    processor = FakeProcessor({"parent": _parent_frame()})

    result = preprocess_registered_tandem_structures(
        processor, ["parent"], config
    )

    assert result.replacements == {"parent": ("parent_A", "parent_B")}
    assert result.entity_parents == {
        "parent_A": "parent",
        "parent_B": "parent",
    }
    domain_a = processor.saved["parent_A"]["frame"]
    domain_b = processor.saved["parent_B"]["frame"]
    assert domain_a[domain_a.group == "ATOM"].auth_seq_id.tolist() == list(range(1, 251))
    assert domain_b[domain_b.group == "ATOM"].auth_seq_id.tolist() == list(range(1, 251))
    assert domain_a[domain_a.group == "ATOM"].parent_auth_seq_id.tolist() == list(range(51, 301))
    assert domain_b[domain_b.group == "ATOM"].parent_auth_seq_id.tolist() == list(range(301, 551))
    assert (domain_a.res_name3l == "RET").sum() == 0
    assert (domain_b.res_name3l == "RET").sum() == 3
    assert result.records[0]["domain_residues"] == 250
    assert result.records[1]["retinal_min_distance"] == pytest.approx(1.0)


def test_length_screen_rejects_configured_parent_below_500(tmp_path):
    config = tmp_path / "domains.json"
    _write_config(config)
    processor = FakeProcessor({"parent": _parent_frame(499)})
    with pytest.raises(ValueError, match="below the 500-residue screen"):
        preprocess_registered_tandem_structures(processor, ["parent"], config)


def test_replacements_expand_datasets_and_pair_mappings():
    replacements = {
        "7pl9": ("TARA_A", "TARA_B"),
        "Tara_RRB_model_0": (
            "Tara_RRB_model_0_A",
            "Tara_RRB_model_0_B",
        ),
    }
    assert apply_structure_replacements(["x", "7PL9", "y"], replacements) == [
        "x",
        "TARA_A",
        "TARA_B",
        "y",
    ]
    assert expand_structure_mapping(
        {"7pl9": "Tara_RRB_model_0"}, replacements
    ) == {
        "TARA_A": "Tara_RRB_model_0_A",
        "TARA_B": "Tara_RRB_model_0_B",
    }


def test_dataset_signature_changes_with_virtual_entities():
    original = dataset_signature(
        {"mo_exp_B": ["7pl9"]}, chain_id="A", retinal_name="RET"
    )
    split = dataset_signature(
        {"mo_exp_B": ["TARA_A", "TARA_B"]},
        chain_id="A",
        retinal_name="RET",
    )
    assert original != split

    changed_boundaries = dataset_signature(
        {"mo_exp_B": ["TARA_A", "TARA_B"]},
        chain_id="A",
        retinal_name="RET",
        preprocessing_signature="updated-config",
    )
    assert split != changed_boundaries


def test_tandem_helix_definitions_are_merged(tmp_path):
    helix_file = tmp_path / "helices.json"
    helix_file.write_text(json.dumps({"existing": {"1": [1, 2]}}))
    records = [
        {
            "structure": "TARA_A",
            "helix_definitions": {"1": [51, 68], "2": [99, 115]},
        }
    ]
    update_tandem_helix_definitions(helix_file, records)
    merged = json.loads(helix_file.read_text())
    assert merged["existing"] == {"1": [1, 2]}
    assert merged["TARA_A"] == {"1": [51, 68], "2": [99, 115]}
