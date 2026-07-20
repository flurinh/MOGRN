"""Tests for fixed-helix transfer onto replacement experimental structures."""

from pathlib import Path

import pandas as pd

from src import helix_analysis


def _frame(sequence: str, start: int = 1) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "atom_name": ["CA"] * len(sequence),
            "auth_chain_id": ["A"] * len(sequence),
            "auth_seq_id": range(start, start + len(sequence)),
            "res_name1l": list(sequence),
        }
    )


def test_sequence_position_map_translates_numbering_without_changing_extent():
    source = _frame("ACDEFGHIK", start=1)
    target = _frame("ACDEFGHIK", start=41)

    assert helix_analysis._sequence_position_map(source, target) == {
        source_position: source_position + 40
        for source_position in range(1, 10)
    }


def test_replacement_helix_transfer_uses_fixed_source_bounds(
    tmp_path: Path, monkeypatch
):
    source = _frame("ACDEFGHIK", start=1)
    target = _frame("ACDEFGHIK", start=41)
    source.to_pickle(tmp_path / "old.pkl")
    alias_file = tmp_path / "aliases.json"
    alias_file.write_text(
        """
        {
          "aliases": {
            "new": {
              "source_structure_id": "old",
              "position_mapping": "sequence_alignment"
            }
          }
        }
        """,
        encoding="utf-8",
    )
    monkeypatch.setattr(helix_analysis, "STRUCTURE_CACHE", tmp_path)
    monkeypatch.setattr(helix_analysis, "HELIX_STRUCTURE_ALIASES", alias_file)

    definitions = {
        "old": {
            str(helix): [helix, helix + 1]
            for helix in range(1, 8)
        }
    }
    structures = {"new": {"df": target}}

    transferred = helix_analysis._add_replacement_helix_definitions(
        structures, definitions
    )

    assert transferred["old"]["1"] == [1, 2]
    assert transferred["new"] == {
        str(helix): [helix + 40, helix + 41]
        for helix in range(1, 8)
    }
