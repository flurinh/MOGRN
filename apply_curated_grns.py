#!/usr/bin/env python3
"""Apply the manually curated ProtOS GRN table to prepared MOGRN structures.

The analysis workflow deliberately produces a raw, uncurated baseline. This
separate application step demonstrates the final GRN system by clearing those
provisional labels and persisting the hand-curated ProtOS annotations.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROTOS_SRC = ROOT / "protos" / "src"
if str(PROTOS_SRC) not in sys.path:
    sys.path.insert(0, str(PROTOS_SRC))

import protos  # noqa: E402

from protos.processing.structure import StructureProcessor  # noqa: E402
from src.curated_grn_storage import (  # noqa: E402
    CURATED_REFERENCE,
    aliases_from_tandem_manifest,
    load_curated_structure_aliases,
    overwrite_stored_structures_with_curated_grns,
)


DEFAULT_DATASETS = ("mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / "data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "opsin_output")
    parser.add_argument("--reference", type=Path, default=CURATED_REFERENCE)
    parser.add_argument("--chain-id", default="A")
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    protos.set_data_path(str(args.data_root))
    processor = StructureProcessor("curated_grn_application")

    structure_ids = []
    for dataset in args.datasets:
        structure_ids.extend(processor.get_dataset_entities(dataset))
    structure_ids = list(dict.fromkeys(structure_ids))
    if not structure_ids:
        raise RuntimeError(
            "No prepared structures found; run python prepare_data.py --rebuild first"
        )

    aliases = load_curated_structure_aliases()
    tandem_manifest_path = args.output_dir / "tandem_structure_preprocessing.json"
    if tandem_manifest_path.is_file():
        tandem_manifest = json.loads(tandem_manifest_path.read_text(encoding="utf-8"))
        aliases.update(aliases_from_tandem_manifest(tandem_manifest))

    summary = overwrite_stored_structures_with_curated_grns(
        processor,
        structure_ids,
        reference_path=args.reference,
        aliases=aliases,
        chain_id=args.chain_id,
        output_dir=args.output_dir,
    )
    result = {
        "stage": "curated_protos_application",
        "input_structures": len(structure_ids),
        "annotated_structures": summary["annotated_structure_count"],
        "skipped_structures": summary["skipped_structure_count"],
        "reference": summary["reference"],
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
