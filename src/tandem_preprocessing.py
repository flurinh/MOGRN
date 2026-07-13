"""Protos-free preprocessing for tandem-rhodopsin protein sequences.

The module screens long sequences for internal repeat evidence and turns a
confident tandem hit into two ordinary sequence records.  It deliberately
stops at sequence preprocessing: downstream annotation and terminal-token
rewriting are separate concerns.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from collections import Counter
from dataclasses import dataclass, replace
from numbers import Real
from pathlib import Path
from statistics import median
from typing import Any, Callable, Iterable, Mapping, Sequence, TextIO
from urllib.parse import quote

from Bio.Align import PairwiseAligner


MIN_TANDEM_LENGTH = 500
MAX_AUTOMATIC_LENGTH = 2500
MANIFEST_SCHEMA = "mogrn.tandem-preprocessing/v1"

# Conservative defaults layered on top of the legacy alignment thresholds.
# They prevent obvious composition-driven repeats and bound automatic work.
MIN_UNIQUE_RESIDUES = 8
MIN_SEQUENCE_ENTROPY = 2.5
MAX_SINGLE_RESIDUE_FRACTION = 0.40
MAX_UNKNOWN_RESIDUE_FRACTION = 0.10
MIN_REPEAT_PERIOD = 180
MAX_REPEAT_PERIOD = 500
NEAR_TIE_SCORE_FRACTION = 0.90
NEAR_TIE_ALIGNED_FRACTION = 0.85
MATERIAL_BOUNDARY_DIFFERENCE = 50
MAX_OPTIMAL_ALIGNMENTS_PER_SPLIT = 32
CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class RepeatThresholds:
    """Minimum evidence required to split a sequence."""

    score: float = 150.0
    aligned_residues: int = 200
    identity: float = 0.35

    def validate(self, *, scan_step: int = 25) -> None:
        """Reject invalid global detector configuration."""

        if isinstance(scan_step, bool) or not isinstance(scan_step, int):
            raise ValueError("scan step must be an integer")
        if scan_step < 1:
            raise ValueError("scan step must be a positive integer")
        if isinstance(self.score, bool) or not isinstance(self.score, Real):
            raise ValueError("minimum repeat score must be numeric")
        if isinstance(self.identity, bool) or not isinstance(self.identity, Real):
            raise ValueError("minimum repeat identity must be numeric")
        score = _finite_number(self.score, "minimum repeat score")
        identity = _finite_number(self.identity, "minimum repeat identity")
        if score <= 0 or score > 2 * MAX_AUTOMATIC_LENGTH:
            raise ValueError(
                "minimum repeat score must be greater than 0 and no greater "
                f"than {2 * MAX_AUTOMATIC_LENGTH}"
            )
        if (
            isinstance(self.aligned_residues, bool)
            or not isinstance(self.aligned_residues, int)
            or not 1 <= self.aligned_residues <= MAX_AUTOMATIC_LENGTH
        ):
            raise ValueError(
                "minimum aligned residues must be an integer between 1 and "
                f"{MAX_AUTOMATIC_LENGTH}"
            )
        if not 0 < identity <= 1:
            raise ValueError(
                "minimum repeat identity must be greater than 0 and at most 1"
            )


@dataclass(frozen=True)
class RepeatEvidence:
    """Sequence-level evidence and proposed repeat-core boundaries."""

    repeat_score: float
    aligned_residues: int
    repeat_identity: float
    core_1_start: int
    core_1_end: int
    core_2_start: int
    core_2_end: int
    test_split: int | None = None
    repeat_period: int | None = None
    ambiguous: bool = False
    ambiguity_reason: str | None = None

    def passes(self, thresholds: RepeatThresholds) -> bool:
        """Return whether all configured evidence thresholds are met."""

        return (
            self.repeat_score >= thresholds.score
            and self.aligned_residues >= thresholds.aligned_residues
            and self.repeat_identity >= thresholds.identity
        )

    def as_manifest(self) -> dict[str, Any]:
        """Return JSON-compatible evidence metadata."""

        return {
            "repeat_score": self.repeat_score,
            "aligned_residues": self.aligned_residues,
            "repeat_identity": self.repeat_identity,
            "test_split": self.test_split,
            "repeat_period": self.repeat_period,
            "proposed_cores": {
                "R1": [self.core_1_start, self.core_1_end],
                "R2": [self.core_2_start, self.core_2_end],
            },
            "ambiguous": self.ambiguous,
            "ambiguity_reason": self.ambiguity_reason,
        }


@dataclass(frozen=True)
class CuratedBoundaries:
    """Explicit, parent-coordinate repeat cores for a known sequence."""

    core_1_start: int
    core_1_end: int
    core_2_start: int
    core_2_end: int
    confidence: str = "curated"

    @classmethod
    def from_intervals(
        cls,
        core_1: Sequence[int],
        core_2: Sequence[int],
        *,
        confidence: str = "curated",
    ) -> "CuratedBoundaries":
        """Build an override from two inclusive 1-based intervals."""

        if (
            isinstance(core_1, (str, bytes))
            or isinstance(core_2, (str, bytes))
            or not isinstance(core_1, Sequence)
            or not isinstance(core_2, Sequence)
            or len(core_1) != 2
            or len(core_2) != 2
        ):
            raise ValueError("each curated core must contain exactly two coordinates")
        if not isinstance(confidence, str) or not confidence.strip():
            raise ValueError("curated confidence must be a nonempty string")
        coordinates = (*core_1, *core_2)
        if any(
            isinstance(coordinate, bool) or not isinstance(coordinate, int)
            for coordinate in coordinates
        ):
            raise ValueError("curated boundary coordinates must be integers")
        return cls(
            core_1_start=core_1[0],
            core_1_end=core_1[1],
            core_2_start=core_2[0],
            core_2_end=core_2[1],
            confidence=confidence,
        )


@dataclass(frozen=True)
class FastaRecord:
    """Minimal FASTA record preserving the original header when possible."""

    record_id: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class PreprocessResult:
    """Output sequences and manifest entry for one parent record."""

    records: tuple[FastaRecord, ...]
    manifest: dict[str, Any]


RepeatDetector = Callable[[str], RepeatEvidence | Mapping[str, Any] | None]


def repeat_aligner() -> PairwiseAligner:
    """Create the local aligner used by the existing tandem-repeat workflow."""

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -0.5
    return aligner


def detect_repeat_evidence(sequence: str, step: int = 25) -> RepeatEvidence | None:
    """Return the strongest internal repeat alignment for an eligible sequence.

    This adapts the evidence calculation from
    ``scripts/detect_annotate_dual_rhodopsins.py`` without any annotation or
    structure-processing dependency.  The hard length screen is repeated here
    so direct detector callers receive the same safe behavior.
    """

    if len(sequence) < MIN_TANDEM_LENGTH or len(sequence) > MAX_AUTOMATIC_LENGTH:
        return None
    RepeatThresholds().validate(scan_step=step)
    sequence = sequence.upper()
    if sequence_complexity_rejection(sequence) is not None:
        return None

    aligner = repeat_aligner()
    candidates: list[RepeatEvidence] = []
    optimal_alignment_overflow = False
    for test_split in range(250, len(sequence) - 249, step):
        alignments = aligner.align(sequence[:test_split], sequence[test_split:])
        for alignment_number, alignment in enumerate(alignments):
            if alignment_number >= MAX_OPTIMAL_ALIGNMENTS_PER_SPLIT:
                optimal_alignment_overflow = True
                break
            evidence = _alignment_evidence(sequence, test_split, alignment)
            if evidence is not None:
                candidates.append(evidence)

    if not candidates:
        return None
    candidates = list({_evidence_placement(item): item for item in candidates}.values())
    best = max(
        candidates,
        key=lambda item: (
            item.repeat_score,
            item.aligned_residues,
            item.repeat_identity,
            -item.core_1_start,
        ),
    )
    if not MIN_REPEAT_PERIOD <= (best.repeat_period or 0) <= MAX_REPEAT_PERIOD:
        return None
    if (
        sequence_complexity_rejection(
            sequence[best.core_1_start - 1 : best.core_1_end]
        )
        is not None
        or sequence_complexity_rejection(
            sequence[best.core_2_start - 1 : best.core_2_end]
        )
        is not None
    ):
        return None

    competing = [
        candidate
        for candidate in candidates
        if candidate is not best
        and _is_near_tied(candidate, best)
        and _materially_different(candidate, best)
    ]
    if optimal_alignment_overflow or competing:
        reasons: list[str] = []
        if competing:
            reasons.append(
                f"{len(competing)} materially different near-tied placement(s)"
            )
        if optimal_alignment_overflow:
            reasons.append(
                "optimal alignment enumeration exceeded the conservative limit"
            )
        return replace(
            best,
            ambiguous=True,
            ambiguity_reason="; ".join(reasons),
        )
    return best


def sequence_complexity_rejection(sequence: str) -> str | None:
    """Return the default low-complexity rejection reason, if any.

    The guard is intentionally conservative and composition-only: at least
    eight canonical amino-acid types, Shannon entropy of 2.5 bits, no residue
    above 40%, and no more than 10% noncanonical residues are required.
    """

    if not sequence:
        return "empty sequence"
    counts = Counter(sequence.upper())
    canonical_counts = {
        residue: count
        for residue, count in counts.items()
        if residue in CANONICAL_AMINO_ACIDS
    }
    canonical_total = sum(canonical_counts.values())
    unknown_fraction = 1 - canonical_total / len(sequence)
    if unknown_fraction > MAX_UNKNOWN_RESIDUE_FRACTION:
        return "more than 10% noncanonical residues"
    if len(canonical_counts) < MIN_UNIQUE_RESIDUES:
        return f"fewer than {MIN_UNIQUE_RESIDUES} canonical residue types"
    maximum_fraction = max(canonical_counts.values()) / canonical_total
    if maximum_fraction > MAX_SINGLE_RESIDUE_FRACTION:
        return "one residue exceeds 40% of the sequence"
    entropy = -sum(
        (count / canonical_total) * math.log2(count / canonical_total)
        for count in canonical_counts.values()
    )
    if entropy < MIN_SEQUENCE_ENTROPY:
        return f"Shannon entropy is below {MIN_SEQUENCE_ENTROPY} bits"
    return None


def _alignment_evidence(
    sequence: str,
    test_split: int,
    alignment: Any,
) -> RepeatEvidence | None:
    blocks_1 = alignment.aligned[0]
    blocks_2 = alignment.aligned[1]
    if not len(blocks_1) or not len(blocks_2):
        return None
    aligned = int(
        sum(
            min(end_1 - start_1, end_2 - start_2)
            for (start_1, end_1), (start_2, end_2) in zip(blocks_1, blocks_2)
        )
    )
    # Preserve the legacy cheap prefilter before residue-level calculations.
    if aligned < 100:
        return None

    matches = 0
    offsets: list[int] = []
    for (start_1, end_1), (start_2, end_2) in zip(blocks_1, blocks_2):
        block_length = min(end_1 - start_1, end_2 - start_2)
        for delta in range(block_length):
            position_1 = int(start_1 + delta)
            position_2 = int(test_split + start_2 + delta)
            matches += sequence[position_1] == sequence[position_2]
            offsets.append(position_2 - position_1)
    return RepeatEvidence(
        repeat_score=float(alignment.score),
        aligned_residues=aligned,
        repeat_identity=matches / aligned,
        test_split=test_split,
        repeat_period=int(round(float(median(offsets)))),
        core_1_start=int(blocks_1[0][0]) + 1,
        core_1_end=int(blocks_1[-1][1]),
        core_2_start=test_split + int(blocks_2[0][0]) + 1,
        core_2_end=test_split + int(blocks_2[-1][1]),
    )


def _evidence_placement(evidence: RepeatEvidence) -> tuple[int, int, int, int]:
    return (
        evidence.core_1_start,
        evidence.core_1_end,
        evidence.core_2_start,
        evidence.core_2_end,
    )


def _is_near_tied(candidate: RepeatEvidence, best: RepeatEvidence) -> bool:
    return (
        candidate.repeat_score >= best.repeat_score * NEAR_TIE_SCORE_FRACTION
        and candidate.aligned_residues
        >= best.aligned_residues * NEAR_TIE_ALIGNED_FRACTION
    )


def _materially_different(candidate: RepeatEvidence, best: RepeatEvidence) -> bool:
    coordinate_difference = max(
        abs(left - right)
        for left, right in zip(
            _evidence_placement(candidate), _evidence_placement(best)
        )
    )
    period_difference = abs(
        (candidate.repeat_period or 0) - (best.repeat_period or 0)
    )
    return (
        coordinate_difference > MATERIAL_BOUNDARY_DIFFERENCE
        or period_difference > MATERIAL_BOUNDARY_DIFFERENCE
    )


def validate_boundaries(
    boundaries: CuratedBoundaries,
    sequence_length: int,
) -> None:
    """Validate inclusive, 1-based, ordered, nonoverlapping core intervals."""

    coordinates = (
        boundaries.core_1_start,
        boundaries.core_1_end,
        boundaries.core_2_start,
        boundaries.core_2_end,
    )
    if any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in coordinates
    ):
        raise ValueError("all boundary coordinates must be integers")
    if not 1 <= boundaries.core_1_start <= boundaries.core_1_end <= sequence_length:
        raise ValueError("R1 core is outside the parent sequence or reversed")
    if not 1 <= boundaries.core_2_start <= boundaries.core_2_end <= sequence_length:
        raise ValueError("R2 core is outside the parent sequence or reversed")
    if boundaries.core_1_end >= boundaries.core_2_start:
        raise ValueError("repeat cores must be ordered and nonoverlapping")
    if (
        not isinstance(boundaries.confidence, str)
        or not boundaries.confidence.strip()
    ):
        raise ValueError("boundary confidence must be a nonempty string")


def ownership_cut(boundaries: CuratedBoundaries) -> int:
    """Return the left residue at the midpoint ownership split (cut|cut+1)."""

    return (boundaries.core_1_end + boundaries.core_2_start) // 2


def preprocess_record(
    parent_id: str,
    sequence: str,
    *,
    description: str = "",
    detector: RepeatDetector | None = None,
    thresholds: RepeatThresholds = RepeatThresholds(),
    curated_boundaries: CuratedBoundaries | None = None,
    scan_step: int = 25,
) -> PreprocessResult:
    """Preprocess one sequence, splitting only a validated confident hit."""

    thresholds.validate(scan_step=scan_step)
    # Whitespace is not a residue, but retain residue letter case so every
    # pass-through sequence is emitted exactly as supplied.
    sequence = "".join(sequence.split())
    base_manifest: dict[str, Any] = {
        "parent_id": parent_id,
        "parent_description": description,
        "parent_length": len(sequence),
        "eligible": len(sequence) >= MIN_TANDEM_LENGTH,
        "minimum_tandem_length": MIN_TANDEM_LENGTH,
        "status": None,
        "boundary_method": None,
        "boundary_confidence": None,
        "evidence": None,
        "segments": [],
        "residual_regions": {},
        "postprocessing": {"terminal_token_rewriting": "not_applied"},
    }

    if len(sequence) < MIN_TANDEM_LENGTH:
        return _pass_through(
            parent_id,
            sequence,
            description,
            base_manifest,
            status="bypassed_too_short",
        )

    if curated_boundaries is None and len(sequence) > MAX_AUTOMATIC_LENGTH:
        return _pass_through(
            parent_id,
            sequence,
            description,
            base_manifest,
            status="pass_through_automatic_length_limit",
        )

    boundary_method: str
    boundary_confidence: str
    if curated_boundaries is not None:
        boundaries = curated_boundaries
        boundary_method = "curated_override"
        boundary_confidence = boundaries.confidence
    else:
        active_detector = detector or (
            lambda candidate: detect_repeat_evidence(candidate, step=scan_step)
        )
        try:
            raw_evidence = active_detector(sequence.upper())
            evidence = _coerce_evidence(raw_evidence, len(sequence))
            if evidence is not None:
                base_manifest["evidence"] = evidence.as_manifest()
                if evidence.ambiguous:
                    return _pass_through(
                        parent_id,
                        sequence,
                        description,
                        base_manifest,
                        status="pass_through_ambiguous_evidence",
                    )
                if not evidence.passes(thresholds):
                    return _pass_through(
                        parent_id,
                        sequence,
                        description,
                        base_manifest,
                        status="pass_through_below_thresholds",
                    )
                boundaries = CuratedBoundaries(
                    evidence.core_1_start,
                    evidence.core_1_end,
                    evidence.core_2_start,
                    evidence.core_2_end,
                    confidence="confident",
                )
        except Exception as error:  # Fail closed at the detector API boundary.
            base_manifest["error"] = f"{type(error).__name__}: {error}"
            return _pass_through(
                parent_id,
                sequence,
                description,
                base_manifest,
                status="pass_through_detector_error",
            )
        if evidence is None:
            return _pass_through(
                parent_id,
                sequence,
                description,
                base_manifest,
                status="pass_through_no_repeat",
            )
        boundary_method = "repeat_alignment"
        boundary_confidence = boundaries.confidence

    try:
        validate_boundaries(boundaries, len(sequence))
    except ValueError as error:
        base_manifest["error"] = str(error)
        base_manifest["boundary_method"] = boundary_method
        base_manifest["boundary_confidence"] = boundary_confidence
        status = (
            "pass_through_invalid_curated_boundaries"
            if curated_boundaries is not None
            else "pass_through_invalid_detected_boundaries"
        )
        return _pass_through(
            parent_id,
            sequence,
            description,
            base_manifest,
            status=status,
        )

    cut = ownership_cut(boundaries)
    intervals = (
        ("R1", boundaries.core_1_start, cut),
        ("R2", cut + 1, boundaries.core_2_end),
    )
    records: list[FastaRecord] = []
    segments: list[dict[str, Any]] = []
    for domain, start, end in intervals:
        virtual_id = stable_virtual_id(parent_id, domain)
        records.append(
            FastaRecord(
                record_id=virtual_id,
                sequence=sequence[start - 1 : end],
                description=(
                    f"parent={parent_id} domain={domain} "
                    f"parent_interval={start}-{end} "
                    f"parent_description={quote(description, safe='')}"
                ),
            )
        )
        segments.append(
            {
                "virtual_id": virtual_id,
                "parent_id": parent_id,
                "parent_description": description,
                "domain": domain,
                "parent_interval": [start, end],
                "length": end - start + 1,
                "boundary_method": boundary_method,
                "boundary_confidence": boundary_confidence,
                "local_to_parent": [
                    {
                        "local_position": local_position,
                        "parent_position": parent_position,
                    }
                    for local_position, parent_position in enumerate(
                        range(start, end + 1), start=1
                    )
                ],
            }
        )

    base_manifest.update(
        {
            "status": (
                "split_curated" if curated_boundaries is not None else "split_detected"
            ),
            "boundary_method": boundary_method,
            "boundary_confidence": boundary_confidence,
            "ownership_cut": {
                "left_parent_position": cut,
                "right_parent_position": cut + 1,
            },
            "core_intervals": {
                "R1": [boundaries.core_1_start, boundaries.core_1_end],
                "R2": [boundaries.core_2_start, boundaries.core_2_end],
            },
            "segments": segments,
            "residual_regions": {
                "n_terminal": _residual_region(
                    sequence, 1, boundaries.core_1_start - 1
                ),
                "c_terminal": _residual_region(
                    sequence, boundaries.core_2_end + 1, len(sequence)
                ),
            },
        }
    )
    return PreprocessResult(tuple(records), base_manifest)


def preprocess_records(
    records: Iterable[FastaRecord],
    *,
    detector: RepeatDetector | None = None,
    thresholds: RepeatThresholds = RepeatThresholds(),
    curated_boundaries: Mapping[str, CuratedBoundaries] | None = None,
    scan_step: int = 25,
) -> tuple[list[FastaRecord], dict[str, Any]]:
    """Preprocess FASTA records and construct the complete run manifest."""

    thresholds.validate(scan_step=scan_step)
    input_records = list(records)
    if curated_boundaries is None:
        curated_boundaries = {}
    elif not isinstance(curated_boundaries, Mapping):
        raise ValueError("curated boundaries must be a mapping keyed by parent ID")
    else:
        curated_boundaries = dict(curated_boundaries)
    invalid_override_ids = [
        parent_id
        for parent_id in curated_boundaries
        if not isinstance(parent_id, str) or not parent_id
    ]
    if invalid_override_ids:
        raise ValueError("curated override IDs must be nonempty strings")
    input_ids = [record.record_id for record in input_records]
    duplicate_inputs = sorted(
        record_id for record_id, count in Counter(input_ids).items() if count > 1
    )
    if duplicate_inputs:
        raise ValueError(
            "duplicate FASTA record ID(s): " + ", ".join(duplicate_inputs)
        )
    unmatched_overrides = sorted(set(curated_boundaries) - set(input_ids))
    if unmatched_overrides:
        raise ValueError(
            "curated override ID(s) not found in FASTA: "
            + ", ".join(unmatched_overrides)
        )
    records_by_id = {record.record_id: record for record in input_records}
    for parent_id, boundaries in curated_boundaries.items():
        if not isinstance(boundaries, CuratedBoundaries):
            raise ValueError(
                f"curated override for {parent_id!r} must be CuratedBoundaries"
            )
        record = records_by_id[parent_id]
        if len(record.sequence) >= MIN_TANDEM_LENGTH:
            try:
                validate_boundaries(boundaries, len(record.sequence))
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"invalid curated override for {parent_id!r}: {error}"
                ) from error

    output_records: list[FastaRecord] = []
    entries: list[dict[str, Any]] = []
    for record in input_records:
        result = preprocess_record(
            record.record_id,
            record.sequence,
            description=record.description,
            detector=detector,
            thresholds=thresholds,
            curated_boundaries=curated_boundaries.get(record.record_id),
            scan_step=scan_step,
        )
        output_records.extend(result.records)
        entries.append(result.manifest)

    output_ids = [record.record_id for record in output_records]
    collisions = sorted(
        record_id for record_id, count in Counter(output_ids).items() if count > 1
    )
    if collisions:
        raise ValueError(
            "output FASTA ID collision(s): " + ", ".join(collisions)
        )

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "coordinate_system": "1-based inclusive parent residue coordinates",
        "detector": {
            "minimum_length": MIN_TANDEM_LENGTH,
            "maximum_automatic_length": MAX_AUTOMATIC_LENGTH,
            "scan_step": scan_step,
            "thresholds": {
                "score": thresholds.score,
                "aligned_residues": thresholds.aligned_residues,
                "identity": thresholds.identity,
            },
            "confidence_guards": {
                "minimum_unique_residues": MIN_UNIQUE_RESIDUES,
                "minimum_shannon_entropy_bits": MIN_SEQUENCE_ENTROPY,
                "maximum_single_residue_fraction": MAX_SINGLE_RESIDUE_FRACTION,
                "maximum_unknown_residue_fraction": MAX_UNKNOWN_RESIDUE_FRACTION,
                "repeat_period_range": [MIN_REPEAT_PERIOD, MAX_REPEAT_PERIOD],
                "near_tie_score_fraction": NEAR_TIE_SCORE_FRACTION,
                "near_tie_aligned_fraction": NEAR_TIE_ALIGNED_FRACTION,
                "material_boundary_difference": MATERIAL_BOUNDARY_DIFFERENCE,
                "maximum_optimal_alignments_per_split": (
                    MAX_OPTIMAL_ALIGNMENTS_PER_SPLIT
                ),
            },
        },
        "records": entries,
    }
    return output_records, manifest


def stable_virtual_id(parent_id: str, domain: str) -> str:
    """Return a deterministic virtual FASTA ID retaining parent and domain."""

    if domain not in {"R1", "R2"}:
        raise ValueError("domain must be R1 or R2")
    return f"{parent_id}__{domain}"


def read_fasta(path: Path) -> list[FastaRecord]:
    """Read ordinary single-line or wrapped FASTA records."""

    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_record_from_fasta(header, chunks))
                header = line[1:].strip()
                chunks = []
                if not header:
                    raise ValueError(f"empty FASTA header on line {line_number}")
            else:
                if header is None:
                    raise ValueError(
                        "FASTA sequence precedes the first header on line "
                        f"{line_number}"
                    )
                chunks.append("".join(line.split()))
    if header is not None:
        records.append(_record_from_fasta(header, chunks))
    if not records:
        raise ValueError("input FASTA contains no records")
    return records


def write_fasta(records: Iterable[FastaRecord], path: Path, width: int = 80) -> None:
    """Write FASTA records, wrapping sequence lines at ``width`` residues."""

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        _write_fasta_stream(records, handle, width=width)


def _write_fasta_stream(
    records: Iterable[FastaRecord],
    handle: TextIO,
    *,
    width: int = 80,
) -> None:
    if width < 1:
        raise ValueError("FASTA line width must be positive")
    for record in records:
        header = record.record_id
        if record.description:
            header = f"{header} {record.description}"
        handle.write(f">{header}\n")
        for offset in range(0, len(record.sequence), width):
            handle.write(f"{record.sequence[offset : offset + width]}\n")


def load_curated_boundaries(path: Path) -> dict[str, CuratedBoundaries]:
    """Load curated overrides from a JSON object keyed by parent FASTA ID."""

    try:
        with path.open(encoding="utf-8") as handle:
            raw = json.load(handle)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"invalid curated-boundary JSON at line {error.lineno}, "
            f"column {error.colno}: {error.msg}"
        ) from error
    if not isinstance(raw, dict):
        raise ValueError("curated-boundary JSON must be an object keyed by parent ID")

    boundaries: dict[str, CuratedBoundaries] = {}
    for parent_id, value in raw.items():
        try:
            if not isinstance(parent_id, str) or not parent_id:
                raise ValueError("parent ID must be a nonempty string")
            if isinstance(value, list):
                if len(value) != 2:
                    raise ValueError("override must contain exactly R1 and R2")
                core_1, core_2 = value
                confidence = "curated"
            elif isinstance(value, dict):
                if "core_1" not in value or "core_2" not in value:
                    raise ValueError("override requires core_1 and core_2")
                core_1 = value["core_1"]
                core_2 = value["core_2"]
                confidence = value.get("confidence", "curated")
            else:
                raise ValueError("override has an unsupported shape")
            boundaries[parent_id] = CuratedBoundaries.from_intervals(
                core_1, core_2, confidence=confidence
            )
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"invalid curated override for {parent_id!r}: {error}"
            ) from error
    return boundaries


def cli_main(argv: Sequence[str] | None = None) -> int:
    """Run the standalone FASTA-to-FASTA preprocessing command."""

    parser = argparse.ArgumentParser(
        description="Split confident tandem-rhodopsin sequences before annotation."
    )
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument(
        "output_fasta",
        type=Path,
        help="expanded FASTA (must share a directory with manifest)",
    )
    parser.add_argument(
        "manifest",
        type=Path,
        help="JSON manifest published transactionally with output FASTA",
    )
    parser.add_argument("--curated-boundaries", type=Path)
    parser.add_argument("--scan-step", type=int, default=25)
    parser.add_argument("--min-repeat-score", type=float, default=150.0)
    parser.add_argument("--min-aligned-residues", type=int, default=200)
    parser.add_argument("--min-repeat-identity", type=float, default=0.35)
    args = parser.parse_args(argv)

    thresholds = RepeatThresholds(
        score=args.min_repeat_score,
        aligned_residues=args.min_aligned_residues,
        identity=args.min_repeat_identity,
    )
    try:
        thresholds.validate(scan_step=args.scan_step)
        _validate_distinct_paths(
            args.input_fasta, args.output_fasta, args.manifest
        )
        overrides = (
            load_curated_boundaries(args.curated_boundaries)
            if args.curated_boundaries
            else {}
        )
        output_records, manifest = preprocess_records(
            read_fasta(args.input_fasta),
            thresholds=thresholds,
            curated_boundaries=overrides,
            scan_step=args.scan_step,
        )
        publish_output_pair(
            output_records,
            manifest,
            args.output_fasta,
            args.manifest,
        )
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        parser.error(str(error))
    return 0


def publish_output_pair(
    records: Iterable[FastaRecord],
    manifest: Mapping[str, Any],
    output_fasta: Path,
    manifest_path: Path,
) -> None:
    """Failure-safely publish a FASTA and its required manifest as one pair.

    Both files are fully serialized and flushed to hidden temporary files in
    their shared destination directory before either final path is touched.
    Existing final files are moved to backups. Synchronous rename or directory
    flush failures remove newly published files and restore both backups.
    """

    destination = output_fasta.expanduser().resolve().parent
    if destination != manifest_path.expanduser().resolve().parent:
        raise ValueError(
            "output FASTA and manifest must share one destination directory "
            "for failure-safe pair publication"
        )
    destination.mkdir(parents=True, exist_ok=True)
    fasta_temp: Path | None = None
    manifest_temp: Path | None = None
    try:
        fasta_temp = _serialize_temp(
            output_fasta,
            lambda handle: _write_fasta_stream(records, handle),
        )
        manifest_temp = _serialize_temp(
            manifest_path,
            lambda handle: _write_manifest_stream(manifest, handle),
        )
        _publish_temp_pair(
            ((fasta_temp, output_fasta), (manifest_temp, manifest_path)),
            destination,
        )
    finally:
        for temporary in (fasta_temp, manifest_temp):
            if temporary is not None:
                _safe_unlink(temporary)


def _serialize_temp(
    destination: Path,
    writer: Callable[[TextIO], None],
) -> Path:
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(
            file_descriptor, "w", encoding="utf-8", newline="\n"
        ) as handle:
            writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        _safe_unlink(temporary)
        raise
    return temporary


def _write_manifest_stream(manifest: Mapping[str, Any], handle: TextIO) -> None:
    json.dump(manifest, handle, indent=2, allow_nan=False)
    handle.write("\n")


def _publish_temp_pair(
    publications: Sequence[tuple[Path, Path]],
    destination_directory: Path,
) -> None:
    backups: dict[Path, Path] = {}
    published: list[Path] = []
    try:
        for _temporary, final_path in publications:
            if os.path.lexists(final_path) and (
                final_path.is_symlink() or not final_path.is_dir()
            ):
                backup = _reserve_backup_path(final_path)
                try:
                    os.replace(final_path, backup)
                except BaseException:
                    _safe_unlink(backup)
                    raise
                backups[final_path] = backup
        _fsync_directory(destination_directory)

        for temporary, final_path in publications:
            os.replace(temporary, final_path)
            published.append(final_path)
        _fsync_directory(destination_directory)
    except BaseException as publication_error:
        rollback_errors: list[str] = []
        for final_path in reversed(published):
            try:
                _safe_unlink(final_path)
            except OSError as error:
                rollback_errors.append(f"remove {final_path}: {error}")
        for final_path, backup in reversed(tuple(backups.items())):
            try:
                os.replace(backup, final_path)
            except OSError as error:
                rollback_errors.append(f"restore {final_path}: {error}")
        try:
            _fsync_directory(destination_directory)
        except OSError as error:
            rollback_errors.append(f"flush rollback directory: {error}")
        if rollback_errors:
            raise RuntimeError(
                "pair publication failed and rollback was incomplete: "
                + "; ".join(rollback_errors)
            ) from publication_error
        raise
    else:
        for backup in backups.values():
            _safe_unlink(backup)
        _fsync_directory(destination_directory)


def _reserve_backup_path(final_path: Path) -> Path:
    file_descriptor, backup_name = tempfile.mkstemp(
        prefix=f".{final_path.name}.",
        suffix=".backup",
        dir=final_path.parent,
    )
    os.close(file_descriptor)
    return Path(backup_name)


def _safe_unlink(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    file_descriptor = os.open(path, flags)
    try:
        os.fsync(file_descriptor)
    finally:
        os.close(file_descriptor)


def _coerce_evidence(
    evidence: RepeatEvidence | Mapping[str, Any] | None,
    sequence_length: int,
) -> RepeatEvidence | None:
    if evidence is None:
        return evidence
    if isinstance(evidence, RepeatEvidence):
        source: Mapping[str, Any] = {
            "repeat_score": evidence.repeat_score,
            "aligned_residues": evidence.aligned_residues,
            "repeat_identity": evidence.repeat_identity,
            "core_1_start": evidence.core_1_start,
            "core_1_end": evidence.core_1_end,
            "core_2_start": evidence.core_2_start,
            "core_2_end": evidence.core_2_end,
            "test_split": evidence.test_split,
            "repeat_period": evidence.repeat_period,
            "ambiguous": evidence.ambiguous,
            "ambiguity_reason": evidence.ambiguity_reason,
        }
    elif isinstance(evidence, Mapping):
        source = evidence
    else:
        raise ValueError("repeat evidence must be a mapping, RepeatEvidence, or None")

    aliases = {
        "core_1_start": "domain_a_aligned_start",
        "core_1_end": "domain_a_aligned_end",
        "core_2_start": "domain_b_aligned_start",
        "core_2_end": "domain_b_aligned_end",
    }
    values: dict[str, Any] = {}
    for field_name in (
        "repeat_score",
        "aligned_residues",
        "repeat_identity",
        "core_1_start",
        "core_1_end",
        "core_2_start",
        "core_2_end",
    ):
        key = field_name if field_name in source else aliases.get(field_name)
        if key is None or key not in source:
            raise ValueError(f"repeat evidence is missing {field_name!r}")
        values[field_name] = source[key]

    ambiguous = source.get("ambiguous", False)
    if not isinstance(ambiguous, bool):
        raise ValueError("repeat evidence ambiguous must be a boolean")
    ambiguity_reason = source.get("ambiguity_reason")
    if ambiguity_reason is not None and not isinstance(ambiguity_reason, str):
        raise ValueError("repeat evidence ambiguity_reason must be a string or null")

    test_split = source.get("test_split")
    if test_split is not None:
        test_split = _bounded_integer(
            test_split,
            "repeat evidence test_split",
            minimum=1,
            maximum=sequence_length - 1,
        )
    repeat_period = source.get("repeat_period")
    if repeat_period is not None:
        repeat_period = _bounded_integer(
            repeat_period,
            "repeat evidence repeat_period",
            minimum=1,
            maximum=sequence_length,
        )

    return RepeatEvidence(
        repeat_score=_bounded_number(
            values["repeat_score"],
            "repeat evidence repeat_score",
            minimum=0,
            maximum=2 * sequence_length,
        ),
        aligned_residues=_bounded_integer(
            values["aligned_residues"],
            "repeat evidence aligned_residues",
            minimum=1,
            maximum=sequence_length,
        ),
        repeat_identity=_bounded_number(
            values["repeat_identity"],
            "repeat evidence repeat_identity",
            minimum=0,
            maximum=1,
        ),
        core_1_start=_bounded_integer(
            values["core_1_start"],
            "repeat evidence core_1_start",
            minimum=1,
            maximum=sequence_length,
        ),
        core_1_end=_bounded_integer(
            values["core_1_end"],
            "repeat evidence core_1_end",
            minimum=1,
            maximum=sequence_length,
        ),
        core_2_start=_bounded_integer(
            values["core_2_start"],
            "repeat evidence core_2_start",
            minimum=1,
            maximum=sequence_length,
        ),
        core_2_end=_bounded_integer(
            values["core_2_end"],
            "repeat evidence core_2_end",
            minimum=1,
            maximum=sequence_length,
        ),
        test_split=test_split,
        repeat_period=repeat_period,
        ambiguous=ambiguous,
        ambiguity_reason=ambiguity_reason,
    )


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (Real, str)):
        raise ValueError(f"{label} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def _bounded_number(
    value: Any,
    label: str,
    *,
    minimum: float,
    maximum: float,
) -> float:
    number = _finite_number(value, label)
    if not minimum <= number <= maximum:
        raise ValueError(f"{label} must be between {minimum} and {maximum}")
    return number


def _bounded_integer(
    value: Any,
    label: str,
    *,
    minimum: int,
    maximum: int,
) -> int:
    number = _finite_number(value, label)
    if not number.is_integer():
        raise ValueError(f"{label} must be an integer")
    integer = int(number)
    if not minimum <= integer <= maximum:
        raise ValueError(f"{label} must be between {minimum} and {maximum}")
    return integer


def _validate_distinct_paths(
    input_path: Path,
    output_path: Path,
    manifest_path: Path,
) -> None:
    labelled_paths = {
        "input FASTA": input_path.expanduser().resolve(),
        "output FASTA": output_path.expanduser().resolve(),
        "manifest": manifest_path.expanduser().resolve(),
    }
    by_path: dict[Path, list[str]] = {}
    for label, path in labelled_paths.items():
        by_path.setdefault(path, []).append(label)
    collisions = [labels for labels in by_path.values() if len(labels) > 1]
    if collisions:
        raise ValueError(
            "input FASTA, output FASTA, and manifest paths must be distinct; "
            + "; ".join("/".join(labels) for labels in collisions)
            + " alias"
        )
    if labelled_paths["output FASTA"].parent != labelled_paths["manifest"].parent:
        raise ValueError(
            "output FASTA and manifest must share one destination directory "
            "for failure-safe pair publication"
        )


def _pass_through(
    parent_id: str,
    sequence: str,
    description: str,
    manifest: dict[str, Any],
    *,
    status: str,
) -> PreprocessResult:
    manifest["status"] = status
    manifest["segments"] = [
        {
            "virtual_id": parent_id,
            "parent_id": parent_id,
            "parent_description": description,
            "domain": None,
            "parent_interval": [1, len(sequence)] if sequence else None,
            "length": len(sequence),
            "boundary_method": None,
            "boundary_confidence": None,
            "local_to_parent": [
                {"local_position": position, "parent_position": position}
                for position in range(1, len(sequence) + 1)
            ],
        }
    ]
    return PreprocessResult((FastaRecord(parent_id, sequence, description),), manifest)


def _residual_region(sequence: str, start: int, end: int) -> dict[str, Any] | None:
    if start > end:
        return None
    return {
        "parent_interval": [start, end],
        "length": end - start + 1,
        "sequence": sequence[start - 1 : end],
    }


def _record_from_fasta(header: str, chunks: Sequence[str]) -> FastaRecord:
    header_fields = header.split(maxsplit=1)
    record_id = header_fields[0]
    description = header_fields[1] if len(header_fields) == 2 else ""
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"FASTA record {record_id!r} has an empty sequence")
    return FastaRecord(record_id, sequence, description)


if __name__ == "__main__":
    raise SystemExit(cli_main())
