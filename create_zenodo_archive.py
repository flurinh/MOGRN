#!/usr/bin/env python3
"""Create a deterministic, checksummed MOGRN data archive for Zenodo."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import stat
import zipfile
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_DATA_FOLDERS = (
    "opsin_output",
    "outputs",
    "new_opsins_outputs",
    "flat_outputs",
    "structures",
    "flat_new_opsins_outputs",
    "property",
    "yaml_configs",
)
EXCLUDED_PARTS = {".git", "__pycache__", ".pytest_cache", ".DS_Store"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".tmp"}
ARCHIVE_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
MANIFEST_NAME = "MANIFEST.sha256"


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def should_exclude(relative_path: Path) -> bool:
    """Return whether a relative archive path is generated scratch state."""

    return bool(EXCLUDED_PARTS.intersection(relative_path.parts)) or (
        relative_path.suffix.lower() in EXCLUDED_SUFFIXES
    )


def iter_archive_files(root: Path, folders: Sequence[str]) -> Iterable[Path]:
    """Yield regular files in deterministic archive-path order."""

    candidates: list[Path] = []
    for folder in folders:
        folder_path = root / folder
        if not folder_path.exists():
            continue
        for path in folder_path.rglob("*"):
            relative = path.relative_to(root)
            if path.is_symlink():
                raise ValueError(f"Refusing to archive symbolic link: {relative}")
            if path.is_file() and not should_exclude(relative):
                candidates.append(relative)
    yield from sorted(candidates, key=lambda path: path.as_posix())


def _write_regular_file(archive: zipfile.ZipFile, source: Path, name: str) -> None:
    """Stream one file into the archive with reproducible metadata."""

    info = zipfile.ZipInfo(name, date_time=ARCHIVE_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = (stat.S_IFREG | 0o644) << 16
    with source.open("rb") as input_handle, archive.open(info, "w") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)


def _write_text(archive: zipfile.ZipFile, name: str, content: str) -> None:
    info = zipfile.ZipInfo(name, date_time=ARCHIVE_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = (stat.S_IFREG | 0o644) << 16
    archive.writestr(info, content.encode("utf-8"))


def create_archive(
    root: Path = Path("."),
    output: Path = Path("zenodo_upload/mogrn_data.zip"),
    folders: Sequence[str] = DEFAULT_DATA_FOLDERS,
) -> tuple[Path, str, int]:
    """Create the archive and return its path, SHA-256, and file count."""

    root = root.resolve()
    output = output if output.is_absolute() else root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    files = list(iter_archive_files(root, folders))
    if not files:
        raise ValueError("No release-data files found in the selected folders")

    manifest_lines = [f"{sha256_file(root / path)}  {path.as_posix()}" for path in files]
    manifest = "\n".join(manifest_lines) + "\n"

    try:
        with zipfile.ZipFile(
            temporary, mode="w", compression=zipfile.ZIP_DEFLATED, allowZip64=True
        ) as archive:
            for path in files:
                _write_regular_file(archive, root / path, path.as_posix())
            _write_text(archive, MANIFEST_NAME, manifest)
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()

    archive_digest = sha256_file(output)
    checksum_path = output.with_suffix(output.suffix + ".sha256")
    checksum_path.write_text(
        f"{archive_digest}  {output.name}\n", encoding="utf-8"
    )
    return output, archive_digest, len(files)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path, default=Path("zenodo_upload/mogrn_data.zip")
    )
    parser.add_argument(
        "--folder",
        action="append",
        dest="folders",
        help="Data folder to include; repeat to replace the default folder set",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output, digest, count = create_archive(
        root=args.root,
        output=args.output,
        folders=tuple(args.folders) if args.folders else DEFAULT_DATA_FOLDERS,
    )
    print(f"Archive: {output}")
    print(f"Files: {count}")
    print(f"SHA-256: {digest}")
    print(f"Checksum: {output.with_suffix(output.suffix + '.sha256')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
