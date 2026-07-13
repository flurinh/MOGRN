#!/usr/bin/env python3
"""Download MOGRN data from Zenodo.

The project data is hosted as a single archive (mogrn_data.zip) on Zenodo.
This script downloads and extracts the archive to restore the data folders.

Usage:
    python download_data.py <record_id>
    python download_data.py <record_id> --overwrite
    python download_data.py <record_id> --list

Data folders included in the archive:
- opsin_output/     (~1.4G) - Main analysis output
- outputs/          (~380M) - Boltz prediction results
- new_opsins_outputs/ (~140M) - Boltz results for new opsins
- flat_outputs/     (~66M)  - Flattened outputs
- structures/       (~32M)  - Structure files
- flat_new_opsins_outputs/ (~6M) - Flattened new opsin outputs
- property/         (~1.5M) - Property data
- yaml_configs/     (~500K) - YAML configuration files
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import urllib.error
import urllib.request
import zipfile
from typing import Any, Iterable

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None  # type: ignore

# Default archive name on Zenodo
DEFAULT_ARCHIVE = "mogrn_data.zip"

ZENODO_API_URL = "https://zenodo.org/api/records/{record_id}"


def _build_request(url: str, token: str | None = None) -> urllib.request.Request:
    """Build a request with optional authorization header."""
    request = urllib.request.Request(url)
    if token:
        request.add_header("Authorization", f"Bearer {token}")
    return request


def fetch_record_metadata(record_id: str, token: str | None = None) -> dict[str, Any]:
    """Fetch metadata for a Zenodo record."""
    url = ZENODO_API_URL.format(record_id=record_id)
    request = _build_request(url, token)
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def iter_record_files(metadata: dict[str, Any]) -> Iterable[dict[str, Any]]:
    """Iterate over files in a Zenodo record."""
    files = metadata.get("files", [])
    if not files:
        raise ValueError("No files declared in Zenodo record metadata.")
    for entry in files:
        if "links" not in entry:
            raise ValueError(f"Missing download links for file entry: {entry.get('key')!r}")
        yield entry


class _NullProgress:
    """Fallback progress indicator when tqdm is not available."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.total = kwargs.get("total")

    def update(self, _amount: int) -> None:
        return

    def close(self) -> None:
        return


def _open_progress(*, total: int | None, desc: str) -> Any:
    """Open a progress bar if tqdm is available."""
    if tqdm is None:
        return _NullProgress(total=total)
    return tqdm(total=total, desc=desc, unit="iB", unit_scale=True, leave=False)


def download_file(
    url: str,
    destination: pathlib.Path,
    token: str | None = None,
    total_size: int | None = None,
    label: str = "",
) -> None:
    """Download a file from a URL to a destination path."""
    request = _build_request(url, token)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(request) as response:
        if total_size is None:
            try:
                total_size = int(response.headers.get("Content-Length", ""))
            except (TypeError, ValueError):
                total_size = None
        progress = _open_progress(total=total_size, desc=label or destination.name)
        try:
            with destination.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    progress.update(len(chunk))
        finally:
            progress.close()


def extract_archive(archive_path: pathlib.Path, output_dir: pathlib.Path) -> None:
    """Extract a zip archive to the output directory."""
    print(f"Extracting {archive_path.name}...")
    with zipfile.ZipFile(archive_path, "r") as zip_ref:
        members = zip_ref.namelist()
        if tqdm is not None:
            for member in tqdm(members, desc="Extracting", unit="files"):
                zip_ref.extract(member, output_dir)
        else:
            zip_ref.extractall(output_dir)
    print(f"Extracted {len(members)} files to {output_dir}")


def get_file_info(metadata: dict[str, Any], filename: str) -> dict[str, Any] | None:
    """Get info for a specific file in the record."""
    for file_info in iter_record_files(metadata):
        if file_info.get("key") == filename:
            return file_info
    return None


def download_and_extract(
    record_id: str,
    output_dir: pathlib.Path,
    token: str | None = None,
    archive_name: str = DEFAULT_ARCHIVE,
    overwrite: bool = False,
    keep_archive: bool = False,
) -> None:
    """Download and extract the data archive from Zenodo."""
    metadata = fetch_record_metadata(record_id, token)
    file_info = get_file_info(metadata, archive_name)

    if file_info is None:
        available = [f.get("key") for f in iter_record_files(metadata)]
        raise ValueError(
            f"Archive '{archive_name}' not found in record. "
            f"Available files: {', '.join(available)}"
        )

    links = file_info.get("links", {})
    download_url = links.get("download") or links.get("self") or links.get("content")
    if not download_url:
        raise ValueError(f"Missing download link for {archive_name}")

    archive_path = output_dir / archive_name

    # Check if data folders already exist
    data_folders = [
        "opsin_output", "outputs", "new_opsins_outputs", "flat_outputs",
        "structures", "flat_new_opsins_outputs", "property", "yaml_configs"
    ]
    existing = [f for f in data_folders if (output_dir / f).exists()]

    if existing and not overwrite:
        print(f"Data folders already exist: {', '.join(existing)}")
        print("Use --overwrite to replace existing data.")
        return

    # Download
    print(f"Downloading {archive_name} from Zenodo record {record_id}...")
    try:
        download_file(
            download_url,
            archive_path,
            token,
            total_size=int(file_info.get("size") or 0) or None,
            label=archive_name,
        )
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Failed to download {archive_name}: {exc}") from exc

    # Extract
    extract_archive(archive_path, output_dir)

    # Clean up archive unless requested to keep
    if not keep_archive:
        archive_path.unlink()
        print(f"Removed archive {archive_path}")

    print("\nDownload complete! Data folders restored.")


def list_record_files(record_id: str, token: str | None = None) -> None:
    """List all files in a Zenodo record."""
    metadata = fetch_record_metadata(record_id, token)
    print(f"\nFiles in Zenodo record {record_id}:")
    print("-" * 60)
    total_size = 0
    for file_info in iter_record_files(metadata):
        filename = file_info.get("key", "unknown")
        size = file_info.get("size", 0)
        total_size += size
        size_mb = size / (1024 * 1024)
        print(f"  {filename:40} {size_mb:>10.2f} MB")
    print("-" * 60)
    print(f"  {'Total:':40} {total_size / (1024 * 1024):>10.2f} MB")


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "record_id",
        help="Zenodo record identifier (e.g. 1234567)",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        type=pathlib.Path,
        help="Directory to extract data into (default: current directory)",
    )
    parser.add_argument(
        "--token",
        default=None,
        help="Zenodo access token for private or embargoed records.",
    )
    parser.add_argument(
        "--archive",
        default=DEFAULT_ARCHIVE,
        help=f"Name of the archive file on Zenodo (default: {DEFAULT_ARCHIVE})",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing data folders.",
    )
    parser.add_argument(
        "--keep-archive",
        action="store_true",
        help="Keep the downloaded archive after extraction.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_files",
        help="List files in the record without downloading.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    args = parse_args(sys.argv[1:] if argv is None else argv)

    try:
        if args.list_files:
            list_record_files(args.record_id, args.token)
            return 0

        download_and_extract(
            args.record_id,
            args.output_dir,
            args.token,
            archive_name=args.archive,
            overwrite=args.overwrite,
            keep_archive=args.keep_archive,
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
