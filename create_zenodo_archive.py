#!/usr/bin/env python3
"""Create a zip archive of all data folders for Zenodo upload."""
import os
import zipfile
from pathlib import Path

# Folders to include in the archive
DATA_FOLDERS = [
    "opsin_output",
    "outputs",
    "new_opsins_outputs",
    "flat_outputs",
    "structures",
    "flat_new_opsins_outputs",
    "property",
    "yaml_configs",
]

# Files/patterns to exclude
EXCLUDE_PATTERNS = {".pyc", "__pycache__", ".DS_Store", ".git"}

def should_exclude(path: str) -> bool:
    """Check if a path should be excluded."""
    parts = Path(path).parts
    return any(exc in parts or path.endswith(exc) for exc in EXCLUDE_PATTERNS)

def create_archive():
    """Create the zip archive."""
    output_dir = Path("zenodo_upload")
    output_dir.mkdir(exist_ok=True)
    archive_path = output_dir / "mogrn_data.zip"

    print(f"Creating archive: {archive_path}")

    total_files = 0
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for folder in DATA_FOLDERS:
            folder_path = Path(folder)
            if not folder_path.exists():
                print(f"  Skipping {folder} (not found)")
                continue

            print(f"  Adding {folder}...")
            folder_files = 0
            for root, dirs, files in os.walk(folder_path):
                # Filter out excluded directories
                dirs[:] = [d for d in dirs if d not in EXCLUDE_PATTERNS]

                for file in files:
                    file_path = Path(root) / file
                    if should_exclude(str(file_path)):
                        continue
                    arcname = str(file_path)
                    zf.write(file_path, arcname)
                    folder_files += 1

            print(f"    Added {folder_files} files")
            total_files += folder_files

    size_mb = archive_path.stat().st_size / (1024 * 1024)
    print(f"\nArchive created: {archive_path}")
    print(f"Total files: {total_files}")
    print(f"Archive size: {size_mb:.1f} MB")

if __name__ == "__main__":
    create_archive()
