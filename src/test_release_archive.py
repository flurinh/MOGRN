"""Tests for deterministic release archives and safe data extraction."""

from __future__ import annotations

import hashlib
import zipfile

import pytest

from create_zenodo_archive import MANIFEST_NAME, create_archive, sha256_file
from download_data import extract_archive, verify_checksum


def test_release_archive_is_deterministic_and_manifested(tmp_path) -> None:
    (tmp_path / "property").mkdir()
    (tmp_path / "structures").mkdir()
    (tmp_path / "property" / "metadata.csv").write_text("name,value\na,1\n")
    (tmp_path / "structures" / "example.cif").write_text("data_example\n")

    first, first_digest, count = create_archive(
        tmp_path, tmp_path / "first.zip", ("structures", "property")
    )
    second, second_digest, _ = create_archive(
        tmp_path, tmp_path / "second.zip", ("structures", "property")
    )

    assert count == 2
    assert first_digest == second_digest
    assert first.read_bytes() == second.read_bytes()
    with zipfile.ZipFile(first) as archive:
        assert archive.namelist() == [
            "property/metadata.csv",
            "structures/example.cif",
            MANIFEST_NAME,
        ]
        manifest = archive.read(MANIFEST_NAME).decode()
    assert hashlib.sha256(b"name,value\na,1\n").hexdigest() in manifest
    assert "property/metadata.csv" in manifest
    assert sha256_file(first) == first_digest


def test_extract_archive_rejects_path_traversal(tmp_path) -> None:
    archive_path = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("../escaped.txt", "bad")

    with pytest.raises(ValueError, match="Unsafe archive member"):
        extract_archive(archive_path, tmp_path / "output")
    assert not (tmp_path / "escaped.txt").exists()


def test_extract_archive_rejects_windows_style_traversal(tmp_path) -> None:
    archive_path = tmp_path / "unsafe-windows.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("..\\escaped.txt", "bad")

    with pytest.raises(ValueError, match="Unsafe archive member"):
        extract_archive(archive_path, tmp_path / "output")


def test_verify_checksum_accepts_match_and_rejects_mismatch(tmp_path) -> None:
    path = tmp_path / "archive.zip"
    path.write_bytes(b"release-data")
    expected = hashlib.md5(path.read_bytes()).hexdigest()  # Zenodo currently reports MD5.

    verify_checksum(path, f"md5:{expected}")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        verify_checksum(path, "md5:00000000000000000000000000000000")
