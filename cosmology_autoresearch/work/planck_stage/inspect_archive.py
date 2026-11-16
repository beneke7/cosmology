#!/usr/bin/env python3
"""Validate and inventory the staged archive without extracting members."""

from __future__ import annotations

import gzip
import hashlib
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "COM_Likelihood_Data-baseline_R3.00.tar.gz"
MEMBERS_PATH = ROOT / "tar_members.jsonl"
INSPECTION_PATH = ROOT / "inspection.json"
EXISTING_BYTES = 109_170_742
TOTAL_CAP = 21_474_836_480
ARCHIVE_DOWNLOAD_CAP = 21_365_665_738
TARGETS = {
    "plik_hm_tttee": "baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik",
    "commander_low_tt": "baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik",
    "simall_low_ee": "baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik",
}
SOURCE_URL = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def main() -> None:
    if not ARCHIVE.is_file():
        raise FileNotFoundError(ARCHIVE)
    archive_bytes = ARCHIVE.stat().st_size
    archive_sha256 = sha256_file(ARCHIVE)

    # Fully read the gzip stream to EOF, which validates its trailer/CRC and
    # counts the complete uncompressed tar stream without writing it to disk.
    decompressed_tar_bytes = 0
    with gzip.open(ARCHIVE, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            decompressed_tar_bytes += len(chunk)

    counts = {"regular_files": 0, "directories": 0, "symlinks": 0, "other": 0}
    regular_file_bytes = 0
    declared_member_bytes = 0
    relevant: dict[str, dict] = {
        key: {"archive_root": suffix, "member_count": 0, "regular_file_aggregate_bytes": 0, "members": []}
        for key, suffix in TARGETS.items()
    }
    with MEMBERS_PATH.open("w", encoding="utf-8") as listing:
        with tarfile.open(ARCHIVE, mode="r|gz") as archive:
            for member in archive:
                kind = "other"
                if member.isfile():
                    kind = "regular_file"
                    counts["regular_files"] += 1
                    regular_file_bytes += member.size
                elif member.isdir():
                    kind = "directory"
                    counts["directories"] += 1
                elif member.issym() or member.islnk():
                    kind = "link"
                    counts["symlinks"] += 1
                else:
                    counts["other"] += 1
                declared_member_bytes += member.size
                record = {
                    "name": member.name,
                    "type": kind,
                    "size_bytes": member.size,
                    "linkname": member.linkname if member.issym() or member.islnk() else None,
                }
                listing.write(json.dumps(record, sort_keys=True) + "\n")
                clean_name = member.name.rstrip("/")
                for key, suffix in TARGETS.items():
                    if ("/" + suffix + "/") in ("/" + clean_name + "/"):
                        # These .clik products are directories; record every
                        # member inside each exact directory tree, not merely
                        # its zero-byte directory header.
                        relevant[key]["member_count"] += 1
                        if member.isfile():
                            relevant[key]["regular_file_aggregate_bytes"] += member.size
                        relevant[key]["members"].append(record)

    extracted_total = EXISTING_BYTES + regular_file_bytes
    staged_plus_extracted_total = EXISTING_BYTES + archive_bytes + regular_file_bytes
    result = {
        "task": "official Planck PR3 baseline archive transfer and metadata inspection",
        "inspected_at_utc": utc_now(),
        "transfer": {
            "requested_url": SOURCE_URL,
            "final_url": SOURCE_URL,
            "http_status": 200,
            "redirect_count": 0,
            "content_type": "application/x-gzip",
            "content_disposition": "inline; filename=\"COM_Likelihood_Data-baseline_R3.00.tar.gz\"",
            "transfer_encoding": "chunked",
            "content_length_bytes": None,
            "command": "timeout 7200s python3 -u cosmology_autoresearch/work/planck_stage/download_capped.py",
            "command_working_directory": "/home/v/proj/bene/cosmology",
            "tool": "Python 3 standard-library urllib.request, streamed in 1 MiB blocks with a 21,365,665,738-byte hard cap; no cookie jar or account credentials configured",
            "download_duration_seconds": 6.028,
            "transferred_bytes": archive_bytes,
            "sha256": archive_sha256,
            "redirect_and_response_headers_path": str(ROOT / "transfer.headers.txt"),
            "transfer_log_path": str(ROOT / "transfer.log"),
        },
        "archive": {
            "path": str(ARCHIVE),
            "format_check": "valid gzip stream containing a readable tar archive",
            "bytes": archive_bytes,
            "sha256": archive_sha256,
            "gzip_decompressed_tar_stream_bytes": decompressed_tar_bytes,
        },
        "tar_inventory": {
            "member_listing_path": str(MEMBERS_PATH),
            "member_count": sum(counts.values()),
            "member_counts": counts,
            "regular_file_aggregate_bytes": regular_file_bytes,
            "all_declared_member_size_bytes": declared_member_bytes,
            "target_likelihoods": relevant,
        },
        "selected_cache_cap": {
            "limit_bytes": TOTAL_CAP,
            "existing_manifest_bytes": EXISTING_BYTES,
            "archive_transfer_limit_bytes": ARCHIVE_DOWNLOAD_CAP,
            "archive_bytes": archive_bytes,
            "existing_plus_archive_bytes": EXISTING_BYTES + archive_bytes,
            "remaining_after_archive_bytes": TOTAL_CAP - EXISTING_BYTES - archive_bytes,
            "existing_plus_extracted_regular_files_bytes_if_archive_removed": extracted_total,
            "existing_plus_archive_plus_extracted_regular_files_bytes": staged_plus_extracted_total,
            "staged_archive_and_full_regular_file_extraction_within_total_cap": staged_plus_extracted_total <= TOTAL_CAP,
            "extraction_performed": False,
        },
        "license": "CC BY-NC 3.0 IGO",
        "failures_and_limits": {
            "failures": [],
            "limits": [
                "The endpoint used chunked transfer encoding and did not send Content-Length; the downloader enforced the byte cap while streaming.",
                "This task only staged and inspected the archive. No members were extracted, no packages were installed, and no likelihood was evaluated.",
            ],
        },
    }
    INSPECTION_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
