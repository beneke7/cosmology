#!/usr/bin/env python3
"""Fetch the exact official PR3 base-plikHM-TTTEEE+lowE chain archive safely."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time
from urllib.parse import urlparse
import urllib.request
import zipfile


URL = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/"
    "cosmoparams/COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip"
)
NAME = URL.rsplit("/", 1)[-1]
OUT_DIR = Path(__file__).resolve().parents[1] / "planck_stage"
ARCHIVE = OUT_DIR / NAME
PART = OUT_DIR / (NAME + ".part")
RECORD = OUT_DIR / "reference_chain_transfer.json"
MAX_BYTES = 200 * 1024 * 1024


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if ARCHIVE.exists() or PART.exists():
        raise RuntimeError(f"refusing to overwrite existing target: {ARCHIVE}")

    request = urllib.request.Request(
        URL,
        headers={
            "User-Agent": "cosmology-autoresearch/Planck-chain-reference",
            "Accept-Encoding": "identity",
        },
    )
    t0 = time.monotonic()
    digest = hashlib.sha256()
    count = 0
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.status != 200:
            raise RuntimeError(f"expected HTTP 200, got {response.status}")
        final_url = response.geturl()
        final = urlparse(final_url)
        if final.scheme != "https" or final.hostname != "irsa.ipac.caltech.edu":
            raise RuntimeError(f"unexpected final URL: {final_url}")
        declared = response.headers.get("Content-Length")
        if declared is not None and int(declared) > MAX_BYTES:
            raise RuntimeError(f"declared size exceeds {MAX_BYTES} byte safety cap")
        with PART.open("xb") as target:
            while True:
                block = response.read(min(1024 * 1024, MAX_BYTES - count + 1))
                if not block:
                    break
                count += len(block)
                if count > MAX_BYTES:
                    raise RuntimeError(f"stream exceeded {MAX_BYTES} byte safety cap")
                digest.update(block)
                target.write(block)
        content_type = response.headers.get_content_type()
        headers = dict(response.headers.items())

    if count == 0:
        raise RuntimeError("empty response body")
    with zipfile.ZipFile(PART) as archive:
        corrupt_member = archive.testzip()
        members = archive.infolist()
        if corrupt_member is not None:
            raise RuntimeError(f"CRC check failed for ZIP member {corrupt_member}")
    PART.replace(ARCHIVE)
    record = {
        "url": URL,
        "final_url": final_url,
        "http_status": 200,
        "content_type": content_type,
        "headers": headers,
        "path": str(ARCHIVE),
        "bytes": count,
        "sha256": digest.hexdigest(),
        "max_bytes": MAX_BYTES,
        "zip_members": len(members),
        "zip_crc_test": "pass",
        "runtime_seconds": time.monotonic() - t0,
        "role": "official posterior comparison / proposal covariance only; not likelihood evidence",
    }
    RECORD.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
