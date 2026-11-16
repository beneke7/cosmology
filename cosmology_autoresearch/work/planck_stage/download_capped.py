#!/usr/bin/env python3
"""Stream the official Planck PR3 likelihood archive under the locked byte cap."""

from __future__ import annotations

import hashlib
import os
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


URL = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz"
MAX_BYTES = 21_365_665_738
EXPECTED_NAME = "COM_Likelihood_Data-baseline_R3.00.tar.gz"
OUT = Path(__file__).resolve().parent
ARCHIVE = OUT / EXPECTED_NAME
PART = OUT / (EXPECTED_NAME + ".part")
HEADERS = OUT / "transfer.headers.txt"
LOG = OUT / "transfer.log"


def log(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    line = f"{stamp} {message}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as stream:
        stream.write(line + "\n")


def fmt_headers(label: str, url: str, status: int, headers) -> str:
    lines = [f"[{label}]", f"URL: {url}", f"HTTP status: {status}"]
    lines.extend(f"{key}: {value}" for key, value in headers.items())
    return "\n".join(lines) + "\n"


class RedirectRecorder(urllib.request.HTTPRedirectHandler):
    def __init__(self):
        super().__init__()
        self.redirects = []

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        self.redirects.append((req.full_url, code, headers, newurl))
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def main() -> int:
    if ARCHIVE.exists() or PART.exists():
        raise RuntimeError("refusing to overwrite an existing archive or partial transfer")
    LOG.write_text("", encoding="utf-8")
    HEADERS.write_text("", encoding="utf-8")
    log(f"start URL={URL}")
    log(f"streaming byte cap={MAX_BYTES}; cookies and account credentials are not configured")

    redirects = RedirectRecorder()
    opener = urllib.request.build_opener(redirects)
    request = urllib.request.Request(
        URL,
        headers={
            "User-Agent": "cosmology-autoresearch/Planck-stage (bounded public archive fetch)",
            "Accept": "application/gzip, application/x-gzip, application/octet-stream;q=0.9, */*;q=0.1",
            "Accept-Encoding": "identity",
        },
        method="GET",
    )
    started = time.monotonic()
    response = opener.open(request, timeout=60)
    with response:
        status = response.status
        final_url = response.geturl()
        headers = response.headers
        header_dump = "\n".join(
            fmt_headers(f"redirect {index + 1}", old_url, code, hdr)
            + f"Redirect target: {new_url}\n"
            for index, (old_url, code, hdr, new_url) in enumerate(redirects.redirects)
        )
        header_dump += fmt_headers("final response", final_url, status, headers)
        HEADERS.write_text(header_dump, encoding="utf-8")
        log(f"HTTP status={status}; final_url={final_url}; redirects={len(redirects.redirects)}")

        if status != 200:
            raise RuntimeError(f"expected HTTP 200, received {status}")
        disposition = headers.get("Content-Disposition", "")
        if EXPECTED_NAME.lower() not in disposition.lower():
            raise RuntimeError(f"unexpected Content-Disposition: {disposition!r}")
        content_type = headers.get_content_type().lower()
        accepted_types = {"application/gzip", "application/x-gzip", "application/octet-stream"}
        if content_type not in accepted_types:
            raise RuntimeError(f"unexpected Content-Type: {content_type!r}")

        declared = headers.get("Content-Length")
        if declared is not None:
            declared_int = int(declared)
            log(f"Content-Length={declared_int}")
            if declared_int > MAX_BYTES:
                raise RuntimeError("declared Content-Length exceeds the streaming byte cap")
        else:
            log(f"Transfer-Encoding={headers.get('Transfer-Encoding', 'not stated')}; Content-Length absent")

        count = 0
        digest = hashlib.sha256()
        with PART.open("xb") as output:
            while True:
                # Read at most one byte past the cap, so an unknown-length stream
                # can be stopped without ever writing beyond the authorized limit.
                chunk = response.read(min(1024 * 1024, MAX_BYTES - count + 1))
                if not chunk:
                    break
                count += len(chunk)
                if count > MAX_BYTES:
                    raise RuntimeError(f"stream exceeded byte cap ({count} > {MAX_BYTES})")
                output.write(chunk)
                digest.update(chunk)
                if count and count % (64 * 1024 * 1024) < len(chunk):
                    log(f"received_bytes={count}")

        with PART.open("rb") as archive_file:
            magic = archive_file.read(3)
        if magic != b"\x1f\x8b\x08":
            raise RuntimeError(f"body is not a gzip stream (first three bytes {magic.hex()})")

        os.replace(PART, ARCHIVE)
        elapsed = time.monotonic() - started
        log(f"complete bytes={count}; sha256={digest.hexdigest()}; elapsed_seconds={elapsed:.3f}")
        log(f"saved archive={ARCHIVE}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        try:
            log(f"failure={type(exc).__name__}: {exc}")
        except Exception:
            print(f"failure={type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        try:
            PART.unlink(missing_ok=True)
        except OSError:
            pass
        raise
