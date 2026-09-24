#!/usr/bin/env python3
"""Fetch a manifest of research context files with bounded, auditable retrieval.

Only ordinary files are downloaded. Archives are never unpacked and downloaded
code is never executed. Successful PDFs may be converted to text by a locally
installed ``pdftotext`` executable.
"""

from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Callable


DEFAULT_ITEM_LIMIT = 25 * 1024 * 1024
DEFAULT_TOTAL_MIB = 200
CHUNK_SIZE = 64 * 1024
REQUEST_TIMEOUT_SECONDS = 30
MAX_ATTEMPTS = 3
ARXIV_INTERVAL_SECONDS = 3.0
RETRY_DELAYS_SECONDS = (0.25, 0.5)
SUCCESS_STATUSES = {"downloaded", "verified_existing"}
KNOWN_SUFFIXES = {
    ".pdf", ".html", ".htm", ".txt", ".csv", ".tsv", ".json", ".xml",
    ".fits", ".fit", ".dat", ".bin", ".zip", ".gz", ".tar", ".npy",
}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class FetchError(Exception):
    """A safe-to-record retrieval or validation failure."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        context: dict[str, Any] | None = None,
        http_status: int | None = None,
        attempts: int | None = None,
    ):
        super().__init__(message)
        self.code = code
        self.context = context or {}
        self.http_status = http_status
        self.attempts = attempts


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _valid_https_url(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise FetchError("invalid_url", f"{label} must be a non-empty HTTPS URL")
    try:
        parsed = urllib.parse.urlsplit(value.strip())
        port = parsed.port
    except ValueError as exc:
        raise FetchError("invalid_url", f"{label} is not a valid HTTPS URL: {exc}") from exc
    if (
        parsed.scheme.lower() != "https" or not parsed.hostname or parsed.username or parsed.password
        or any(char.isspace() for char in parsed.netloc) or (port is not None and not 1 <= port <= 65535)
    ):
        raise FetchError("invalid_url", f"{label} must be an HTTPS URL without embedded credentials")
    return value.strip()


def _source_url(source: dict[str, Any]) -> tuple[str, str]:
    raw_url = source.get("url")
    url = _valid_https_url(raw_url, "url")
    download_url = source.get("download_url")
    if download_url is None:
        return url, url
    return url, _valid_https_url(download_url, "download_url")


def _default_filename(source: dict[str, Any], request_url: str) -> str:
    source_id = source.get("id")
    if not isinstance(source_id, str) or not source_id.strip():
        raise FetchError("invalid_id", "id must be a non-empty string")
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", source_id.strip()).strip("._-")
    if not slug:
        raise FetchError("invalid_id", "id has no safe filename characters")
    suffix = Path(urllib.parse.urlsplit(request_url).path).suffix.lower()
    if suffix not in KNOWN_SUFFIXES:
        kind = source.get("kind")
        suffix = {"paper": ".pdf", "documentation": ".html", "data": ".dat", "code": ".txt"}.get(kind, ".dat")
    return slug + suffix


def safe_relative_filename(source: dict[str, Any], request_url: str) -> str:
    """Return a POSIX relative name, rejecting traversal and platform tricks."""
    raw = source.get("filename")
    if raw is None:
        raw = _default_filename(source, request_url)
    if not isinstance(raw, str) or not raw or "\\" in raw or "\x00" in raw:
        raise FetchError("unsafe_filename", "filename must be a safe relative path")
    # PurePosixPath catches POSIX absolutes and parent components. Reject a
    # Windows drive prefix too, even though it is not special on Unix.
    path = PurePosixPath(raw)
    if path.is_absolute() or re.match(r"^[A-Za-z]:", raw) or any(part in ("", ".", "..") for part in raw.split("/")):
        raise FetchError("unsafe_filename", "filename must stay within the output directory")
    if not path.parts or path.as_posix() != raw:
        raise FetchError("unsafe_filename", "filename must be a normalized relative path")
    return path.as_posix()


def _validate_source(source: Any) -> tuple[dict[str, Any], str, str, str]:
    if not isinstance(source, dict):
        raise FetchError("invalid_source", "each source must be a JSON object")
    source_id = source.get("id")
    if not isinstance(source_id, str) or not source_id.strip():
        raise FetchError("invalid_id", "id must be a non-empty string")
    kind = source.get("kind")
    if kind not in {"paper", "data", "documentation", "code"}:
        raise FetchError("invalid_kind", "kind must be paper, data, documentation, or code")
    priority = source.get("priority")
    if priority not in {"core", "optional"}:
        raise FetchError("invalid_priority", "priority must be core or optional")
    required = source.get("required", False)
    if not isinstance(required, bool):
        raise FetchError("invalid_required", "required must be a boolean")
    fetch = source.get("fetch", True)
    if not isinstance(fetch, bool):
        raise FetchError("invalid_fetch", "fetch must be a boolean")
    origin_url, request_url = _source_url(source)
    filename = safe_relative_filename(source, request_url)
    expected_sha = source.get("sha256")
    if expected_sha is not None and (not isinstance(expected_sha, str) or not SHA256_RE.fullmatch(expected_sha)):
        raise FetchError("invalid_sha256", "sha256 must contain exactly 64 hexadecimal characters")
    item_limit = source.get("max_bytes", DEFAULT_ITEM_LIMIT)
    if isinstance(item_limit, bool) or not isinstance(item_limit, int) or item_limit <= 0:
        raise FetchError("invalid_max_bytes", "max_bytes must be a positive integer")
    return source, origin_url, request_url, filename


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            manifest = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise FetchError("manifest_read_error", f"cannot read manifest: {exc}") from exc
    if not isinstance(manifest, dict) or not isinstance(manifest.get("sources"), list):
        raise FetchError("invalid_manifest", "manifest must be a JSON object with a sources array")
    return manifest


def _safe_target(output_root: Path, relative_name: str) -> Path:
    root = output_root.resolve()
    target = root.joinpath(*PurePosixPath(relative_name).parts)
    # Do not follow a symlink in an existing path component, even if it points
    # back inside the output tree. This keeps writes predictable and local.
    current = root
    for part in PurePosixPath(relative_name).parts[:-1]:
        current = current / part
        if current.is_symlink():
            raise FetchError("unsafe_filename", "filename traverses a symbolic link")
    if target.is_symlink():
        raise FetchError("unsafe_filename", "target is a symbolic link")
    try:
        target.resolve(strict=False).relative_to(root)
    except ValueError as exc:
        raise FetchError("unsafe_filename", "filename resolves outside the output directory") from exc
    return target


def _hash_file(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        while True:
            block = stream.read(CHUNK_SIZE)
            if not block:
                break
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _lock_record_is_verified(previous: Any, filename: str, digest: str, request_url: str) -> bool:
    return (
        isinstance(previous, dict)
        and previous.get("status") in SUCCESS_STATUSES
        and previous.get("path") == filename
        and previous.get("file_sha256") == digest
        and previous.get("requested_url") == request_url
    )


class _ArxivPacer:
    def __init__(self) -> None:
        self.last_request_start: float | None = None

    def wait_if_needed(self, url: str) -> None:
        host = (urllib.parse.urlsplit(url).hostname or "").lower().rstrip(".")
        if host not in {"arxiv.org", "www.arxiv.org", "export.arxiv.org"}:
            return
        now = time.monotonic()
        if self.last_request_start is not None:
            wait_for = ARXIV_INTERVAL_SECONDS - (now - self.last_request_start)
            if wait_for > 0:
                time.sleep(wait_for)
        self.last_request_start = time.monotonic()


def _retryable_http_status(status: int) -> bool:
    return status == 429 or 500 <= status <= 599


def _parse_content_length(headers: Any) -> int | None:
    try:
        value = headers.get("Content-Length")
    except AttributeError:
        return None
    if value is None:
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 0 else None


def _content_type(headers: Any) -> str | None:
    try:
        value = headers.get("Content-Type")
    except AttributeError:
        return None
    if value is None:
        return None
    return str(value).split(";", 1)[0].strip().lower() or None


def _download_once(
    request_url: str,
    temp_path: Path,
    item_limit: int,
    total_remaining: int,
    pacer: _ArxivPacer,
) -> dict[str, Any]:
    req = urllib.request.Request(
        request_url,
        headers={"User-Agent": "cosmology-autoresearch-context/1.0"},
        method="GET",
    )
    pacer.wait_if_needed(request_url)
    with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SECONDS) as response:
        status = response.getcode()
        headers = response.headers
        length = _parse_content_length(headers)
        content_type = _content_type(headers)
        response_url = response.geturl()
        response_context = {
            "http_status": status,
            "content_type": content_type,
            "resolved_url": response_url,
        }
        if urllib.parse.urlsplit(response_url).scheme.lower() != "https":
            raise FetchError("https_downgrade", "server redirected the HTTPS request to a non-HTTPS URL", context=response_context)
        if length is not None and length > item_limit:
            raise FetchError("item_size_limit", f"Content-Length {length} exceeds item limit {item_limit}", context={**response_context, "content_length": length})
        if length is not None and length > total_remaining:
            raise FetchError("total_size_limit", f"Content-Length {length} exceeds remaining total limit {total_remaining}", context={**response_context, "content_length": length})
        digest = hashlib.sha256()
        count = 0
        with temp_path.open("wb") as out:
            while True:
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                count += len(chunk)
                if count > item_limit:
                    raise FetchError("item_size_limit", f"download exceeds item limit {item_limit}", context={**response_context, "bytes_observed": count})
                if count > total_remaining:
                    raise FetchError("total_size_limit", f"download exceeds remaining total limit {total_remaining}", context={**response_context, "bytes_observed": count})
                digest.update(chunk)
                out.write(chunk)
            if length is not None and count != length:
                raise FetchError(
                    "content_length_mismatch",
                    f"received {count} bytes but Content-Length declared {length}",
                    context={
                        **response_context,
                        "bytes": count,
                        "content_length": length,
                        "file_sha256": digest.hexdigest(),
                    },
                )
            out.flush()
            os.fsync(out.fileno())
    return {
        "http_status": status,
        "content_type": content_type,
        "resolved_url": response_url,
        "bytes": count,
        "file_sha256": digest.hexdigest(),
    }


def _download_with_retries(
    request_url: str,
    temp_path: Path,
    item_limit: int,
    total_remaining: int,
    pacer: _ArxivPacer,
) -> tuple[dict[str, Any], int, int | None, str | None]:
    last_error: Exception | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        http_status: int | None = None
        try:
            if temp_path.exists():
                temp_path.unlink()
            result = _download_once(request_url, temp_path, item_limit, total_remaining, pacer)
            return result, attempt, result.get("http_status"), None
        except urllib.error.HTTPError as exc:
            http_status = exc.code
            error_context = {
                "http_status": exc.code,
                "content_type": _content_type(exc.headers),
                "resolved_url": exc.geturl(),
            }
            try:
                exc.close()
            except Exception:
                pass
            last_error = exc
            retryable = _retryable_http_status(exc.code)
            message = f"HTTP {exc.code}: {exc.reason}"
        except FetchError as exc:
            # Validation and configured byte-limit failures are deterministic;
            # preserve their response metadata but do not retry them.
            if exc.attempts is None:
                exc.attempts = attempt
            raise
        except (urllib.error.URLError, TimeoutError, OSError, http.client.HTTPException) as exc:
            last_error = exc
            retryable = True
            message = str(exc)
        if not retryable or attempt >= MAX_ATTEMPTS:
            code = "http_error" if http_status is not None else "network_error"
            raise FetchError(
                code, message, context=error_context if http_status is not None else None,
                http_status=http_status, attempts=attempt,
            ) from last_error
        time.sleep(RETRY_DELAYS_SECONDS[min(attempt - 1, len(RETRY_DELAYS_SECONDS) - 1)])
    raise FetchError("network_error", str(last_error or "retrieval failed"), attempts=MAX_ATTEMPTS)


def _is_expected_pdf(filename: str, request_url: str, content_type: str | None) -> bool:
    return (
        Path(filename).suffix.lower() == ".pdf"
        or Path(urllib.parse.urlsplit(request_url).path).suffix.lower() == ".pdf"
        or content_type == "application/pdf"
    )


def _validate_download(
    source: dict[str, Any], filename: str, request_url: str, result: dict[str, Any]
) -> None:
    content_type = result.get("content_type")
    if _is_expected_pdf(filename, request_url, content_type):
        try:
            with Path(result["temp_path"]).open("rb") as stream:
                signature = stream.read(5)
        except OSError as exc:
            raise FetchError("pdf_validation_error", f"cannot read downloaded PDF: {exc}") from exc
        if signature != b"%PDF-":
            raise FetchError(
                "invalid_pdf",
                "expected a PDF (%PDF- header), received different content",
                context={key: result[key] for key in ("http_status", "content_type", "resolved_url", "bytes", "file_sha256")},
            )
    expected_sha = source.get("sha256")
    if expected_sha is not None and result["file_sha256"].lower() != expected_sha.lower():
        raise FetchError(
            "sha256_mismatch",
            "downloaded SHA-256 does not match the manifest",
            context={key: result[key] for key in ("http_status", "content_type", "resolved_url", "bytes", "file_sha256")},
        )


def _extract_pdf_text(pdf_path: Path, output_root: Path, replace_existing: bool = False) -> dict[str, Any]:
    executable = shutil.which("pdftotext")
    if executable is None:
        return {"extractor_status": "missing_pdftotext", "text_path": None, "extractor_error": None}
    text_path = pdf_path.with_suffix(".txt")
    if text_path.is_symlink():
        return {"extractor_status": "failed", "text_path": None, "extractor_error": "refusing to follow a symbolic link at the text companion path"}
    if text_path.exists() and not text_path.is_file():
        return {"extractor_status": "failed", "text_path": None, "extractor_error": "text companion path exists and is not a regular file"}
    if text_path.exists() and not replace_existing:
        return {"extractor_status": "text_companion_exists", "text_path": text_path.relative_to(output_root).as_posix(), "extractor_error": None}
    fd, temp_name = tempfile.mkstemp(prefix=f".{text_path.name}.", suffix=".part", dir=text_path.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        proc = subprocess.run(
            [executable, "-layout", str(pdf_path), str(temp_path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            timeout=60,
            check=False,
        )
        if proc.returncode != 0:
            error = proc.stderr.decode("utf-8", "replace")[:500]
            return {"extractor_status": "failed", "text_path": None, "extractor_error": error or f"exit {proc.returncode}"}
        os.replace(temp_path, text_path)
        return {"extractor_status": "completed", "text_path": text_path.relative_to(output_root).as_posix(), "extractor_error": None}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"extractor_status": "failed", "text_path": None, "extractor_error": str(exc)[:500]}
    finally:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass


def _atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_path, path)
    finally:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass


def _read_previous_lock(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            lock = json.load(stream)
    except FileNotFoundError:
        return {}
    except (OSError, json.JSONDecodeError):
        return {}
    records = lock.get("sources") if isinstance(lock, dict) else None
    if not isinstance(records, list):
        return {}
    return {str(item.get("id")): item for item in records if isinstance(item, dict) and item.get("id") is not None}


def _record_base(source: Any, index: int, priority: str) -> dict[str, Any]:
    if isinstance(source, dict):
        return {
            "id": source.get("id"),
            "kind": source.get("kind"),
            "priority": source.get("priority"),
            "required": source.get("required") is True,
            "selected": source.get("fetch", True) is True and (
                priority == "all" or source.get("priority") == priority
            ),
            "note": source.get("note"),
            "manifest_index": index,
            "status": "pending",
        }
    return {"id": None, "manifest_index": index, "status": "invalid_source"}


def fetch_manifest(
    manifest: dict[str, Any],
    output_dir: Path,
    priority: str,
    max_total_bytes: int,
    refresh: bool = False,
    dry_run: bool = False,
    progress: Callable[[int, int, dict[str, Any], str], None] | None = None,
) -> tuple[dict[str, Any], int]:
    output_root = output_dir.resolve()
    lock_path = output_root / "manifest.lock.json"
    previous = _read_previous_lock(lock_path) if not dry_run else {}
    source_rows = manifest["sources"]
    records: list[dict[str, Any]] = [_record_base(s, i, priority) for i, s in enumerate(source_rows)]
    prepared: dict[int, tuple[dict[str, Any], str, str, str]] = {}
    targets: dict[str, list[int]] = {}
    ids: dict[str, list[int]] = {}

    for i, raw in enumerate(source_rows):
        record = records[i]
        try:
            source, origin_url, request_url, filename = _validate_source(raw)
            prepared[i] = (source, origin_url, request_url, filename)
            ids.setdefault(source["id"], []).append(i)
            record.update({
                "id": source["id"], "kind": source["kind"], "priority": source["priority"],
                "required": source.get("required", False),
                "selected": source.get("fetch", True) and (priority == "all" or source["priority"] == priority),
                "fetch": source.get("fetch", True), "origin_url": origin_url,
                "requested_url": request_url, "path": filename, "manifest_sha256": source.get("sha256"),
                "max_bytes": source.get("max_bytes", DEFAULT_ITEM_LIMIT),
            })
            if source.get("fetch", True) and (priority == "all" or source["priority"] == priority):
                targets.setdefault(filename, []).append(i)
            else:
                record["status"] = "not_selected"
        except FetchError as exc:
            record.update({"status": "failed", "error_code": exc.code, "error": str(exc)})

    # Refuse ambiguous paths before writing anything for any colliding source.
    collision_indices = {i for indices in targets.values() if len(indices) > 1 for i in indices}
    for i in collision_indices:
        records[i].update({"status": "failed", "error_code": "destination_collision", "error": "multiple selected sources target the same filename"})
    duplicate_id_indices = {i for indices in ids.values() if len(indices) > 1 for i in indices}
    for i in duplicate_id_indices:
        records[i].update({"status": "failed", "error_code": "duplicate_id", "error": "source IDs must be unique"})

    selected = {
        i for i, item in prepared.items()
        if item[0].get("fetch", True) and (priority == "all" or item[0]["priority"] == priority)
    }
    if dry_run:
        for i in selected - collision_indices - duplicate_id_indices:
            records[i]["status"] = "planned"
        lock = {
            "schema_version": 1, "generated_at": utc_now(), "priority": priority,
            "total_limit_bytes": max_total_bytes, "dry_run": True, "sources": records,
        }
        dry_failures = sum(
            1 for row in records
            if row.get("selected") is True and row.get("required") is True and row.get("status") == "failed"
        )
        return lock, dry_failures

    output_root.mkdir(parents=True, exist_ok=True)
    remaining = max_total_bytes
    pacer = _ArxivPacer()

    def current_lock() -> dict[str, Any]:
        return {
            "schema_version": 1,
            "generated_at": utc_now(),
            "priority": priority,
            "total_limit_bytes": max_total_bytes,
            "sources": records,
        }

    # Persist validated selections and all not-selected entries before the
    # first retrieval, then checkpoint each source as soon as it is complete.
    _atomic_write_json(lock_path, current_lock())
    for i in range(len(source_rows)):
        if i not in selected or i in collision_indices or i in duplicate_id_indices:
            if i in selected:
                if progress:
                    progress(i + 1, len(source_rows), records[i], "finished")
                _atomic_write_json(lock_path, current_lock())
            continue
        source, origin_url, request_url, filename = prepared[i]
        record = records[i]
        if progress:
            progress(i + 1, len(source_rows), record, "starting")
        record["started_at"] = utc_now()
        try:
            target = _safe_target(output_root, filename)
            target.parent.mkdir(parents=True, exist_ok=True)
            # Parent directories may have been replaced by symlinks after the
            # first containment check; validate again immediately before use.
            target = _safe_target(output_root, filename)
            expected_sha = source.get("sha256")
            item_limit = source.get("max_bytes", DEFAULT_ITEM_LIMIT)
            if target.exists():
                if not target.is_file():
                    raise FetchError("existing_target_not_file", "existing target is not a regular file")
                local_hash, local_size = _hash_file(target)
                if expected_sha is not None and local_hash.lower() != expected_sha.lower():
                    raise FetchError("existing_sha256_mismatch", "existing file SHA-256 does not match the manifest; refusing to replace it")
                record.update({"existing_bytes": local_size, "existing_sha256": local_hash})
                was_verified = expected_sha is not None or _lock_record_is_verified(previous.get(str(source["id"])), filename, local_hash, request_url)
                if was_verified and not refresh:
                    if local_size > item_limit:
                        raise FetchError("item_size_limit", f"existing file size {local_size} exceeds item limit {item_limit}")
                    if local_size > remaining:
                        raise FetchError("total_size_limit", f"existing file size {local_size} exceeds remaining total limit {remaining}")
                    remaining -= local_size
                    old = previous.get(str(source["id"]), {})
                    record.update({
                        "status": "verified_existing", "resolved_url": old.get("resolved_url"),
                        "http_status": old.get("http_status"), "content_type": old.get("content_type"),
                        "bytes": local_size, "file_sha256": local_hash, "attempts": 0,
                        "retrieved_at": old.get("retrieved_at"), "verified_at": utc_now(),
                    })
                    if target.suffix.lower() == ".pdf":
                        record.update(_extract_pdf_text(target, output_root))
                    continue
                if not refresh and not was_verified:
                    raise FetchError("existing_unverified_refused", "existing file has no matching manifest hash or prior verified lock record")

            fd, temp_name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".part", dir=target.parent)
            os.close(fd)
            temp_path = Path(temp_name)
            try:
                result, attempts, http_status, error = _download_with_retries(
                    request_url, temp_path, item_limit, remaining, pacer
                )
                result["temp_path"] = str(temp_path)
                record.update({
                    "resolved_url": result["resolved_url"], "http_status": http_status,
                    "content_type": result["content_type"], "bytes": result["bytes"],
                    "file_sha256": result["file_sha256"], "attempts": attempts,
                })
                _validate_download(source, filename, request_url, result)
                os.replace(temp_path, target)
            finally:
                try:
                    temp_path.unlink()
                except FileNotFoundError:
                    pass
            remaining -= result["bytes"]
            record.update({
                "status": "downloaded", "resolved_url": result["resolved_url"],
                "http_status": http_status, "content_type": result["content_type"],
                "bytes": result["bytes"], "file_sha256": result["file_sha256"],
                "attempts": attempts, "retrieved_at": utc_now(), "error": error,
            })
            if _is_expected_pdf(filename, request_url, result["content_type"]):
                record.update(_extract_pdf_text(target, output_root, replace_existing=refresh))
        except FetchError as exc:
            record.update(exc.context)
            record.update({"status": "failed", "error_code": exc.code, "error": str(exc)})
            if exc.code in {"http_error", "network_error"}:
                record["attempts"] = exc.attempts
                record["http_status"] = exc.http_status
            elif exc.attempts is not None:
                record["attempts"] = exc.attempts
            record["failed_at"] = utc_now()
        except (OSError, ValueError) as exc:
            record.update({"status": "failed", "error_code": "local_io_error", "error": str(exc), "failed_at": utc_now()})
        finally:
            if progress:
                progress(i + 1, len(source_rows), record, "finished")
            _atomic_write_json(lock_path, current_lock())

    failed_required = sum(
        1 for row in records
        if row.get("selected") is True and row.get("required") is True and row.get("status") == "failed"
    )
    lock = current_lock()
    _atomic_write_json(lock_path, lock)
    return lock, failed_required


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("sources.json"), help="source manifest JSON (default: sources.json)")
    parser.add_argument("--output", type=Path, default=Path("context"), help="output directory (default: context)")
    parser.add_argument("--priority", choices=("core", "all"), default="core", help="select core sources or all fetchable sources")
    parser.add_argument("--max-total-mib", type=float, default=DEFAULT_TOTAL_MIB, help=f"maximum total bytes retrieved (default: {DEFAULT_TOTAL_MIB} MiB)")
    parser.add_argument("--refresh", action="store_true", help="re-fetch and atomically replace existing files after verification")
    parser.add_argument("--dry-run", action="store_true", help="validate and print the selection without network or filesystem writes")
    return parser


def _cli_progress(index: int, total: int, record: dict[str, Any], phase: str) -> None:
    source_id = record.get("id") or f"source-{index}"
    if phase == "starting":
        message = f"[{index}/{total}] fetching {source_id}"
    else:
        details = record.get("error") or f"{record.get('bytes', 0)} bytes"
        message = f"[{index}/{total}] {source_id}: {record.get('status')} ({details})"
    print(message, file=sys.stderr, flush=True)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not math.isfinite(args.max_total_mib) or args.max_total_mib <= 0 or int(args.max_total_mib * 1024 * 1024) <= 0:
        parser.error("--max-total-mib must be positive")
    max_total_bytes = int(args.max_total_mib * 1024 * 1024)
    try:
        manifest = _load_manifest(args.manifest)
        lock, failed_required = fetch_manifest(
            manifest,
            args.output,
            args.priority,
            max_total_bytes,
            refresh=args.refresh,
            dry_run=args.dry_run,
            progress=None if args.dry_run else _cli_progress,
        )
    except FetchError as exc:
        print(f"fetch_context: {exc.code}: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"fetch_context: local_io_error: {exc}", file=sys.stderr)
        return 2
    if args.dry_run:
        print(json.dumps(lock, indent=2, sort_keys=True))
    else:
        counts: dict[str, int] = {}
        for item in lock["sources"]:
            status = str(item.get("status", "unknown"))
            counts[status] = counts.get(status, 0) + 1
        print(json.dumps({"lock": str(args.output / "manifest.lock.json"), "statuses": counts}, sort_keys=True))
    return 1 if failed_required else 0


if __name__ == "__main__":
    raise SystemExit(main())
