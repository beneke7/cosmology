"""Offline tests for context bootstrap data validation and provenance."""

from __future__ import annotations

import hashlib
import json
import os
import copy
import sys
import tempfile
import unittest
from email.message import Message
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fetch_context


class FakeResponse:
    def __init__(self, body: bytes, *, content_type: str = "application/octet-stream", url: str = "https://example.test/data", status: int = 200, content_length: int | None = None):
        self.body = body
        self.offset = 0
        self.status = status
        self.url = url
        self.headers = Message()
        self.headers["Content-Type"] = content_type
        self.headers["Content-Length"] = str(len(body) if content_length is None else content_length)

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        return False

    def getcode(self):
        return self.status

    def geturl(self):
        return self.url

    def read(self, amount: int = -1) -> bytes:
        if amount < 0:
            amount = len(self.body) - self.offset
        block = self.body[self.offset : self.offset + amount]
        self.offset += len(block)
        return block


def one_source(**overrides):
    source = {
        "id": "desi-mean",
        "kind": "data",
        "url": "https://example.test/releases/mean.csv",
        "filename": "data/desi_mean.csv",
        "priority": "core",
        "required": True,
        "fetch": True,
    }
    source.update(overrides)
    return {"sources": [source]}


class FetchContextTests(unittest.TestCase):
    def test_download_records_integrity_and_http_provenance(self):
        body = b"z,mean\n0.1,0.3\n"
        manifest = one_source(sha256=hashlib.sha256(body).hexdigest())
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            output = root / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(
                body, content_type="text/csv", url="https://cdn.example.test/releases/mean.csv"
            )) as open_url:
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024)
            open_url.assert_called_once()
            self.assertEqual(failures, 0)
            self.assertEqual((output / "data/desi_mean.csv").read_bytes(), body)
            saved = json.loads((output / "manifest.lock.json").read_text(encoding="utf-8"))
            row = saved["sources"][0]
            self.assertEqual(row["status"], "downloaded")
            self.assertEqual(row["resolved_url"], "https://cdn.example.test/releases/mean.csv")
            self.assertEqual(row["http_status"], 200)
            self.assertEqual(row["content_type"], "text/csv")
            self.assertEqual(row["bytes"], len(body))
            self.assertEqual(row["file_sha256"], hashlib.sha256(body).hexdigest())
            self.assertEqual(lock["sources"][0]["requested_url"], manifest["sources"][0]["url"])

    def test_lock_checkpoints_completed_sources_before_the_next_download(self):
        first = one_source()["sources"][0]
        second = dict(first, id="second", url="https://example.test/second.csv", filename="data/second.csv")
        manifest = {"sources": [first, second]}
        real_write = fetch_context._atomic_write_json
        snapshots = []

        def capture_checkpoint(path, value):
            snapshots.append(copy.deepcopy(value))
            real_write(path, value)

        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            responses = [FakeResponse(b"first"), FakeResponse(b"second")]
            with mock.patch.object(fetch_context.urllib.request, "urlopen", side_effect=responses), mock.patch.object(
                fetch_context, "_atomic_write_json", side_effect=capture_checkpoint
            ):
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024)
            self.assertEqual(failures, 0)
            self.assertGreaterEqual(len(snapshots), 4)  # initial, after each source, final
            after_first = snapshots[1]["sources"]
            self.assertEqual(after_first[0]["status"], "downloaded")
            self.assertEqual(after_first[1]["status"], "pending")
            self.assertEqual(lock["sources"][1]["status"], "downloaded")

    def test_pdf_html_challenge_is_rejected_and_cleaned(self):
        manifest = one_source(
            id="paper", kind="paper", url="https://example.test/paper.pdf",
            filename="papers/paper.pdf", sha256=None,
        )
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(
                b"<html>challenge</html>", content_type="text/html", url="https://cdn.example.test/paper.pdf"
            )):
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024)
            row = lock["sources"][0]
            self.assertEqual(failures, 1)
            self.assertEqual(row["status"], "failed")
            self.assertEqual(row["error_code"], "invalid_pdf")
            self.assertEqual(row["content_type"], "text/html")
            self.assertEqual(row["http_status"], 200)
            self.assertEqual(row["file_sha256"], hashlib.sha256(b"<html>challenge</html>").hexdigest())
            self.assertFalse((output / "papers/paper.pdf").exists())
            self.assertFalse(list(output.rglob("*.part")))

    def test_valid_pdf_is_kept_and_missing_extractor_is_recorded(self):
        body = b"%PDF-1.4\nsmall test fixture\n"
        manifest = one_source(
            id="paper", kind="paper", url="https://example.test/paper.pdf",
            filename="papers/paper.pdf", sha256=hashlib.sha256(body).hexdigest(),
        )
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(
                body, content_type="application/pdf", url="https://cdn.example.test/paper.pdf"
            )), mock.patch.object(fetch_context.shutil, "which", return_value=None):
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024)
            row = lock["sources"][0]
            self.assertEqual(failures, 0)
            self.assertEqual((output / "papers/paper.pdf").read_bytes(), body)
            self.assertEqual(row["extractor_status"], "missing_pdftotext")
            self.assertIsNone(row["text_path"])

    def test_pdf_extractor_rejects_text_companion_symlink(self):
        if not hasattr(os, "symlink"):
            self.skipTest("symbolic links are unavailable")
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            papers = output / "papers"
            papers.mkdir(parents=True)
            pdf_path = papers / "paper.pdf"
            pdf_path.write_bytes(b"%PDF-1.4\nfixture\n")
            external = Path(temp) / "external.txt"
            external.write_text("preserve", encoding="utf-8")
            try:
                (papers / "paper.txt").symlink_to(external)
            except OSError as exc:
                self.skipTest(f"symbolic links are unavailable: {exc}")
            with mock.patch.object(fetch_context.shutil, "which", return_value="/usr/bin/pdftotext"):
                result = fetch_context._extract_pdf_text(pdf_path, output, replace_existing=True)
            self.assertEqual(result["extractor_status"], "failed")
            self.assertIn("symbolic link", result["extractor_error"])
            self.assertEqual(external.read_text(encoding="utf-8"), "preserve")

    def test_dry_run_is_offline_and_does_not_write(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", side_effect=AssertionError("network used")):
                lock, failures = fetch_context.fetch_manifest(one_source(), output, "core", 1000, dry_run=True)
            self.assertEqual(failures, 0)
            self.assertEqual(lock["sources"][0]["status"], "planned")
            self.assertFalse(output.exists())

    def test_unverified_existing_file_is_preserved_and_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            target = output / "data/desi_mean.csv"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"local user content")
            with mock.patch.object(fetch_context.urllib.request, "urlopen", side_effect=AssertionError("must not replace")):
                lock, failures = fetch_context.fetch_manifest(one_source(), output, "core", 1024)
            self.assertEqual(failures, 1)
            self.assertEqual(lock["sources"][0]["error_code"], "existing_unverified_refused")
            self.assertEqual(target.read_bytes(), b"local user content")

    def test_manifest_hash_mismatch_during_refresh_keeps_verified_file(self):
        original = b"z,mean\n0.1,0.3\n"
        wrong = b"z,mean\n0.1,9.9\n"
        manifest = one_source(sha256=hashlib.sha256(original).hexdigest())
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            target = output / "data/desi_mean.csv"
            target.parent.mkdir(parents=True)
            target.write_bytes(original)
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(wrong, content_type="text/csv")):
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024, refresh=True)
            self.assertEqual(failures, 1)
            self.assertEqual(lock["sources"][0]["error_code"], "sha256_mismatch")
            self.assertEqual(target.read_bytes(), original)
            self.assertFalse(list(output.rglob("*.part")))

    def test_relative_subdirectories_are_allowed_but_traversal_is_rejected(self):
        good = one_source(filename="data/releases/mean.csv")
        self.assertEqual(fetch_context._validate_source(good["sources"][0])[3], "data/releases/mean.csv")
        bad = one_source(filename="data/../../outside.csv")
        with tempfile.TemporaryDirectory() as temp:
            lock, failures = fetch_context.fetch_manifest(bad, Path(temp) / "context", "core", 1024)
            self.assertEqual(failures, 1)
            self.assertEqual(lock["sources"][0]["error_code"], "unsafe_filename")

    def test_size_limit_is_a_recorded_failure_with_response_metadata(self):
        manifest = one_source(max_bytes=3)
        body = b"too long"
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(body, content_type="text/csv")):
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024)
            row = lock["sources"][0]
            self.assertEqual(failures, 1)
            self.assertEqual(row["error_code"], "item_size_limit")
            self.assertEqual(row["http_status"], 200)
            self.assertEqual(row["content_type"], "text/csv")
            self.assertEqual(row["resolved_url"], "https://example.test/data")
            self.assertFalse(list(output.rglob("*.part")))

    def test_truncated_response_is_not_accepted_as_a_successful_data_file(self):
        manifest = one_source(sha256=None)
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(
                b"short", content_type="text/csv", content_length=10
            )):
                lock, failures = fetch_context.fetch_manifest(manifest, output, "core", 1024)
            row = lock["sources"][0]
            self.assertEqual(failures, 1)
            self.assertEqual(row["error_code"], "content_length_mismatch")
            self.assertEqual(row["bytes"], 5)
            self.assertEqual(row["content_length"], 10)
            self.assertFalse((output / "data/desi_mean.csv").exists())

    def test_prior_lock_is_not_reused_after_source_url_changes(self):
        body = b"sample"
        original = one_source(sha256=None)
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "context"
            with mock.patch.object(fetch_context.urllib.request, "urlopen", return_value=FakeResponse(body)):
                fetch_context.fetch_manifest(original, output, "core", 1024)
            changed = one_source(url="https://other.example.test/data.csv", sha256=None)
            with mock.patch.object(fetch_context.urllib.request, "urlopen", side_effect=AssertionError("must not overwrite")):
                lock, failures = fetch_context.fetch_manifest(changed, output, "core", 1024)
            self.assertEqual(failures, 1)
            self.assertEqual(lock["sources"][0]["error_code"], "existing_unverified_refused")
            self.assertEqual((output / "data/desi_mean.csv").read_bytes(), body)


if __name__ == "__main__":
    unittest.main()
