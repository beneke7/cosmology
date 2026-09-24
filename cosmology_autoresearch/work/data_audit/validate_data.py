#!/usr/bin/env python3
"""Reproducible provenance and numerical smoke checks for the data audit."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from background_bao import BAOData  # noqa: E402


EXPECTED_SHA256 = {
    "desi_dr2_mean.txt": "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585",
    "desi_dr2_cov.txt": "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509",
    "des_dovekie_hd.csv": "2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf",
    "des_dovekie_stat_sys.npz": "ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b",
}

EXPECTED_BAO_ROWS = [
    (0.295, "DV_over_rs"),
    (0.510, "DM_over_rs"),
    (0.510, "DH_over_rs"),
    (0.706, "DM_over_rs"),
    (0.706, "DH_over_rs"),
    (0.934, "DM_over_rs"),
    (0.934, "DH_over_rs"),
    (1.321, "DM_over_rs"),
    (1.321, "DH_over_rs"),
    (1.484, "DM_over_rs"),
    (1.484, "DH_over_rs"),
    (2.330, "DH_over_rs"),
    (2.330, "DM_over_rs"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_snana_hubble_diagram(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the release's SNANA-style whitespace table without changing order."""
    columns: list[str] | None = None
    records: list[dict[str, str]] = []
    for line_no, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "VARNAMES:":
            columns = fields[1:]
            continue
        if fields[0] == "SN:":
            if columns is None or len(fields[1:]) != len(columns):
                raise ValueError(f"{path}:{line_no}: malformed SN row")
            records.append(dict(zip(columns, fields[1:], strict=True)))
    if columns is None:
        raise ValueError(f"{path}: no VARNAMES header")
    return columns, records


def main() -> None:
    data_dir = ROOT / "context" / "data"
    hashes = {name: sha256(data_dir / name) for name in EXPECTED_SHA256}
    assert hashes == EXPECTED_SHA256, hashes

    bao = BAOData.from_files(
        data_dir / "desi_dr2_mean.txt", data_dir / "desi_dr2_cov.txt"
    )
    observed_rows = list(zip(bao.z.tolist(), bao.observable.tolist(), strict=True))
    assert observed_rows == EXPECTED_BAO_ROWS, observed_rows
    bao_eigenvalues = np.linalg.eigvalsh(bao.covariance)
    np.linalg.cholesky(bao.covariance)
    bao_condition = float(bao_eigenvalues[-1] / bao_eigenvalues[0])

    mean_text = (data_dir / "desi_dr2_mean.txt").read_text(encoding="utf-8")
    cov = np.loadtxt(data_dir / "desi_dr2_cov.txt", dtype=np.float64)
    correlation = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    off_diagonal = correlation - np.diag(np.diag(correlation))
    bao_blocks = [(0,), (1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12)]
    block_number = np.empty(len(bao.z), dtype=int)
    for block_id, rows in enumerate(bao_blocks):
        block_number[list(rows)] = block_id
    cross_block_mask = block_number[:, None] != block_number[None, :]
    pair_correlations = [
        float(correlation[rows[0], rows[1]]) for rows in bao_blocks if len(rows) == 2
    ]

    hd_path = data_dir / "des_dovekie_hd.csv"
    hd_columns, hd_rows = read_snana_hubble_diagram(hd_path)
    required = {"CID", "IDSURVEY", "zHD", "zHEL", "MU", "MUERR"}
    assert required.issubset(hd_columns), hd_columns
    assert len({row["CID"] for row in hd_rows}) == len(hd_rows)
    for col in ("zHD", "zHEL", "MU", "MUERR"):
        vals = np.asarray([float(row[col]) for row in hd_rows], dtype=np.float64)
        assert np.all(np.isfinite(vals)), col
        if col == "MUERR":
            assert np.all(vals > 0.0)
    zhd = np.asarray([float(row["zHD"]) for row in hd_rows])
    zhel = np.asarray([float(row["zHEL"]) for row in hd_rows])
    mu = np.asarray([float(row["MU"]) for row in hd_rows])
    muerr = np.asarray([float(row["MUERR"]) for row in hd_rows])
    idsurvey: dict[str, int] = {}
    for row in hd_rows:
        sid = row["IDSURVEY"]
        idsurvey[sid] = idsurvey.get(sid, 0) + 1

    npz_path = data_dir / "des_dovekie_stat_sys.npz"
    with np.load(npz_path, allow_pickle=False) as archive:
        n = int(archive[archive.files[0]][0])
        packed_precision = archive[archive.files[1]]
        archive_keys = list(archive.files)
    assert n == len(hd_rows), (n, len(hd_rows))
    assert packed_precision.dtype == np.float32
    assert packed_precision.size == n * (n + 1) // 2
    precision = np.zeros((n, n), dtype=np.float64)
    upper = np.triu_indices(n)
    precision[upper] = packed_precision.astype(np.float64)
    lower = np.tril_indices(n, -1)
    precision[lower] = precision.T[lower]
    assert np.all(np.isfinite(precision))
    assert np.array_equal(precision, precision.T)
    np.linalg.cholesky(precision)
    p_eigenvalues = eigvalsh(
        precision,
        subset_by_index=[0, n - 1],
        check_finite=False,
        driver="evr",
    )
    precision_condition = float(p_eigenvalues[-1] / p_eigenvalues[0])

    code = (ROOT / "work" / "data_audit" / "Dovekie_cosmosis_likelihood.py").read_text(
        encoding="utf-8"
    )
    official_reader = "Table.read(filename, format='ascii.csv')" in code
    has_commas = "," in hd_path.read_text(encoding="utf-8")

    result = {
        "status": "independently_checked",
        "sha256": hashes,
        "desi_dr2": {
            "row_count": len(bao.z),
            "row_order": [
                {"z": float(z), "observable": str(observable), "value": float(value)}
                for z, observable, value in zip(
                    bao.z, bao.observable, bao.value, strict=True
                )
            ],
            "covariance_shape": list(bao.covariance.shape),
            "covariance_exactly_symmetric": bool(
                np.array_equal(bao.covariance, bao.covariance.T)
            ),
            "covariance_cholesky": "pass",
            "covariance_eigenvalue_min": float(bao_eigenvalues[0]),
            "covariance_eigenvalue_max": float(bao_eigenvalues[-1]),
            "covariance_condition_2": bao_condition,
            "max_abs_paired_correlation": float(np.max(np.abs(off_diagonal))),
            "exactly_block_diagonal_by_tracer": bool(
                np.all(cov[cross_block_mask] == 0.0)
            ),
            "max_abs_cross_tracer_covariance": float(
                np.max(np.abs(cov[cross_block_mask]))
            ),
            "within_pair_correlations_in_row_order": pair_correlations,
            "unit_labels_in_mean_file": sorted(set(bao.observable.tolist())),
            "mean_uses_rs_label": "_over_rs" in mean_text,
        },
        "des_dovekie": {
            "hubble_diagram_columns": hd_columns,
            "row_count": len(hd_rows),
            "unique_cid_count": len({row["CID"] for row in hd_rows}),
            "zHD_positive_rows": int(np.count_nonzero(zhd > 0)),
            "zHD_range": [float(zhd.min()), float(zhd.max())],
            "zHEL_range": [float(zhel.min()), float(zhel.max())],
            "MU_range_mag": [float(mu.min()), float(mu.max())],
            "MUERR_range_mag": [float(muerr.min()), float(muerr.max())],
            "IDSURVEY_counts": dict(sorted(idsurvey.items())),
            "archive_keys": archive_keys,
            "packed_matrix_shape": [n, n],
            "packed_entries": int(packed_precision.size),
            "packed_dtype": str(packed_precision.dtype),
            "stored_object": "inverse STAT+SYS covariance",
            "inverse_covariance_cholesky": "pass",
            "inverse_covariance_eigenvalue_min": float(p_eigenvalues[0]),
            "inverse_covariance_eigenvalue_max": float(p_eigenvalues[-1]),
            "covariance_condition_2_from_inverse": precision_condition,
            "likelihood_reader_is_ascii_csv": official_reader,
            "hubble_diagram_contains_csv_commas": has_commas,
            "reader_format_mismatch": bool(official_reader and not has_commas),
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
