"""Arithmetic and input-contract tests for the background BAO starter."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad

from background_bao import (
    BAOData,
    GaussianBAOLikelihood,
    comoving_integral,
    e2_cpl,
    fit_model,
    predict_bao,
)


class BackgroundNumericsTests(unittest.TestCase):
    def test_cpl_lcdm_limit_and_eds_integral(self):
        z = np.asarray([0.0, 0.05, 0.4, 1.5, 3.0], dtype=np.float64)
        np.testing.assert_allclose(
            e2_cpl(z, 0.3, -1.0, 0.0), 0.3 * (1.0 + z) ** 3 + 0.7, rtol=2e-15
        )
        np.testing.assert_allclose(
            comoving_integral(z[1:], 1.0),
            2.0 * (1.0 - 1.0 / np.sqrt(1.0 + z[1:])),
            rtol=3e-13,
            atol=3e-14,
        )
        scalar = predict_bao(0.5, "DM_over_rs", 30.0, 0.3)
        vector = predict_bao([0.5], ["DM_over_rs"], 30.0, 0.3)[0]
        self.assertAlmostEqual(float(scalar), float(vector), places=14)

    def test_vector_gauss_legendre_matches_independent_quad(self):
        z = np.asarray([0.02, 0.2, 0.7, 1.8, 3.0], dtype=np.float64)
        actual = comoving_integral(z, 0.29, -0.87, 0.31)
        expected = np.asarray([
            quad(lambda x: 1.0 / np.sqrt(float(e2_cpl(x, 0.29, -0.87, 0.31))), 0.0, zi,
                 epsabs=1e-13, epsrel=1e-13)[0]
            for zi in z
        ])
        np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=3e-14)

    def test_parser_preserves_observable_row_order_and_cholesky_likelihood(self):
        with tempfile.TemporaryDirectory() as directory:
            mean = Path(directory) / "mean.txt"
            cov = Path(directory) / "cov.txt"
            mean.write_text(
                "# [z] [value at z] [quantity]\n"
                "0.3 7.1 DV_over_rs\n"
                "0.5 20.2 DH_over_rs\n"
                "0.5 13.0 DM_over_rs\n",
                encoding="utf-8",
            )
            cov.write_text("0.04 0.002 0\n0.002 0.09 -0.01\n0 -0.01 0.16\n", encoding="utf-8")
            data = BAOData.from_files(mean, cov)
        self.assertEqual(data.observable.tolist(), ["DV_over_rs", "DH_over_rs", "DM_over_rs"])
        self.assertEqual(data.z.dtype, np.dtype(np.float64))
        prediction = predict_bao(data.z, data.observable, 30.0, 0.3)
        residual = data.value - prediction
        chol = np.linalg.cholesky(data.covariance)
        whitened = np.linalg.solve(chol, residual)
        expected_chi2 = float(whitened @ whitened)
        self.assertAlmostEqual(GaussianBAOLikelihood(data).chi2(30.0, 0.3), expected_chi2, places=12)

    def test_parser_rejects_undocumented_observable_and_bad_covariance(self):
        with tempfile.TemporaryDirectory() as directory:
            mean = Path(directory) / "mean.txt"
            cov = Path(directory) / "cov.txt"
            mean.write_text("0.5 12.0 DM_over_rd\n", encoding="utf-8")
            cov.write_text("1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsupported observable"):
                BAOData.from_files(mean, cov)
            with self.assertRaisesRegex(ValueError, "symmetric"):
                BAOData([0.5, 0.6], [1.0, 2.0], ["DM_over_rs", "DH_over_rs"], [[1.0, 0.1], [0.0, 1.0]])
            with self.assertRaisesRegex(ValueError, "positive definite"):
                BAOData([0.5, 0.6], [1.0, 2.0], ["DM_over_rs", "DH_over_rs"], [[1.0, 2.0], [2.0, 1.0]])

    def test_nested_synthetic_lambda_is_recovered_by_cpl_fit(self):
        z = np.asarray([0.2, 0.3, 0.5, 0.5, 0.7, 0.7, 0.9, 1.1, 1.3, 1.3,
                        1.5, 1.8, 2.0, 2.33, 2.33, 2.7], dtype=np.float64)
        observable = np.asarray(["DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs",
                                 "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs",
                                 "DM_over_rs", "DH_over_rs", "DV_over_rs", "DM_over_rs",
                                 "DH_over_rs", "DH_over_rs", "DM_over_rs", "DV_over_rs"])
        values = predict_bao(z, observable, 29.5, 0.31, -1.0, 0.0)
        data = BAOData(z, values, observable, np.diag(np.full(len(z), 0.01**2)))
        recovered = fit_model(data, "cpl", starts=8)
        self.assertLess(recovered["chi2"], 1e-6)
        np.testing.assert_allclose(
            [recovered["alpha"], recovered["Omega_m"], recovered["w0"], recovered["wa"]],
            [29.5, 0.31, -1.0, 0.0], rtol=0.0, atol=2e-4,
        )


if __name__ == "__main__":
    unittest.main()
