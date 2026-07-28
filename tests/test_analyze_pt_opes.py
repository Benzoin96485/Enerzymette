"""Tests for validation/analyze_pt_opes.py.

Optional integration coverage uses local data under the gitignored
``validation/`` directory when present; otherwise those cases are skipped.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ase import Atoms

from validation.analyze_pt_opes import assess_run, compute_generic_metrics, write_csv_summary

_REPO_ROOT = Path(__file__).resolve().parents[1]
_OPTIONAL_RUN_DIR = _REPO_ROOT / "validation" / "cmtr1_pt_opes_scan1a_test"


class AnalyzePtOpesTests(unittest.TestCase):
    def test_generic_metrics_detects_transfer(self) -> None:
        donor = 0
        proton = 1
        acceptor = 2
        scope = {
            "index_donor": donor,
            "transfer_protons": [proton],
            "acceptors": [acceptor],
        }
        bound = Atoms(
            symbols=["O", "H", "O"],
            positions=[[0, 0, 0], [0.96, 0, 0], [3, 0, 0]],
        )
        transferred = Atoms(
            symbols=["O", "H", "O"],
            positions=[[0, 0, 0], [2.9, 0.1, 0], [2.95, 0, 0]],
        )
        metrics = compute_generic_metrics([bound, transferred], scope)
        self.assertTrue(metrics["generic_transfer_observed"])
        self.assertEqual(metrics["n_transfer_frames"], 1)

    @unittest.skipUnless(
        (_OPTIONAL_RUN_DIR / "plumed.traj.xyz").exists(),
        "validation trajectory missing",
    )
    def test_assess_existing_validation_run(self) -> None:
        result = assess_run(_OPTIONAL_RUN_DIR, system_profile="cmtr1")
        self.assertTrue(result["traj_exists"])
        self.assertTrue(result["scope_exists"])
        self.assertIn("generic", result)
        self.assertIn("system_specific", result)

    def test_write_csv_summary(self) -> None:
        rows = [
            {
                "run_dir": "run_a",
                "sufficient": False,
                "generic": {
                    "generic_transfer_observed": False,
                    "structure_drift": True,
                    "min_d_donor_h": 0.9,
                    "max_d_donor_h": 4.0,
                    "min_d_acc_h": 2.0,
                    "max_d_acc_h": 3.0,
                    "n_transfer_frames": 0,
                },
                "system_specific": {},
                "failure_modes": ["structure_drift"],
            }
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / "summary.csv"
            write_csv_summary(rows, csv_path)
            content = csv_path.read_text()
            self.assertIn("structure_drift", content)


if __name__ == "__main__":
    unittest.main()
