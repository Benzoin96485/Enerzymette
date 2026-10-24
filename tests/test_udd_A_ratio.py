"""Tests for the UDD strength calibration by A_ratio in active learning."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from enerzymette.altoolkit.launcher import active_learning_launcher
from enerzymette.cli import get_parser
from tests.test_multi_initial_al import _write_simulation_yaml, _write_xyz


def _make_launcher(tmpdir: str, A_ratio: float = 0.36) -> active_learning_launcher:
    xyz_path = os.path.join(tmpdir, "init.xyz")
    sim_path = os.path.join(tmpdir, "simulation.yaml")
    _write_xyz(xyz_path)
    _write_simulation_yaml(sim_path, xyz_path)
    return active_learning_launcher(
        output_path=tmpdir,
        tmp_path=tmpdir,
        calculator_patch_key="uma_calculator",
        plumed_patch_key="sammt",
        simulation_config_path=sim_path,
        extraction_config_path=os.path.join(tmpdir, "extraction.yaml"),
        annotation_config_path=os.path.join(tmpdir, "annotation.yaml"),
        training_config_path=os.path.join(tmpdir, "training.yaml"),
        n_iterations=1,
        A_ratio=A_ratio,
    )


def _write_prediction(tmpdir: str, with_E_var_grad: bool = True) -> str:
    rng = np.random.default_rng(0)
    n_frames, n_atoms = 3, 4
    data = {
        "Fa": [rng.normal(scale=1e-2, size=(n_atoms, 3)) for _ in range(n_frames)],
        "E_var": list(rng.uniform(1e-4, 2e-4, size=n_frames)),
    }
    if with_E_var_grad:
        data["E_var_grad"] = [rng.normal(scale=1e-5, size=(n_atoms, 3)) for _ in range(n_frames)]
    prediction_dir = os.path.join(tmpdir, "prediction", "processed_dataset_0")
    os.makedirs(prediction_dir)
    pd.DataFrame(data).to_pickle(os.path.join(prediction_dir, "model-prediction.pkl"))
    return os.path.join(tmpdir, "prediction", "processed_dataset_*", "model-prediction.pkl")


class UDDARatioTests(unittest.TestCase):
    def test_get_udd_A_scales_bias_force_to_A_ratio(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            launcher = _make_launcher(tmpdir, A_ratio=0.36)
            glob_pattern = _write_prediction(tmpdir)
            with mock.patch.object(
                launcher, "_get_prediction_path", return_value=(None, None, None, glob_pattern)
            ):
                B, na = 1.5e-5, 4
                A = launcher._get_udd_A(1, B, na)
            prediction = pd.read_pickle(glob_pattern.replace("*", "0"))
            scale = 1 / (2 * na * B)
            F_bias = [A * scale * g for g in prediction["E_var_grad"]]
            ratio = (
                np.mean(np.abs(np.concatenate(F_bias)))
                / np.mean(np.abs(np.concatenate(list(prediction["Fa"]))))
            )
            self.assertAlmostEqual(ratio, 0.36)

    def test_get_udd_A_without_E_var_grad_keeps_config_A(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            launcher = _make_launcher(tmpdir)
            glob_pattern = _write_prediction(tmpdir, with_E_var_grad=False)
            with mock.patch.object(
                launcher, "_get_prediction_path", return_value=(None, None, None, glob_pattern)
            ):
                self.assertIsNone(launcher._get_udd_A(1, 1.5e-5, 4))

    def test_prediction_model_params_updater_appends_E_var_grad_layer(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            launcher = _make_launcher(tmpdir)
            layers = [
                {"name": "ShallowEnsembleReduce", "params": {"train_only": True, "var": ["E"]}},
                {"name": "Force"},
                {"name": "ShallowEnsembleReduce", "params": {"eval_only": True, "var": ["E", "Fa"]}},
            ]
            model_param = {"layers": [dict(layer) for layer in layers]}
            launcher._prediction_model_params_updater(1, model_param)
            self.assertEqual(model_param["layers"][:-1], layers)
            self.assertEqual(model_param["layers"][-1], {"name": "EnergyVarianceGradient"})

    def test_cli_A_ratio_default_and_override(self) -> None:
        with mock.patch.object(sys, "argv", ["enerzymette", "enerzyme_active_learning"]):
            self.assertEqual(get_parser().A_ratio, 0.36)
        with mock.patch.object(sys, "argv", ["enerzymette", "enerzyme_active_learning", "-ar", "0"]):
            self.assertEqual(get_parser().A_ratio, 0.0)


if __name__ == "__main__":
    unittest.main()
