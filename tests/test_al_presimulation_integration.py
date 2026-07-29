"""AL launcher integration smoke: restrained pre-sim + pool writeback policy.

Exercises the full ``active_learning_launcher.launch()`` loop with a fake
``enerzyme`` CLI (no host paths, Slurm, or GPU models required). Verifies:

1. When the simulation method is ``standard_steered_md`` and
   ``n_presimulation_steps_per_iteration > 0``, the pre-simulation config is
   switched to ``standard_restrained_md``.
2. Structure-pool XYZ files are updated by pre-simulation only — never by the
   steered-MD end frame.
"""

from __future__ import annotations

import json
import os
import pickle
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import List, Optional
from unittest import mock

import ase.io
import numpy as np
import yaml
from ase import Atoms

from enerzymette.altoolkit.launcher import active_learning_launcher

# Distinct displacements so we can tell seed / pre-sim / steered frames apart.
_SEED = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]])
_PRESIM = _SEED + np.array([[0.1, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.0, 0.0]])
_STEERED = _SEED + np.array([[0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])


class _EVar:
    def mean(self):
        return self

    def item(self) -> float:
        return 1.0e-3


def _atoms(positions) -> Atoms:
    return Atoms("CCC", positions=np.asarray(positions, dtype=float))


def _write_xyz(path: str, positions) -> None:
    ase.io.write(path, _atoms(positions), format="extxyz")


def _write_traj(path: str, positions) -> None:
    """Write ≥2 frames; ``collect_trajectory`` skips the first frame."""
    frame = _atoms(positions)
    ase.io.write(path, [frame, frame], format="extxyz")


def _write_minimal_yamls(root: Path, structure_file: str) -> dict:
    sim = {
        "System": {"structure_file": structure_file, "charge": 0, "multiplicity": 1},
        "Simulation": {
            "task": "plumed",
            "idx_start_from": 0,
            "integrate": {
                "integrator": "Langevin",
                "n_step": 4,
                "time_step": 0.5,
                "friction": 0.01,
                "temperature_in_K": 300,
            },
            "plumed_config_generator": {
                "name": "BondReactionConfigGenerator",
                "method": "standard_steered_md",
            },
            "sampling": {
                "params": {
                    "plumed_config": {
                        "dump_interval": 1,
                        "lower_bound": -2.0,
                        "upper_bound": 2.0,
                        "forming_bond": {"atom1": {"index": 1}, "atom2": {"index": 2}},
                        "breaking_bond": {"atom1": {"index": 0}, "atom2": {"index": 1}},
                        "max_bond_length": 3.0,
                    }
                }
            },
            "uncertainty_calculator": {
                "name": "UDD",
                "params": {"A": 0.08, "B": 1.0e-6},
            },
        },
    }
    extract = {
        "Datahub": {
            "data_format": "pickle",
            "data_path": None,
            "features": {"N": None, "Q": "Q", "Ra": "Ra", "Za": "Za"},
            "neighbor_list": "full",
            "preload": True,
            "targets": None,
        },
        "Extractor": {
            "fragment_per_frame": 1,
            "fragment_radius": 3.0,
            "local_uncertainty_radius": 4.0,
            "n_centers": 1,
            "reference_mol_path": None,
        },
        "Metric": None,
        "Trainer": {"inference_batch_size": 1},
    }
    annotate = {"Supplier": {"path": None}}
    train = {"Trainer": {"batch_size": 1, "max_epochs": 1, "num_workers": 0}}
    model = {
        "Modelhub": {
            "internal_FFs": {
                "FF01": {
                    "active": True,
                    "architecture": "FakeNet",
                    "suffix": "",
                    "layers": [
                        {"name": "Embedding", "params": {}},
                        {"name": "ShallowEnsembleReduce", "params": {"train_only": True}},
                        {"name": "Force", "params": {}},
                    ],
                },
            }
        },
        "Datahub": {
            "data_path": None,
            "features": {"N": None, "Za": "Za", "Ra": "Ra", "Q": "Q"},
            "targets": {"E": "E"},
            "transforms": {},
        },
        "Trainer": {},
    }

    paths = {
        "simulation": root / "simulation.yaml",
        "extraction": root / "extraction.yaml",
        "annotation": root / "annotation.yaml",
        "training": root / "training.yaml",
        "model": root / "model.yaml",
    }
    for key, cfg in (
        ("simulation", sim),
        ("extraction", extract),
        ("annotation", annotate),
        ("training", train),
        ("model", model),
    ):
        with open(paths[key], "w") as handle:
            yaml.dump(cfg, handle, default_flow_style=False)
    return {k: str(v) for k, v in paths.items()}


class _FakeEnerzyme:
    """Minimal enerzyme stand-in used through ``subprocess.Popen``."""

    def __init__(self) -> None:
        self.presim_methods: List[Optional[str]] = []
        self.steered_methods: List[Optional[str]] = []

    def __call__(self, args, *a, stdout=None, stderr=None, **kw):
        argv = list(args)
        if not argv or Path(argv[0]).name != "enerzyme":
            raise AssertionError(f"unexpected command: {argv!r}")
        cmd = argv[1]
        if cmd == "simulate":
            self._simulate(argv)
        elif cmd == "predict":
            self._predict(argv)
        elif cmd == "extract":
            self._extract(argv)
        elif cmd == "annotate":
            self._annotate(argv)
        elif cmd == "train":
            self._train(argv)
        else:
            raise AssertionError(f"unsupported enerzyme subcommand: {cmd}")
        return mock.Mock(wait=lambda: 0, returncode=0)

    def _flag(self, argv: List[str], flag: str) -> str:
        return argv[argv.index(flag) + 1]

    def _simulate(self, argv: List[str]) -> None:
        config_path = self._flag(argv, "-c")
        out_dir = self._flag(argv, "-o")
        with open(config_path) as handle:
            cfg = yaml.safe_load(handle)
        method = (
            cfg.get("Simulation", {})
            .get("plumed_config_generator", {})
            .get("method")
        )
        # Pre-sim is written as presimulation.yaml; steered as simulation.yaml.
        if Path(config_path).name == "presimulation.yaml":
            self.presim_methods.append(method)
            # Restrained plumed writes plumed.traj.xyz; launcher renames it.
            _write_traj(os.path.join(out_dir, "plumed.traj.xyz"), _PRESIM)
        else:
            self.steered_methods.append(method)
            _write_traj(os.path.join(out_dir, "plumed.traj.xyz"), _STEERED)

    def _predict(self, argv: List[str]) -> None:
        out_dir = self._flag(argv, "-o")
        # Match launcher glob: processed_dataset_*/{model}-prediction.pkl
        model_cfg = self._flag(argv, "-mc")
        with open(model_cfg) as handle:
            suffix = yaml.safe_load(handle)["Modelhub"]["internal_FFs"]["FF01"]["suffix"]
        # Model config written for prediction uses the iteration suffix already.
        # Glob pattern uses _get_model_str(i); easiest is one matching file via
        # scanning the out_dir name ..._prediction → model stem.
        stem = Path(out_dir).name.replace("_prediction", "")
        dest_dir = Path(out_dir) / "processed_dataset_0"
        dest_dir.mkdir(parents=True, exist_ok=True)
        with open(dest_dir / f"{stem}-prediction.pkl", "wb") as handle:
            pickle.dump({"E_var": _EVar()}, handle)
        del suffix  # kept for clarity / future assertions

    def _extract(self, argv: List[str]) -> None:
        out_dir = Path(self._flag(argv, "-o"))
        stem = out_dir.name.replace("_extraction", "")
        (out_dir / f"{stem}_fragments.sdf").write_text("fake\n")

    def _annotate(self, argv: List[str]) -> None:
        out_root = Path(self._flag(argv, "-o"))
        # annotation.yaml lives under extraction dir; output is {model}_fragments
        # Launcher sets -o to AL output_path and expects
        # ``{output}/{model_str(i)}_fragments/fragments.pkl``.
        # Discover the iteration model stem from Supplier path in annotation cfg.
        cfg_path = self._flag(argv, "-c")
        with open(cfg_path) as handle:
            cfg = yaml.safe_load(handle)
        supplier = Path(cfg["Supplier"]["path"])
        stem = supplier.name.replace("_fragments.sdf", "")
        ann_dir = out_root / f"{stem}_fragments"
        ann_dir.mkdir(parents=True, exist_ok=True)
        datapoints = [
            {"Za": np.array([6, 6, 6]), "Ra": _PRESIM.copy(), "Q": 0, "S": 0, "E": 0.0},
            {"Za": np.array([6, 6, 6]), "Ra": _STEERED.copy(), "Q": 0, "S": 0, "E": 0.1},
        ]
        with open(ann_dir / "fragments.pkl", "wb") as handle:
            pickle.dump(datapoints, handle)

    def _train(self, argv: List[str]) -> None:
        out_dir = Path(self._flag(argv, "-o"))
        # training path is ``{model_str(i+1)}_training``; enerzyme writes model dir
        # named ``{model_str(i+1)}`` inside it.
        model_name = out_dir.name.replace("_training", "")
        (out_dir / model_name).mkdir(parents=True, exist_ok=True)
        (out_dir / model_name / "ok").write_text("1")


class ActiveLearningPresimulationIntegrationTests(unittest.TestCase):
    def _make_workspace(self, tmp: str):
        root = Path(tmp)
        seed_xyz = root / "seed.xyz"
        _write_xyz(str(seed_xyz), _SEED)
        paths = _write_minimal_yamls(root, str(seed_xyz))
        pretrain = root / "pretrain"
        pretrain.mkdir(parents=True)
        # Launcher reads ``{pretrain}/config.yaml`` when -p is set.
        shutil.copy(paths["model"], pretrain / "config.yaml")
        (pretrain / "FF01-FakeNet").mkdir()
        (pretrain / "FF01-FakeNet" / "ok").write_text("1")
        out = root / "out"
        tmpdir = root / "tmp"
        out.mkdir()
        tmpdir.mkdir()
        return paths, str(pretrain), str(out), str(tmpdir)

    def _launch(self, paths, pretrain, out, tmp, *, n_iterations: int, n_presim: int):
        fake = _FakeEnerzyme()
        # Avoid importing fairchem via the uma calculator patch on CPU login nodes.
        calc_stub = Path(tmp) / "calc_stub.py"
        calc_stub.write_text("# stub calculator patch for AL smoke\n")
        with mock.patch(
            "enerzymette.altoolkit.launcher.get_calculator_patch",
            return_value=str(calc_stub),
        ):
            launcher = active_learning_launcher(
                output_path=out,
                tmp_path=tmp,
                calculator_patch_key="uma",
                plumed_patch_key="bond_reaction",
                simulation_config_path=paths["simulation"],
                extraction_config_path=paths["extraction"],
                annotation_config_path=paths["annotation"],
                training_config_path=paths["training"],
                model_config_path=paths["model"],
                pretrain_path=pretrain,
                n_iterations=n_iterations,
                n_presimulation_steps_per_iteration=n_presim,
                training_ratio=0.5,
            )
        with mock.patch.object(subprocess, "Popen", side_effect=fake):
            launcher.launch()
        return launcher, fake

    def test_full_loop_presim_updates_pool_not_steered(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            paths, pretrain, out, tmpdir = self._make_workspace(tmp)
            launcher, fake = self._launch(
                paths, pretrain, out, tmpdir, n_iterations=2, n_presim=8
            )

            pool_path = launcher.structure_pool_state["entries"][0]["path"]
            pool_xyz = ase.io.read(pool_path, index=-1)

            # Iteration 0: first use → no pre-sim. Steered must not rewrite pool.
            self.assertEqual(fake.presim_methods, ["standard_restrained_md"])
            self.assertEqual(
                fake.steered_methods,
                ["standard_steered_md", "standard_steered_md"],
            )
            # After iter1 pre-sim, pool must equal pre-sim coords (not steered).
            np.testing.assert_allclose(pool_xyz.positions, _PRESIM, atol=1e-8)
            self.assertFalse(np.allclose(pool_xyz.positions, _STEERED))

            # Pre-sim yaml on disk also records the switched method.
            presim_yaml = Path(out) / "FF01-FakeNet-1_md" / "presimulation.yaml"
            self.assertTrue(presim_yaml.is_file())
            with open(presim_yaml) as handle:
                presim_cfg = yaml.safe_load(handle)
            self.assertEqual(
                presim_cfg["Simulation"]["plumed_config_generator"]["method"],
                "standard_restrained_md",
            )
            self.assertEqual(presim_cfg["Simulation"]["task"], "plumed")

            # Full loop produced next model and completed flags.
            self.assertTrue((Path(out) / "FF01-FakeNet-1").exists())
            self.assertTrue((Path(out) / "FF01-FakeNet-2").exists())
            self.assertTrue(
                (Path(out) / "FF01-FakeNet-0_md" / "simulation_completed").exists()
            )
            self.assertTrue(
                (Path(out) / "FF01-FakeNet-1_md" / "presimulation_completed").exists()
            )

            with open(Path(out) / "structure_pool.json") as handle:
                state = json.load(handle)
            self.assertTrue(state["entries"][0]["run"])

    def test_zero_presim_leaves_pool_at_seed_after_steered(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            paths, pretrain, out, tmpdir = self._make_workspace(tmp)
            launcher, fake = self._launch(
                paths, pretrain, out, tmpdir, n_iterations=2, n_presim=0
            )
            self.assertEqual(fake.presim_methods, [])
            pool_xyz = ase.io.read(
                launcher.structure_pool_state["entries"][0]["path"], index=-1
            )
            np.testing.assert_allclose(pool_xyz.positions, _SEED, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
