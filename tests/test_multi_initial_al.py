"""Tests for multi-system initial structure pool in active learning."""

from __future__ import annotations

import inspect
import json
import os
import tempfile
import unittest
from unittest import mock

import ase.io
import yaml
from ase import Atoms

from enerzymette.altoolkit.launcher import active_learning_launcher
from enerzymette.altoolkit.structures_manifest import load_systems_manifest


def _write_xyz(path: str, n_atoms: int = 3) -> None:
    atoms = Atoms("H" * n_atoms, positions=[[0, 0, i] for i in range(n_atoms)])
    ase.io.write(path, atoms, format="extxyz")


def _write_simulation_yaml(
    path: str,
    structure_file: str,
    reference_pdb_file: str | None = None,
    proton_transfer: bool = False,
) -> None:
    plumed_config = {
        "dump_interval": 20,
        "lower_bound": -2,
        "upper_bound": 2,
        "substrate": "X",
        "nucleophile": "Y",
    }
    if reference_pdb_file is not None:
        plumed_config["reference_pdb_file"] = reference_pdb_file
    if proton_transfer:
        plumed_config["proton_transfer"] = True

    config = {
        "System": {
            "structure_file": structure_file,
            "charge": 0,
            "multiplicity": 1,
        },
        "Simulation": {
            "task": "plumed",
            "idx_start_from": 1,
            "integrate": {
                "integrator": "Langevin",
                "n_step": 10,
                "time_step": 0.5,
                "friction": 0.01,
                "temperature_in_K": 300,
            },
            "plumed_config_generator": {"name": "get_sammt_config"},
            "sampling": {
                "params": {
                    "plumed_config": plumed_config,
                }
            },
            "uncertainty_calculator": {
                "name": "UDD",
                "params": {"A": 0.08, "B": 1.0e-6},
            },
        },
    }
    with open(path, "w") as handle:
        yaml.dump(config, handle, default_flow_style=False)


class StructuresManifestTests(unittest.TestCase):
    def test_load_manifest_without_explicit_reference_xyz(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "ref.pdb")
            xyz_path = os.path.join(tmpdir, "init.xyz")
            sim_path = os.path.join(tmpdir, "simulation.yaml")
            manifest_path = os.path.join(tmpdir, "systems.yaml")

            open(pdb_path, "w").close()
            _write_xyz(xyz_path, n_atoms=4)
            _write_simulation_yaml(sim_path, xyz_path, reference_pdb_file=pdb_path)
            with open(manifest_path, "w") as handle:
                yaml.dump(
                    {
                        "systems": [
                            {
                                "name": "sys_a",
                                "reference_pdb": pdb_path,
                                "simulation_config": sim_path,
                            }
                        ]
                    },
                    handle,
                )

            entries = load_systems_manifest(manifest_path)
            self.assertEqual(len(entries), 1)
            self.assertEqual(entries[0].name, "sys_a")
            self.assertIsNone(entries[0].reference_xyz)
            self.assertFalse(entries[0].reference_xyz_explicit)

    def test_load_manifest_with_explicit_reference_xyz(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "ref.pdb")
            xyz_path = os.path.join(tmpdir, "init.xyz")
            sim_path = os.path.join(tmpdir, "simulation.yaml")
            manifest_path = os.path.join(tmpdir, "systems.yaml")

            open(pdb_path, "w").close()
            _write_xyz(xyz_path, n_atoms=4)
            _write_simulation_yaml(sim_path, xyz_path, reference_pdb_file=pdb_path)
            with open(manifest_path, "w") as handle:
                yaml.dump(
                    {
                        "systems": [
                            {
                                "name": "sys_a",
                                "reference_pdb": pdb_path,
                                "reference_xyz": xyz_path,
                                "simulation_config": sim_path,
                            }
                        ]
                    },
                    handle,
                )

            entries = load_systems_manifest(manifest_path)
            self.assertEqual(entries[0].reference_xyz, os.path.abspath(xyz_path))
            self.assertTrue(entries[0].reference_xyz_explicit)

    def test_proton_transfer_raises_not_implemented(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "ref.pdb")
            xyz_path = os.path.join(tmpdir, "init.xyz")
            sim_path = os.path.join(tmpdir, "simulation.yaml")
            manifest_path = os.path.join(tmpdir, "systems.yaml")

            open(pdb_path, "w").close()
            _write_xyz(xyz_path)
            _write_simulation_yaml(
                sim_path, xyz_path, reference_pdb_file=pdb_path, proton_transfer=True
            )
            with open(manifest_path, "w") as handle:
                yaml.dump(
                    {
                        "systems": [
                            {
                                "name": "sys_a",
                                "reference_pdb": pdb_path,
                                "simulation_config": sim_path,
                            }
                        ]
                    },
                    handle,
                )

            with self.assertRaises(NotImplementedError):
                load_systems_manifest(manifest_path)


class ActiveLearningLauncherMultiSystemTests(unittest.TestCase):
    def _make_launcher(self, tmpdir: str, manifest_path: str) -> active_learning_launcher:
        sim_path = os.path.join(tmpdir, "fallback_simulation.yaml")
        _write_xyz(os.path.join(tmpdir, "fallback.xyz"))
        _write_simulation_yaml(sim_path, os.path.join(tmpdir, "fallback.xyz"))
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
            initial_structures_config_path=manifest_path,
        )

    def _write_two_system_manifest(self, tmpdir: str) -> str:
        pdb_a = os.path.join(tmpdir, "a.pdb")
        pdb_b = os.path.join(tmpdir, "b.pdb")
        xyz_a = os.path.join(tmpdir, "a.xyz")
        xyz_b = os.path.join(tmpdir, "b.xyz")
        sim_a = os.path.join(tmpdir, "sim_a.yaml")
        sim_b = os.path.join(tmpdir, "sim_b.yaml")
        manifest_path = os.path.join(tmpdir, "systems.yaml")

        open(pdb_a, "w").close()
        open(pdb_b, "w").close()
        _write_xyz(xyz_a, n_atoms=3)
        _write_xyz(xyz_b, n_atoms=5)
        _write_simulation_yaml(sim_a, xyz_a, reference_pdb_file=pdb_a)
        _write_simulation_yaml(sim_b, xyz_b, reference_pdb_file=pdb_b)
        with open(manifest_path, "w") as handle:
            yaml.dump(
                {
                    "systems": [
                        {
                            "name": "sys_a",
                            "reference_pdb": pdb_a,
                            "reference_xyz": xyz_a,
                            "simulation_config": sim_a,
                        },
                        {
                            "name": "sys_b",
                            "reference_pdb": pdb_b,
                            "reference_xyz": xyz_b,
                            "simulation_config": sim_b,
                        },
                    ]
                },
                handle,
            )
        return manifest_path

    def test_topology_paths_are_namespaced_per_system(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            launcher = self._make_launcher(tmpdir, manifest_path)
            mol_a = launcher._cluster_mol_path_for_system("sys_a")
            mol_b = launcher._cluster_mol_path_for_system("sys_b")
            self.assertEqual(
                mol_a,
                os.path.join(tmpdir, "topology", "sys_a", "cluster.mol"),
            )
            self.assertEqual(
                mol_b,
                os.path.join(tmpdir, "topology", "sys_b", "cluster.mol"),
            )
            self.assertNotEqual(mol_a, mol_b)

    def test_init_structure_pool_from_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            launcher = self._make_launcher(tmpdir, manifest_path)
            launcher._init_structure_pool()

            state_path = os.path.join(tmpdir, "structure_pool.json")
            self.assertTrue(os.path.exists(state_path))
            with open(state_path, "r") as handle:
                state = json.load(handle)
            self.assertEqual(len(state["entries"]), 2)
            self.assertEqual(state["entries"][0]["name"], "sys_a")
            self.assertEqual(state["entries"][1]["name"], "sys_b")
            self.assertTrue(os.path.exists(state["entries"][0]["path"]))
            self.assertTrue(os.path.exists(state["entries"][1]["path"]))
            self.assertEqual(
                state["entries"][0]["cluster_mol_path"],
                os.path.join(tmpdir, "topology", "sys_a", "cluster.mol"),
            )

    def test_pool_rotation_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            launcher = self._make_launcher(tmpdir, manifest_path)
            launcher.structure_pool_state = {
                "entries": [
                    {"path": "a.xyz", "run": False},
                    {"path": "b.xyz", "run": False},
                ]
            }
            idx0, entry0 = launcher._get_structure_pool_entry(0)
            idx1, entry1 = launcher._get_structure_pool_entry(1)
            idx2, entry2 = launcher._get_structure_pool_entry(2)
            self.assertEqual((idx0, entry0["path"]), (0, "a.xyz"))
            self.assertEqual((idx1, entry1["path"]), (1, "b.xyz"))
            self.assertEqual((idx2, entry2["path"]), (0, "a.xyz"))

    def test_load_simulation_config_for_entry_injects_reference_pdb(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            launcher = self._make_launcher(tmpdir, manifest_path)
            entry = {
                "simulation_config": os.path.join(tmpdir, "sim_a.yaml"),
                "reference_pdb": os.path.join(tmpdir, "a.pdb"),
            }
            config = launcher._load_simulation_config_for_entry(entry)
            plumed_config = config["Simulation"]["sampling"]["params"]["plumed_config"]
            self.assertEqual(
                plumed_config["reference_pdb_file"],
                os.path.abspath(entry["reference_pdb"]),
            )

    def test_global_reference_flags_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            sim_path = os.path.join(tmpdir, "fallback_simulation.yaml")
            with self.assertRaises(ValueError):
                active_learning_launcher(
                    output_path=tmpdir,
                    tmp_path=tmpdir,
                    calculator_patch_key="uma_calculator",
                    plumed_patch_key="sammt",
                    simulation_config_path=sim_path,
                    extraction_config_path=os.path.join(tmpdir, "extraction.yaml"),
                    annotation_config_path=os.path.join(tmpdir, "annotation.yaml"),
                    training_config_path=os.path.join(tmpdir, "training.yaml"),
                    n_iterations=1,
                    initial_structures_config_path=manifest_path,
                    reference_pdb_path=os.path.join(tmpdir, "a.pdb"),
                )

    def test_initial_scan_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            sim_path = os.path.join(tmpdir, "fallback_simulation.yaml")
            with self.assertRaises(NotImplementedError):
                active_learning_launcher(
                    output_path=tmpdir,
                    tmp_path=tmpdir,
                    calculator_patch_key="uma_calculator",
                    plumed_patch_key="sammt",
                    simulation_config_path=sim_path,
                    extraction_config_path=os.path.join(tmpdir, "extraction.yaml"),
                    annotation_config_path=os.path.join(tmpdir, "annotation.yaml"),
                    training_config_path=os.path.join(tmpdir, "training.yaml"),
                    n_iterations=1,
                    initial_structures_config_path=manifest_path,
                    initial_scan=True,
                )

    def test_get_multi_system_topology_runs_bond_per_system(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            launcher = self._make_launcher(tmpdir, manifest_path)

            def _fake_bond(pdb_path, img_path, mol_path, reference_sdf=None):
                os.makedirs(os.path.dirname(mol_path), exist_ok=True)
                open(img_path, "w").close()
                open(mol_path, "w").close()

            with mock.patch.object(
                launcher, "_run_enerzyme_bond", side_effect=_fake_bond
            ) as mock_bond, mock.patch(
                "rdkit.Chem.MolFromMolFile", return_value=mock.Mock()
            ), mock.patch(
                "rdkit.Chem.GetFormalCharge", return_value=0
            ), mock.patch(
                "enerzymette.altoolkit.launcher.get_indices",
                return_value=[],
            ):
                launcher._get_multi_system_topology()
                self.assertEqual(mock_bond.call_count, 2)
                calls = mock_bond.call_args_list
                self.assertEqual(
                    calls[0].args[1],
                    os.path.join(tmpdir, "topology", "sys_a", "cluster.png"),
                )
                self.assertEqual(
                    calls[0].args[2],
                    os.path.join(tmpdir, "topology", "sys_a", "cluster.mol"),
                )
                self.assertEqual(
                    calls[1].args[1],
                    os.path.join(tmpdir, "topology", "sys_b", "cluster.png"),
                )

    def test_multi_system_reuse_reads_pool_entry_not_previous_iteration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = self._write_two_system_manifest(tmpdir)
            launcher = self._make_launcher(tmpdir, manifest_path)
            pool_dir = os.path.join(tmpdir, "structure_pool")
            os.makedirs(pool_dir, exist_ok=True)
            entry_path = os.path.join(pool_dir, "001.xyz")
            initial_structure_path = os.path.join(tmpdir, "initial_structure.xyz")
            _write_xyz(entry_path, n_atoms=7)

            launcher.structure_pool_state = {
                "entries": [
                    {"path": os.path.join(pool_dir, "000.xyz"), "run": True},
                    {"path": entry_path, "run": True},
                ]
            }
            launcher.n_presimulation_steps_per_iteration = 0

            with mock.patch.object(launcher, "_get_simulation_path") as mock_paths:
                launcher._prepare_iteration_initial_structure(
                    1,
                    initial_structure_path,
                    tmpdir,
                    os.path.join(tmpdir, "presim.yaml"),
                    os.path.join(tmpdir, "md.traj.xyz"),
                    os.path.join(tmpdir, "presim_done"),
                    os.path.join(tmpdir, "model_config.yaml"),
                    {},
                )
                mock_paths.assert_not_called()

            copied = ase.io.read(initial_structure_path, index=-1)
            self.assertEqual(len(copied), 7)

    def test_steered_md_does_not_update_structure_pool(self) -> None:
        """Pool writeback after steered MD was removed on purpose."""
        source = inspect.getsource(active_learning_launcher._simulation_step)
        self.assertNotIn('updated pool entry', source)
        self.assertNotIn('from steered MD', source)
        # The only ase.io.write(entry["path"], ...) lives in pre-simulation prep.
        prep = inspect.getsource(
            active_learning_launcher._prepare_iteration_initial_structure
        )
        self.assertIn('ase.io.write(entry["path"]', prep)
        self.assertNotIn('ase.io.write(entry["path"]', source)


if __name__ == "__main__":
    unittest.main()
