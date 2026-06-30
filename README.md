# Enerzymette
A collection of small scripts that are useful for Enerzyme. Install it in the cloned git repo by
`pip install -e .`

## IDPP
Transition path interpolation between reactants and products with IDPP method implemented in [ASE](https://wiki.fysik.dtu.dk/ase/ase/neb.html#ase.mep.neb.idpp_interpolate).

(1) Smidstrup, S.; Pedersen, A.; Stokbro, K.; Jónsson, H. Improved Initial Guess for Minimum Energy Path Calculations. The Journal of Chemical Physics 2014, 140 (21), 214106. https://doi.org/10.1063/1.4878664.

Usage:
```
enerzymette idpp -r <reactant xyz path> -p <product xyz path> -o <output xyz path> -n <number of images> -c <terachem input file that contains fix atom constraints>
```

If multiple frames are stored in the reactant or the product xyz file (for example, an optimization trajectory), the last frame will be used.

## Terachem timing
Calculate total wall time from TeraChem output files, handling multiple iteration tables and identifying unfinished calculations.

The script analyzes TeraChem output files to extract timing information from iteration tables, providing detailed statistics including total wall time, average time per iteration, and identification of unfinished calculations.

Usage:
```
enerzymette terachem_timing -f <terachem output file path>
```

## Request Terachem calculation for ORCA optimizer

ORCA can be used as an external optimizer that calls energy and gradient calculation in another quantum chemistry software. See:
https://www.faccts.de/docs/orca/6.1/tutorials/workflows/extopt.html

Here this script is combining the advance of ORCA's optimizer and the speed of Terachem's GPU-accelerated quantum chemistry. The script generates a TeraChem input and run a TeraChem single-point gradient calculation using geometry and settings provided by an ORCA optimizer (via an `.extinp.tmp` file) by reading it, applying a TeraChem input template, runnning the calculation, and writing the results in ORCA-compatible format (via an `.engrad` file). In subsequence optimization cycles, the script cleans the scratch directory and uses the molecular orbital for the scf initial guess in the next cycle automatically.

Usage:
```
enerzymette orca_terachem_request -i <orca extinp.tmp file> -t <terachem input template>
```

- The `-i`/`--input` argument specifies the ORCA external input file (typically ending with `.extinp.tmp`).
- The `-t`/`--template` argument specifies a TeraChem input template file, which provides the calculation settings.

For this purpose, besides enerzymette, you need to prepare

1. Your orca input file including
    ```
    ! ExtOpt # The simple input keyword "ExtOpt" means ORCA is used as an external optimizer.
             # Other simple input keywords for the optimization task, like `Opt` or `NEB-CI`, should be also present as usual.

    %method  # replace the method block for quantum chemistry with this
    ProgExt "/full/path/to/your/enerzymette/wrapper"
    end
    ```
2. A template Terachem input file that stores the quantum chemistry parameters, like basis set, functional, solvent model, scf convergence criterions, etc.
3. All relevant file to your quantum chemistry parameters in the correct place, for example, when the `pcm_radii` is set to `read`, the `pcm_radii_file` should be set accordingly as well.
4. Your enerzymette wrapper `.sh` file that contains
    ```
    enerzymette orca_terachem_request -i $1 -t <your_template_terachem_input_file>
    ```
5. Setting the permission of your enerzymette wrapper by
    ```
    chmod +x <your_enerzymette_wrapper>
    ```
    Then `<your_enerzymette_wrapper> <orca_generated_extinp_tmp_file>` will request the Terachem calculation and collect the results as ORCA requires 
6. Make sure the `terachem` is an executable in your environment.

Then run the ORCA job as usual.

## Enerzyme Scan Launcher

Automates NNP-driven flexible bond scans. For each elementary reaction it calls `enerzyme simulate` to optimize the reactant, scan a bond-distance CV, optimize the product, and analyze the energy profile. System settings (charge, spin, frozen atoms, scan bond) are read from a TeraChem reference input (`-q`). If an intermediate minimum appears left of the transition state, the launcher chains another scan from that structure. Final results go to `rate_determining_ts/` (TS geometry and `results.json`).

Usage:
```
enerzymette enerzyme_scan \
    -r <reactant.xyz> \
    -o <output_dir> \
    -m <model_dir> \
    -q <scan.in> \
    -n 25
```

- `-r` — initial reactant structure (typically a DFT-optimized reactant).
- `-q` — TeraChem input with `constraint_freeze` and `constraint_scan` (bond indices); parsed by `terachem/io.py`.
- `-n` — number of scan steps along the bond CV.

**Supporting utilities:**

- `scantoolkit/launcher.py` — `EnerzymeScanLauncher` (opt → scan → opt loop and reaction chaining).
- `scantoolkit/io.py` — `update_terachem_scan` CLI to refresh coordinates in a scan input after structure update.
- `altoolkit/get_index.py` — select backbone/Cα atom indices from PDB for constraint setup.
- `mep_util.py` — CI, intermediate, and product-frame detection on the scan path.

Also supports PLUMED CV scans (`-pp`, `-psc`), YAML `scan_config` as an alternative to TeraChem input, and `altoolkit/launcher.py` for active-learning steered MD.

## Enerzyme NEB Launcher

Automates NNP-driven NEB: ORCA optimizes the image chain while Enerzyme supplies energies and gradients via an external optimizer (same ExtOpt pattern as [orca_terachem_request](#request-terachem-calculation-for-orca-optimizer)). The launcher starts `enerzyme listen`, writes ORCA NEB input plus a wrapper calling `enerzyme request`, monitors stdout for intermediate minima, and chains sub-reactions when needed. Requires `ORCA_PATH`.

Usage:
```
enerzymette enerzyme_neb \
    -r <reactant.xyz> \
    -p <product.xyz> \
    -o <output_dir> \
    -m <model_dir> \
    -q <reference.in> \
    -c <server.yaml> \
    -n 8 -b 5001 -i stdout
```

- `-q` — Reference TeraChem input file for charge, spin, and `constraint_freeze` atoms.
- `-c` — Enerzyme server config (`cuda`, `dtype`, `neighbor_list`, …).
- `-n` — number of NEB images (including endpoints).
- `-b` — server port; `-i stdout` interrupts ORCA when an intermediate minimum is detected.

**Supporting utilities:**

- `nebtoolkit/launcher.py` — `EnerzymeNEBLauncher` (server lifecycle, ORCA subprocess, restart/backtrack, chaining).
- `nebtoolkit/io.py` — `write_orca_neb_in` for ORCA NEB input generation; MEP parsing from ORCA output.
- `nebtoolkit/analysis.py`, `nebtoolkit/network.py` — convergence/CI checks and port selection.
- `altoolkit/get_index.py` — backbone atom indices for frozen constraints in NEB input.
- `mep_util.py` — intermediate and rate-determining-step detection (shared with the scan launcher).

Also configurable: spring constants (`--min_spring_constant` / `--max_spring_constant`), ORCA optimizer (`--optimization_method`), and `-i none` to disable early interruption.


