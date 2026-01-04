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


## Enerzyme NEB Launcher


