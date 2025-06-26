# Enerzymette
A collection of small scripts that are useful for Enerzyme

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
