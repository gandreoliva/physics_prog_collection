# Lane-Emden equation

Numerical integration of the Lane-Endem equation.

## Requirements

Library `shtawa` (https://github.com/gandreoliva/shtawa)

## Workflow
* Numerical integration and explanation of the theory: `lane_endem.f90`
    * Edit polytropic index
    * Compile with `make`
    * Run with `./lane_endem.bin > structure.dat`
* The results can be plotted with the notebook `analysis.ipynb` (edit data file name)