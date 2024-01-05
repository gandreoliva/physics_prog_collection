# Neutron star interior

Solves the Tolman-Oppenheimer-Volkoff (TOV) equation in a fixed grid, assuming a polytrope. The TOV equation is derived analytically with Maxima in `tov.wxm` (open with wxmaxima or run directly from maxima in the command line).

## Workflow
For dydactical purposes, there are two source files that solve the same equations:
    * `tov_structure.f90` outputs the internal structure of the neutron star (density, pressure, etc.)
    * `tov_metric.f90` outputs also the metric components in the grid

* Compile each program with `make structure` or `make metric`.
* Run with `./tov_structure.bin < parameters/inputfile.inp > structure.dat` and `./tov_metric.bin < parameters/inputfile.inp > metric.dat`, where `inputfile.inp` is one of the files in the `parameters/` directory.
* The output can be plotted and analyzed with the Jupyterlab notebook `analysis.ipynb`.