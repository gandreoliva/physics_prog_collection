# Simple 1D hydro code - stellar oscillations

1D hydrodynamics grid code based on a staggered grid (ZEUS-like), and application to the problem of stellar oscillations. Originally developed for the Advanced Lab. in Astrophysics, University of Tübingen.

## Workflow

For testing the hydro code with the propagation of an isothermal shock:
* Compile with `gfortran -o hydro_shock_test hydro_shock_test.f90` (see internal documentation and the book Bodenheimer, Laughlin, Rozyczka and Yorke (2007) "Numerical Methods in Astrophysics", Taylor and Francis. Sections 6.3.2 and 6.5.)
* Create the directory `data/` for data output
* Run as `./hydro_shock_test`.
* The output can be plotted with the notebook `hydro_plot.ipynb`

For simulating stellar oscillations:
* First, generate an equilibrium configuration with `python lane_emden.py` (see internal documdentation)
* Second, introduce perturbations that become oscillations with the hydro code. Compile with `gfortran -o hydro hydro.f90`
* Create the directory `data/` for data output
* Run as `./hydro`
* The output can be plotted in a similar way than what is done in `hydro_plot.ipynb`
