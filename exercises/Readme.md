# Exercises in computational physics

TODO: each exercise is independent of each other and should be in their own directory. Warning: there are some incomplete or outdated exercises.

List of exercises:

* `newtonian_dynamics*` : oscillator with damping. (Warning: solution with the Euler method!)
* `analytical_mech_dpend*` : double pendulum. Workflow:
    * `maxima -b  analytical_mech_dpend.mac` to produce the equations of motion from the (analytical) Lagrangian formalism (depends on Maxima)
    * solve numerically with the .f90 file
    * visualize with the .py file (Warning: depends on vpython!)
* `poisson_freeflow*` : solution of the Poisson equation for free flow and a simple obstacle (the obstacle files are in `poisson_obstacles/`). An in-depth documentation is in `poisson_doc/poisson_doc.tex`. There is a python script for plotting the result of the Fortran programs.
* `electrst_conduct_*` : it should solve the Poisson equation in polar coordinates. Warning: not properly tested!
* `classical_1d_wave*` : solution of a 1d wave equation (f90 program) and plotting script for the results (py)
* `angular_momentum_sum.py` : simple script for generating the table sum of two quantum mechanical angular momentum states