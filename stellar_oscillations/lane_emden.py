from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# Running order
# (( Lane-Emden eqn. )) -->  hydro code  --> plotting scripts
"""
Solver for the Lane-Emden equation (hydrostatic stellar structure).
This is used as an initial condition for the hydrodynamical stellar structure
code.
"""



header=r"""
--------------------------------------------------------------------------
 Lane-Emden equation   *  . ·                > Stellar equilibrium config.
  ::: Adv Lab Astron Astroph - CPT - Uni Tübingen - 2020
==========================================================================
"""

print(header)

# ..........................................................................
# 1. Numerical solution (RK2)
# --------------------------------------------------------------------------

r"""
System of equations:
    dw/dz = xi
    dxi/dz = - 2*xi/z - w**n

In order to write this in the form Y' = F(X,Y(X)), let's define:
    X  (scalar) := z
    Y  (array)  := ( w     , xi             )
    Y' (array)  := ( dw/dz , dxi/dz         )
    F  (array)  := ( xi    , -2*xi/z - w**n )

Heun method (2nd order Runge-Kutta):

    K1 = F(X,   Y)
    K2 = F(X + dX,   Y + dX*K1 )
    Y = Y + dX/2 * (K1 + K2)  <-- new Y
    X = X + dX                <-- new X

"""

## Polytropic index|
n = 1.5


def F(X,Y):
    z = X
    w = Y[0]
    xi = Y[1]
    return np.array([xi, -2*xi/z - w**n ])

## Initial conditions (at z = 0)
z = 0
w = 1
xi = 0

dz = 0.0001
wmin = 1e-6 # minimum value allowed for the adimensional density
# wmin = 1e-20


# Convert to the variables of the num. solution
# (this is done for dydactical reasons only!)
X = z
Y = np.array([w, xi])
dX = dz

# Arrays that collect the solution
ws = []
zs = []
xis = []

# while the adim. density is higher than 'zero'...
# (at zero density, we hit the radius of the star)

while Y[0] > wmin:
    X = X + dX
    K1 = F(X, Y)
    K2 = F(X+dX, Y+dX*K1)
    Y = Y + dX/2 * (K1 + K2)

    # store the results
    ws.append(Y[0])
    zs.append(X)


zs = np.array(zs)
ws = np.array(ws)

print(f"Lane-Emden equation solved for n = {n}.")

# 8< - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Interpolation and output of the results (for the hydro code)
# --------------------------------------------------------------------------
Rstar = zs[-1]
interp_func = interpolate.interp1d(zs, ws)

## Parameters from the hydro simulation|
rmin = 0
rmax = 1.5*Rstar
nr = 100

dr = (rmax - rmin)/nr

# (i)  We want rb[i] = ((i-1) + 0.5)*dr, for i from 1 to nr
#      but arange doesn't include the endpoints, so we add 1 to the right limit
hydro_grid = (np.arange(-1,nr+1)+0.5)*dr
# w_grid: array filled with the density floor value
w_grid = np.zeros(hydro_grid.shape) + 2*wmin

# value for the left ghost cell
w_grid[0] = interp_func(hydro_grid[1])

# w_grid: set the values of the density with the solution found
for i in range(hydro_grid.shape[0]):
    if hydro_grid[i] < Rstar and hydro_grid[i] > 0:
        w_grid[i] = interp_func(hydro_grid[i])


with open("equil_model.dat",'w') as outfile:
    print(float(n),file=outfile)
    print(Rstar,file=outfile)
    for w_i in w_grid:
        print(w_i**n,file=outfile)

print(f"Output for a grid with nr = {nr}, rmin = {rmin:.2g}, rmax = {rmax:.2g}.")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - >8



# 8< - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. Output of results for the Chandrasekhar limit
# --------------------------------------------------------------------------
print("--- Chandrasekhar limit ---")
print("   (only works for n = 3)")
print(f" z1 = {zs[-1]}")
print(f" xi1 = {Y[1]}")
k0 = zs[-1]**2*Y[1]
print(f" z1**2*xi1 = {k0} ")
from scipy.constants import G, h, c, m_p
M_sun = 1.98847e30 # kg
mu_e = 2
MCh = np.sqrt(6)/(32*np.pi) * (h*c/G)**(3/2) * (2/mu_e)**2 * np.abs(k0)/(m_p**2)
print(f" Chandrasekhar mass = {MCh/M_sun:.4g} M_sun")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - >8


# 8< - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Output of results in physical units
# --------------------------------------------------------------------------
print("--- Physical units, Sun-like star ---")
from scipy.constants import G, h, c, m_p
M_sun = 1.98847e30 # kg
R_sun = 6.95700e8 # m
mean_dens_sun = M_sun/(4/3*np.pi*R_sun**3)

A = zs[-1]/R_sun
central_density = abs(zs[-1]/Y[1])*mean_dens_sun

K_gas = (4*np.pi*G)/((n+1)*A**2)*central_density**(1-1/n)
print(f" K_gas = {K_gas:.3g} [SI units]")
print(f" stellar radius = {(zs[-1]/A)/R_sun:.4g} R_sun")
print(f" central density = {central_density:.3g} kg/m^3 = {central_density/1000:.3g} g/cm^3")
print(f" central density / avg density = {central_density/mean_dens_sun:.4g}")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - >8


# ..........................................................................
# 5. Plotting
# --------------------------------------------------------------------------
plot_mode = "grid"

if plot_mode == "code":
    plt.xlabel("z [code units]")
    plt.ylabel("w [code units]")
    plt.plot(zs, ws)

# depends on part 3
if plot_mode == "grid":
    plt.xlabel("z [code units]")
    plt.ylabel("w [code units]")
    plt.plot(hydro_grid, w_grid)

# depends on part 4
if plot_mode == "physical":
    plt.rc('mathtext', fontset='cm') # allows LaTeX
    plt.xlabel(r"r [$R_\odot$]")
    plt.ylabel(r"$\rho$ [$\bar \rho_\odot$]")
    plt.plot((zs/A)/R_sun, (central_density*ws**n)/mean_dens_sun)


print(" ")

plt.show()
