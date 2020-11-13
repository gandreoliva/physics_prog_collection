import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sys
import scipy.interpolate

if len(sys.argv) < 2:
    print("Usage: python poisson_plot_plt.py datafile")
    sys.exit(0)

nr = 30
nph = 200

potential = np.loadtxt(sys.argv[1]).T
r = np.linspace(0.1,1,nr+1,endpoint=True)
ph = np.linspace(0,2*np.pi,nph+1,endpoint=True)

r_grid, ph_grid = np.meshgrid(r,ph)
x, y = r_grid*np.cos(ph_grid), r_grid*np.sin(ph_grid)

dr = np.diff(r)[0]
dph = np.diff(ph)[0]

# gradient in polar coords ((i): potential is transposed)
Er = (np.diff(potential, axis=1)/dr)[:-1,:]
Eph = (np.diff(potential,axis=0)/(r_grid[:-2,:-1]*dph))[:,:-1]

Ex = Er*np.cos(ph_grid[:-2,:-2]) - Eph*np.sin(ph_grid[:-2,:-2])
Ey = Er*np.sin(ph_grid[:-2,:-2]) + Eph*np.cos(ph_grid[:-2,:-2])

# interpolation of electric field in polar coords to a regular cartesian grid
xreg = np.linspace(-1,1,80)
yreg = np.linspace(-1,1,80)
xreg_grid, yreg_grid = np.meshgrid(xreg,yreg)
dist_reg_grid = np.sqrt(xreg_grid**2+yreg_grid**2)

Exreg = scipy.interpolate.griddata((x[:-2,:-2].flatten(), y[:-2,:-2].flatten()), Ex.flatten(), (xreg_grid, yreg_grid), method='linear')
Eyreg = scipy.interpolate.griddata((x[:-2,:-2].flatten(), y[:-2,:-2].flatten()), Ey.flatten(), (xreg_grid, yreg_grid), method='linear', )

# get rid of the inner boundary
Exreg[dist_reg_grid < 0.1] = np.nan
Eyreg[dist_reg_grid < 0.1] = np.nan
Enorm = np.sqrt(Exreg**2 + Eyreg**2)


plt.pcolormesh(x,y, potential)
plt.colorbar()
# plt.quiver(x[:-2,:-2],y[:-2,:-2], (Ex/Enorm), (Ey/Enorm)) # polar
streamplot_lw = 0.5+5*Enorm/np.amax(Enorm[~np.isnan(Enorm)])
plt.streamplot(xreg,yreg, Exreg, Eyreg, color='k', density=3, linewidth=streamplot_lw)
# plt.quiver(xreg_grid,yreg_grid,Exreg/Enorm,Eyreg/Enorm, scale=50, alpha=0.2) # cart
plt.show()
