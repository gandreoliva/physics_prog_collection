import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sys

if len(sys.argv) < 2:
    print("Usage: python poisson_plot_plt.py datafile")
    sys.exit(0)

data = np.loadtxt(sys.argv[1])
n = data.shape[0]+1
x = np.linspace(0,1,n,endpoint=True)
y = np.linspace(0,1,n,endpoint=True)

cellsize = np.diff(x)[0]

vy, vx = np.gradient(data)
vy[np.sqrt(vx**2+vy**2)>0.5] = 0
vx[np.sqrt(vx**2+vy**2)>0.5] = 0

plt.pcolormesh(x,y, data, vmin=0, vmax=8,cmap='PuBu_r')
plt.quiver(x[:-1]+cellsize/2, y[:-1]+cellsize/2, -vx, -vy)
plt.show()
