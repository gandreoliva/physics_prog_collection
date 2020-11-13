import numpy as np
import matplotlib.pyplot as plt

"Animation of the results of classical_1d_wave.f90"

data_type = np.dtype([('dt',float),('dx',float),\
                    ('nt',int),('nx',int)])

params = np.genfromtxt('out2.dat',dtype=data_type)
dt = params['dt']
dx = params['dx']
nt = params['nt']+1
nx = params['nx']+1

data = np.fromfile("out1.dat", dtype=np.float32)
data = data.reshape(nt,nx)
# Rows: j (t), columns: i (x)

x = np.arange(0,nx*dx,dx)


for j in range(0,nt,1):

    plt.plot(x,data[j,:])

    # Animation settings
    plt.pause(0.0005) # controls frames per second
    plt.clf() # clears figure for next frame
    plt.cla() # clears artist for next frame
    plt.xlim(0,0.3) # scene x range
    plt.ylim(-0.05,0.05) # scene y range
    #plt.axis("off")

plt.show()
