import numpy as np
import matplotlib.pyplot as plt

"Plot phase space for the sim. newtonian_dynamics.f90"

data_type = np.dtype([("t", np.float32),\
                    ("x", np.float32),\
                    ("y", np.float32),\
                    ("z", np.float32),\
                    ("vx", np.float32),\
                    ("vy", np.float32),\
                    ("vz", np.float32)])
data = np.fromfile("out1.dat", dtype=data_type)

plt.scatter(data["x"],data["vx"])

plt.show()
