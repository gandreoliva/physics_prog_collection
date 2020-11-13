import numpy as np
from vpython import *

"Animation of the trajectory of a particle (newtonian_dynamics.f90)"

data_type = np.dtype([("t", np.float32),\
                    ("theta", np.float32),\
                    ("phi", np.float32)])
data = np.fromfile("out1.dat", dtype=data_type)

ball1 = sphere(radius=0.1, make_trail=True, trail_color=color.red)
ball2 = sphere(radius=0.1, make_trail=True, trail_color=color.orange)
rope1 = cylinder(radius=0.02)
rope2 = cylinder(radius=0.02)

scene.range = 3

l = 1
m = 0.4

for i in range(0,data["t"].shape[0],2):
    rate(50)
    theta = data["theta"][i]
    phi = data["phi"][i]

    ball1.pos = vec(l*sin(theta),-l*cos(theta),0)
    ball2.pos = vec(l*sin(theta)+l*sin(phi),-l*cos(theta)-l*cos(phi),0)
    rope1.pos = vec(0,0,0)
    rope1.axis = ball1.pos
    rope2.pos = ball1.pos
    rope2.axis = ball2.pos - ball1.pos
