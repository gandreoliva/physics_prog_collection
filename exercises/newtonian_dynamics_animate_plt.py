import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

"Animation of the trajectory of a particle (newtonian_dynamics.f90)"

data_type = np.dtype([("t", np.float32),\
                    ("x", np.float32),\
                    ("y", np.float32),\
                    ("z", np.float32),\
                    ("vx", np.float32),\
                    ("vy", np.float32),\
                    ("vz", np.float32)])
data = np.fromfile("out1.dat", dtype=data_type)


bim = plt.imread("assets/man.png") # background image
oi = OffsetImage(plt.imread("assets/ball.png"),zoom=0.5) # image for points

for i in range(0,data["t"].shape[0],2):
    # Background still image
    plt.imshow(bim, origin="upper", extent=[-2,-1.3,0,2])

    plt.scatter([data["x"][i]],[data["y"][i]])

    # Subtitutes points by images
    ab = AnnotationBbox(oi, (data["x"][i], data["y"][i]), frameon=False)
    plt.gca().add_artist(ab)


    # Animation settings
    plt.pause(0.001) # controls frames per second
    plt.clf() # clears figure for next frame
    plt.cla() # clears artist for next frame
    plt.xlim(-2,2) # scene x range
    plt.ylim(-2,2) # scene y range
    #plt.axis("off")

plt.show()
