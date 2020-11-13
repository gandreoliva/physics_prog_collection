import matplotlib.pyplot as plt
import numpy as np
import sys

# Lab procedure
# Lane-Emden eqn. -->  hydro code  --> (( plotting scripts ))
r"""
Displays an animation of the density as a function of position from the hydro code
Usage:
    python animate_rh.py
(Click the matplotlib window and press Q to quit during the animation)
"""

out_dir = "data/"

t = np.loadtxt(f"{out_dir}/t.dat")
rh = np.loadtxt(f"{out_dir}/rh.dat")
ra = np.loadtxt(f"{out_dir}/ra.dat")
rb = np.loadtxt(f"{out_dir}/rb.dat")

iterations = t.shape[0]

print("Press 'q' to quit")
def press(event):
    if event.key == 'q':
        quit()

rh_min, rh_max = np.amin(rh[1:-1]),np.amax(rh[1:-1])
rb_min, rb_max = np.amin(rb[:]),np.amax(rb[:])

print(f"rh_max = {rh_max},  rh_min = {rh_min}")


for nt in np.arange(0,iterations-1,10):
    plt.cla()
    plt.clf()
    plt.xlim(rb_min,rb_max)
    plt.ylim(0,2)
    plt.xlabel('r [code units]')
    plt.ylabel('density [code units]')

    plt.plot(rb[:],rh[nt,:],'.-',label=f"{nt}")
    plt.pause(0.01)
    plt.gcf().canvas.mpl_connect('key_press_event', press)


quit()
