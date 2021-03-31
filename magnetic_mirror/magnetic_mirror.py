from vpython import *
import numpy as np
import matplotlib.pyplot as plt
import random

N = 20
kT = 10 # eV
L = 1 # tamaño inicial de la caja
rL0 = 0.001 # radio de Larmor de referencia
# Tipo de campo magnético
Btype = "dipole" # dipole, helmholtz, uniform, no

scene.range = 2*L
scene.width = 1200
scene.height = 600

# constantes
m_p = 1.67262192369e-27 # masa del protón, kg
qe = 1.602176634e-19 # carga del electrón, C
kB = 1.380649e-23 # constante de Boltzmann, SI
eV = qe # electronvolt, J

# escalas y unidades
B0 = sqrt(m_p*kT*eV)/(qe*rL0)
v0 = sqrt(kT*eV/m_p) # most probable speed
print(f"B0 = {B0}")
print(f"v0 = {v0}")


def randvec(L,samples=10):
    x = L*random.randint(-samples,samples)/(samples)
    y = L*random.randint(-samples,samples)/(samples)
    z = L*random.randint(-samples,samples)/(samples)
    return vec(x,y,z)


def get_maxwellian_speed():
    # http://astro.physics.ncsu.edu/urca/course_files/Lesson20/index.html
    vmin = 0.5 # units of v0
    vmax = 1.5
    # Generate a random speed
    vrand = vmin + (vmax-vmin)*random.random()
    # Maxwell-Boltzmann distribution, normalized to v0
    #    (Probab.) = f_MB(v) dv
    P_MB = sqrt(2/pi)*vrand**2*exp(-vrand/2)
    # Generate a random probability
    Prand = random.random()
    if Prand > P_MB:
        return get_maxwellian_speed()
    else:
        return vrand # units of v0

def set_maxwellian_distrib(particles,kT):
    for i,particle in enumerate(particles):
        dir = randvec(1,samples=10)
        dir = dir/mag(dir)
        particle.vel = dir*get_maxwellian_speed()*v0



refbox = box(width=2*L,length=2*L,height=2*L,opacity=0.1)
ions = [sphere(cpos = randvec(L),color=color.cyan, mass=m_p,
            make_trail = True, trail_radius=0, trail_type="curve",interval=4,
            q = qe, radius=0.03, lpos=[]) for i in range(N//2)]
# electrons have a big mass to be visible
electrons = [sphere(cpos = randvec(L),color=color.red, mass=m_p/5,
            make_trail = True, trail_radius=0, trail_type="curve",interval=4,
            q = -qe, radius=0.03, lpos=[]) for i in range(N//2)]
particles = ions + electrons



set_maxwellian_distrib(particles,kT)



if Btype == "dipole":
    # dipolos puntuales, mostrar posición
    d_dip = 8 # media separación
    d1 = sphere(radius=0.05,color=color.yellow,
        pos=vec(0,0,d_dip),dir=vec(0,0,-1))
    d2 = sphere(radius=0.05,color=color.yellow,
        pos=vec(0,0,-d_dip),dir=vec(0,0,1))



dt = (L/v0)*0.0005
tf = 3e4*dt
t = 0
times = []

print("Calculando...")
l1 = label( pos=vec(0,0,1), text="Calculando...")

while t < tf:
    t = t+dt
    times.append(t)
    for particle in particles:

        if mag(particle.vel) > 50*v0:
            print("Warning: dt too large, ignoring particle.")
            continue

        if Btype == "uniform":
            B = 0.01*B0*vec(0,0,1) # campo uniforme

        if Btype == "helmholtz":
            # Campo variable en z (Bobina de Helmholtz)
            z_ref = 1 # z donde B = B0
            d_dip = 8 # separación de bobinas
            Bz = B0*(z_ref)**3/(particle.cpos.z + d_dip)**3 + B0*(z_ref)**3/(particle.cpos.z - d_dip)**3
            B = vec(0,0,Bz)

        if Btype == "dipole":
            # Dipolos puntuales
            zref = 1
            r_pd1 = particle.cpos - d1.pos
            r_pd2 = particle.cpos - d2.pos
            e_r_pd1 = r_pd1/mag(r_pd1)
            e_r_pd2 = r_pd2/mag(r_pd2)

            B1 = B0*zref**3*( 3*e_r_pd1*(d1.dir.dot(e_r_pd1)) - d1.dir )/mag(r_pd1)**3
            B2 = B0*zref**3*( 3*e_r_pd2*(d1.dir.dot(e_r_pd2)) - d2.dir )/mag(r_pd2)**3
            B = B1 + B2

        if Btype == "no":
            particle.acc = vec(0,0,0) # expansión libre

        if Btype != "no":
            # Fuerza de Lorentz, Runge-Kutta 2
            acc1 = particle.q*B.cross(particle.vel)/particle.mass
            dv1 = 1/2*dt*acc1
            particle.acc = particle.q*B.cross(particle.vel+dv1)/particle.mass

        particle.vel = particle.vel + particle.acc*dt
        particle.cpos = particle.cpos + particle.vel*dt
        particle.lpos.append(particle.cpos)

print("Visualizando...")
l1.visible=False

for i in range(0,len(times),20):
    rate(30)
    for particle in particles:
        particle.pos = particle.lpos[i]
