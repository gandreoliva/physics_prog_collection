program double_pendulum
  implicit none
  integer :: file1
  real :: theta, phi, t, dt
  real :: thetad, phid
  real, parameter :: l = 1
  real, parameter :: m = 0.4
  real, parameter :: tmax = 15
  real, parameter :: g = 9.8
  real, parameter :: pi = 3.141592
  real, dimension(1:2) :: elerhs


  open(newunit=file1, file="out1.dat", form="unformatted", access="stream", status="replace")

  theta = pi/2+pi/3
  phi = pi/4
  thetad = 0
  phid = 0

  t = 0
  dt = 0.01

  do while (t<tmax)
    include "out_analytical_mech_dpend_ele.f90"

    thetad = thetad + elerhs(1)*dt
    phid = phid + elerhs(2)*dt
    theta = theta + thetad*dt
    phi = phi + phid*dt

    write(file1) t, theta, phi
    t = t + dt
  end do


end program double_pendulum
