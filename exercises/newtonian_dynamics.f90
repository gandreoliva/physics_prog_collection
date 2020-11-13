! Newtonian dynamics

program newtonian_dynamics
  integer :: file1
  real, dimension(1:3) :: pos, vel, acc, force
  real :: mass, t, dt
  real, parameter :: tmax = 5

  open(newunit=file1, file="out1.dat", form="unformatted", access="stream", status="replace")

  pos = [1,0,0] ! m
  vel = [0,0,0] ! m/s

  dt = 0.01 ! s
  t = 0
  mass = 1.3 ! kg

  do while (t<tmax)
    ! Harmonic oscillator
    force = -50*pos - 0.5*vel
    acc = force/mass

    vel = vel + acc*dt
    pos = pos + vel*dt

    write(file1) t, pos, vel

    t = t + dt
  end do

end program newtonian_dynamics
