! Classical 1D wave
! -----------------
! The wave equation f_{tt} = c^2 f_{xx} with finite differences is
! ( f(x,t+dt) - 2f(x,t) + f(x,t-dt) )/dt = c^2 ( f(x+dx,t) - 2f(x,t) + f(x-dx,t) )/dx
!
! We solve for f(x,t+dt) and define alpha^2 = c^2*dt/dx and get
! f(x,t+dt) = 2*(1-alpha^2)*f(x,t) + alpha^2*(f(x+dx,t)+f(x-dx,t)) - f(x,t-dt)
!
! We will use the conditions:
! f(0,t) = 0 (boundary condition, fixed end)
! f(l,t) = 0 (boundary condition, fixed end)
! f(x,0) = initial shape of wave
! ( f(x,t+dt) - f(x,t) )/dt = initial velocity of wave
!
! General procedure:
!   We initialize a matrix f(i,j) that will contain all timesteps and x cells.
!   Then, we set the boundary conditions for all t, and the initial condition
!   for all x, accompanied by the condition for the first timestep, for which
!   we use the initial velocity of the wave. Finally, we use the wave equation,
!   with the obvious care to iterate only in the "inner cells" of the
!   computational grid, that is, excluding the boundaries (2nd order differentials)

program classical_1d_wave
  implicit none
  integer, parameter :: nx = 60
  integer, parameter :: nt = 500
  real, dimension(0:nx,0:nt) :: f
  real :: dt, dx, x, c, alpha2, l, tmax
  integer :: i, j
  real, parameter :: tension = 60 ! N
  real, parameter :: dens = 8.6e-4 ! kg/m

  integer :: file1,file2

  open(newunit=file1, file="out1.dat", form="unformatted", access="stream", status="replace")
  open(newunit=file2, file="out2.dat", form="formatted", status="replace")

  l = 0.3  ! m
  tmax = 5e-6 ! s
  dt = tmax/nt
  dx = l/nx
  c = sqrt(tension/dens)
  alpha2 = c**2*dt/dx

  ! Check for stability. alpha2 should be < 1.
  if (alpha2 > 1) print*, "Warning: alpha2 > 1. Unstable."

  ! Fixed ends
  do j = 1, nt
    f(0,j) = 0
    f(nx,j) = 0
  end do
  ! See comment (1).

  ! This should call the initial shape of the string
  f(0,0) = 0
  f(nx,0) = 0

  do i = 0, nx-1
    x = i*dx
    ! Initial shape of string
    f(i,0) = 0.02*exp(-(80*(x-0.1))**2)
    ! Condition with initial vel.: f(i,1) = f(i,0) + init_vel*dt
    f(i,1) = f(i,0) + 0.02*100*exp(-(80*(x-0.1))**2)*dt
  end do
  ! See comments (2), (3).

  do j = 1, nt-1
    do i = 1, nx-1
      ! Wave equation (recursion relation)
      f(i,j+1) = 2*(1-alpha2)*f(i,j) + alpha2*(f(i+1,j)+f(i-1,j)) - f(i,j-1)
      ! See comment (4)
    end do
  end do

  write(file1) f
  write(file2, "(2f10.6,2i8)") dt, dx, nt, nx


end program classical_1d_wave

! Other comments
! --------------
! (1) open ends are set with f_x(end,t) = 0, but warning:
! that kind of boundary also has reflexions (string tied to a column
! but free to move vertically)
!
! (2) the initial shape of the string could be a piecewise function
! by using ifs.
!
! (3) Interesting cases:
! - Traveling pulse: f(i,0) = 0.02*exp(-(80*(x-0.1))**2),
!     init_vel = 0.02*100*exp(-(80*(x-0.1))**2)
! - Standing wave f(i,0) = -2*x*(x-l), init_vel = 0
!
! (4) one could propagate the wave through different materials
! by simply defining a function alpha2(i) instead of a constant
