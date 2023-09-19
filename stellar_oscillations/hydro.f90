!!! Stellar oscillations - 1D hydro+gravity code
!!! compilation: gfortran -o hydro hydro.f90
!!! run: ./hydro
!!!
!!! Running order
!!!     Lane-Emden eqn. --> (( hydro code )) --> plotting scripts

module units
  !! Handles unit conversions and sets constants
  !! --------------------------------------------------------
  !! Example of use in `program`:
  !!    use units
  !!    n_gas = 1.5d0 ! <--- Set n_gas
  !!    unit_density = 162d3 ! <--- Set central density [ kg/m^3 ]
  !!    call set_code_constants
  !!    call set_units
  !! After that, any quantity (e.g. density, pressure, etc.) can be computed as
  !!  (quantity in SI units) = (unit_quantity)*(adimensional quantity)
  use, intrinsic :: iso_fortran_env, dp => real64
  implicit none
  real(dp), parameter :: grav_constant = 6.67430d-11 ! m^3*kg^-1*s^-2
  real(dp), parameter :: solar_mass = 1.98847d30 ! kg
  real(dp), parameter :: solar_radius = 6.95700d8 ! m
  real(dp), parameter :: pi = 3.141592653589793
  real(dp) :: unit_time, unit_length, unit_velocity
  real(dp) :: unit_density, unit_pressure, k_gas, k_gas_si, n_gas
  save
contains
    subroutine set_code_constants!(n_gas)
        !! Sets the constants in code units
        !! ------------------------------
        !! Requires the previous setup of n_gas
        k_gas = 4*pi/(n_gas + 1)
    end subroutine

    subroutine set_units!(central_density, n_gas)
      !! Computes the code units given the central density of the star
      !! and the polytropic index n
      !! ------------------------------
      !! Requires the previous setup of:
      !!    unit_density in kg/m^3, n_gas, k_gas_si
      unit_time = 1/sqrt(grav_constant*unit_density)
      unit_length = sqrt((n_gas+1)*k_gas_si*unit_density**(1/n_gas-1)/(4*pi*grav_constant))
      unit_velocity = unit_length/unit_time
      unit_pressure = k_gas_si*unit_density**(1+1/n_gas)
    end subroutine
end module

module hydro_solver
  !! 1D hydrodynamics solver with self-gravity
  !! --------------------------------------------------------
  !! Simple finite differences with donor cell upwind scheme and operator splitting. See
  !! Bodenheimer, Laughlin, Rozyczka and Yorke (2007) "Numerical Methods in Astrophysics",
  !! Taylor and Francis. Sections 6.3.2 and 6.5.
  !! --------------------------------------------------------
  !! /!\ All variables in this module are in code units only!
  use, intrinsic :: iso_fortran_env, dp => real64
  use units, only: n_gas, k_gas
  implicit none
  !! Parameters
  integer, parameter :: nr = 100 ! number of cells
  real(dp), parameter, private :: pi = 3.141592653589793
  !! primitive variables of the problem
  real(dp), dimension(0:nr+1) :: rh ! (adimensional) density (def. at centers)
  real(dp), dimension(0:nr+1) :: p ! (adimensional) pressure (def. at centers)
  real(dp), dimension(0:nr+2) :: u ! (adimensional) velocity (def. at walls)
  real(dp), dimension(0:nr+2) :: ra ! cell walls
  real(dp), dimension(0:nr+1) :: rb ! cell centers
  !! variables of the numerical solution
  real(dp) :: cfl ! Courant-Friedrich-Levy number
  real(dp) :: t ! time
  real(dp) :: dt ! difference of time
  real(dp) :: dr ! difference of radial position: regular grid
  real(dp), dimension(1:nr) :: dvolb ! volumes of all (centered/space) cells
  real(dp), dimension(2:nr) :: dvola ! volumes of all (wall/momentum) cells
  real(dp), dimension(1:nr+1) :: surfa ! surface area of all walls
  real(dp), dimension(1:nr) :: surfb ! surface area of all centers
  real(dp), dimension(1:nr+1), private :: w ! momentum
  real(dp), dimension(1:nr+1), private :: avg_rh ! Average density
  real(dp), dimension(1:nr), private :: avg_u ! Average velocity
  real(dp), dimension(1:nr), private :: cs ! sound speed
  !! fluxes
  real(dp), dimension(1:nr+1), private :: m_flux ! mass flux
  real(dp), dimension(1:nr), private :: w_flux ! momentum flux
  integer :: status
  save
contains
  subroutine set_grid!(dr)
    !! Builds the grids, and computes the areas and volumes of the cells
    !! ------------------------------
    !! Requires the previous setup of dr
    integer :: i

    !! creation of walls
    do i=0,nr+2
      ra(i) = (i-1)*dr
      !! expl. of (i-1): the left ghost cell is in a negative position
    end do
    !! creation of cell centers
    do i=0,nr+1
      rb(i) = ((i-1)+1/2d0)*dr
    end do
    !! calculation of volumes of the cells
    do i=1,nr
      dvolb(i) = 4*pi/3*( ra(i+1)**3 - ra(i)**3 )
    end do
    do i=2,nr
      dvola(i) = 4*pi/3*( rb(i)**3 - rb(i-1)**3 )
    end do
    !! calculation of the areas of the cells
    surfa(1:nr+1) = 4*pi*ra(1:nr+1)**2
    surfb(1:nr) = 4*pi*rb(1:nr)**2
  end subroutine

  subroutine apply_boundary_cond
    !! no velocity at boundaries => no outflow
    u(1) = 0
    u(nr+1) = 0
    !! ghost cells: velocity reflected
    u(0) = -u(2)
    u(nr+2) = -u(nr)
    !! density: reflective boundaries: d rh/ dr = 0
    rh(0) = rh(1)
    rh(nr+1) = rh(nr)
  end subroutine

  subroutine apply_restrictions
    integer :: i

    do i = 0, nr+1
      if ( rh(i) <= 0 ) then
        rh(i) = 1d-6 !! density floor
        print*, " /!\ Density floor reached"
        status = -1
      end if
      if ( rh(i) > huge(0d0) ) stop "*** (X) rh became infinity!"
      if ( isnan(rh(i)) ) stop "*** (X) rh became NaN!"
    end do
  end subroutine

  subroutine check_cfl
    real(dp) :: recom_dt
    !! Check CFL criterion
    cfl = 0.75
    cs = sqrt( (1+1/n_gas)*p(1:nr)/rh(1:nr) )
    recom_dt = cfl*minval(dr/(cs+abs(u(1:nr))))

    if ( dt > recom_dt ) then
      print*, " /!\ CFL criterion not satisfied! You should choose a lower dt!"
      print*, " t = ", t, " dt =", dt, " recom_dt = ", recom_dt
      status = -1
    end if
  end subroutine

  !! Evolution steps
  subroutine advect_density
    integer :: i
    real(dp) :: upwind_rh
    !! a. Density advection
    !! a.1. compute the mass flux
    do i = 1,nr+1
      !! decide what value of density to use (determine upwind density)
      if ( u(i) > 0 ) then
        upwind_rh = rh(i-1)
      else
        upwind_rh = rh(i)
      end if
      !! set the mass flux
      m_flux(i) = upwind_rh*u(i)
    end do
    !! a.2. evolve density (continuity equation)
    !!    ---> we obtain rh(t + dt , r)
    do i = 1,nr
      rh(i) = rh(i) - dt/dvolb(i) * ( surfa(i+1)*m_flux(i+1) - surfa(i)*m_flux(i) )
    end do
  end subroutine

  subroutine advect_momentum
    integer :: i
    real(dp) :: upwind_w

    !! b. Momentum advection
    !! b.1. compute the current momentum
    do i = 1,nr+1
      avg_rh(i) = 1/2d0*( rh(i) + rh(i-1) )
      w(i) = u(i)* avg_rh(i)
    end do

    !! b.2. compute the momentum flux
    do i = 1,nr
      !! decide what value of momentum to use (determine upwind momentum)
      avg_u(i) = 1/2d0*( u(i) + u(i+1) )
      if ( avg_u(i) > 0 ) then
        upwind_w = w(i)
      else
        upwind_w = w(i+1)
      end if
      !! set the momentum flux
      w_flux(i) = upwind_w*avg_u(i)
    end do

    !! b.3. evolve momentum (momentum equation, only advection terms)
    !!     ---> we obtain w( t + first_step, r )
    !! N.b.: u(i) is defined at the right hand side of the cell, which is
    !! a different convention than Bodenheimer et al 2007 Sect 6.5.
    do i = 2,nr
      w(i) = w(i) - dt/dvola(i) * ( surfb(i)*w_flux(i) - surfb(i-1)*w_flux(i-1) )
      u(i) = w(i)/avg_rh(i)
    end do
  end subroutine

  subroutine sources_step
    integer :: i

    !! compute pressure
    p = k_gas*rh**(1 + 1/n_gas)

    !! evolve velocity (momentum equation, force terms)
    !! ---> we obtain u(t + dt)
    do i = 2, nr
      u(i) = u(i) + dt*( f_ext(i) - 1/avg_rh(i)*(p(i) - p(i-1))/(rb(i) - rb(i-1)) )
    end do
  end subroutine

!! External specific forces implemented

  function f_ext(i)
    !! Total external force. Here one can turn gravity on and off, for example,
    !! or define an external constant force.
    integer :: i
    real(dp) :: f_ext
    ! f_ext = f_grav(i)
    f_ext = 0d0
  end function

  function f_grav(i)
    !! Computes the gravitational force for cell i
    integer :: i,j
    real(dp) :: f_grav, m
    !! calculation of the enclosed mass
    m = 0d0
    do j = 1,i-1
      m = m + rh(j)*dvolb(j)
    end do
    !! gravitational specific force (code stellar units)
    f_grav = -m/(ra(i)**2)
  end function

end module


program stellar_oscillations
  !! Main program (the execution starts here)
  !! --------------------------------------
  use, intrinsic :: iso_fortran_env, dp => real64
  use hydro_solver
  ! use stellar_units
  implicit none
  integer :: nt, i, ntmax, output_frequency
  integer :: rh_outfile, u_outfile, t_outfile, ra_outfile, rb_outfile, inputfile
  character(*), parameter :: out_dir = "data/"
  !! Controls what's printed on the screen: 0: step info, 1: detailed, 2: minimal
  integer, parameter :: verbose = 2
  real(dp), dimension(0:nr+1) :: init_p, init_cs
  real(dp) :: rmin, rmax, rstar

  status = 0 ! Tells us if there are warnings or errors

  print*, "--------------------------------------------------------------------------"
  print*, " Stellar oscillations   *  . ·           > 1D hydro+gravity code"
  print*, "=========================================================================="

  !! Parameters
  !! ([un]comment to select: reading parameters from file or prescription)

  !! Parameters: FROM FILE
  open(newunit=inputfile,  file="equil_model.dat", action="read")
  read(inputfile,*) n_gas
  read(inputfile,*) rstar

  !! Parameters: ANALYTICAL
  ! rstar = pi
  ! n_gas = 1d0

  !! (end of selection by [un]commenting)

  !! Setting of units and grid
  rmax = 1.5*rstar ! maximum r position in the grid [code units]
  rmin = 0d0 ! minimum r position in the grid [code units]
  !! The number of cells in r, `nr`, is adjusted in the module `hydro_solver`
  dr = (rmax - rmin)/nr

  call set_code_constants
  call set_grid


  !! Initial conditions
  !! ([un]comment to select: reading initial condition from file or prescription)

  !! initial density: ANALYTICAL
  ! rh = 0.001d0
  ! rh(0) = 1d0
  ! do i = 1, nr
  !   if (rb(i) < 3.14 ) then
  !     rh(i) = sin(rb(i))/rb(i)
  !   end if
  ! end do

  !! initial density: FROM FILE
  do i = 0, nr+1
    read(inputfile,*) rh(i)
  end do

  !! (end of selection by [un]commenting)


  !! initial velocity
  u = 0.d0

  !! velocity perturbations
  init_p = k_gas*rh**(1 + 1/n_gas)
  init_cs = sqrt((1 +1/n_gas)*init_p/rh)
  do i = 0, nr+1
    if ( rb(i) <= rstar) then
      u(i) = 0.2d0 * init_cs(i) * sin( pi*rb(i)/rstar )
    end if
  end do

  !! time control
  dt = 0.001
  ntmax = 3000
  output_frequency = 100

  !! data files
  open(newunit=rh_outfile,  file=out_dir//"rh.dat", status="replace", action="write")
  open(newunit=u_outfile,   file=out_dir//"u.dat", status="replace", action="write")
  open(newunit=t_outfile,   file=out_dir//"t.dat", status="replace", action="write")
  open(newunit=ra_outfile,file=out_dir//"ra.dat", status="replace", action="write")
  open(newunit=rb_outfile,file=out_dir//"rb.dat", status="replace", action="write")

  !! output grid
  write(ra_outfile,*) ra
  write(rb_outfile,*) rb

  do nt = 0, ntmax
    if (verbose == 0) write(*,"(a,i6,a,e12.2)") "step", nt, " | t =",t
    if (verbose == 1) print*, ":::::::: step ", nt, "::::::::", " t =",t

    if (verbose == 1) print*, "advecting density..."
    call advect_density

    if (verbose == 1) print*, "advecting momentum..."
    call advect_momentum

    if (verbose == 1) print*, "forces step..."
    call sources_step

    if (verbose == 1) print*, "applying boundary conditions..."
    call apply_boundary_cond

    if (verbose == 1) print*, "applying restrictions..."
    call apply_restrictions

    if (verbose == 1) print*, "checking time step..."
    t = t + dt
    call check_cfl

    if (verbose == 2) then
      write(*,'(a,a,f5.1,a,$)') char(13), "evolving in time... ", (nt)/(1.0*abs(0)+abs(ntmax))*100.0, "%"
    end if

    !! output to files
    if (mod(nt,output_frequency) == 0) then
      write(rh_outfile,*) rh
      write(u_outfile,*) u
      write(t_outfile,*) t
    end if
  end do

  close(rh_outfile)
  close(u_outfile)
  close(t_outfile)
  close(ra_outfile)
  close(rb_outfile)

  if (status /= 0) print*, "  /!\ The program finished with warnings."
  if (verbose == 2) print*, " "

end program
