program tov
  ! Integrates TOV equations and outputs mass and the pressure and density structure
  implicit none
  real, parameter :: pi = 3.141592
  real, parameter :: G = 6.6743e-11, c=2.99792e8, km=1e3, msun=1.988e30

  integer, parameter :: n_r = 100 ! Number of cells in 'r' - 1
  real :: r_end, dr, k_pol, gamma_pol
  integer :: i
  real, dimension(0:n_r) :: p,r,rho,m

  ! Boundary conditions and parameters
  r_end = 40 ! km
  read(*,*) rho(0) ! g/cm**3
  rho(0) = rho(0)*1e3 * G/c**2*km**2 ! km**-2, see notes 1,3
  read(*,*) k_pol ! see note 2
  read(*,*) gamma_pol


  ! Grid generation
  dr = r_end/n_r
  r(:) = dr/10
  do i = 1,n_r
    r(i) = r(i-1) + dr
  end do

  p(0) = k_pol*rho(0)**gamma_pol
  m(0) = 0

  print "(4a18)", "#r (km)   ", "m (Msun)   ", "p (Pa)   ", "rho (kg/m^3)   "

  do i = 1,n_r
    ! Explicit integration: mass
    !  (m(i) - m(i-1))/dr = 4*pi*r(i-1)**2*rho(i-1)
    m(i) = m(i-1) + 4*pi*r(i-1)**2*rho(i-1) * dr

    ! Explicit integration: TOV equation
    ! (p(i) - p(i-1))/dr = rhs_tov(r(i-1),p(i-1),rho(i-1),m(i-1))
    p(i) = p(i-1) + rhs_tov(r(i-1),p(i-1),rho(i-1),m(i-1)) * dr
    rho(i) = (p(i)/k_pol)**(1/gamma_pol)
    print*, r(i), m(i)*km*c**2/(G*msun), p(i)*c**4/(G*km**2), rho(i)*c**2/(G*km**2)

  end do

contains
  function rhs_tov(r,p,rho,m)
    real :: r,p,rho,m,rhs_tov
    rhs_tov = -((4*p*pi*r**3+m)*(rho+p))/(r*(r-2*m))
  end function

end program


! Notes:

! 1. Conversion SI <=> geometrized
!  according to [Wald 1984 General Relativity, UChicago: p470], for x with
!  D[x,SI] = L**n * T**m * M**p
!  => x[g] = c**m * (G/c**2)**p * x[SI], and D[x,g] = L**(n+m+p)

! 2. Polytropic EOS
!  See parameters/ directory
!  The code can be run with those files as input from the command line.
!  * For a non relativistic neutron gas, P = 1.9*hbar**2*mn**(-8/3)*rho**(5/3)
!     ---> gamma_pol = 5/3, k_pol = 7.349
!     (bad values but first theoretical approach; order-of-magnitude correct only)
!  * BU0 model for a polytrope of moderate stiffness: K[g] ~ 100 km**2 for Gamma = 2

! 3. We are assuming relativistic energy density = rho; a rest mass distiction
!   must be taken into account.
