program tov
  ! Integrates TOV equations and outputs mass, pressure and density
  ! gfortran tov_fields.f90
  implicit none
  real, parameter :: pi = 3.141592
  real, parameter :: G = 6.6743e-11, c=2.99792e8, km=1e3, msun=1.988e30

  integer, parameter :: n_r = 50 ! Number of cells in 'r' - 1
  real :: r_end, dr, k_pol, gamma_pol
  integer :: i
  real, dimension(0:n_r) :: p,r,rho,m
  real :: phi, lambda

  ! Boundary conditions and parameters
  r_end = 30 ! km
  rho(0) = 3e15*1e3 * G/c**2*km**2 ! km**-2, see notes 1,3
  k_pol = 7.25 ! see note 2
  gamma_pol = 5/3.

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
!  * For a relativistic gas, P = 0.8*hbar*c/(mn)**(4/3) * rho**(4/3)
!     => K[SI] = 0.8*hbar*c/(mn)**(4/3) = 1.23e10
!     Dimensions: D[hbar]*D[c]/D[mn] = M**(-1/3)*L**3*T**(-2)
!     => K[g] = c**(-2) * (G/c**2)**(-1/3) * K[SI] = 151 m**(2/3)
!  * In [arXiv:1406.3775v1]: K[g] = 100 km**2 for Gamma = 2
!  * For a non relativistic gas, P = 1.9*hbar**2*mn**(-8/3)*rho**(5/3)
!     => K[SI] = 5341
!     Dimensions: D[P]/D[rho]**(5/3) = L**4*M**(-2/3)*T**-2
!     => D[K,g] = L**(4/3), and K[g] = c**-2*(G/c**2)**(-2/3) * K[SI]
!     => K[g] = 72 470 m**(4/3) = 7.25 km**(4/3)

! 3. We are assuming relativistic energy density = rho; a rest mass distiction
!   must be taken into account.
