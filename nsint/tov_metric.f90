program tov
  ! Integrates the TOV equations and outputs metric potentials
  implicit none
  real, parameter :: pi = 3.141592
  real, parameter :: G = 6.6743e-11, c=2.99792e8, km=1e3, msun=1.988e30

  integer, parameter :: n_r = 100
   ! Number of cells in 'r' - 1
  real :: r_end, dr, k_pol, gamma_pol
  integer :: i,ii
  real, dimension(0:n_r) :: p,r,rho,m

  logical :: reached_surface ! stops integration once the surface is reached
  real, dimension(0:n_r) :: pphi ! metric potential Phi
  integer :: isurf ! Index where the surface is located

  ! Boundary conditions and parameters
  r_end = 40 ! km
  read(*,*) rho(0) ! g/cm**3
  rho(0) = rho(0)*1e3 * G/c**2*km**2 ! km**-2, see notes 1,3
  read(*,*) k_pol ! see note 2
  read(*,*) gamma_pol
  reached_surface = .false.

  ! Grid generation
  dr = r_end/n_r
  r(:) = dr/10
  do i = 1,n_r
    r(i) = r(i-1) + dr
  end do

  p(0) = k_pol*rho(0)**gamma_pol
  m(0) = 0


  print "(4a15)", "#r (km)    ", "Phi      ", "g_tt       ",  "g_rr       "

  do i = 1,n_r
    ! Explicit integration: mass
    !  (m(i) - m(i-1))/dr = 4*pi*r(i-1)**2*rho(i-1)
    m(i) = m(i-1) + 4*pi*r(i-1)**2*rho(i-1) * dr

    ! Explicit integration: TOV equation
    ! (p(i) - p(i-1))/dr = rhs_tov(r(i-1),p(i-1),rho(i-1),m(i-1))
    p(i) = p(i-1) + rhs_tov(r(i-1),p(i-1),rho(i-1),m(i-1)) * dr
    rho(i) = (p(i)/k_pol)**(1/gamma_pol)

    if (p(i) < 0) then
      reached_surface = .true.
      isurf = i
    end if

    if (reached_surface .eqv. .true.) then
      p(i) = 0
      m(i) = m(i-1)
      rho(i) = 0
    end if

  end do


  ! See note 4
  pphi(n_r) = 0.5*log(1 - 2*m(n_r)/r(n_r))

  do ii = 1,n_r
    i = n_r - ii
    ! Explicit integration: metric potential Phi
    ! backwards in r
    ! (pphi(i+1) - pphi(i))/dr = (4*p*pi*r**3+m)/(r**2-2*m*r)
    pphi(i) = pphi(i+1) - (4*p(i)*pi*r(i)**3+m(i))/(r(i)**2-2*m(i)*r(i))
  end do


  ! Print results
  do i = 0,n_r
    print "(4g15.6)", r(i), pphi(i), exp(2*pphi(i)), (1-2*m(i)/r(i))**(-1)
  end do


contains
  function rhs_tov(r,p,rho,m)
    real :: r,p,rho,m,rhs_tov
    rhs_tov = -((4*p*pi*r**3+m)*(rho+p))/(r*(r-2*m))
  end function

end program


! NOTE:

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

! 4. Backwards integration of the metric potential Phi. The boundary condition
!   is that at large distances outside of the neutron star, the potential
!   must tend to Schwarzschild, since if we integrate anallytically the
!   equation for Phi with p=0, m(r) = M, we get
!     Phi(r>R) = 1/2*ln(1-2*M/r)
!  ++ Warning:the numerical method is very bad (Euler: balance the number of
!   steps [exponential error!] vs the stepsize). For a serious application,
!   another numerical method must be chosen.
