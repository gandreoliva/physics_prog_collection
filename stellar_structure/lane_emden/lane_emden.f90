!! Requires the module 'shtawa' for solving the ODE

!! Stellar structure
!! =================
!! Solves the Lane-Endem equation for the hydrostatic structure of a star, given a polytropic index.
!! The solution can be plotted using the Jupyter notebook.
!!
!! Theory
!! ------
!! For a static star (rho: density), it must be true that
!!    (i) mass (dM) is conserved within a shell (radius r, thickness dr):
!!          dM = 4*pi*r**2*rho*dr ;
!!    (ii) the sum of forces (gravity vs internal pressure P) is zero:
!!          dP/dr = -G*M_r*rho/r**2 .
!! Combining both equations and adimensionalizing using the transformation
!!      rho = rho_c*th**n,  r = a*xi,   a = sqrt((n+1)*K*rho_c**(1/n-1)/(4*pi*G))
!! one arrives at the Lane-Endem equation,
!!      1/xi**2 * d/dxi (xi**2*dth/dxi) = -th**n,
!! which has the following boundary conditions:
!!      rho(r=0) = rho_c  ==> th(xi=0) = 1
!!      and th'(xi=0) = 0.
!! The mass of the star is
!!      M = int_0^R 4*pi*r**2*rho*dr = 4*pi*a**3*rho_c * int_0^xi1 xi**2*th**n*dxi
!!
!! Implementation
!! --------------
!! The Lane-Endem equation can be re-written as a system of first-order ODEs:
!!      dth/dxi = -om/xi**2 (i)
!!      dom/dxi = th**n * xi**2  (ii).
!! The boundary conditions become initial conditions:
!!      th(0) = 1, th'(0) = 0.
!! Comparing (ii) to the enclosed mass in the star,
!!      M = 4*pi*a**3*rho_c * om ,
!! which means that om(0) = 0 (enclosed mass at the center is zero).
!! Warning: even though the initial conditions can be physically justified for xi = r = 0,
!! (i) is not defined there. That means that the r.h.s. must not be evaluated at xi=0.

!! Notation warning
!! ----------------
!! Compared to Kippenhahn+ 2012 "Stellar Structure and Evolution" Springer, Ch 19,
!! the following substitutions are required:
!!   xi --> z
!!   th  --> w

program lane_emden
    use iso_fortran_env, only: output_unit, wp => real64
    use shtawa, only: rk4
    implicit none

    real(wp), dimension(1:2) :: state
    real(wp) :: xi, th, om, n, dxi
    integer :: niter

    !! initial conditions
    th = 1
    om = 0
    xi = 0

    !! polytropic index
    n = 1.5
    state = [th, om]
    
    dxi = 0.001
    niter = 10000

    call rk4(rhs, state, xi, dxi*niter, niter,&
        &outfile=output_unit, stopping_cond=has_hit_surface, &
        &outfreq=niter/1000)


contains
    function rhs(xi,state)
        real(wp), intent(in) :: xi
		real(wp), dimension(:), intent(in) :: state
		real(wp), dimension(size(state)) :: rhs

        associate (th => state(1), om => state(2),&
            & dth_dxi => rhs(1), dom_dxi => rhs(2) )
            
            if (xi > 0d0) then
                dth_dxi = -om/xi**2
                dom_dxi = th**n * xi**2
            else !! initial cond. (rhs otherwise is undefined)
                dth_dxi = 0d0
                dom_dxi = 0d0
            end if

        end associate

    end function

    function has_hit_surface(xi,state)
		real(wp), intent(in) :: xi
		real(wp), dimension(:), intent(in) :: state
		logical :: has_hit_surface

        has_hit_surface = .false.

        associate (th => state(1) )
            
            !! if the adimensional density is almost zero or negative, stop
            if (state(1) < 1d-6 ) then
                has_hit_surface = .true.
            end if

        end associate

    end function

end program