! Static star
! ===========
! Qualitatively solves the stellar structure equations with the shooting method
! and manual adjustements to satisfy the central boundary conditions.
! Adapted from Carroll & Ostlie (2014). For teaching purposes only.

module constants
    use iso_fortran_env, only: wp => real64
    implicit none
    real(wp), parameter :: pi = 3.141592653589793d0
    real(wp), parameter :: G = 6.67430d-8 ! dyn cm^2 g^-2
    real(wp), parameter :: c_light = 2.99792458e+10 ! cm/s
    real(wp), parameter :: gamma_ad = 5/3.d0 ! monoatomic gas
    real(wp), parameter :: sigma_B = 5.670374419d-5 ! erg cm^-2 s^-1 K^-4
    real(wp), parameter :: k_B = 1.38d-16 ! erg K^-1
    real(wp), parameter :: a_rad = 4*sigma_B/c_light
    real(wp), parameter :: solar_luminosity = 3.828d33 ! erg/s
    real(wp), parameter :: solar_mass = 1.988416d33 ! g
    real(wp), parameter :: solar_metallicity = 0.0134
    real(wp), parameter :: solar_H_frac = 0.7381
    real(wp), parameter :: m_proton = 1.672621925d-24 ! g
    save
end module

module microphysics
    ! Computation of microphysics operating on scalar quantities (not arrays)
    use iso_fortran_env, only: wp => real64
    use constants
    implicit none
    real(wp), public :: molec_weight
    real(wp), private :: X, Y, Z
    save
contains
    subroutine initalize_microphysics(H_frac, He_frac, metallicity)
        real(wp), intent(in) :: H_frac, He_frac, metallicity
        X = H_frac
        Y = He_frac
        Z = metallicity
        molec_weight = 1/(2*X + 0.75*Y + 0.5*Z)
    end subroutine

    function get_density_from_eos(pressure,temperat) result(density)
        real(wp) :: density, pressure, temperat
        density = molec_weight*m_proton/k_B * pressure/temperat
    end function

    function get_opacity(density,temperat) result(opacity)
        real(wp) :: density, temperat
        real(wp) :: opacity
        real(wp) :: k_boundfree, k_freefree, k_elecscatt, gaunt_to_guillotine, gff

        ! bound-free opacity (~> Kramer's law) [photoionization]
        gaunt_to_guillotine = 2.82d0*(density*(1+X))**0.2d0
        k_boundfree = 4.34d25/gaunt_to_guillotine * Z*(1+X) * density/temperat**3.5d0
        ! free-free opacity (~> Kramer's law) [Bremsstrahlung]
        gff = 1d0 ! Gaunt factor
        k_freefree = 3.68d22*gff*(1-Z)*(1+X)*density/temperat**3.5d0
        ! electron scattering [Thomson scattering]
        k_elecscatt = 0.2*(1+X)
        opacity = k_boundfree + k_freefree + k_elecscatt
    end function

    function get_nuclear_energy_generated(density,temperat) result(nuclear_energy_generation)
        real(wp) :: density, temperat
        real(wp) :: T6
        real(wp) :: fx, f_pp, psi_pp, g_pp, energen_pp
        real(wp) :: X_CNO, g_CNO, energen_CNO
        real(wp) :: nuclear_energy_generation ! erg g^-1 s^-1

        T6 = temperat*1d-6

        ! pp-chain
        fx = 0.133*X*sqrt(density*(3+X))/T6**1.5d0
        f_pp = 1 + fx*X
        psi_pp = 1 + 1.412d8 * (1d0/X - 1) * exp(-49.98*T6**(-1/3d0))
        g_pp = 1 + 0.0123d0*T6**(1/3d0) + 0.0109d0*T6**(2/3d0) + 0.0009*T6
        energen_pp = 2.38d6*density*X**2*f_pp*psi_pp*g_pp*T6**(-2/3d0)*exp(-33.8d0*T6**(-1/3d0))

        ! CNO cycle
        X_CNO = Z/2
        g_CNO = 1 + 0.0027d0*T6**(-1/3d0) - 0.00778d0*T6**(2/3d0) - 0.00149d0*T6
        energen_CNO = 8.67d27*g_CNO*density*X*X_CNO*T6**(-2/3d0)*exp(-152.28d0*T6**(-1/3d0))

        nuclear_energy_generation = energen_pp + energen_CNO
    end function
end module


program static_star
    use iso_fortran_env, only: wp => real64
    use constants
    use microphysics
    implicit none
    integer, parameter :: imax = 100
    real(wp), dimension(1:imax) :: r, prs, mass, lum
    real(wp), dimension(1:imax) :: k_opac, rho, temp, enuc
    real(wp) :: stellar_mass, stellar_luminosity, effective_temperature
    real(wp) :: stellar_radius, dr, dlT_dlP, dP_dr
    real(wp) :: H_frac, He_frac, metallicity
    integer :: i
    character(len=1) :: transport_type ! = {r: radiative, c: convective}

    ! input quantities: Sun
    stellar_mass = solar_mass
    stellar_luminosity = solar_luminosity
    effective_temperature = 5772 ! K
    H_frac = solar_H_frac
    metallicity = solar_metallicity
    He_frac = 1 - H_frac - metallicity

    ! input quantities: Meissa A
    ! stellar_mass = 34*solar_mass
    ! stellar_luminosity = 200000*solar_luminosity
    ! effective_temperature = 36000 ! K
    ! H_frac = solar_H_frac
    ! metallicity = solar_metallicity
    ! He_frac = 1 - H_frac - metallicity


    ! external boundary conditions
    prs(imax) = 0 ! pressure
    mass(imax) = stellar_mass ! enclosed mass
    lum(imax) = stellar_luminosity ! shell luminosity
    temp(imax) = effective_temperature ! temperature
    rho(imax) = 0 ! density
    stellar_radius = sqrt(stellar_luminosity/(4*pi*sigma_B*effective_temperature**4))
    r(imax) = stellar_radius
    dr = -stellar_radius/imax
    transport_type = 'r'

    call initalize_microphysics(H_frac, He_frac, metallicity)

    print*, "r(i)/Rsun     mass(i)/Msun    temp(i)       rho(i)       prs(i)      lum(i)/Lsun    conv|rad"


    i = imax


    do while (i > imax-20)
        call structure_near_surface(i)
        print"(6es14.5,' ',a)", r(i)/stellar_radius, mass(i)/solar_mass, temp(i), rho(i), prs(i), lum(i)/solar_luminosity, transport_type
        i = i - 1
    end do

    do while (i > 0 )
        r(i-1) = r(i) + dr

        ! microphysics
        rho(i) = get_density_from_eos(prs(i),temp(i))
        k_opac(i) = get_opacity(rho(i),temp(i))
        enuc(i) = get_nuclear_energy_generated(rho(i),temp(i))

        ! *** structure equations ***
        !   (1) hydrostatic equilibrium (cons. of momentum)
        prs(i-1) = prs(i) - G*mass(i)*rho(i)*dr/r(i)**2

        !   (2) conservation of mass
        mass(i-1) = mass(i) + 4*pi*r(i)**2 * rho(i) * dr

        !   (3) conservation of energy
        lum(i-1) = lum(i) + 4*pi*r(i)**2 * rho(i) * enuc(i) * dr

        ! Energy transport
        if (transport_type == 'c') then
            !   (4.1) energy transfer in convective zones
            dP_dr = ( prs(i-1) - prs(i) ) / dr
            temp(i-1) = temp(i) + (1-1d0/gamma_ad) * temp(i)/prs(i) * dP_dr * dr
        else
            !   (4.2) energy transfer in radiative zones
            temp(i-1) = temp(i) - (3/(16*pi*a_rad*c_light)) * &
                &   k_opac(i) * rho(i) * lum(i) / (r(i)**2 * temp(i)**3) * dr
        end if

        dlT_dlP = log(temp(i-1)/temp(i))/log(prs(i-1)/prs(i))
        ! test to decide whether the next shell is radiative or convective
        if (dlT_dlP > (gamma_ad - 1)/gamma_ad) then
            transport_type = 'c' ! covective
        else
            transport_type = 'r' ! radiative
        end if

        print"(6es14.5,' ',a)", r(i)/stellar_radius, mass(i)/solar_mass, temp(i), rho(i), prs(i), lum(i)/solar_luminosity, transport_type
        i = i - 1
    end do

contains
    subroutine structure_near_surface(i)
        ! Near the surface, we need to expand the structure equations (otherwise, e.g. rho = 0 because prs = 0)
        integer, intent(in) :: i
        real(wp) :: gaunt_to_guillotine, a_bf, a_ff, a_tot, ad_pt_const

        ! Calculation of opacity factors using Kramer's law form
        if (i == imax) then
            gaunt_to_guillotine = 0.01d0 ! boundary value
            ad_pt_const = 0.3d0 ! adiabatic pressure constant, arbitrary at the boundary
        else
            gaunt_to_guillotine = 2.82d0*(rho(i)*(1+H_frac))**0.2d0
            ad_pt_const = prs(i)/temp(i)**(gamma_ad/(gamma_ad-1))
        end if
        ! coefficient for bound-free opacity
        a_bf = 4.34d25/gaunt_to_guillotine * metallicity*(1+H_frac)
        ! coefficient for free-free opacity
        a_ff = 3.68d22*1*(1-metallicity)*(1+H_frac)
        a_tot = a_bf + a_ff

        r(i-1) = r(i) + dr
        mass(i-1) = mass(i)
        lum(i-1) = lum(i)

        if (transport_type == 'r') then
            ! We assume near the surface lum(i) = L, mass(i) = M. From the structure eqns (with rad),
            ! one builds dP/dT = (dP/dr)/(dT/dr). Using only bound-free and free-free opacity (Kramer's law)
            ! and the equation of state, one eliminates the dep. on rho and can separate variables.
            ! Then, one can integrate and find P(T). Finally, with the structure eqn on dT/dr one
            ! can integrate to get T(r). Similar procedure for convection.
            temp(i-1) = (molec_weight*m_proton/k_B) * G*mass(i-1)/4.25 &
                &       * (1/r(i-1) - 1/stellar_radius)
            prs(i-1) = sqrt((1/4.25d0)*(16*pi*a_rad*c_light/3)*(k_B/(a_tot*molec_weight*m_proton))) &
                &       * sqrt(G*mass(i-1)/lum(i-1)) * temp(i-1)**4.25d0
        else
            temp(i-1) = (gamma_ad - 1)/gamma_ad * G*mass(i-1) * (molec_weight*m_proton/k_B) &
                &       * (1/r(i-1) - 1/stellar_radius)
            prs(i-1) = ad_pt_const * temp(i-1)**(gamma_ad/(gamma_ad-1))
        end if

        rho(i-1) = get_density_from_eos(prs(i-1),temp(i-1))
        k_opac(i-1) = get_opacity(rho(i-1),temp(i-1))
        enuc(i-1) = get_nuclear_energy_generated(rho(i-1),temp(i-1))

        ! determine whether the next layer will be convective or radiative
        if (i < imax) then
            dlT_dlP = log(temp(i-1)/temp(i))/log(prs(i-1)/prs(i))
            if (dlT_dlP > (gamma_ad - 1)/gamma_ad) then
                transport_type = 'c' ! covective
            else
                transport_type = 'r' ! radiative
            end if
        end if

    end subroutine

end program