Polytropic EOS
==============

Structure of the files
----------------------
central density, in g/cm**3 (real)
polytrope k (real)
polytrope Gamma (real)

Note: in the literature about polytropes, Gamma = 1 + 1/n

List of polytropes
------------------

* nrel_neutrongas:
    For a non relativistic neutron gas, P = 1.9*hbar**2*mn**(-8/3)*rho**(5/3)
    => K[SI] = 5377
    Dimensions: D[P]/D[rho]**(5/3) = L**4*M**(-2/3)*T**-2
    => D[K,g] = L**(4/3), and K[g] = c**-2*(G/c**2)**(-2/3) * K[SI]
    => K[g] = 72 470 m**(4/3) = 7.349 km**(4/3)
    ---> gamma_pol = 5/3, k_pol = 7.349
    (bad approach: order-of-magnitude correct only, but small and light neutron star)

 * rel_neutrongas:
    For an ultra-relativistic neutron gas, P = 0.8*hbar*c/(mn)**(4/3) * rho**(4/3)
    => K[SI] = 0.8*hbar*c/(mn)**(4/3) = 2e9
    Dimensions: D[hbar]*D[c]/D[mn] = M**(-1/3)*L**3*T**(-2)
    => K[g] = c**(-2) * (G/c**2)**(-1/3) * K[SI] = 1.515 km**(2/3)
    ---> gamma_pol = 4/3, k_pol = 1.515
    (bad approach: way too large neutron star! It's here just for completeness)

 * BU0: empirical
    K[g] ~ 100 km**2 for Gamma = 2
    This corresponds to a polytropic EOS that yields average mass and radius neutron stars;
    EOS of moderate stiffness. See, e.g., Stergioulas et al 2004 MNRAS 352 1089.