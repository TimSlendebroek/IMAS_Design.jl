function alpha_t(dd::IMAS.dd) # from https://doi.org/10.1088/1741-4326/ad89da

    R0 = dd.equilibrium.vacuum_toroidal_field.r0
    B_T = dd.equilibrium.vacuum_toroidal_field.b0[]
    a_minor = dd.equilibrium.time_slice[].boundary.minor_radius
    kappa = dd.equilibrium.time_slice[].boundary.elongation
    delta = dd.equilibrium.time_slice[].boundary.triangularity
    n_e_edge = dd.core_profiles.profiles_1d[].electrons.density_thermal[end]
    Z_eff = dd.core_profiles.profiles_1d[].zeff[end]
    T_e_edge = dd.core_profiles.profiles_1d[].electrons.temperature[end]

    A = R0 / a_minor
    B_theta = IMAS.mks.μ0 *  dd.equilibrium.time_slice[].global_quantities.ip / (2π *  a_minor)

    k_eff = sqrt((1 + kappa^2 * (1 + 2*delta^2 - 1.2*delta^3)) / 2)

    q_cyl = (B_T / B_theta) * (k_eff / A)

    αt = 3.13e-18 * R0 * q_cyl^2 * (n_e_edge * Z_eff) / (T_e_edge^2)

    return αt
end