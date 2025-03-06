import properties
import numpy as np
import conversions
import initialize as init
import parameters as params
from scipy.linalg import solve_banded
from scipy.ndimage import gaussian_filter

# 3/2/2025

def setup_diffusivity(N, x, D_layers, layer_boundaries):

    # Assign diffusivities to spatial points based on layer boundaries
    D = np.zeros(N)
    for j in range(0, len(layer_boundaries)):
        for i in range(N):
            if j == 0 and x[i] < layer_boundaries[0]: D[i] = D_layers[0]
            if j > 0 and (layer_boundaries[j-1] < x[i] < layer_boundaries[j]):
                D[i] = D_layers[j]
    return D


def CNFV(D_layers, bottom, layer_thicknesses, model_dt, steps, C, int_flux, first_sat_thick):

    # Parameters
    N, L = params.diff_n_dz, sum(layer_thicknesses) # Number of layers (control volumes), Length of domain
    dx, dt = L / N, model_dt / steps # Spatial step size, Time step size

    layer_boundaries = np.cumsum(layer_thicknesses)

    # Spatial grid
    x = np.linspace(dx / 2, L - dx / 2, N) # Cell-centered grid

    first_sat_idx = 0
    D = setup_diffusivity(N, x, D_layers, layer_boundaries)  # Assign diffusivity to each spatial point
    if len(np.where(D != 0.0)[0]) > 0: first_sat_idx = np.where(D != 0.0)[0][0]
    second_sat_idx = first_sat_idx + int(first_sat_idx + first_sat_thick / dx)

    # Diffusivities at the interfaces
    d_interfaces = np.zeros(N + 1)
    for i in range(1, N): d_interfaces[i] = D[i]
    d_interfaces[0], d_interfaces[-1] = D[0], D[-1]

    r = (dt / (2 * dx**2)) * d_interfaces  # Crank-Nicolson Coefficients

    N = len(np.where(x <= bottom)[0]) + 1  # Number of diffusion layers (saturated layers)

    # Construct the A (implicit) matrix
    lower, main, upper = np.zeros(N), np.zeros(N), np.zeros(N)  # Sub, Main, and Super diagonals
    for i in range(1, N - 1): lower[i], upper[i], main[i] = -r[i - 1], -r[i], 1 + r[i] + r[i - 1]
    main[0] = main[-1] = 1  # Neumann boundary conditions (zero-flux at edges)
    upper[0] = lower[-1] = -1
    A_banded = np.vstack([np.hstack([0, upper[:-1]]), main, np.hstack([lower[1:], 0])])

    # Construct the B (explicit) matrix
    B = np.zeros((N, N))
    for i in range(N):
        B[i, i] = 1 - r[i + 1] - r[i]
        if i > 0:
            B[i, i-1] = r[i]
        if i < N - 1:
            B[i, i + 1] = r[i + 1]
    B[0, 0] = B[-1, -1] = 1  # Neumann condition
    B[0, 1] = B[-1, -2] = -1

    diffgrid_mass_1 = np.sum(C) * dx

    for n in range(steps):  # Time-stepping loop

        #C[first_sat_idx:second_sat_idx] += ((dt * int_flux) / dx) / (second_sat_idx - first_sat_idx)

        RHS = B @ C[0:N]  # Right-hand side for the linear system

        C_new = solve_banded((1, 1), A_banded, RHS)  # Solve the linear system A * C_new = RHS

        C = np.append(C_new, C[N:])

        diffgrid_mass_n = np.sum(C) * dx  # Optional: Print mass to verify conservation

    diffgrid_mass_2 = np.sum(C) * dx

    return C

def diffusion(species, chd, dtd, sfd, forcing_t, dt):

    # CR-solver
    thicks, depths = init.define_layers()
    iceline = depths[np.where(forcing_t['tem'] <= 0.0)[0][0]-1]
    watline = depths[np.where(forcing_t['sat'] < 100.0)[0][-1]] if len(np.where(forcing_t['sat'] < 100.0)[0]) > 0 else depths[-1]
    unsaturated_ids = np.where(forcing_t['sat'] < 100.0)[0]

    U1 = chd[species]['conc']

    int_flux = 0.0
    first_sat_thick = 0.0
    if iceline != params.total_depth:

        if len(unsaturated_ids) > 0 and np.nanmax(forcing_t['efd'][species]) != 0.0:

            # Calculate the water-to-air interface flux, and modify the boundary layers
            first_sat_layer = unsaturated_ids[-1] + 1
            second_sat_layer = first_sat_layer + 1
            first_sat_depth, second_sat_depth = depths[first_sat_layer], depths[second_sat_layer]
            first_sat_thick = second_sat_depth - first_sat_depth

            # First guess of an interface flux
            int_flux = interface_flux(species, chd, forcing_t, unsaturated_ids[-1] + 1)
            # If the flux will result in removal > availability in the first saturated layer, cap it
            if int_flux < 0.0 and abs((int_flux * dt) / thicks[first_sat_layer]) > U1[first_sat_layer]:
                int_flux = -(U1[first_sat_layer] * thicks[first_sat_layer])/dt

            # Modify the first saturated layer
            if params.interface_flux: U1[first_sat_layer] += (int_flux * dt) / thicks[first_sat_layer]
            sfd[species] = int_flux

        # Convert exponential model grid to linear diffusion grid
        U1b = conversions.modelgrid2diffgrid(chd[species]['conc'])

        U2b = CNFV(forcing_t['efd'][species], iceline, thicks, dt, params.diff_n_dt, U1b, int_flux, first_sat_thick)

        # Convert diffusion to model grid
        U2 = conversions.diffgrid2modelgrid(U2b)

        # Ensure change = 0 if diffusivity = 0 in every layer
        if np.nanmax(forcing_t['efd'][species]) == 0.0: U2 = U1

        dtd[species]['prof'] = (U2 - U1) / dt  # Record the rates of change

    else:

        # Ensure change = 0 if diffusivity = 0 in every layer
        if np.nanmax(forcing_t['efd'][species]) == 0.0: U2 = U1

        dtd[species]['prof'] = (U1 - U1) / dt  # Record the rates of change

    return(dtd)


def interface_flux(species, chd, forcing_t, layer):
    # Calculates flux across the water-air interface (mol/m2/day)

    tran_vels = forcing_t['tvel'][species]
    equi_conc = forcing_t['eqc'][species]

    intflux = 0.0 if forcing_t['tem'][layer] <= 0.0 \
        else -1.0 * tran_vels[layer] * (chd[species]['conc'][layer] - equi_conc[layer])

    return(intflux)

