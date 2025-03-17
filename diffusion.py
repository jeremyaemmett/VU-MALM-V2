import numpy as np
import conversions
import initialize as init
import parameters as params
from scipy.linalg import solve_banded

# 3/17/2025

def setup_diffusivity(N, x, D_layers, layer_boundaries):

    """ Assign diffusivities to spatial points based on layer boundaries
    :param N:
    :param x:
    :param D_layers:
    :param layer_boundaries:
    :return:
    """

    D = np.zeros(N)
    for j in range(0, len(layer_boundaries)):
        for i in range(N):
            if j == 0 and x[i] < layer_boundaries[0]: D[i] = D_layers[0]
            if j > 0 and (layer_boundaries[j-1] < x[i] < layer_boundaries[j]):
                D[i] = D_layers[j]
    return D


def CNFV(D_layers, top, bottom, thicks, model_dt, steps, C, species):

    """ Solve diffusion using a Crank-Nicolson finite volume method
    :param D_layers:
    :param top:
    :param bottom:
    :param thicks:
    :param model_dt:
    :param steps:
    :param C:
    :param species:
    :return:
    """

    # Parameters
    N, L = params.diff_n_dz, sum(thicks) # Number of layers (control volumes), Length of domain
    dx, dt = L / N, model_dt / steps # Spatial step size, Time step size

    # Spatial grid
    x = np.linspace(dx / 2, L - dx / 2, N) # Cell-centered grid

    # Assign diffusivity to each spatial point
    layer_boundaries = np.cumsum(thicks)
    D = setup_diffusivity(N, x, D_layers, layer_boundaries)

    # Number of diffusion layers (saturated layers)
    indices = np.where((x > top) & (x <= bottom))[0]
    N = len(indices)

    # Solve diffusion if there's more than 2 un-saturated layers
    if N > 2:

        # Diffusivities at the interfaces
        d_interfaces = np.zeros(N + 1)
        for i in range(1, N): d_interfaces[i] = D[i + indices[0]]
        d_interfaces[0], d_interfaces[-1] = D[indices[0]], D[indices[-1]]

        r = (dt / (2 * dx**2)) * d_interfaces  # Crank-Nicolson Coefficients

        # Construct the A (implicit) matrix
        lower, main, upper = np.zeros(N), np.zeros(N), np.zeros(N)  # Sub, Main, and Super diagonals
        for i in range(1, N - 1): lower[i], upper[i], main[i] = -r[i - 1], -r[i], 1 + r[i] + r[i - 1]
        # Neumann boundary condition
        main[0] = main[-1] = 1
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
        # Neumann boundary condition
        B[0, 0] = B[-1, -1] = 1
        B[0, 1] = B[-1, -2] = -1

        diffgrid_mass_1 = np.sum(C[indices]) * dx  # Optional: Integrated mass 1 for conservation verification

        for n in range(steps):  # Time-stepping loop

            RHS = B @ C[indices]

            C_new = solve_banded((1, 1), A_banded, RHS)  # Solve the linear system A * C_new = RHS

            C[indices] = C_new

            diffgrid_mass_n = np.sum(C) * dx  # Optional: Integrated mass 2 for conservation verification

        diffgrid_mass_2 = np.sum(C[indices]) * dx
        C[indices] *= diffgrid_mass_1/diffgrid_mass_2

    return C

def transport(species, chd, dtd, ifd, forcing_t, dt, level_change):

    """ Predict diffusive transport fluxes
    :param species:
    :param chd:
    :param dtd:
    :param ifd:
    :param forcing_t:
    :param dt:
    :param level_change:
    :return:
    """

    # CR-solver
    thicks, depths = init.define_layers()
    iceline = depths[np.where(forcing_t['tem'] <= 0.0)[0][0]]
    watline = depths[np.where(forcing_t['sat'] < 100.0)[0][-1]] if len(np.where(forcing_t['sat'] < 100.0)[0]) > 0 else depths[-1] # Don't change this from depths[np.where(forcing_t['sat'] < 100.0)[0][-1]]!
    unsaturated_ids = np.where(forcing_t['sat'] < 100.0)[0]

    U1 = chd[species]['conc']

    # If the entire soil column is thawed and there's saturated layers above the iceline,
    # run the diffusion scheme in the saturated layers
    int_flux, first_sat_thick = 0.0, 0.0
    if iceline != params.total_depth and iceline > watline:

        if len(unsaturated_ids) > 0 and np.nanmax(forcing_t['efd'][species]) != 0.0:

            # Note some important indices and depths
            first_sat_layer = unsaturated_ids[-1] + 1

            # Calculate a water-air interface flux. Defaults to zero if there's no gradient at the interface.
            int_flux = 0.0
            if U1[first_sat_layer] > 1e-10 or U1[first_sat_layer-1] > 1e-10:
                # First guess of an interface flux
                int_flux = interface_flux(species, chd, forcing_t, unsaturated_ids[-1] + 1)
                # If the flux will result in removal > availability in the first saturated layer, cap it
                if int_flux < 0.0 and abs((int_flux * dt) / thicks[first_sat_layer]) > U1[first_sat_layer]:
                    int_flux = -(U1[first_sat_layer] * thicks[first_sat_layer])/dt
            ifd[species] = int_flux

            # Modify the first saturated layer and overlying sub-saturated layers
            if params.interface_flux:
                U1[first_sat_layer] += (int_flux * dt) / thicks[first_sat_layer]

        # Parameters
        N, L = params.diff_n_dz, sum(thicks)  # Number of layers (control volumes), Length of domain
        dx = L / N

        # Convert exponential model grid to linear diffusion grid
        # if species == 'ch4': print('ch4 before conversion', np.dot(chd['ch4']['conc'], thicks))
        U1b = conversions.modelgrid2diffgrid(chd[species]['conc'], [watline, iceline])
        # if species == 'ch4': print('ch4 after U1 conversion', np.sum(U1b) * dx)

        U2b = CNFV(forcing_t['efd'][species], watline, iceline, thicks, dt, params.diff_n_dt, U1b, species)

        # Convert diffusion to model grid
        # if species == 'ch4': print('lin 2: ', np.sum(U2b) * dx)
        U2 = conversions.diffgrid2modelgrid(U2b, [watline, iceline])
        # if species == 'ch4': print('exp 2: ', np.dot(U2, thicks))

        # Modify the overlying unsaturated layers
        if len(unsaturated_ids) > 0 and np.nanmax(forcing_t['efd'][species]) != 0.0:
            if params.interface_flux:
                U2[0:first_sat_layer] -= (int_flux * dt) / depths[first_sat_layer-1]
                col_tot_dist_conc = np.dot(U2[0:first_sat_layer], thicks[0:first_sat_layer]) / depths[first_sat_layer - 1]
                U2[0:first_sat_layer] = col_tot_dist_conc

        # Ensure change = 0 if diffusivity = 0 in every layer
        if np.nanmax(forcing_t['efd'][species]) == 0.0: U2 = U1

        dtd[species]['prof'] = (U2 - U1) / dt  # Record the rates of change

    else:

        # Ensure change = 0 if diffusivity = 0 in every layer
        if np.nanmax(forcing_t['efd'][species]) == 0.0: U2 = U1

        dtd[species]['prof'] = (U1 - U1) / dt  # Record the rates of change

    dtd[species]['prof'][abs(dtd[species]['prof']) < 1e-10] = 0.0  # Suppress numerical noise around the minimum value
    #dtd[species]['prof'][forcing_t['tem'] <= 0.0] = 0.0  # Ensure zero change in frozen layers

    return(dtd, ifd)


def interface_flux(species, chd, forcing_t, layer):

    """ Calculates the flux across the water-air interface (mol/m2/day)
    :param species:
    :param chd:
    :param forcing_t:
    :param layer:
    :return:
    """

    tran_vels = forcing_t['tvel'][species]

    intflux = 0.0 if forcing_t['tem'][layer] <= 0.0 or params.interface_flux == False \
        else -1.0 * tran_vels[layer] * (chd[species]['conc'][layer] - chd[species]['conc'][layer-1])

    return(intflux)

