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


def CNFV(D_layers, top, bottom, thicks, model_dt, steps, C, int_flux, first_sat_thick, species):

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
        #print(d_interfaces)
        #print(D)
        #print(indices[0], D[i + indices[0]])
        #print(d_interfaces)

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

        diffgrid_mass_1 = np.sum(C[indices]) * dx
        C1 = np.copy(C)
        for n in range(steps):  # Time-stepping loop

            RHS = B @ C[indices]

            C_new = solve_banded((1, 1), A_banded, RHS)  # Solve the linear system A * C_new = RHS

            #CA, CB, CC = C[0:indices[0]], C_new, C[(indices[-1]+1):]
            #C = np.append(CA, np.append(CB, CC))
            C[indices] = C_new

            diffgrid_mass_n = np.sum(C) * dx  # Optional: Print mass to verify conservation

        diffgrid_mass_2 = np.sum(C[indices]) * dx
        C[indices] *= diffgrid_mass_1/diffgrid_mass_2
        C2 = np.copy(C)
        ###if species == 'ch4': print('CH4 before/after CNFV: ', np.sum(C1) * dx, np.sum(C2) * dx)
        #    print(np.sum(C1) * dx)
        #    print(np.sum(C2) * dx)
        #    print(np.sum(C1) * dx - np.sum(C2) * dx)

        #C *= (diffgrid_mass_1/diffgrid_mass_2)

        #if species == 'ch4': print('1, 2, 1-2: ', diffgrid_mass_1, diffgrid_mass_2, diffgrid_mass_1 - diffgrid_mass_2)

    return C

def transport(species, chd, dtd, sfd, forcing_t, dt):

    # CR-solver
    thicks, depths = init.define_layers()
    iceline = depths[np.where(forcing_t['tem'] <= 0.0)[0][0]-1]
    watline = depths[np.where(forcing_t['sat'] < 100.0)[0][-1]] if len(np.where(forcing_t['sat'] < 100.0)[0]) > 0 else depths[-1]
    unsaturated_ids = np.where(forcing_t['sat'] < 100.0)[0]

    U1 = chd[species]['conc']

    int_flux = 0.0
    first_sat_thick = 0.0
    if iceline != params.total_depth:  # If the entire soil column is NOT frozen

        if len(unsaturated_ids) > 0 and np.nanmax(forcing_t['efd'][species]) != 0.0:

            # Calculate the water-to-air interface flux, and modify the boundary layers
            first_sat_layer = unsaturated_ids[-1] + 1
            second_sat_layer = first_sat_layer + 1
            first_sat_depth, second_sat_depth = depths[first_sat_layer], depths[second_sat_layer]
            first_sat_thick = second_sat_depth - first_sat_depth

            # Calculate a water-air interface flux. Defaults to zero if there's no gradient at the interface.
            int_flux = 0.0
            if U1[first_sat_layer] > 1e-10 or U1[first_sat_layer-1] > 1e-10:
                # First guess of an interface flux
                int_flux = interface_flux(species, chd, forcing_t, unsaturated_ids[-1] + 1)
                # If the flux will result in removal > availability in the first saturated layer, cap it
                if int_flux < 0.0 and abs((int_flux * dt) / thicks[first_sat_layer]) > U1[first_sat_layer]:
                    int_flux = -(U1[first_sat_layer] * thicks[first_sat_layer])/dt

            # Modify the first saturated layer and overlying sub-saturated layers
            # Modify the overlying unsaturated layers
            if params.interface_flux:
                U1[first_sat_layer] += (int_flux * dt) / thicks[first_sat_layer]
                U1[0:first_sat_layer] -= (int_flux * dt) / depths[first_sat_layer-1]
                col_tot_dist_conc = np.dot(U1[0:first_sat_layer], thicks[0:first_sat_layer]) / depths[first_sat_layer - 1]
                U1[0:first_sat_layer] = col_tot_dist_conc
            sfd[species] = int_flux
        ###if species == 'ch4': print('ch4 before conversion', np.dot(chd['ch4']['conc'], thicks))

        # Parameters
        N, L = params.diff_n_dz, sum(thicks)  # Number of layers (control volumes), Length of domain
        dx = L / N

        # Convert exponential model grid to linear diffusion grid
        U1b = conversions.modelgrid2diffgrid(chd[species]['conc'])
        ###if species == 'ch4': print('ch4 after U1 conversion', np.sum(U1b) * dx)

        U2b = CNFV(forcing_t['efd'][species], watline, iceline, thicks, dt, params.diff_n_dt, U1b, int_flux, first_sat_thick, species)

        # Convert diffusion to model grid
        U2 = conversions.diffgrid2modelgrid(U2b)
        ###if species == 'ch4': print('ch4 after U2 conversion', np.dot(U2, thicks))

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
        else -1.0 * tran_vels[layer] * (chd[species]['conc'][layer] - chd[species]['conc'][layer-1])

    return(intflux)

