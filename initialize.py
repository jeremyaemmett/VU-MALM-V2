""" Variable definitions and initial values, informed by the user's input in parameters.py """

import parameters as params
import numpy as np
import conversions
import system

# 3/2/2025

def define_layers():

    layer_thicknesses = []
    current_depth = 0
    thickness = params.dz0

    while current_depth + thickness < params.total_depth:
        # Round to the nearest cm
        rounded_thickness = round(thickness * 100) / 100.0
        layer_thicknesses.append(rounded_thickness)
        current_depth += rounded_thickness

        # Increase thickness exponentially
        thickness *= params.growth_rate

    # Add the last layer to reach the total depth, ensuring it's thicker then the previous
    last_layer_thickness = params.total_depth - current_depth

    last_layer_thickness = round(last_layer_thickness * 100) / 100.0
    layer_thicknesses.append(last_layer_thickness) # List of layer thicknesses (m)
    layer_thicknesses = np.array([np.round(layer_thicknesses[i], 2) for i in range(0, len(layer_thicknesses))])

    depths = np.cumsum(layer_thicknesses) # Layer bottom depths (m)
    depths = np.array([np.round(depths[i], 2) for i in range(0, len(depths))])

    return layer_thicknesses, depths


def define_atmfillin():

    unique_chemicals = system.list_user_chemicals()

    layer_thicknesses, depths = define_layers()
    E = np.zeros(len(depths))  # Array of zeros, used repeatedly to initialize profiles

    atmfillin_dictionary = {}
    for i in range(0, len(unique_chemicals)): atmfillin_dictionary[unique_chemicals[i]] = E

    return atmfillin_dictionary


def define_diffusion():

    unique_chemicals = system.list_user_chemicals()

    layer_thicknesses, depths = define_layers()
    E = np.zeros(len(depths))  # Array of zeros, used repeatedly to initialize profiles

    diffusion_dictionary = {}
    for i in range(0, len(unique_chemicals)): diffusion_dictionary[unique_chemicals[i]] = {'prof': E, 'surf': 0.0}

    return diffusion_dictionary


def define_plantdiff():

    unique_chemicals = system.list_user_chemicals()

    layer_thicknesses, depths = define_layers()
    E = np.zeros(len(depths))  # Array of zeros, used repeatedly to initialize profiles

    plantdiff_dictionary = {}
    for i in range(0, len(unique_chemicals)): plantdiff_dictionary[unique_chemicals[i]] = {'prof': E, 'surf': 0.0}

    return plantdiff_dictionary


def define_chemistry():

    unique_chemicals = system.list_user_chemicals()
    reactions_involving_chemicals = system.list_reactions_involving_chemicals()

    layer_thicknesses, depths = define_layers()

    # Initial depth-uniform species concentrations
    initc = 1e-10

    # Initial concentration
    E = np.zeros(len(depths))  # Array of zeros, used repeatedly to initialize profiles
    L = sum(layer_thicknesses)  # Length of the domain
    dx = sum(layer_thicknesses) / params.diff_n_dz  # Spatial step size
    x = np.linspace(dx / 2, L - dx / 2, params.diff_n_dz)  # Cell-centered grid

    chem_dict = {}
    for i in range(0, len(unique_chemicals)):
        chemical = unique_chemicals[i]
        chem_dict[chemical] = {'conc': E + initc, 'netp': E, 'prod': E, 'cons': E}
        if params.gaussian_profile:
            chem_dict[chemical]['conc'] = conversions.diffgrid2modelgrid(10.0 * np.exp(-((x - L / 4.0) ** 2) / 0.01))
        cons_dict, prod_dict = {}, {}
        for m in range(0, len(reactions_involving_chemicals[chemical]['cons'])):
            cons_dict[reactions_involving_chemicals[chemical]['cons'][m]] = E
        for n in range(0, len(reactions_involving_chemicals[chemical]['prod'])):
            prod_dict[reactions_involving_chemicals[chemical]['prod'][n]] = E
        chem_dict[chemical]['cons2'] = cons_dict
        chem_dict[chemical]['prod2'] = prod_dict

    # Michaelis-Menten term decomposition for each pathway
    mmkt_dict = {'vmax': {'acmg': E, 'aesp': E, 'ansp': E, 'hoag': E, 'h2mg': E, 'mtox': E},
                 'subs': {'acmg': E, 'aesp': E, 'ansp': E, 'hoag': E, 'h2mg': E, 'mtox': E},
                 'refr': {'acmg': E, 'aesp': E, 'ansp': E, 'hoag': E, 'h2mg': E, 'mtox': E}}

    return chem_dict, mmkt_dict


def define_srffluxes():

    unique_chemicals = system.list_user_chemicals()

    srffluxes_dictionary = {}
    for i in range(0, len(unique_chemicals)): srffluxes_dictionary[unique_chemicals[i]] = 0.0

    return srffluxes_dictionary


def define_intfluxes():

    unique_chemicals = system.list_user_chemicals()

    srffluxes_dictionary = {}
    for i in range(0, len(unique_chemicals)): srffluxes_dictionary[unique_chemicals[i]] = 0.0

    return srffluxes_dictionary


def variables():

    chemistry_dictionary, mmrdecomp_dictionary = define_chemistry()
    atmfillin_dictionary = define_atmfillin()
    diffusion_dictionary = define_diffusion()
    plantdiff_dictionary = define_plantdiff()
    srffluxes_dictionary = define_srffluxes()
    intfluxes_dictionary = define_intfluxes()
    dict_dict = {'chd': chemistry_dictionary, 'mmd': mmrdecomp_dictionary,
                 'atd': atmfillin_dictionary, 'dtd': diffusion_dictionary, 'ptd': plantdiff_dictionary,
                 'ifd': intfluxes_dictionary, 'sfd': srffluxes_dictionary}

    return chemistry_dictionary, mmrdecomp_dictionary, atmfillin_dictionary, \
           diffusion_dictionary, plantdiff_dictionary, intfluxes_dictionary, \
           srffluxes_dictionary, dict_dict
