""" Physical properties of the model layers, e.g. effective diffusivities, temperature factors, equilibrium values """

import parameters as params
import numpy as np
import system

# 2/28/2025

def diffusivity(layer_temps):

    temps_k = np.array(layer_temps) + 273.15
    sec_to_day, th, phi = 8.64e4, 273.15, 298.0 # sec-to-day conversion, reference temperatures

    chems = system.list_user_chemicals()

    adfs = {} # Air diffusivities
    for i in range(0, len(chems)):
        chem = chems[i]
        adfs[chem] = sec_to_day * params.f_d_a * params.cprops[chem]['DA1'] * (temps_k/th)**params.cprops[chem]['DA2']
        adfs[chem][temps_k < 0.0] = 0.0

    wdfs = {} # Water diffusivities
    for i in range(0, len(chems)):
        chem = chems[i]
        wdfs[chem] = sec_to_day * params.f_d_w * params.cprops[chem]['DW1'] * (temps_k/phi)
        if chem == 'co2': wdfs[chem] = sec_to_day * params.f_d_w * params.cprops[chem]['DW1'] * np.exp(-2032.6/temps_k)
        wdfs[chem][temps_k < 0.0] = 0.0

    diffs_dictionary = {'air': adfs, 'wat': wdfs}

    return diffs_dictionary


def effective_diffusivity(satpct, layer_temps):

    eff_diffs, eff_diffs_full = {}, {}
    diffs = diffusivity(layer_temps)
    ww, aw = 0.01 * satpct, 1.0 - 0.01 * satpct  # Weights representing relative water and air volumes, for averaging
    for s in list(diffs['air'].keys()):

        eff_diffs[s] = (ww * diffs['wat'][s] + aw * diffs['air'][s]) / (ww + aw)
        eff_diffs_full[s] = (ww * diffs['wat'][s] + aw * diffs['air'][s]) / (ww + aw)

        # Confine diffusion calculations to the saturated region - easiest way is to apply zero-diffusivity masks here
        eff_diffs[s][layer_temps <= 0.0] = 0.0
        #eff_diffs[s] = diffs['wat'][s]
        #eff_diffs[s][satpct < 100.0] = 0.0

        #eff_diffs[s][layer_temps <= 10.0] = (layer_temps[layer_temps <= 10.0] / 10.0) * eff_diffs[s][layer_temps <= 10.0]

    return eff_diffs, eff_diffs_full


def temperature_factor_q10(t, q10):

    # Output: 'f_tX' = soil temperature factors ()
    # Input: temperature 't' (C); temperature sensitivity 'q10' (1/10K)

    f_t = q10 ** ((t - params.ft_max_temp) / 10.0) # t - 30.0
    # print('f_t: ',f_t)
    f_t[t < 0.0] = 0.0
    f_t[t > params.ft_max_temp] = 1.0 # > 30.0

    return f_t


def temperature_factors_1_and_2(t):

    # Output: 'f_t1, f_t2' = soil temperature factors ()
    # Input: temperature 't' (C)

    t_min1, t_max1, t_opt1 = params.t_min1, params.t_max1, params.t_opt1
    t_min2, t_max2, t_opt2 = params.t_min2, params.t_max2, params.t_max2

    f_t1 = ((t - t_min1) * (t - t_max1)) / ((t - t_min1) * (t - t_max1) - ((t - t_opt1) ** 2.0))
    f_t1[t < t_min1] = 0.0
    f_t1[t > t_max1] = 0.0
    f_t2 = ((t - t_min2) * (t - t_max2)) / ((t - t_min2) * (t - t_max2) - ((t - t_opt2) ** 2.0))
    f_t2[t < t_min2] = 0.0
    f_t2[t > t_max2] = 0.0

    return f_t1, f_t2


def q10_temperature_factors(t):

    reacts = list(params.reactions.keys())
    q10_temperature_factors = {}
    for i in range(0, len(reacts)):
        react, q10 = reacts[i], params.reactions[reacts[i]]['q10']['val']
        f_t = temperature_factor_q10(t, q10)
        q10_temperature_factors[react] = f_t

    return q10_temperature_factors


def equilibrium_concentration(temps):

    th = 273.15  # Reference temp (K)
    t2 = np.copy(temps)
    t2 += 273.15 # Convert C to K

    chems = system.list_user_chemicals()
    henrys = {}
    for i in range(0, len(chems)):
        henrys[chems[i]] = params.cprops[chems[i]]['H1'] * np.exp(params.cprops[chems[i]]['H2'] * ((1.0/t2)-(1.0/th)))

    partial_pressures = partial_pressure()

    keys = list(henrys.keys())
    equilibrium_concentrations = {}
    for i in range(0, len(keys)):
        equilibrium_concentrations[keys[i]] = partial_pressures[keys[i]] * henrys[keys[i]]

    return(equilibrium_concentrations)


def partial_pressure():

    chems = system.list_user_chemicals()
    partial_pressures = {}
    for i in range(0, len(chems)): partial_pressures[chems[i]] = (params.cprops[chems[i]]['airvol'] / 1e2) * 101.325e3

    return(partial_pressures)


def dissolved_organic_carbon(org_cpls, temps, sats, thks):

    return params.k_cpool * (org_cpls / thks) * temperature_factor_q10(temps, params.doc_prod_q10) * 0.01 * sats


def schmidt_number(temps):

    chems = system.list_user_chemicals()
    schmidts = {}
    for i in range(0, len(chems)):
        chem = chems[i]
        schmidts[chem] = params.cprops[chem]['S1'] - params.cprops[chem]['S2'] * temps + \
                         params.cprops[chem]['S3'] * (temps ** 2) - params.cprops[chem]['S4'] * (temps ** 3)

    return(schmidts)


def transfer_velocity(schm_nums):

    u_10, n = 0.0, -0.5
    phi_600 = 2.07 + 0.215 * u_10 ** 1.7

    # Calculate transfer velocities, converting phi [cm/hr] to phi [m/day] with a 0.24x factor
    # Force a 'zero schmidt number = zero transfer velocity' condition to avoid division by zero in the 'nth' power
    chems = system.list_user_chemicals()
    phis = {}
    for i in range(0, len(chems)):
        chem = chems[i]
        if np.max(schm_nums[chem]) == 0.0:
            phis[chem] = 0.0 * schm_nums[chem]
        else:
            phis[chem] = 0.24 * phi_600 * (schm_nums[chem] / 600.0) ** n

    return(phis)