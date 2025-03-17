""" Physical and dimensional conversions """

from datetime import datetime, date, timedelta
import parameters as params
import initialize as init
import pandas as pd
import numpy as np

# 3/2/2025

def key2plotparams(key):

    if key == 'tem': label, cmap, clr = 'Temperature [C]', 'rainbow', 'black'
    if key == 'ft1': label, cmap, clr = 'Temp factor 1', 'rainbow', 'black'
    if key == 'ft2': label, cmap, clr = 'Temp factor 2', 'rainbow', 'black'
    if key == 'wat': label, cmap, clr = 'Liquid Water [%]', 'Blues', 'blue'
    if key == 'ice': label, cmap, clr = 'Ice [%]', 'Blues', 'blue'
    if key == 'org': label, cmap, clr = 'Organic [%]', 'Oranges', 'orange'
    if key == 'min': label, cmap, clr = 'Mineral [%]', 'Wistia', 'brown'
    if key == 'air': label, cmap, clr = 'Air [%]', 'gray', 'gray'
    if key == 'sat': label, cmap, clr = 'Pore Saturation [%]', 'rainbow', 'black'
    if key == 'doc': label, cmap, clr = 'DOC', 'Oranges', 'orange'
    if key == 'roo': label, cmap, clr = 'Root Density [kg/m3]', 'rainbow', 'black'
    if key == 'rxa': label, cmap, clr = 'Root X-Sectional Area [m2/m3]', 'rainbow', 'black'

    if key == 'ansp': label, cmap, clr = 'ANSP', 'rainbow', 'black'
    if key == 'aesp': label, cmap, clr = 'AESP', 'rainbow', 'black'
    if key == 'hoag': label, cmap, clr = 'HOAG', 'rainbow', 'black'
    if key == 'h2mg': label, cmap, clr = 'H2MG', 'rainbow', 'black'
    if key == 'acmg': label, cmap, clr = 'ACGM', 'rainbow', 'black'
    if key == 'mtox': label, cmap, clr = 'MTOX', 'rainbow', 'black'
    if key == 'aero': label, cmap, clr = 'AERO', 'rainbow', 'black'

    if key == 'ch4': label, cmap, clr = 'CH4', 'Purples', 'purple'
    if key == 'o2' : label, cmap, clr = 'O2', 'Greens', 'green'
    if key == 'co2': label, cmap, clr = 'CO2', 'Blues', 'blue'
    if key == 'h2' : label, cmap, clr = 'H2', 'Oranges', 'orange'
    if key == 'ace': label, cmap, clr = 'Acetate', 'Reds', 'red'
    if key == 'c'  : label, cmap, clr = 'C', 'Wistia', 'brown'
    #if key == 'doc': label, cmap = 'DOC', 'Wistia'

    return(label, cmap, clr)


def doy2datetime(year, doy):

    return datetime(year, 1, 1) + timedelta(doy - 1)


def datetime2doy(date_time):

    day_of_year = date_time.timetuple().tm_yday

    return day_of_year


def interp2grid(depths, values, newgrid):

    ynew = np.interp(newgrid, depths, values)

    return newgrid, ynew


def diffgrid2modelgrid(lin_grid, cons_depths = []):

    exp_thicks, exp_depths = init.define_layers()

    lin_thicks = params.total_depth / params.diff_n_dz
    lin_depths = np.linspace(lin_thicks, params.total_depth, params.diff_n_dz)

    # Interpolate the gas concentrations onto the exponentially spaced grid
    exp_grid = np.interp(exp_depths, lin_depths, lin_grid)

    if len(cons_depths) == 0:  # Correct the entire column

        # Find the scaling factor needed for mass conservation
        lin_mass = np.sum(lin_grid) * lin_thicks
        exp_mass = np.dot(exp_grid, exp_thicks)
        scaling_factor = lin_mass / exp_mass

        # Scale the concentrations on the exponential grid to conserve mass
        exp_grid *= scaling_factor

    else:  # Correct only within the specified depth range

        # Find the depth indices where mass conservation should be applied
        exp_idx_1 = np.where(exp_depths == cons_depths[0])[0][0]
        exp_idx_2 = np.where(exp_depths == cons_depths[1])[0][0]
        lin_idx_1 = np.where(abs(lin_depths - cons_depths[0]) == np.min(abs(lin_depths - cons_depths[0])))[0][0]
        lin_idx_2 = np.where(abs(lin_depths - cons_depths[1]) == np.min(abs(lin_depths - cons_depths[1])))[0][0]

        lin_mass_a = np.sum(lin_grid[0:lin_idx_1]) * lin_thicks
        lin_mass_b = np.sum(lin_grid[lin_idx_1:lin_idx_2]) * lin_thicks
        lin_mass_c = np.sum(lin_grid[lin_idx_2:]) * lin_thicks

        exp_mass_a = np.dot(exp_grid[0:exp_idx_1], exp_thicks[0:exp_idx_1])
        exp_mass_b = np.dot(exp_grid[exp_idx_1:exp_idx_2], exp_thicks[exp_idx_1:exp_idx_2])
        exp_mass_c = np.dot(exp_grid[exp_idx_2:], exp_thicks[exp_idx_2:])

        # Find the scale factor needed for mass conservation
        lin_mass = np.sum(lin_grid[lin_idx_1:]) * lin_thicks
        exp_mass = np.dot(exp_grid[exp_idx_1:], exp_thicks[exp_idx_1:])
        scaling_factor = lin_mass / exp_mass

        # Find the scale factor needed for mass conservation
        sclf_a, sclf_b, sclf_c = lin_mass_a / exp_mass_a, lin_mass_b / exp_mass_b, lin_mass_c / exp_mass_c

        # Scale the concentrations on the exponential grid to conserve mass
        exp_grid[0:exp_idx_1] *= sclf_a
        exp_grid[exp_idx_1:exp_idx_2] *= sclf_b
        exp_grid[exp_idx_2:] *= sclf_c

    return exp_grid


def map_depth_values(exp_depth_grid, high_res_depth_grid):
    # Initialize an array to store mapped values
    mapped_values = np.zeros_like(high_res_depth_grid)

    # For each index in the higher resolution grid, find the corresponding depth in the exponential grid
    exp_index = 0
    for i, depth in enumerate(high_res_depth_grid):
        # While the depth is greater than the current exponential depth
        while exp_index < len(exp_depth_grid) - 1 and depth > exp_depth_grid[exp_index + 1]:
            exp_index += 1

        # Assign the value of the exponential grid at the corresponding depth range
        mapped_values[i] = exp_depth_grid[exp_index]

    return mapped_values


def modelgrid2diffgrid(exp_grid, cons_depths = []):

    exp_thicks, exp_depths = init.define_layers()

    lin_thicks = params.total_depth / params.diff_n_dz
    lin_depths = np.linspace(lin_thicks, params.total_depth, params.diff_n_dz)
    lin_depths = [np.round(lin_depths[i], 2) for i in range(0, params.diff_n_dz)]

    # Interpolate the gas concentrations onto the linearly spaced grid
    lin_grid = np.interp(lin_depths, exp_depths, exp_grid)

    #print('exp before: ', np.dot(exp_grid[:], exp_thicks[:]))

    if len(cons_depths) == 0:  # Correct the entire column

        # Mass conservation by scaling the concentrations
        lin_mass = np.sum(lin_grid) * lin_thicks
        exp_mass = np.dot(exp_grid, exp_thicks)

        # Scale the concentrations on the exponential grid to conserve mass
        scaling_factor = exp_mass / lin_mass
        lin_grid *= scaling_factor

    else:  # Correct only within the specified depth range

        # Find the depth indices where mass conservation should be applied
        exp_idx_1 = np.where(exp_depths == cons_depths[0])[0][0]
        exp_idx_2 = np.where(exp_depths == cons_depths[1])[0][0]
        lin_idx_1 = np.where(abs(lin_depths - cons_depths[0]) == np.min(abs(lin_depths - cons_depths[0])))[0][0]
        lin_idx_2 = np.where(abs(lin_depths - cons_depths[1]) == np.min(abs(lin_depths - cons_depths[1])))[0][0]

        lin_mass_a = np.sum(lin_grid[0:lin_idx_1]) * lin_thicks
        lin_mass_b = np.sum(lin_grid[lin_idx_1:lin_idx_2]) * lin_thicks
        lin_mass_c = np.sum(lin_grid[lin_idx_2:]) * lin_thicks

        exp_mass_a = np.dot(exp_grid[0:exp_idx_1], exp_thicks[0:exp_idx_1])
        exp_mass_b = np.dot(exp_grid[exp_idx_1:exp_idx_2], exp_thicks[exp_idx_1:exp_idx_2])
        exp_mass_c = np.dot(exp_grid[exp_idx_2:], exp_thicks[exp_idx_2:])

        # Find the scale factor needed for mass conservation
        sclf_a, sclf_b, sclf_c = exp_mass_a / lin_mass_a, exp_mass_b / lin_mass_b, exp_mass_c / lin_mass_c

        # Scale the concentrations on the exponential grid to conserve mass
        lin_grid[0:lin_idx_1] *= sclf_a
        lin_grid[lin_idx_1:lin_idx_2] *= sclf_b
        lin_grid[lin_idx_2:] *= sclf_c

    #print('lin after: ', np.sum(lin_grid[:]) * lin_thicks)

    return lin_grid


def harmonic_mean(a, b):

    return 2 * a * b / (a + b)


def volpct2kgm3(volpct, density):

    return (volpct / 100.0) * density


def unitm3s2unitm2(unitm3s, thicks):

    return np.dot(unitm3s, thicks)


def ddmmyyyy2doy(ddmmyyyy):
    doy = datetime(int(ddmmyyyy[4:]), int(ddmmyyyy[2:4]), int(ddmmyyyy[0:2])).timetuple().tm_yday
    return(doy)


def datedprofiles2grid(days, profiles, n_times):

    # Transpose the values array so that each row corresponds to the values at a specific index across days
    transp_profs = np.array(profiles).T  # Each row will represent a different value index over the days

    interp_profs, target_day = [], 0.0
    for p in range(0, n_times):
        interp_profs.append(np.array([np.interp(target_day, days, transp_profs[i])
                                      for i in range(transp_profs.shape[0])]))
        target_day += 365.0 / n_times

    return np.transpose(np.array(interp_profs))


def temporally_extrude_profile(profile, ntimes):

    return np.transpose(np.array([profile for n in range(0, ntimes)]))


def vertically_extrude_profile(profile, ndepths):

    return np.array([profile for n in range(0, ndepths)])


def get_var_name(var):
    for name, value in globals().items():
        if value is var:
            return name


def resampleDTvals(datetimes, values, resampling_interval):

    datetime_strings = datetimes  # Sample input data: datetime strings and corresponding values
    dates = pd.to_datetime(datetime_strings)  # Convert datetime strings to pandas datetime objects

    df = pd.DataFrame({'date': dates, 'value': values})  # Create a pandas DataFrame
    df.set_index('date', inplace=True)  # Set the date column as the index

    resampled_means = df.resample(resampling_interval).mean()  # Resample the data to weekly frequency
    resampled_dates = resampled_means.index
    resampled_values = resampled_means['value'].values

    return(resampled_dates, resampled_values)
