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


def diffgrid2modelgrid(diffgrid, *args, **kwargs):

    cons_ids = kwargs.get('cons_ids', None)  # Optional index range in which to force mass conservation
    thicks, depths = init.define_layers()

    diff_dz = params.total_depth / params.diff_n_dz
    diff_depths = np.linspace(diff_dz, params.total_depth, params.diff_n_dz)

    # Interpolate the gas concentrations onto the exponentially spaced grid
    modelgrid = np.interp(depths, diff_depths, diffgrid)

    # Mass conservation by scaling the concentrations
    linear_mass = np.trapz(diffgrid, diff_depths)  # Integral over the linear grid
    exp_mass = np.trapz(modelgrid, depths)  # Integral over the exponential grid

    # Mass conservation by scaling the concentrations
    linear_mass = np.sum(diffgrid) * diff_dz
    exp_mass = np.dot(modelgrid, thicks)

    # Scale the concentrations on the exponential grid to conserve mass
    scaling_factor = linear_mass / exp_mass
    modelgrid *= scaling_factor

    return modelgrid


def modelgrid2diffgrid(modelgrid, *args, **kwargs):

    cons_ids = kwargs.get('cons_ids', None)  # Optional index range in which to force mass conservation
    thicks, depths = init.define_layers()

    diff_dz = params.total_depth / params.diff_n_dz
    diff_depths = np.linspace(diff_dz, params.total_depth, params.diff_n_dz)
    diff_depths = [np.round(diff_depths[i], 2) for i in range(0, params.diff_n_dz)]

    # Interpolate the gas concentrations onto the linearly spaced grid
    diffgrid = np.interp(diff_depths, depths, modelgrid)

    # Mass conservation by scaling the concentrations
    diffgrid_mass = np.trapz(diffgrid, diff_depths)  # Integral over the linear grid
    diffgrid_mass = np.sum(diffgrid) * diff_dz
    modlgrid_mass = np.trapz(modelgrid, depths)  # Integral over the exponential grid
    modlgrid_mass = np.dot(modelgrid, thicks)

    # Scale the concentrations on the exponential grid to conserve mass
    scaling_factor = modlgrid_mass / diffgrid_mass
    diffgrid *= scaling_factor

    return diffgrid


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
