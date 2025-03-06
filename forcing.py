""" Forcing data processing for VU-MALM ingestion """

from scipy.ndimage import gaussian_filter
import parameters as params
import conversions as conv
import numpy as np
import properties
import initialize
import system
import csv
import os

# 3/2/2025


def list_directory_files(directory, extension):
    for file in os.listdir(directory):
        if file.endswith(extension):
            print(os.path.join(directory, file))


def read_csv_array(filename):
    with open(filename, newline='') as csvfile:
        data = list(csv.reader(csvfile))

    return data


def forcing_data():

    # Forcing data [cg_depths, cg_times]

    forcing_file_directory = params.main_directory + 'data/' + params.site + '/forcing/' + params.forcing_file + '/'

    tem_data = np.array(read_csv_array(forcing_file_directory + 'T.txt'), dtype='float64')
    wat_data = 100.0 * np.array(read_csv_array(forcing_file_directory + 'water.txt'), dtype='float64')
    ice_data = 100.0 * np.array(read_csv_array(forcing_file_directory + 'ice.txt'), dtype='float64')
    org_data = 100.0 * np.array(read_csv_array(forcing_file_directory + 'soilOrganic.txt'), dtype='float64')
    min_data = 100.0 * np.array(read_csv_array(forcing_file_directory + 'soilMineral.txt'), dtype='float64')
    air_data = 100.0 - (wat_data + ice_data) - org_data - min_data
    tim_data = (np.array(read_csv_array(forcing_file_directory + 'TIMESTAMP.txt'), dtype='float64'))[0]
    dep_data = np.array(read_csv_array(forcing_file_directory + 'z.txt'), dtype='float64')

    # Re-shaping and normalization
    org_data = np.full(np.shape(tem_data), org_data)
    min_data = np.full(np.shape(tem_data), min_data)
    air_data = np.full(np.shape(tem_data), air_data)
    tim_data = tim_data - tim_data[0]
    dep_data = dep_data - dep_data[0]
    dep_data = ((np.transpose(dep_data))[0])[0:-1]

    # Interpolate the forcing data to the microbe model's depth grid [n_layers x n_forcing_time-steps]

    # Number of forcing time records
    nt_fc = np.shape(tem_data)[1]

    thicks, depths = initialize.define_layers()

    tem_data_interp = np.zeros((len(depths), nt_fc))
    wat_data_interp = np.zeros((len(depths), nt_fc))
    ice_data_interp = np.zeros((len(depths), nt_fc))
    org_data_interp = np.zeros((len(depths), nt_fc))
    min_data_interp = np.zeros((len(depths), nt_fc))
    air_data_interp = np.zeros((len(depths), nt_fc))

    # Number of forcing depth layers needed to encompass the model's layer depths, for interpolation
    ir = len(dep_data[dep_data >= np.min(-depths)]) + 1

    # Make quantity profiles at each time-step, interpolated to the model's depth grid
    for i in range(0, nt_fc):
        xnew, tem_data_interp[:, i] = conv.interp2grid(-dep_data[0:ir], tem_data[0:ir, i], depths)
        xnew, wat_data_interp[:, i] = conv.interp2grid(-dep_data[0:ir], wat_data[0:ir, i], depths)
        xnew, ice_data_interp[:, i] = conv.interp2grid(-dep_data[0:ir], ice_data[0:ir, i], depths)
        xnew, org_data_interp[:, i] = conv.interp2grid(-dep_data[0:ir], org_data[0:ir, i], depths)
        xnew, min_data_interp[:, i] = conv.interp2grid(-dep_data[0:ir], min_data[0:ir, i], depths)
        xnew, air_data_interp[:, i] = conv.interp2grid(-dep_data[0:ir], air_data[0:ir, i], depths)

    # Derived predetermined fields

    # Model depths, thicknesses
    dep_data, thk_data = conv.temporally_extrude_profile(depths, nt_fc), conv.temporally_extrude_profile(thicks, nt_fc)

    # Column-integrated organic carbon concentrations
    n_t = np.shape(org_data_interp)[1]
    org_kgm3 = np.transpose(np.array([conv.volpct2kgm3(org_data_interp[:, n], params.org_dens) for n in range(0, n_t)]))
    org_cpls = np.array([conv.unitm3s2unitm2(org_kgm3[:, n], thk_data[:, n])         for n in range(0, n_t)])
    cpl_data = conv.vertically_extrude_profile(org_cpls, len(depths))

    # Saturation %'s
    sat_data_interp = 100.0 * (wat_data_interp) / (wat_data_interp + air_data_interp)
    sat_data_interp[sat_data_interp > 90.0] = 100.0
    # Remove any subsaturated layer artifacts below the water table
    for i in range(0, np.shape(sat_data_interp)[1]):
        sat_ids = np.where(sat_data_interp[:, i] == 100.0)[0]
        discontinuity_exists = np.max(np.diff(sat_ids)) > 1
        if discontinuity_exists: sat_data_interp[sat_ids[0]:sat_ids[-1], i] = 100.0

    # Temperature factors 1 and 2
    f_t1, f_t2 = properties.temperature_factors_1_and_2(tem_data_interp)

    # q10 Temperature factors
    ftqs = properties.q10_temperature_factors(tem_data_interp)

    # Dissolved organic carbon concentrations
    doc_data = properties.dissolved_organic_carbon(cpl_data, tem_data_interp, sat_data_interp, thk_data)

    # Effective diffusivities
    eff_diffs, eff_diffs_full = properties.effective_diffusivity(sat_data_interp, tem_data_interp)

    # Equilibrium concentrations
    equ_concs = properties.equilibrium_concentration(tem_data_interp)

    # Schmidt numbers
    schm_nums = properties.schmidt_number(tem_data_interp)

    # Transfer velocities
    tran_vels = properties.transfer_velocity(schm_nums)

    # Root densities
    days, profiles = system.datedprofiles2arrays('roo', depths)
    roo_densy = conv.datedprofiles2grid(days, profiles, nt_fc)

    # Root xsectional areas
    roo_xsect = 0.05 * roo_densy

    forcing_data_dictionary = {'tem': tem_data_interp, 'ft1': f_t1, 'ft2': f_t2, 'wat': wat_data_interp,
                               'ice': ice_data_interp, 'org': org_data_interp, 'min': min_data_interp,
                               'air': air_data_interp, 'sat': sat_data_interp, 'doc': doc_data, 'roo': roo_densy,
                               'rxa': roo_xsect, 'tim': tim_data, 'dep': depths, 'efd': eff_diffs,
                               'efdf': eff_diffs_full, 'eqc': equ_concs, 'ftq': ftqs, 'schm': schm_nums,
                               'tvel': tran_vels}

    return forcing_data_dictionary


def forcing_data_t(fdd, t):

    profiles = {}

    # Closest time index in the forcing data dictionary to the model datetime
    closest_index = np.where(abs(fdd['tim'] - t) == min(abs(fdd['tim'] - t)))[0][0]

    for key in fdd:
        key_dimensions = len(np.shape(fdd[key]))
        if key_dimensions == 0:
            profiles[key] = {}
            for sub_key in list((fdd[key]).keys()):
                profiles[key][sub_key] = fdd[key][sub_key][:, closest_index]
        if key_dimensions == 2:  # Depth x Time data
            profiles[key] = fdd[key][:, closest_index]
        if key_dimensions == 3:  # Type x Depth x Time data
            profiles[key] = fdd[key][:, :, closest_index]

    profiles['tim'] = fdd['tim'][closest_index]
    profiles['dep'] = fdd['dep']

    return profiles
