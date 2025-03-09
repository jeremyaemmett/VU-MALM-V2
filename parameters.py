""" User-specified parameters """

from datetime import datetime

# 3/2/2025

### Main directory
main_directory = 'C:/Users/Jeremy/Desktop/VUMALM/'  # Main
###

# Data directory
site = 'palsa_low'  # Site
###

### Forcing file
forcing_file = 'palsa_low_Drying0p0025_z00p1_rootDepth0p1'  # Forcing file
###

### Main process switches
microbe_flag = True  # microbe chemistry
diffusion_flag = True  # vertical gas transport by diffusion
plant_flag = True  # vertical gas transport by plants
###

### Secondary process switches
instant_diffusion = False
interface_flux = True
###

### Operational switches
jump2plots = False
gaussian_profile = True
###

### Grid parameters
dz0 = 0.02  # initial layer thickness (m)
growth_rate = 1.05  # factor by which layer thickness grows exponentially ()
total_depth = 1.0  # total depth of the layer column (m)
dt = 0.05  # main time-step (days)
diff_n_dt = 10  # number of diffusion sub-time-steps within dt
diff_n_dz = 100  # number of diffusion layers (uniform thickness) to total_depth
write_dt = 0.25  # output write interval in days
start_datetime = datetime(2022, 1, 1, 0, 0, 0)  # Starting datetime
end_datetime = datetime(2022, 12, 31, 0, 0, 0)  # Ending datetime
###

### Model parameters
gc_per_molc = 12.01  #
org_dens = 175.0  # kg/m3 # Organic density
k_cpool = 0.02  # DOC:C Ratio
k_ch4_prod = 0.5  # CH4 prod ratio
doc_prod_q10 = 2.5  # Temp sens. of DOC prod - good oxygen scavenging < 0.2; 1-5 allowed
f_d_w = 0.8  # Reduction factor for diff in water
f_d_a = 0.8  # Reduction factor for diff in air
tau = 1.5  # Root tort.
a_m_a = 0.085  # Root ending area per dry biomass
s_l_a = 20.0  # Specific leaf area
sigma = 0.8  # Peat porosity
k = 0.0001  # Time constant
t_min1 = 0.0  # 0.0 Min soil temp for hom.ace.genesis
t_max1 = 20.0  # Maximum...
t_opt1 = 10.0  # Optimum...
t_min2 = 0.0  # 20.0 Min temp for H2 m.genesis
t_max2 = 50.0  # Maximum...
t_opt2 = 35.0  # Optimum...
ft_max_temp = 30.0  # Max temp for ace m.genesis, m.trophy, and DOC decomp
ph_min = 3.0  # Min pH for pH factor
ph_max = 9.0  # Maximum...
ph_opt = 6.2  # Optimum...
###

### Chemical reactions

mcbs_dict = None
subs_list, cons_dict, prod_dict = ['doc'],       {},                      {'ace': 1.0, 'co2': 0.5, 'h2': 1.0 / 6.0}
vmax_dict = {'lab': 'v_doc_prod_ace_max',                'val': 0.5,          'rand': [0.3, 0.7, 0]}
q10_dict =  {'lab': 'ace_prod_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_doc_prod_ace'],                  'val': [10.0],       'rand': [[5, 15, 0]]}
ansp_dict = {'name': 'ansp', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

mcbs_dict = None
subs_list, cons_dict, prod_dict = ['doc', 'o2'], {'doc': 1.0, 'o2': 1.0}, {'ace': 1.0, 'co2': 0.5}
vmax_dict = {'lab': 'v_doc_prod_ace_max',                'val': 0.5,          'rand': [0.3, 0.7, 0]}
q10_dict =  {'lab': 'ace_prod_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_doc_prod_ace', 'k_ace_prod_o2'], 'val': [10.0, 0.04], 'rand': [[5, 15, 0], [0.01, 0.1, 0]]}
aesp_dict = {'name': 'aesp', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

mcbs_dict = {'lab': 'homa', 'init': 1.0, 'grow': 0.2, 'mort': 0.06, 'gcpc': 1e-13}
subs_list, cons_dict, prod_dict = ['co2', 'h2'], {'co2': 2.0, 'h2': 4.0}, {'ace': 1.0}
vmax_dict = {'lab': 'v_h2_prod_ace_max',                 'val': 0.15,         'rand': [0.01, 0.3, 0]}
q10_dict =  {'lab': 'ace_prod_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_co2_prod_ace', 'k_h2_prod_ace'], 'val': [0.05, 0.01], 'rand': [[0.01, 0.1, 0], [0.01, 0.1, 0]]}
hoag_dict = {'name': 'hoag', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

mcbs_dict = {'lab': 'h2gn', 'init': 1.0, 'grow': 0.2, 'mort': 0.06, 'gcpc': 1e-13}
subs_list, cons_dict, prod_dict = ['co2', 'h2'], {'co2': 1.0, 'h2': 4.0}, {'ch4': 1.0}
vmax_dict = {'lab': 'v_h2_prod_ch4_max',                 'val': 0.15,         'rand': [0.01, 0.3, 0]}
q10_dict =  {'lab': 'ch4_prod_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_co2_prod_ch4', 'k_h2_prod_ch4'], 'val': [0.05, 0.01], 'rand': [[0.01, 0.1, 0], [0.01, 0.1, 0]]}
h2mg_dict = {'name': 'h2mg', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

mcbs_dict = {'lab': 'acgn', 'init': 1.0, 'grow': 0.3, 'mort': 0.06, 'gcpc': 1e-13}
subs_list, cons_dict, prod_dict = ['ace'],       {'ace': 1.0},            {'ch4': 1.0, 'co2': 1.0}
vmax_dict = {'lab': 'v_ace_cons_max',                    'val': 0.5,          'rand': [0.3, 0.7, 0]}
q10_dict =  {'lab': 'ch4_prod_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_ace_prod_ch4'],                  'val': [0.05],       'rand': [[0.01, 0.1, 0]]}
acmg_dict = {'name': 'acmg', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

mcbs_dict = {'lab': 'mtph', 'init': 1.0, 'grow': 0.4, 'mort': 0.06, 'gcpc': 1e-13}
subs_list, cons_dict, prod_dict = ['ch4', 'o2'], {'ch4': 1.0, 'o2': 2.0}, {'co2': 1.0}
vmax_dict = {'lab': 'v_ch4_oxid_max',                    'val': 0.5,          'rand': [0.3, 0.7, 0]}
q10_dict =  {'lab': 'ch4_oxid_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_ch4_oxid_ch4', 'k_ch4_oxid_o2'], 'val': [0.05, 0.02], 'rand': [[0.01, 0.1, 0], [0.01, 0.1, 0]]}
mtox_dict = {'name': 'mtox', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

mcbs_dict = None
subs_list, cons_dict, prod_dict = ['doc', 'o2'], {'doc': 1.0, 'o2': 2.0}, {}
vmax_dict = {'lab': 'k_aer',                             'val': 0.2,          'rand': [0.3, 0.7, 0]}
q10_dict =  {'lab': 'doc_prod_q10',                      'val': 2.5,          'rand': [1.0, 5.0, 0]}
k_dict =    {'lab': ['k_aer_doc, k_aer_o2'],             'val': [10.0, 0.22], 'rand': [[5.0, 15.0, 0], [0.01, 0.3, 0]]}
aero_dict = {'name': 'aero', 'vmax': vmax_dict, 'q10': q10_dict, 'k': k_dict, 'subs': subs_list,
             'cons': cons_dict, 'prod': prod_dict, 'mcbs': mcbs_dict}

reactions = {'ansp': ansp_dict, 'aesp': aesp_dict, 'hoag': hoag_dict, 'h2mg': h2mg_dict, 'acmg': acmg_dict,
             'mtox': mtox_dict, 'aero': aero_dict}
###

### Chemical properties
cprops = {'ch4': {'H1': 1.3e-3, 'H2': 1700.0, 'S1': 1898.0, 'S2': 110.10, 'S3': 2.834000, 'S4': 0.027910,
                  'DA1': 1.90e-5, 'DA2': 1.820, 'DW1': 1.50e-9, 'airvol': 0.000179, 'colmap': 'Purples'},
          'o2':  {'H1': 1.3e-3, 'H2': 1500.0, 'S1': 1800.6, 'S2': 120.10, 'S3': 3.781800, 'S4': 0.047608,
                  'DA1': 1.80e-5, 'DA2': 1.820, 'DW1': 2.40e-9, 'airvol': 20.94600, 'colmap': 'Greens'},
          'co2': {'H1': 3.4e-2, 'H2': 2400.0, 'S1': 1911.0, 'S2': 113.70, 'S3': 2.967000, 'S4': 0.029430,
                  'DA1': 1.47e-5, 'DA2': 1.792, 'DW1': 1.81e-6, 'airvol': 0.041200, 'colmap': 'Blues'},
          'h2':  {'H1': 7.8e-4, 'H2': 530.00, 'S1': 629.95, 'S2': 34.691, 'S3': 0.868100, 'S4': 0.008400,
                  'DA1': 6.68e-5, 'DA2': 1.820, 'DW1': 5.11e-9, 'airvol': 0.000050, 'colmap': 'Oranges'},
          'ace': {'H1': 0.0000, 'H2': 0.0000, 'S1': 0.0000, 'S2': 0.0000, 'S3': 0.000000, 'S4': 0.000000,
                  'DA1': 0.00000, 'DA2': 0.000, 'DW1': 0.00000, 'airvol': 0.000000, 'colmap': 'Reds'},
          'c':   {'H1': 0.0000, 'H2': 0.0000, 'S1': 0.0000, 'S2': 0.0000, 'S3': 0.000000, 'S4': 0.000000,
                  'DA1': 0.00000, 'DA2': 0.00000, 'DW1': 0.00000, 'airvol': 0.000000, 'colmap': 'Wistia'},
          'doc': {'H1': 0.0000, 'H2': 0.0000, 'S1': 0.0000, 'S2': 0.0000, 'S3': 0.000000, 'S4': 0.000000,
                  'DA1': 0.00000, 'DA2': 0.00000, 'DW1': 0.00000, 'airvol': 0.000000, 'colmap': 'Wistia'}}
