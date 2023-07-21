'''
the clean2cmb you need

reads best fit from BBCompSep and writes .yml with those params
'''

# the covariance stays the same in all computations
# read in best fit from chi2.npz file
# call predicted spectral to create cellmodel
# clean cell model best fit
# clean sampled cell
# compute chi2 from the two 

import numpy as np
import yaml
from utils.params import PATH_DICT


def write_cell_model(chi2_path):

    '''
    creates dictionary from chi2.npz file (params, best fit values)
    '''

    chi2_file = np.load(chi2_path + 'chi2.npz')
    dict_chi2 = dict(zip(chi2_file['names'], chi2_file['params']))

    return dict_chi2


def dictmodel_fromparams(chi2_path):

    '''
    writes down .yml file to call bbcompsep and generate cell_model with best fit params
    (contained in chi2_path)
    '''

    dict_chi2 = write_cell_model(chi2_path)

    dict_cmbmodel = {'cmb_templates':[PATH_DICT['bbpipe_path'] + 'examples/data/camb_lens_nobb.dat',
                                   PATH_DICT['bbpipe_path'] + 'examples/data/camb_lens_r1.dat'],
                'params': {'r_tensor': ['r_tensor', 'tophat', [-1, float(dict_chi2['r_tensor']), 1]],
                            'A_lens': ['A_lens', 'tophat', [0.0, float(dict_chi2['A_lens']), 1.0, 2.0]]}}

    dict_comp1 = {'name': 'Dust',
                    'sed': 'Dust',
                    'cl': {'BB': 'ClPowerLaw'},
                    'sed_parameters': {'beta_d': ['beta_d', 'Gaussian', [ float(dict_chi2['beta_d']), 0.5]],
                                        'temp_d': ['temp', 'fixed', [19.6]],
                                        'nu0_d': ['nu0', 'fixed', [353.0]]},
                    'cl_parameters': {'BB': {'amp_d_bb': ['amp', 'tophat', [0.0, float(dict_chi2['amp_d_bb']), 'inf']],
                                            'alpha_d_bb': ['alpha', 'tophat', [-1.0, float(dict_chi2['alpha_d_bb']), 0.0]],
                                            'l0_d_bb': ['ell0', 'fixed', [80.0]]}},
                    'decorr': {'decorr_amp_s': ['decorr_amp', 'fixed', [1.0]],
                               'decorr_nu01_s': ['decorr_nu01', 'fixed', [353.0]],
                               'decorr_nu02_s': ['decorr_nu02', 'fixed', [217.0]]}
                }

    dict_comp2 = {'name': 'Synchrotron',
                    'sed': 'Synchrotron',
                    'cl': { 'BB': 'ClPowerLaw'},
                    'sed_parameters': { 'beta_s': ['beta_pl', 'Gaussian', [float(dict_chi2['beta_s']), 0.6]],
                                        'nu0_s' : ['nu0', 'fixed', [23.]]},
                    'cl_parameters': {'BB': {   'amp_s_bb': ['amp', 'tophat', [0., float(dict_chi2['amp_s_bb']), "inf"]],
                                                'alpha_s_bb': ['alpha', 'tophat', [-1., float(dict_chi2['alpha_s_bb']), 0.]],
                                                'l0_s_bb': ['ell0', 'fixed', [80.]] }},
                    'decorr': {'decorr_amp_s': ['decorr_amp', 'fixed', [1.0]],
                               'decorr_nu01_s': ['decorr_nu01', 'fixed', [40.0]],
                               'decorr_nu02_s': ['decorr_nu02', 'fixed', [23.0]]}
    }

    dict_comp1['cross'] = { 'epsilon_ds': ['component_2', 'tophat', [-1., float(dict_chi2['epsilon_ds']), 1.]]}

    dict_fgmodel = {'component_1': dict_comp1, 'component_2': dict_comp2}

    dict_config = {'global': {'nside': 256, 'compute_dell': False},
                   'modules': 'bbpower',
                   'launcher': 'local',
                   'stages': [{'name': 'BBCompSep', 'nprocess': 1}],
                   'config': '/global/common/software/act/iabril/python/DustNGwBBP/decorr/bbcompsep_decorr.yml',
                   'resume': False,
                   'BBCompSep': {'sampler': 'predicted_spectra', 'predict_at_minimum':False, 'predict_to_sacc': True,
                                 'likelihood_type': 'chi2', 'pol_channels': ['B'], 'l_min': 30, 'l_max': 300,
                                 'cmb_model': dict_cmbmodel, 'fg_model': dict_fgmodel}
    }
    with open(chi2_path + 'bestfitcell.yml', 'w', encoding = 'utf-8') as file:
        yaml.dump(dict_config, file, default_flow_style=True)


