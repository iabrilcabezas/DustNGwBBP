'''
dict_bbpw
    stores and builds default dictionaries to write BBCompSep config file
'''
from copy import deepcopy
from utils.params import PATH_DICT

DICT_COMP1 =   {'name': 'Dust',
                    'sed': 'Dust',
                    'cl': {'BB': 'ClPowerLaw'},
                    'sed_parameters': {'beta_d': ['beta_d', 'Gaussian', [1.59, 0.11]],
                                        'temp_d': ['temp', 'fixed', [19.6]],
                                        'nu0_d': ['nu0', 'fixed', [353.0]]},
                    'cl_parameters': {'BB': {'amp_d_bb': ['amp', 'tophat', [0.0, 5.0, 'inf']],
                                            'alpha_d_bb': ['alpha', 'tophat', [-1.0, -0.2, 0.0]],
                                            'l0_d_bb': ['ell0', 'fixed', [80.0]]}}}
DICT_COMP2 = {  'name': 'Synchrotron',
                    'sed': 'Synchrotron',
                    'cl': { 'BB': 'ClPowerLaw'},
                    'sed_parameters': { 'beta_s': ['beta_pl', 'Gaussian', [-3.0, 0.3]],
                                        'nu0_s' : ['nu0', 'fixed', [23.]]},
                    'cl_parameters': {'BB': {   'amp_s_bb': ['amp', 'tophat', [0., 2., "inf"]],
                                                'alpha_s_bb': ['alpha', 'tophat', [-1., -0.4, 0.]],
                                                'l0_s_bb': ['ell0', 'fixed', [80.]] }}}
DICT_COMP1_MOMENTS = {   'amp_d_beta': ['amp_beta', 'tophat', [-10., 0., 10.]],
                                    'gamma_d_beta': ['gamma_beta', 'tophat', [-6., -3., -2.]]}
DICT_COMP2_MOMENTS = {   'amp_s_beta': ['amp_beta', 'tophat', [-10., 0., 10.]],
                                    'gamma_s_beta': ['gamma_beta', 'tophat', [-6., -3., -2.]]}
DICT_COMP1_CROSS = { 'epsilon_ds': ['component_2', 'tophat', [-1., 0., 1.]]}

DICT_CMBMODEL = {'cmb_templates': [PATH_DICT['bbpipe_path'] + 'examples/data/camb_lens_nobb.dat',
                                   PATH_DICT['bbpipe_path'] + 'examples/data/camb_lens_r1.dat'],
                'params': {'r_tensor': ['r_tensor', 'tophat', [-0.1, 0.0, 0.1]],
                            'A_lens': ['A_lens', 'tophat', [0.0, 1.0, 2.0]]}}

DICT_COMP1_DECORR = { 'decorr_amp_d' : ['decorr_amp', 'tophat', [0.9, 1.0, 1.1]],
                      'decorr_nu0_d' : ['decorr_nu0', 'fixed', [353.0]]}
DICT_COMP2_DECORR = { 'decorr_amp_s' : ['decorr_amp', 'tophat', [0.9, 1.0, 1.1]],
                      'decorr_nu0_s' : ['decorr_nu0', 'fixed', [23.0]]}


def get_dict_fgmodel(name_comp, cro, mom, decorr):
    '''
    Returns dictionary of foreground model used by BBCompSep according to arguments

    ** Parameters **
    name_comp: str ('dc0' or 'dcs')
        Components (dust, or dust + synchrotron) of model
    cr: (cross) bool
        Include cross term between dust and synchrotron
    mom: (moments) bool
        Use component separation method
    '''

    dict_comp1 = deepcopy(DICT_COMP1)
    dict_comp2 = deepcopy(DICT_COMP2)

    if cro:
        dict_comp1['cross'] = DICT_COMP1_CROSS

    if mom:
        dict_comp1['moments'] = DICT_COMP1_MOMENTS
        dict_comp2['moments'] = DICT_COMP2_MOMENTS

    if decorr:
        dict_comp1['decorr']  = DICT_COMP1_DECORR
        dict_comp2['decorr']  = DICT_COMP2_DECORR

    dict_fgmodel = {'component_1' : dict_comp1}

    if mom:
        dict_fgmodel['use_moments'] = True
        dict_fgmodel['moments_lmax'] = 299 # 300 gives a bug if lmax of data is 300

    if name_comp == 'dcs':
        dict_fgmodel['component_2'] = dict_comp2

    return dict_fgmodel
