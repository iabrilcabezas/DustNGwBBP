from utils.params import PATH_DICT

dict_comp1 =   {'name': 'Dust',
                    'sed': 'Dust',
                    'cl': {'BB': 'ClPowerLaw'},
                    'sed_parameters': {'beta_d': ['beta_d', 'Gaussian', [1.59, 0.11]],
                                        'temp_d': ['temp', 'fixed', [19.6]],
                                        'nu0_d': ['nu0', 'fixed', [353.0]]},
                    'cl_parameters': {'BB': {'amp_d_bb': ['amp', 'tophat', [0.0, 5.0, 'inf']],
                                            'alpha_d_bb': ['alpha', 'tophat', [-1.0, -0.2, 0.0]],
                                            'l0_d_bb': ['ell0', 'fixed', [80.0]]}}}
dict_comp2 = {  'name': 'Synchrotron',
                    'sed': 'Synchrotron',
                    'cl': { 'BB': 'ClPowerLaw'},
                    'sed_parameters': { 'beta_s': ['beta_pl', 'Gaussian', [-3.0, 0.3]],
                                        'nu0_s' : ['nu0', 'fixed', [23.]]},
                    'cl_parameters': {'BB': {   'amp_s_bb': ['amp', 'tophat', [0., 2., "inf"]],
                                                'alpha_s_bb': ['alpha', 'tophat', [-1., -0.4, 0.]],
                                                'l0_s_bb': ['ell0', 'fixed', [80.]] }}}
dict_comp1_moments = {   'amp_d_beta': ['amp_beta', 'tophat', [-10., 0., 10.]],
                                    'gamma_d_beta': ['gamma_beta', 'tophat', [-6., -3., -2.]]}
dict_comp2_moments = {   'amp_s_beta': ['amp_beta', 'tophat', [-10., 0., 10.]],
                                    'gamma_s_beta': ['gamma_beta', 'tophat', [-6., -3., -2.]]}
dict_comp1_cross = { 'epsilon_ds': ['component_2', 'tophat', [-1., 0., 1.]]}

dict_cmbmodel = {'cmb_templates': [PATH_DICT['bbpipe_path'] + 'examples/data/camb_lens_nobb.dat',
                                   PATH_DICT['bbpipe_path'] + 'examples/data/camb_lens_r1.dat'],
                'params': {'r_tensor': ['r_tensor', 'tophat', [-0.1, 0.0, 0.1]],
                            'A_lens': ['A_lens', 'tophat', [0.0, 1.0, 2.0]]}}


def get_dict_fgmodel(name_comp, cross = False, moments = False):

    '''
    
    '''

    if cross:
        dict_comp1['cross'] = dict_comp1_cross

    if moments:
        dict_comp1['moments'] = dict_comp1_moments
        dict_comp2['moments'] = dict_comp2_moments

    dict_fgmodel = {'component_1' : dict_comp1}

    if moments:
        dict_fgmodel['use_moments'] = True
        dict_fgmodel['moments_lmax'] = 300

    if name_comp == 'dcs':
        dict_fgmodel['component_2'] = dict_comp2

    return dict_fgmodel
