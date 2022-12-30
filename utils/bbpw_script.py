'''
bbpw_script
    writes config .yml file required by BBPower
'''

import yaml
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP, POLARIZATION_cov
from utils.sed import get_band_names
from utils.params_compsep import get_configcompsep

band_names = get_band_names()

assert POLARIZATION_cov == 'B', 'script only for B polarization'

def write_config_yml_script(type_cov, nside_bbpw, dict_compsep,
                             cross = False, moments = False, n_iters = int(1e3)):

    '''
    Writes .yml file used by BBPower in BBCompSep (no plotting!)
    Only interested in B- polarization channel

    * Parameters *
    type_cov: 'w' or 'wt'
        type of covariance
    nside_bbpw: int
        map resolution that bbpower uses
    dict_compsep: dict
        lmin : 30
        lmax : 300,
        bands: 'all'
    '''

    lmin_bbcomp, lmax_bbcomp, bands_bbcomp, name_configcompsep = get_configcompsep(dict_compsep)

    name_inputs = '_'.join([NAME_RUN, 'Cl'])
    name_config = '_'.join([NAME_RUN, 'Cl'])

    if cross:
        assert NAME_COMP == 'dcs', 'cross between 1 and 2 specified but there is no 2'
        name_config += '_C'
    else:
        name_config += '_0'

    dict_comp1 =   {'name': 'Dust',
                    'sed': 'Dust',
                    'cl': {'BB': 'ClPowerLaw'},
                    'sed_parameters': {'beta_d': ['beta_d', 'Gaussian', [1.59, 0.11]],
                                        'temp_d': ['temp', 'fixed', [19.6]],
                                        'nu0_d': ['nu0', 'fixed', [353.0]]},
                    'cl_parameters': {'BB': {'amp_d_bb': ['amp', 'tophat', [0.0, 5.0, 'inf']],
                                            'alpha_d_bb': ['alpha', 'tophat', [-1.0, -0.2, 0.0]],
                                            'l0_d_bb': ['ell0', 'fixed', [80.0]]}}}

    if cross:
        dict_comp1['cross'] = { 'epsilon_ds': ['component_2', 'tophat', [-1., 0., 1.]]}

    dict_comp2 = {  'name': 'Synchrotron',
                    'sed': 'Synchrotron',
                    'cl': { 'BB': 'ClPowerLaw'},
                    'sed_parameters': { 'beta_s': ['beta_pl', 'Gaussian', [-3.0, 0.3]],
                                        'nu0_s' : ['nu0', 'fixed', [23.]]},
                    'cl_parameters': {'BB': {   'amp_s_bb': ['amp', 'tophat', [0., 2., "inf"]],
                                                'alpha_s_bb': ['alpha', 'tophat', [-1., -0.4, 0.]],
                                                'l0_s_bb': ['ell0', 'fixed', [80.]] }}}

    if moments:
        dict_comp1['moments'] = {   'amp_d_beta': ['amp_beta', 'tophat', [-10., 0., 10.]],
                                    'gamma_d_beta': ['gamma_beta', 'tophat', [-6., -3., -2.]]}
        dict_comp2['moments'] = {   'amp_s_beta': ['amp_beta', 'tophat', [-10., 0., 10.]],
                                    'gamma_s_beta': ['gamma_beta', 'tophat', [-6., -3., -2.]]}

    dict_fgmodel = {'component_1' : dict_comp1}

    if moments:
        dict_fgmodel['use_moments'] = True
        dict_fgmodel['moments_lmax'] = 300
        name_config += '_M'
    else:
        name_config += '_0'

    if NAME_COMP == 'dcs':
        dict_fgmodel['component_2'] = dict_comp2

    nwalk = 24
    if (moments & cross):
        nwalk = 36

    name_inputs += '_' + type_cov
    name_config += '_' + type_cov
    path_inputs = PATH_DICT['output_path'] +  name_inputs

    name_config = name_config + '_' + name_configcompsep
    path_config = PATH_DICT['output_path'] + 'config_files/' + name_config + '.yml'

    dict_config = {'global': {'nside': nside_bbpw, 'compute_dell': False},
                   'modules': 'bbpower',
                   'launcher': 'local',
                   'stages': [{'name': 'BBCompSep', 'nprocess': 1}], #{'name': 'BBPlotter', 'nprocess': 1}],
                   'inputs': {'cells_coadded': path_inputs + '_tot.fits',
                              'cells_fiducial': path_inputs + '_fid.fits',
                              'cells_noise' : path_inputs  + '_noi.fits'},
                   'config': path_config,
                   'resume': False,
                   'output_dir': PATH_DICT['output_path'] + 'outputs/',
                   'log_dir': PATH_DICT['output_path'] + 'outputs/',
                   'pipeline_log': PATH_DICT['output_path'] + 'outputs/log.txt',
                   'BBCompSep': {'sampler': 'emcee',
                                 'nwalkers': nwalk,
                                 'n_iters': n_iters,
                                 'likelihood_type': 'h&l',
                                 'pol_channels': ['B'],
                                 'l_min': lmin_bbcomp,
                                 'l_max': lmax_bbcomp,
                                 'bands': bands_bbcomp,
                                 'cmb_model': {'cmb_templates': ['/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower/examples/data/camb_lens_nobb.dat',
                                                                 '/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower/examples/data/camb_lens_r1.dat'],
                                               'params': {'r_tensor': ['r_tensor', 'tophat', [-0.1, 0.0, 0.1]],
                                                          'A_lens': ['A_lens', 'tophat', [0.0, 1.0, 2.0]]}},
                                 'fg_model': dict_fgmodel}
                                 }

    with open(path_config, 'w', encoding = 'utf-8') as file:
        yaml.dump(dict_config, file, default_flow_style=True)
