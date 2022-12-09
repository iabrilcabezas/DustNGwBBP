'''
bbpw_script
    writes config .yml file required by BBPower
'''

import yaml
from utils.params import PATH_DICT, NAME_RUN


def write_config_yml_script(type_cov, nside_bbpw):

    '''
    Writes .yml file used by BBPower in BBCompSep (no plotting!)

    * Parameters *
    type_cov: 'w' or 'wt'
        type of covariance
    nside_bbpw: int
        map resolution that bbpower uses

    '''

    name_inputs = '_'.join([NAME_RUN, 'all', 'Cl', type_cov])
    path_inputs = PATH_DICT['output_path'] + name_inputs
    path_config = PATH_DICT['output_path'] + 'config_files/' + name_inputs + '.yml'

    dict_config = {'global': {'nside': nside_bbpw, 'compute_dell': False},
                   'modules': 'bbpower',
                   'launcher': 'local',
                   'stages': [{'name': 'BBCompSep', 'nprocess': 1}], #{'name': 'BBPlotter', 'nprocess': 1}],
                   'inputs': {'cells_coadded': path_inputs + '_tot.fits',
                              'cells_fiducial': path_inputs + '_fid.fits',
                              'cells_noise': path_inputs + '_noi.fits'},
                   'config': path_config,
                   'resume': False,
                   'output_dir': PATH_DICT['output_path'] + 'outputs/',
                   'log_dir': PATH_DICT['output_path'] + 'outputs/',
                   'pipeline_log': PATH_DICT['output_path'] + 'outputs/log.txt',
                   'BBCompSep': {'sampler': 'emcee',
                                 'nwalkers': 24,
                                 'n_iters': 1000,
                                 'likelihood_type': 'h&l',
                                 'pol_channels': ['B'],
                                 'l_min': 30,
                                 'l_max': 120,
                                 'cmb_model': {'cmb_templates': ['/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower/examples/data/camb_lens_nobb.dat',
                                                                 '/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower/examples/data/camb_lens_r1.dat'],
                                               'params': {'r_tensor': ['r_tensor', 'tophat', [-0.1, 0.0, 0.1]],
                                                          'A_lens': ['A_lens', 'tophat', [0.0, 1.0, 2.0]]}},
                                 'fg_model': {'component_1': {'name': 'Dust',
                                                              'sed': 'Dust',
                                                              'cl': {'BB': 'ClPowerLaw'},
                                                              'sed_parameters': {'beta_d': ['beta_d', 'Gaussian', [1.59, 0.11]],
                                                                                 'temp_d': ['temp', 'fixed', [19.6]],
                                                                                 'nu0_d': ['nu0', 'fixed', [353.0]]},
                                                              'cl_parameters': {'BB': {'amp_d_bb': ['amp', 'tophat', [0.0, 5.0, 'inf']],
                                                                                       'alpha_d_bb': ['alpha', 'tophat', [-1.0, -0.2, 0.0]],
                                                                                       'l0_d_bb': ['ell0', 'fixed', [80.0]]}}}}},
                #    'BBPlotter': {'lmax_plot': 120,
                #                  'plot_coadded_total': False,
                #                  'plot_noise': False,
                #                  'plot_nulls': False,
                #                  'plot_likelihood': True}
                                 }

    with open(path_config, 'w') as file:
        yaml.dump(dict_config, file, default_flow_style=True)
