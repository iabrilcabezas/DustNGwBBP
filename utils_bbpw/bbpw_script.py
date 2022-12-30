'''
bbpw_script
    writes config .yml file required by BBPower
'''

import yaml
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP, POLARIZATION_cov
from utils.sed import get_band_names
from utils_bbpw.params import get_configcompsep
from utils_bbpw.dict_bbpw import get_dict_fgmodel, dict_cmbmodel

band_names = get_band_names()

assert POLARIZATION_cov == 'B', 'script only for B polarization'

def write_config_yml_script(type_cov, params_bbpw, dict_compsep,
                            cross = False, moments = False):

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
    params_bbpw: dict
        nside_bbpw: 256
        niters: int(1e3)
        nwalk: 24
    '''

    dict_bbcomp = dict(get_configcompsep(dict_compsep))

    name_inputs = '_'.join([NAME_RUN, 'Cl'])
    name_config = '_'.join([NAME_RUN, 'Cl'])

    if cross:
        assert NAME_COMP == 'dcs', 'cross between 1 and 2 specified but there is no 2'
        name_config += '_C'
    else:
        name_config += '_0'

    if moments:
        name_config += '_M'
    else:
        name_config += '_0'

    if moments & cross:
        params_bbpw['nwalk'] = 36

    name_inputs += '_' + type_cov
    path_inputs = PATH_DICT['output_path'] +  name_inputs

    name_config += '_' + type_cov + '_' + dict_bbcomp['name_config']
    path_config = PATH_DICT['output_path'] + 'config_files/' + name_config + '.yml'

    dict_fgmodel = get_dict_fgmodel(NAME_COMP, cross, moments)

    dict_config = {'global': {'nside': params_bbpw['nside'], 'compute_dell': False},
                   'modules': 'bbpower',
                   'launcher': 'local',
                   'stages': [{'name': 'BBCompSep', 'nprocess': 1}],
                   'inputs': {'cells_coadded': path_inputs + '_tot.fits',
                              'cells_fiducial': path_inputs + '_fid.fits',
                              'cells_noise' : path_inputs  + '_noi.fits'},
                   'config': path_config,
                   'resume': False,
                   'output_dir': PATH_DICT['output_path'] + 'outputs/',
                   'log_dir': PATH_DICT['output_path'] + 'outputs/',
                   'pipeline_log': PATH_DICT['output_path'] + 'outputs/log.txt',
                   'BBCompSep': {'sampler': 'emcee',
                                 'nwalkers': params_bbpw['nwalk'],
                                 'n_iters': params_bbpw['niter'],
                                 'likelihood_type': 'h&l',
                                 'pol_channels': ['B'],
                                 'l_min': dict_bbcomp['lmin'],
                                 'l_max': dict_bbcomp['lmax'],
                                 'bands': dict_bbcomp['bands'],
                                 'cmb_model': dict_cmbmodel,
                                 'fg_model': dict_fgmodel}
                                 }

    with open(path_config, 'w', encoding = 'utf-8') as file:
        yaml.dump(dict_config, file, default_flow_style=True)