'''
params
'''
import yaml


def get_namerun(kwargs_dict):

    '''
    Returns string that contains all arguments used in the code run

    ** Parameters **
    kwargs_dict: dict
            'width_l'
            'apo_deg'
            'smooth_deg'

    ** Returns **
    str
    '''

    return '_'.join(map(str, kwargs_dict.values()))

def pathnames(machine):

    '''
    Returns path to data depending on machine

    ** Parameters **
    machine: str
            machine where code runs ('cori' or 'perl')

    **Returns**
    dict
        stores path to data
        'planck_path': Planck data
    '''
    config_load = Config(yaml.load(open('config.yaml'), yaml.FullLoader))

    dict_path = {}

    bicep_bk15 = config_load.path_param.bicep_bk15
    bbpipe_path = config_load.path_param.bbpipe
    output_path = config_load.path_param.output
    camb_cmb_lens_nobb_path = config_load.path_param.camb_cmb_lens_nobb

    if machine == 'cori':

        dict_path['planck_data'] = config_load.path_param.planck_path_cori

    if machine == 'perl':

        dict_path['planck_data'] = config_load.path_param.planck_path_perl

    dict_path['BK15_data']   = bicep_bk15
    dict_path['bbpipe_path'] = bbpipe_path
    dict_path['output_path'] = output_path
    dict_path['camb_cmb_lens_nobb'] = camb_cmb_lens_nobb_path

    return dict_path

class GlobalConfig:
    '''
    doc
    '''
    def __init__(self, param):
        self.machine    = param['machine']
        self.experiment = param['experiment']
        self.nside      = param['nside']

class BpwConfig:
    '''
    a
    '''
    def __init__(self, param):
        self.lmin = param['lmin']
        self.nbands = param['nbands']
        self.dell = param['dell']

class PolConfig:
    '''
    a
    '''
    def __init__(self, param):
        self.pol_cov = param['pol_cov']

class PathConfig:
    '''
    a
    '''
    def __init__(self, param):
        self.output = param['output_path']
        self.bicep_bk15 = param['bicep_BK15']
        self.bbpipe  = param['bbpipe_path']
        self.planck_path_cori = param['planck_path_cori']
        self.planck_path_perl = param['planck_path_perl']
        self.camb_cmb_lens_nobb = param['camb_cmb_lens_nobb']

class MaskConfig:
    '''
    a
    '''
    def __init__(self, param):
        self.mask_type  = param['mask_type']
        self.apo_deg    = param['apo_deg']
        self.smooth_deg = param['smooth_deg']
        self.dell_nmt   = param['dell_nmt']


class CosmoConfig:

    '''
    a
    '''
    def __init__(self, param):
        self.CMB = param['CMB']
        self.dust = param['dust']
        self.model = param['model']

class BandConfig:

    '''
    bandnames
    '''
    
    def __init__(self, param):
        self.so = param['SO']
        self.bicep = param['BICEP']

class Config:
    '''
    a
    '''
    def __init__(self, param):
        self.global_param = GlobalConfig(param['global'])
        self.bpw_param    = BpwConfig(param['bandpowers'])
        self.pol_param    = PolConfig(param['polarization'])
        self.path_param   = PathConfig(param['paths'])
        self.mask_param   = MaskConfig(param['mask'])
        self.cosmo_param  = CosmoConfig(param['cosmology'])
        self.band_names   = BandConfig(param['band_names'])


config = Config(yaml.load(open('config.yaml'), yaml.FullLoader))

NSIDE   = config.global_param.nside
MACHINE = config.global_param.machine
EXPERIMENT = config.global_param.experiment
MTYPE   = config.mask_param.mask_type
DELL_NMT = config.mask_param.dell_nmt
LMIN = config.bpw_param.lmin
DELL = config.bpw_param.dell
NBANDS = config.bpw_param.nbands
POLARIZATION_cov = config.pol_param.pol_cov

band_names_config = config.band_names

cosmo_params = config.cosmo_param
cmb_params = cosmo_params.CMB
dust_params = cosmo_params.dust
model_params = cosmo_params.model

A_dust_BB = dust_params['A_dust_BB']
EB_dust = dust_params['EB_dust']
alpha_dust_EE = dust_params['alpha_dust_EE']
alpha_dust_BB = dust_params['alpha_dust_BB']
beta_dust = dust_params['beta_dust']
temp_dust = dust_params['temp_dust']
nu0_dust = dust_params['nu0_dust']
Alens = cmb_params['Alens']
lnorm_PL = model_params['lnorm_PL']

path_dict = dict(pathnames(MACHINE))

namerun_dict = {**config.global_param.__dict__, **config.bpw_param.__dict__, \
                **config.mask_param.__dict__, **config.pol_param.__dict__}
name_run  = get_namerun(namerun_dict)
