'''
params
'''
import yaml

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

    name_run = '_'.join(map(str, kwargs_dict.values()))

    return name_run

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
    config = Config(yaml.load(open('config.yaml'), yaml.FullLoader))


    path_dict = {}

    bicep_BK15 = '/global/cfs/cdirs/act/data/iabril/BICEP/'
    bbpipe_path = '/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower/'
    output_path = config.path_param.output

    if machine == 'cori':

    #   figure_path = '/global/cscratch1/sd/iabril/plots/'
        planck_path = '/global/cscratch1/sd/iabril/PlanckData/'
        path_dict['planck_data'] = planck_path

    if machine == 'perl':

    #   figure_path = '/pscratch/sd/i/iabril/plots/'
        planck_path = '/pscratch/sd/i/iabril/data/PlanckData/'
        path_dict['planck_data'] = planck_path

    path_dict['BK15_data']   = bicep_BK15
    path_dict['bbpipe_path'] = bbpipe_path
    path_dict['output_path'] = output_path

    return path_dict
