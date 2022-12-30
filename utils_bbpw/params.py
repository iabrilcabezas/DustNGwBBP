'''
params
    config file read in and definition of variables for BBCompSep

'''
import yaml
from utils.params import EXPERIMENT

def get_niterbands():

    if EXPERIMENT == 'so':
        n_iterbands = config_bbpw.biters.so
    else:
        n_iterbands = config_bbpw.biters.other

    return n_iterbands

class Config_bbpw:

    '''
    1st level class for config file

    global_param: global parameters
    bpw_param: bandpower parameters
    '''

    def __init__(self, param):
        self.global_param = GlobalConfig_bbpw(param['params'])
        self.bpw_param    = BpwConfig_bbpw(param['compsep'])
        self.biters       = Biters_bbpw(param['bands_iter'])

class Biters_bbpw:

    def __init__(self, param):
        self.so = param['so']
        self.other = param['other']

class GlobalConfig_bbpw:

    '''
    2nd level class of config file
    contains info on global params for BBCompSep
    niter: int
        number iterations MCMC
    nwalk: int
        number of walkers MCMC
    nside: int
        resolution for analysis
    '''

    def __init__(self, param):
        self.niter = param['niter']
        self.nside = param['nside']
        self.nwalk = param['nwalk']

class BpwConfig_bbpw:


    '''
    2nd level class of config file
    contains info on global params for BBCompSep
    lmin: int
        min ell for analysis
    lmax: int
        max ell for analysis
    # bands: list
    #     bands used for analysis
    '''

    def __init__(self, param):
        self.lmin  = param['lmin']
        self.lmax  = param['lmax']
    #    self.bands = param['bands']

with open('config_bbpw.yaml', 'r', encoding = 'utf-8') as config_file:
    config_bbpw = Config_bbpw(yaml.load(config_file, yaml.FullLoader))

dict_params_bbpw = config_bbpw.global_param.__dict__
dict_ells_bbpw  = config_bbpw.bpw_param.__dict__
niterbands = get_niterbands()
