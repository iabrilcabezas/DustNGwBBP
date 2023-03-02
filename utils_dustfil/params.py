'''
params
    imports global config parameters for DustFilament sims analysis
'''

import yaml

class PathConfig:

    '''
    2nd level class of config file
    contains path info

    base: str
        location of all DustFilaments sims
    sim: str
        specifies current run
    base_name: str
        initial characters of all output sims names
    out: str
        folder on base + sim to store output
    mid_name: str
        middle characters of output sims names
    end_name: str
        end characters of output sims names
    '''

    def __init__(self, param):
        self.base = param['base']
        self.sim  = param['sim']
        self.out  = param['output']
        self.base_name  = param['base_name']
        self.mid_name   = param['mid_name']
        self.end_name   = param['end_name']

class GlobalConfig:

    '''
    2nd class

    freq: float
        normalization SED frequency
    nside: int
        resolution of sims
    '''

    def __init__(self, param):
        self.nside = param['nside']
        self.freq  = param['freq']
        self.nsims = param['nsims']
        self.mask  = param['mask']
        self.lmin  = param['lmin']
        self.lmax  = param['lmax']
        self.dell  = param['dell']
        self.tcov  = param['type_cov']

class NormConfig:

    '''
    2nd class

    alphad: float
        alpha scaling for dust
    ampd: float
        amplitude of power dust spectrum
    '''

    def __init__(self, param):
        self.ampdd   = param['ampd']
        self.alphad  = param['alphad']

class ConfigDF:


    '''
    1st level class of config file
    
    paths: path info        
    glob: global params
    mask: region of interest (SO patch)
    '''

    def __init__(self, param):
        self.paths  = PathConfig(param['paths'])
        self.glob = GlobalConfig(param['global'])
        self.norm = NormConfig(param['norm'])


with open('config_dustfil.yaml', 'r', encoding = 'utf-8') as config_file:
    config_dfil = ConfigDF(yaml.load(config_file, yaml.FullLoader))

DF_FREQ = config_dfil.glob.freq
DF_MASK = config_dfil.glob.mask
DF_NSIMS = config_dfil.glob.nsims
DF_NSIDE = config_dfil.glob.nside
DF_LMIN = config_dfil.glob.lmin
DF_LMAX = config_dfil.glob.lmax
DF_DELL = config_dfil.glob.dell
TYPE_DCOV = config_dfil.glob.tcov

# construct full path
DF_BASE_PATH = config_dfil.paths.base + config_dfil.paths.sim + config_dfil.paths.base_name
DF_END_NAME = config_dfil.paths.mid_name + f'{DF_FREQ}' + config_dfil.paths.end_name
DF_OUTPUT_PATH = config_dfil.paths.base + config_dfil.paths.sim + config_dfil.paths.out

DF_ALPHA = config_dfil.norm.alphad
DF_AMP   = config_dfil.norm.ampdd

DF_NAME_RUN = f'{DF_FREQ}_{DF_MASK}_{DF_NSIDE}_{abs(DF_ALPHA)}_{DF_AMP}_{DF_LMIN}_{DF_LMAX}_{DF_DELL}_{TYPE_DCOV}'
print(DF_NAME_RUN)
