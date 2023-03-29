'''
params
    config file read in and definition of variables for BBCompSep

'''
import yaml
from utils.params import EXPERIMENT, LMIN, NBANDS, DELL
from utils.sed import get_band_names
from utils.bandpowers import get_ell_arrays

band_names = get_band_names()

SO_bands_equivalent = {'LF1': 'band1', 'LF2':'band2', 'MF1':'band3',
                       'MF2': 'band4', 'UHF1': 'band5', 'UHF2':'band6'}

LMAX = get_ell_arrays(LMIN, DELL, NBANDS)[0]


def get_dictwnamecompsep(dict_bbcomp):

    '''
    Updates dictionary to include name to identify parameters in
    dictionary containing band parameters for BBCompSep analysis

    ** Parameters **
    dict_bbcomp: dict
        lmin : 30
        lmax : 300,
        bands: 'all' or ['UHF1', 'UHF2'], etc.

    ** Returns **
    dict_bbcomp: dict
        name_config
    '''

    bands_bbcomp = dict_bbcomp['bands']
    lmin_bbcomp  = dict_bbcomp['lmin']
    lmax_bbcomp  = dict_bbcomp['lmax']

    if bands_bbcomp != 'all':
        assert all(bb in band_names for bb in bands_bbcomp), 'bands are not in the instrument bands'

    if bands_bbcomp != 'all':
        name_bands =  '_'.join(bands_bbcomp)
    else:
        name_bands = bands_bbcomp

    name_configcompsep = '_'.join([str(lmin_bbcomp), str(lmax_bbcomp), name_bands])

    dict_bbcomp['name_config'] = name_configcompsep

    # update name of bands to match tracer:
    if bands_bbcomp != 'all':
        dict_bbcomp['bands_tr'] = [SO_bands_equivalent.get(b) for b in bands_bbcomp]
    else:
        dict_bbcomp['bands_tr'] = 'all'

    return dict_bbcomp

def get_niterbands():

    '''
    Return list containing bands on which to run BBCompSep analysis
    Depends on EXPERIMENT only
    '''

    if EXPERIMENT == 'so':
        n_iterbands = config_bbpw.biters.so
    else:
        n_iterbands = config_bbpw.biters.other

    return n_iterbands

class ConfigBB:

    '''
    1st level class for config file

    global_emcee: global emcee paramters
    ranges: bandpower parameters
    '''

    def __init__(self, param):
        self.global_emcee = GlobalConfigBB(param['global_emcee'])
        self.ranges       = BpwConfigBB(param['ranges'])
        self.biters       = BitersBB(param['bands_iter'])

class BitersBB:

    '''
    2nd level class of config file
    contains info on bands for BBCompSep to iterate
    so: list
        band combinations in SO experiment
    other: list
        band combinations in other experiments (['all'])
    '''

    def __init__(self, param):
        self.so = param['so']
        self.other = param['other']

class GlobalConfigBB:

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

class BpwConfigBB:


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
    config_bbpw = ConfigBB(yaml.load(config_file, yaml.FullLoader))

dict_params_bbpw = config_bbpw.global_emcee.__dict__
dict_ells_bbpw  = config_bbpw.ranges.__dict__
niterbands = get_niterbands()

lmax = config_bbpw.ranges.lmax
lmin = config_bbpw.ranges.lmin

assert (lmax == LMAX) & (lmin == LMIN), 'bbpw will use different ell range from binned cov'

