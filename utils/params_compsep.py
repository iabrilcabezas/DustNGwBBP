from utils.sed import get_band_names

band_names = get_band_names()

def get_configcompsep(dict_bbcomp):

    '''
    dict_bbcomp: dict
        lmin : 30
        lmax : 300,
        bands: 'all' or ['UHF1', 'UHF2'], etc.
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

    return (lmin_bbcomp, lmax_bbcomp, bands_bbcomp, name_configcompsep)

# bbcomp:
#     lmin: 30
#     lmax: 300
#     bands:  'all' #['UHF1','UHF2']
#     self.bbcomp       = BBCompConfig(param['bbcomp'])
#     class BBCompConfig:
#     def __init__(self, param):
#         self.lmin = param['lmin']
#         self.lmax = param['lmax']
#         self.bands = param['bands']
# LMIN_BBCOMP = config.bbcomp.lmin
# LMAX_BBCOMP = config.bbcomp.lmax
# BANDS_BBCOMP = config.bbcomp.bands
# if BANDS_BBCOMP != 'all':
#     name_bands =  '_'.join(BANDS_BBCOMP)
# else:
#     name_bands = BANDS_BBCOMP
# name_configcompsep = '_'.join([str(LMIN_BBCOMP), str(LMAX_BBCOMP), name_bands])
# print(name_configcompsep)
# from utils.params import LMAX_BBCOMP, LMIN_BBCOMP, BANDS_BBCOMP, name_configcompsep
