'''
chi2 and PTE calculation
'''

import numpy as np
from pte.sfid_class_pte import SfClass
from pte.sfid_class_pte import ptechi2_gvsng
from utils.sed import get_band_names
from utils.params import LMAX_BBCOMP, LMIN_BBCOMP, BANDS_BBCOMP
from utils.params import name_configcompsep
from utils.params import PATH_DICT, NAME_RUN

band_names = get_band_names()
NSIMS = int(1e4)

Sf_all = SfClass(bands = 'all',lmin_bbp =  LMIN_BBCOMP, lmax_bbp = LMAX_BBCOMP)

if BANDS_BBCOMP != 'all':
    assert all(bb in band_names for bb in BANDS_BBCOMP), 'bands are not in the instrument bands'

    band_dict = Sf_all.name_band2trac()
    bands_sf = [band_dict.get(key) for key in BANDS_BBCOMP]

    Sf = SfClass(bands = bands_sf ,lmin_bbp =  LMIN_BBCOMP, lmax_bbp = LMAX_BBCOMP)

else:
    Sf = Sf_all

chi_g_array, chi_ng_array, p_g_array, p_ng_array =  ptechi2_gvsng(NSIMS, Sf)

name_chi2pte = PATH_DICT['output_path'] + 'results_pte/' + NAME_RUN + '_' + name_configcompsep

np.savetxt(name_chi2pte + '_chi2_g.txt',  chi_g_array)
np.savetxt(name_chi2pte + '_chi2_ng.txt', chi_ng_array)
np.savetxt(name_chi2pte + '_pval_g.txt', p_g_array)
np.savetxt(name_chi2pte + '_pval_ng.txt', p_ng_array)
