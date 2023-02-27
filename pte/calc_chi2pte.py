'''
function to do chi2 and PTE calculation
'''
from copy import deepcopy
import numpy as np
from pte.sfid_class_pte import SfClass
from pte.sfid_class_pte import ptechi2_gvsng
from utils.sed import get_band_names
from utils.params import PATH_DICT, NAME_RUN
from utils_bbpw.params import get_dictwnamecompsep

band_names = get_band_names()

def get_chi2andpte(type_cov, dict_compsep, nsims= int(1e4)):

    '''
    Writes chi2 and PTEs calculated with G and NG covariances for NAME_RUN fiducial Cell model

    ** Parameters **
    dict_compsep: dict
        dict containing ell range and bands of analysis
    nsims: int
        number of simulations to run
    '''

    # include name of band+ell specification in dictionary
    dict_bbcomp = deepcopy(dict(get_dictwnamecompsep(dict_compsep)))

    # Sfclass object with all bands
    s_fid_all = SfClass(type_cov = type_cov, bands = 'all', \
                        lmin_bbp =  dict_bbcomp['lmin'], lmax_bbp = dict_bbcomp['lmax'])

    if dict_bbcomp['bands'] != 'all':
        assert all(bb in band_names for bb in dict_bbcomp['bands']), \
                'bands are not in the instrument bands'

        band_dict = s_fid_all.name_band2trac()
        bands_sf = [band_dict.get(key) for key in dict_bbcomp['bands']]
        # Sfclass object with user-specified bands (if != all)
        s_fid = SfClass(type_cov = type_cov, bands = bands_sf , \
                    lmin_bbp = dict_bbcomp['lmin'], lmax_bbp = dict_bbcomp['lmax'])

    else:
        s_fid = s_fid_all

    # compute chi2 and pvalues
    chi_g_array, chi_ng_array, p_g_array, p_ng_array =  ptechi2_gvsng(nsims, s_fid)

    # name of run
    name_chi2pte = PATH_DICT['output_path'] + 'results_pte/' + \
                    NAME_RUN + '_' + dict_bbcomp['name_config'] + '_' + type_cov

    # save to file
    np.savetxt(name_chi2pte + '_chi2_g.txt',  chi_g_array)
    np.savetxt(name_chi2pte + '_chi2_ng.txt', chi_ng_array)
    np.savetxt(name_chi2pte + '_pval_g.txt', p_g_array)
    np.savetxt(name_chi2pte + '_pval_ng.txt', p_ng_array)
