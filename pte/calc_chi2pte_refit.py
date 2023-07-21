'''
function to do chi2 and PTE calculation
'''
from copy import deepcopy
import numpy as np
from pte.sfid_class_pte_refit import SfClass_refit
from pte.sfid_class_pte_refit import get_chi2
from utils.sed import get_band_names
from utils_bbpw.params import get_dictwnamecompsep

band_names = get_band_names()

def get_chi2andpte(dict_compsep, ni, no):

    '''
    Writes chi2 and PTEs calculated with G and NG covariances for NAME_RUN fiducial Cell model

    ** Parameters **
    dict_compsep: dict
        dict containing ell range and bands of analysis
    ns: int
        number of simulations to run
    '''

    # include name of band+ell specification in dictionary
    dict_bbcomp = deepcopy(dict(get_dictwnamecompsep(dict_compsep)))

    chi2_array = np.zeros(no)

    for cc,i in enumerate(range(ni, ni+no)):
        if i%100 == 0:
            print(i)

        # Sfclass object with all bands
        s_fid_all = SfClass_refit(nosim = i, bands = 'all', \
                            lmin_bbp =  dict_bbcomp['lmin'], lmax_bbp = dict_bbcomp['lmax'])

        if dict_bbcomp['bands'] != 'all':
            assert all(bb in band_names for bb in dict_bbcomp['bands']), \
                    'bands are not in the instrument bands'

            band_dict = s_fid_all.name_band2trac()
            bands_sf = [band_dict.get(key) for key in dict_bbcomp['bands']]
            # Sfclass object with user-specified bands (if != all)
            s_fid = SfClass_refit(nosim = i, bands = bands_sf , \
                        lmin_bbp = dict_bbcomp['lmin'], lmax_bbp = dict_bbcomp['lmax'])

        else:
            s_fid = s_fid_all

        # compute chi2 and pvalues
        chi2_array[cc], _ = get_chi2(s_fid)

    # print(dof)
    # print(chi2_array)
    # name of run
    name_chi2 = '/global/cfs/cdirs/act/data/iabril/BBPower/230525_sample/sims/' + dict_bbcomp['name_config'] + f'_{ni}_{no}'

    # save to file
    np.savetxt(name_chi2 + '_chi2_g.txt',  chi2_array)
