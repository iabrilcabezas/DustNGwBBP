'''
compute_couplingmatrix
    computes couplingmatrix with namaster
'''

import pymaster as nmt
import numpy as np
from utils.namaster import get_bandpowers, get_template_wmask, get_couplingmatrix
from utils.params import NSIDE, MTYPE
from utils.params import name_couplingmatrix_w, name_couplingmatrix_wt

def compute_couplingmatrix(**kwargs_dict):

    '''
    Writes to file the coupling matrix (with and without template) for a given mask

    Coupling matrix is squared, of length 3 * NSIDE - 1

    ** Parameters **
    machine: str
            machine where code is run ('perl' or 'cori')
    nside: int
            target resolution of maps
    mtype: str
            mask type ('bicep' or 'full')
    kwargs_dict: dict
            'width_l'
            'apo_deg'
            'smooth_deg'
    '''

    # get bandpowers
    bin_bpw, _ = get_bandpowers(NSIDE, kwargs_dict['dell_nmt'])
    # get template and mask:
    template, w_mask, w2omega = get_template_wmask(NSIDE, MTYPE, **kwargs_dict)

    f_d2 = nmt.NmtField(w_mask * w_mask, None, spin = 0)
    f_dtilde2 = nmt.NmtField((w_mask * template)**2, None, spin = 0)

    mw2_matrix =  get_couplingmatrix(f_d2, f_d2, bin_bpw)
    mwt2_matrix = get_couplingmatrix(f_dtilde2, f_dtilde2, bin_bpw)

    np.savetxt(name_couplingmatrix_w , mw2_matrix)
    np.savetxt(name_couplingmatrix_wt, mwt2_matrix)

    return w2omega
