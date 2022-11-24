'''
coupling matrix

'''

import pymaster as nmt
import numpy as np
import yaml
from utils.namaster import get_bandpowers, get_template_wmask, get_couplingmatrix
from utils.params import pathnames, get_namerun, Config

config = Config(yaml.load(open('config.yaml'), yaml.FullLoader))

NSIDE   = config.global_param.nside
MACHINE = config.global_param.machine
MTYPE   = config.mask_param.mask_type
DELL_NMT = config.bpw_param.dell_nmt

path_dict = dict(pathnames(MACHINE))

def compute_couplingmatrix(**kwargs_dict):

    '''
    Writes to file the coupling matrix (with and without template) for a given mask

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

    name_run  = get_namerun(kwargs_dict)
    # get bandpowers
    bin_bpw, _ = get_bandpowers(NSIDE, kwargs_dict['dell_nmt'])
    # get template and mask:
    template, w_mask = get_template_wmask(MACHINE, NSIDE, MTYPE, **kwargs_dict)

    f_d2 = nmt.NmtField(w_mask * w_mask, None, spin = 0)
    f_dtilde2 = nmt.NmtField((w_mask * template)**2, None, spin = 0)

    mw2_matrix =  get_couplingmatrix(f_d2, f_d2, bin_bpw)
    mwt2_matrix = get_couplingmatrix(f_dtilde2, f_dtilde2, bin_bpw)

    np.savetxt(path_dict['output_path'] + '_'.join([MTYPE, name_run]) + '_couplingM_w.txt',
                 mw2_matrix)
    np.savetxt(path_dict['output_path'] + '_'.join([MTYPE, name_run]) + '_couplingM_wt.txt',
                 mwt2_matrix)
