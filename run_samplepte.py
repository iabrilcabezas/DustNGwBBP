'''
    runs bestfit pte (generates saccs). for further analysis with bbpower
'''

import os
from pte.cell_sample import SfClassPte, get_random_cell
from utils.params import PATH_DICT
from utils.params import LMIN, DELL, NBANDS
from utils.bandpowers import get_ell_arrays

LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

WEIGHT = 'Cl'
TYPE_COV = 'dfwt'

sf_0 = SfClassPte(type_cov = TYPE_COV, bands = 'all', \
                        lmin_bbp = LMIN, lmax_bbp = LMAX)

for i in range(5001, 10000):
    print(i)

    newpath = PATH_DICT['output_path'] + f'{i}/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        os.makedirs(newpath + '/w/')
        os.makedirs(newpath + f'/{TYPE_COV}/')
        os.makedirs(newpath + f'/{TYPE_COV}/LF/')
        os.makedirs(newpath + f'/{TYPE_COV}/MF/')
        os.makedirs(newpath + f'/{TYPE_COV}/UHF/')

    get_random_cell(i, sf_0, TYPE_COV)
