'''
    runs bestfit pte (generates saccs). for further analysis with bbpower
'''

import os
from pte.bestfitpte_cell import generatecell_bestfitpte
from utils.params import PATH_DICT

WEIGHT = 'Cl'
TYPE_COV = 'dfwt'

for i in range(10):
    print(i)

    newpath = PATH_DICT['output_path'] + f'{i}/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        os.makedirs(newpath + '/w/')
        os.makedirs(newpath + f'/{TYPE_COV}/')

    generatecell_bestfitpte(WEIGHT, TYPE_COV, i)
