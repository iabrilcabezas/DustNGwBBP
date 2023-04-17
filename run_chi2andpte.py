'''
chi2 and PTE calculation
    studies the effect of band selection in analysis
'''

from pte.calc_chi2pte import get_chi2andpte
from utils_bbpw.params import dict_ells_bbpw, niterbands

for type_cov in ['wt' , 'dfwt' ]: # 'df00', 'dfwtm']:
    print(type_cov)
    #iterate over all band combinations
    for biter in niterbands:
        dict_ells_bbpw['bands'] = biter
        get_chi2andpte(type_cov, dict_ells_bbpw, nsims= int(1e5))
