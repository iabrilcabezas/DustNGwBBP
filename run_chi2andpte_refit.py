'''
chi2 and PTE calculation
    studies the effect of band selection in analysis
'''

from pte.calc_chi2pte_refit import get_chi2andpte
from utils_bbpw.params import dict_ells_bbpw, niterbands

#iterate over all band combinations
for biter in niterbands:
    dict_ells_bbpw['bands'] = biter
    get_chi2andpte(dict_ells_bbpw, ni = int(1e4), no = int(1e4)) # ns= int(1e4))
