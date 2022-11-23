'''
ell stuff
'''

import numpy as np

def get_ell_arrays(lmin, dell, nbands):

    '''
    ell arrays from parameters
    '''

    lmax = lmin + dell * nbands
    larr_all = np.arange(lmax+1)
    lbands = np.linspace(lmin,lmax,nbands+1,dtype=int)
    leff = 0.5*(lbands[1:]+lbands[:-1])

    return (lmax, larr_all, lbands, leff)

def 
