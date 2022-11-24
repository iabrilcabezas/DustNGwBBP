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

def dell2cell_lmax(lmax):

    '''
    standard conversion Dell to Cell of ell array 0, ... , lmax

    **Parameters**
    lmax: int
    '''

    ell = np.arange(lmax + 1)

    cl2dl=ell*(ell+1)/(2*np.pi)
    dl2cl=np.zeros_like(cl2dl)
    dl2cl[1:] = 1/(cl2dl[1:])

    return dl2cl
