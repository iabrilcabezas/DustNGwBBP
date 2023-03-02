'''
bandpowers
    defines ell binning and Cl to Dl transformations

'''

import numpy as np

def get_ell_arrays(lmin, dell, nbands):

    '''
    Retunrns ell array products

    ** Parameters **
    lmin: int
        starting ell
    dell: int
        spacing of ell bands
    nbands: int
        number of ell bands

    ** Returns **
    lmax: int
        maximum ell
    larr_all: np.array()
        ell array from ell = 0 to ell = lmax (step = 1)
    lbands: np.array()
        ell edges of each band
    leff: np.array()
        effective ell of each band
    '''

    lmax = lmin + dell * nbands
    larr_all = np.arange(lmax+1)
    lbands = np.linspace(lmin,lmax,nbands+1,dtype=int)
    leff = [np.mean(np.arange(x,y)) for x,y in zip(lbands[:-1], lbands[1:])] # 0.5*(lbands[1:]+lbands[:-1])

    return (lmax, larr_all, lbands, leff)

def dell2cell_lmax(lmax):

    '''
    Returns conversion factor from D_ell to C_ell for ell = 0, 1, ..., lmax
    D_ell = ell * (ell + 1) / 2pi * C_ell

    **Parameters**
    lmax: int
        Maximum ell of array
    '''

    if not isinstance(lmax, int):
        print("non-integer input given to function dell2cell_lmax; num="+str(lmax))
        return 0

    ell = np.arange(lmax + 1)

    cl2dl=ell*(ell+1)/(2*np.pi)
    dl2cl=np.zeros_like(cl2dl)
    # ell = 0 is ill-defined
    dl2cl[1:] = 1/(cl2dl[1:])

    return dl2cl
