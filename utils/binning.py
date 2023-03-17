'''
binning
    binning and slicing of 1- and 2-D arrays
'''

import numpy as np
from utils.params import POLARIZATION_cov, NSIDE, LMIN, DELL, NBANDS
from utils.sed import get_band_names
from utils.bandpowers import get_ell_arrays

LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)
band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncombs = len(indices_tr[0])

# Only work with 1D and 2D arrays
ARR_TYPE = [1,2]
ERROR_MSG_ARR_TYPE = 'method only works for 1D and 2D arrays!'


def cut_array(array, ell_cl, lmin, lmax):

    '''
    Cuts array in terms of ell range using array indexing

    Parameters
    ----------
    array: np.array()
            array to be cut. must be square matrix if 2D
    ell_cl: np.array()
            ell each value of array corresponds to
    lmin: float
            min ell
    lmax: float
            max ell

    Returns
    ----------
    np.array()
            array indexed into new shape
    '''

    imin = np.where(ell_cl == lmin)[0][0]
    imax = np.where(ell_cl == lmax)[0][0]

     # determine if 1D or 2D:
    shape = len(array.shape)
    assert shape in ARR_TYPE, ERROR_MSG_ARR_TYPE

    if shape == 1: # 1D

        return array[imin:imax]

    if shape == 2: # 2D

        assert array.shape[0] == array.shape[1], "2D array is not a square matrix!"

        return array[imin:imax, imin:imax]

    return None

def rebin(array, new_shape):

    '''
    Rebin 1D or 2D array

    Applies equal weights using reshape() method

    Parameters
    ----------
    array: np.array()
            array to rebin
    new_shape: list
            target shape of array

    Returns
    ----------
    np.array()
            array rebinned to new shape
    '''

    # determine if 1D or 2D:
    old_shape = len(array.shape)

    assert old_shape in ARR_TYPE, ERROR_MSG_ARR_TYPE

    if old_shape == 1: # 1D

        # from https://stackoverflow.com/questions/21921178/binning-a-numpy-array
        shape = array.shape[0] // new_shape[0]

        return array.reshape(-1, shape).mean(axis = 1)

    if  old_shape == 2: # 2D

        # from https://scipython.com/blog/binning-a-2d-array-in-numpy/
        shape = (new_shape[0], array.shape[0] // new_shape[0],
                new_shape[1], array.shape[1] // new_shape[1])

        return array.reshape(shape).mean(-1).mean(1)

    return None

def bin_covs(matrix):

    '''
    Bins unbinned covariance matrix according to NBANDS, LMIN AND LMAX range

    * Parameters*

    matrix: np.array()
        cov matrix to be binned

    * Returns *
    np.array()
    
    '''
    assert matrix.shape == (ncombs, 3 * NSIDE, ncombs, 3* NSIDE),'is this your unbinned cov matrix?'

    # save binned:
    cov_bin = np.zeros([ncombs, NBANDS, ncombs, NBANDS])
    for i_tr in range(ncombs):
        for j_tr in range(ncombs):
            cov_bin[i_tr,:, j_tr,:]  = rebin(cut_array(matrix[i_tr,:, j_tr,:], \
                                                np.arange(3 * NSIDE), LMIN, LMAX), \
                                                [NBANDS, NBANDS])

    return cov_bin
