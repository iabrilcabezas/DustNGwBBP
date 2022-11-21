import numpy as np

# Only work with 1D and 2D arrays
arr_type = [1,2]
error_msg_arr_type = 'method only works for 1D and 2D arrays!'


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
    assert shape in arr_type, error_msg_arr_type

    if shape == 1: # 1D
        
        return array[imin:imax]

    elif shape == 2: # 2D

        assert array.shape[0] == array.shape[1], "2D array is not a square matrix!"

        return array[imin:imax, imin:imax]

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

    assert old_shape in arr_type, error_msg_arr_type
    
    if old_shape == 1: # 1D
        
        # from https://stackoverflow.com/questions/21921178/binning-a-numpy-array
        shape = array.shape[0] // new_shape[0]

        return array.reshape(-1, shape).mean(axis = 1)
   
    elif old_shape == 2: # 2D

        # from https://scipython.com/blog/binning-a-2d-array-in-numpy/
        shape = (new_shape[0], array.shape[0] // new_shape[0],
                new_shape[1], array.shape[1] // new_shape[1])
        
        return array.reshape(shape).mean(-1).mean(1)
