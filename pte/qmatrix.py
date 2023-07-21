'''
qmatrix

functions required to obtained cleaned CMB spectra

'''

import numpy as np
import sacc
from utils.sed import get_band_names
from utils.params import POLARIZATION_cov, NBANDS
band_names = get_band_names()

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncross = len(indices_tr[0])

def clean_cells(Q, input_path, output_path):

    '''
    clean cells 
    '''
    
    cell_sampled = sacc.Sacc.load_fits(input_path).mean
    cell_sampled = cell_sampled.reshape((ncross, NBANDS))
    
    # reshape into cell(nu, nuprime) shape
    cell_nunup_ell = np.zeros((NBANDS, nfreqs, nfreqs))
    for i in range(NBANDS):
        # https://stackoverflow.com/questions/17527693/transform-the-upper-lower-triangular-part-of-a-symmetric-matrix-2d-array-into
        X = np.zeros((nfreqs, nfreqs))
        X[np.triu_indices(X.shape[0], k = 0)] = cell_sampled[:,i]
        X = X + X.T - np.diag(np.diag(X))
        cell_nunup_ell[i] = X
        
    # clean:
    clean_cell = np.zeros(NBANDS) + np.nan
    for i in range(NBANDS):
        clean_cell[i] = np.matmul( (Q[i,:]).transpose(), np.matmul(cell_nunup_ell[i], Q[i,:]) )
    
    np.savetxt(output_path, clean_cell)

