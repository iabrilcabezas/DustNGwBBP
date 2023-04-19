'''
clean cmb where we calculate new pte values
'''

import numpy as np
from astropy.io import fits
from utils.sed import get_convolved_seds, get_band_names
from utils.bandpowers import get_ell_arrays
from dustngwbbp.compute_cl import get_windows, import_bandpasses
import utils.noise_calc as nc
from utils.params import LMIN, DELL, NBANDS, POLARIZATION_cov, PATH_DICT, NAME_RUN
from pte.sfid_class_pte import SfClass

band_names = get_band_names()
LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps = nfreqs * nmodes
ncross = (nmaps * (nmaps + 1)) // 2

indices_tr = np.triu_indices(nfreqs)
ncombs = len(indices_tr[0])
assert ncross == ncombs

ncomp = 3
type_cov = 'dfwt'

counter = np.arange(len(indices_tr[0]))
dicttr = {}
for i, counter_val in enumerate(counter):
    dicttr[(indices_tr[0][i], indices_tr[1][i])] = counter_val

fsky = nc.get_fsky()
weight = 'Cl'

bpss = import_bandpasses()
S = get_convolved_seds(band_names, bpss)

sens=2
knee=1
ylf=1

nell=np.zeros([nfreqs,LMAX+1])
_,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,LMAX+1,1, atm_noise = True)

windows = get_windows(weight)
n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)
noise_array = n_bpw

invnoise_ell = np.zeros((NBANDS, nfreqs, nfreqs)) + np.nan

for i in range(NBANDS):
    invnoise_ell[i,:,:]  = np.linalg.inv( np.diag( noise_array[:,i] ) )  # nfreq x nfreq

A_ell = np.zeros(( NBANDS, ncomp, ncomp)) + np.nan
B_ell = np.zeros(( NBANDS, ncomp, nfreqs)) + np.nan
for i in range(NBANDS):
    A_ell[i,:,:] = np.matmul( S, np.matmul(invnoise_ell[i], S.transpose()) )
    B_ell[i,:,:] = np.matmul( S, invnoise_ell[i] )

Q_ell = np.zeros(( NBANDS, nfreqs )) + np.nan
for i in range(NBANDS):
    for j in range(nfreqs):
        Q_ell[i,j] = sum(A_ell[i,0,:] * B_ell[i,:,j])

for no_sim in range(1000):

    sf_0 = SfClass(type_cov = type_cov, bands = 'all', lmin_bbp = LMIN, lmax_bbp = LMAX)
    # obtain cell array and covariances
    cell_array, bbcovar_g , bbcovar_ng, _, _ = sf_0.get_cellandcov()
    bbcovar_ng = bbcovar_ng.reshape([ncross, NBANDS, ncross, NBANDS])
    bbcovar_g = bbcovar_g.reshape([ncross, NBANDS, ncross, NBANDS])

    # run nsims simulations of cell arrays
    random_cells = np.random.default_rng().multivariate_normal(cell_array, bbcovar_ng)

    random_cells = random_cells.reshape((NBANDS, ncross), order = 'F')
    # cell(nu, nuprime):
    cell_nunup_ell = np.zeros((NBANDS, nfreqs, nfreqs))
    for i in range(NBANDS):
        # https://stackoverflow.com/questions/17527693/transform-the-upper-lower-triangular-part-of-a-symmetric-matrix-2d-array-into
        X = np.zeros((nfreqs, nfreqs))
        X[np.triu_indices(X.shape[0], k = 0)] = random_cells[i,:]
        X = X + X.T - np.diag(np.diag(X))
        cell_nunup_ell[i] = X

    clean_cell = np.zeros(NBANDS) + np.nan
    for i in range(NBANDS):
        clean_cell[i] = np.matmul( (Q_ell[i,:]).transpose(), np.matmul(cell_nunup_ell[i], Q_ell[i,:]) )

    cov_cleaned_g = np.zeros((NBANDS, NBANDS))
    cov_cleaned_ng = np.zeros((NBANDS, NBANDS))

    for l1 in range(NBANDS): 
        for l2 in range(NBANDS):
            sumag = 0
            sumang = 0
            for n1 in range(nfreqs):
                for n2 in range(nfreqs):
                    cros1 = (n1,n2)
                    cros1_val = dicttr[tuple(sorted(cros1))]
                    for n3 in range(nfreqs):
                        for n4 in range(nfreqs):
                            cros2 = (n3,n4)
                            cros2_val = dicttr[tuple(sorted(cros2))]
                            qval = (Q_ell[l1, n1] * Q_ell[l1, n2] * Q_ell[l2, n3] * Q_ell[l2, n4]) 
                            sumag  += qval * bbcovar_g[ cros1_val, l1, cros2_val, l2]
                            sumang += qval * bbcovar_ng[cros1_val, l1, cros2_val, l2]
            cov_cleaned_g[l1, l2] = sumag
            cov_cleaned_ng[l1,l2] = sumang

    # save clean cell and covariances with no_sim

    fits_clean = fits.PrimaryHDU(clean_cell)
    fits_clean.writeto(PATH_DICT['output_path'] + \
                       '_'.join([NAME_RUN, 'clean_cell', f'{no_sim}']) +  '.fits', overwrite = True)
    fits_clean_cov_g = fits.PrimaryHDU(cov_cleaned_g)
    fits_clean_cov_g.writeto(PATH_DICT['output_path'] + \
                             '_'.join([NAME_RUN, 'cov_clean_g', f'{no_sim}']) + \
                             '.fits', overwrite = True)
    fits_clean_cov_ng = fits.PrimaryHDU(cov_cleaned_ng)
    fits_clean_cov_ng.writeto(PATH_DICT['output_path'] + \
                              '_'.join([NAME_RUN, 'cov_clean_ng', f'{no_sim}']) +  \
                              '.fits', overwrite = True)
