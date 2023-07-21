'''
add_cov2decorr
    adds gaussian and non-gaussian cov matrices to data with decorrelation
'''

import numpy as np
import copy
from astropy.io import fits
import sacc
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP, LMIN, DELL
from utils.params import POLARIZATION_cov, NBANDS, EXPERIMENT
from utils.sed import get_band_names
import utils.noise_calc as nc
from utils.bandpowers import get_ell_arrays
from dustngwbbp.compute_cl import get_windows, add_tracers, add_powerspectra

LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)


band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
pols = [POLARIZATION_cov]
nmaps=nmodes*nfreqs
ncross = (nmaps * (nmaps + 1))//2

def add_cov2decorr(type_cov, weight):

    '''
    reads in cells with decorrelation, and adds covariance to saccs
    saves gaussian (w) and non gaussian (type_cov) covariance
    '''

    # load fiducial power spectrum for NAME_RUN parameters
    name_sf = PATH_DICT['output_path'] + 'cells_model.fits'

    s_fg = sacc.Sacc.load_fits(name_sf)
    s_fng = sacc.Sacc.load_fits(name_sf)

    # no. of bandpowers after ell range selection
    tr_names = sorted(list(s_fg.tracers.keys()))

    ell_b = s_fg.get_ell_cl('cl_' + 2 * pols[0].lower(), \
                                            tr_names[0], tr_names[0])[0]
    n_bpws = len(ell_b)


    # populate with gaussian and NG covariance
    cov_ng = fits.open(PATH_DICT['input_path'] + NAME_RUN + '_' + \
                            '_'.join([NAME_COMP,'Cov', 'bin']) + f'_{type_cov}.fits')[0].data
    cov_ng = cov_ng.reshape([ncross* n_bpws,ncross * n_bpws ])
    s_fng.add_covariance(cov_ng)
    # repeat:
    cov_g = fits.open(PATH_DICT['input_path'] + NAME_RUN + '_' + \
                            '_'.join([NAME_COMP,'Cov', 'bin']) + '_w.fits')[0].data
    cov_g = cov_g.reshape([ncross* n_bpws,ncross * n_bpws ])
    s_fg.add_covariance(cov_g)

    # save:
    s_fg.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, 'w']) + \
                    '_tot.fits', overwrite = True)
    s_fng.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, type_cov]) + \
                    '_tot.fits', overwrite = True)

def add_covandnoise2decorr():
    '''
    reads in cells with decorrelation, and adds covariance and noise! to saccs
    saves gaussian (w) and non gaussian (type_cov) covariance
    '''

    # load fiducial power spectrum for NAME_RUN parameters
    name_sf = PATH_DICT['output_path'] + 'cells_model.fits'

    cell_model = sacc.Sacc.load_fits(name_sf)
    cell_model_matrix = cell_model.mean.reshape((ncross, NBANDS))

    # cell(nu, nuprime):
    cell_nunup_ell = np.zeros((nfreqs, nfreqs, NBANDS))
    for i in range(NBANDS):
        # https://stackoverflow.com/questions/17527693/transform-the-upper-lower-triangular-part-of-a-symmetric-matrix-2d-array-into
        X = np.zeros((nfreqs, nfreqs))
        X[np.triu_indices(X.shape[0], k = 0)] = cell_model_matrix[:,i]
        X = X + X.T - np.diag(np.diag(X))
        cell_nunup_ell[:,:, i] = X

    weight = 'Cl'
    windows = get_windows(weight)
    fsky = nc.get_fsky()

    if EXPERIMENT == 'so':

        # N_ell
        sens=2
        knee=1
        ylf=1
        nell=np.zeros([nfreqs,LMAX+1])
        _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,LMAX+1,1, atm_noise = True)

    else:
        print('only implemented for SO EXPERIMENT')
        return None

    noise_array=np.sum(nell[:,None,:]*windows[None,:,:],axis=2) # bin

    sn_array = copy.deepcopy(cell_nunup_ell)
    for ib in range(nfreqs):
        sn_array[ib,ib,:] += noise_array[ib]

    s_g = add_tracers(LARR_ALL)
    s_ng = add_tracers(LARR_ALL)
    s_g = add_powerspectra(s_g, sn_array, LEFF, True, 'Cl')
    s_ng = add_powerspectra(s_ng, sn_array, LEFF, True, 'Cl')

    # populate with gaussian and NG covariance
    cov_ng = fits.open(PATH_DICT['input_path'] + NAME_RUN + '_' + \
                            '_'.join([NAME_COMP,'Cov', 'bin']) + '_dfwt.fits')[0].data
    cov_ng = cov_ng.reshape([ncross* NBANDS,ncross * NBANDS ])
    s_ng.add_covariance(cov_ng)
    # repeat:
    cov_g = fits.open(PATH_DICT['input_path'] + NAME_RUN + '_' + \
                            '_'.join([NAME_COMP,'Cov', 'bin']) + '_w.fits')[0].data
    cov_g = cov_g.reshape([ncross* NBANDS,ncross * NBANDS ])
    s_g.add_covariance(cov_g)

    # save:
    s_g.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, 'w']) + \
                    '_tot.fits', overwrite = True)
    s_ng.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, 'dfwt']) + \
                    '_tot.fits', overwrite = True)
