'''
sample cell from NG and generate sacc with this and {NG, G} covariances
'''

import numpy as np
from pte.sfid_class_pte import SfClass
from utils.params import PATH_DICT, NAME_RUN
from utils.params import LMIN, DELL, NBANDS, POLARIZATION_cov
from utils.bandpowers import get_ell_arrays
from utils.sed import get_band_names
from dustngwbbp.compute_cl import add_tracers, add_powerspectra_1d

LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps = nfreqs * nmodes
ncross = (nmaps * (nmaps + 1)) // 2

def generatecell_bestfitpte(weight, type_cov, no_sim):

    '''
    reads in fiducial model, and generates new cells from NG covariance matrix
    saves to sacc files, containing {G, NG} cov matrices
    '''

    sf_0 = SfClass(type_cov = type_cov, bands = 'all', \
                        lmin_bbp = LMIN, lmax_bbp = LMAX)
        # obtain cell array and covariances
    cell_array, bbcovar_g , bbcovar_ng, _, _ = sf_0.get_cellandcov()

    # run nsims simulations of cell arrays
    random_cells = np.random.default_rng().multivariate_normal(cell_array, bbcovar_ng)
    # reshape:
    random_cells = random_cells.reshape((NBANDS, ncross), order = 'F')

    # 2 new sacc files by adding tracers:
    sd_g = add_tracers(LARR_ALL)
    sd_ng = add_tracers(LARR_ALL)
    # add power spectra
    sd_g  = add_powerspectra_1d(sd_g,  random_cells, LEFF, True, weight)
    sd_ng = add_powerspectra_1d(sd_ng, random_cells, LEFF, True, weight)

    # add covariances:
    sd_g.add_covariance(bbcovar_g)
    sd_ng.add_covariance(bbcovar_ng)

    sd_g.save_fits(PATH_DICT['output_path'] + f'{no_sim}/w/' + \
                   '_'.join([NAME_RUN, weight, 'w']) + \
                   '_tot.fits', overwrite = True)
    sd_ng.save_fits(PATH_DICT['output_path'] + f'{no_sim}/{type_cov}/' + \
                   '_'.join([NAME_RUN, weight, type_cov]) + \
                   '_tot.fits', overwrite = True)
