'''
compute_cov
    computes covariance matrix as expressed in the Appendix
    computes effective covariance matrix for final power spectrum
'''

import numpy as np
import sacc
from astropy.io import fits
from utils.params import NBANDS, NSIDE, LMIN, DELL, POLARIZATION_cov, PATH_DICT, NAME_RUN
from utils.sed import get_band_names
from utils.bandpowers import get_ell_arrays

band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)

lmax, _, _, _ = get_ell_arrays(LMIN, DELL, NBANDS)

def compute_cov(ctype):

    '''
    Compute covariance (with and without modulating template) of ctype power spectrum
    Covariances are unbinned, from ell = 0 to ell = 3 * NSIDE - 1

    This is because coupling matrix is given in that shape. and binning of covariance is only
    linear once the final covariance has been calculated (cannot bin elements first)

    ** Parameters **

    ctype: 'dust' or 'all'
        type of power spectrum
    '''

    assert len(POLARIZATION_cov == 1), 'covariance formula is for spin 0 field only!'
    name_pol_cov = 'cl_' + 2 * POLARIZATION_cov.lower()

    # READ IN CL
    s_d = sacc.Sacc.load_fits(PATH_DICT['output_path'] + \
                                '_'.join([NAME_RUN, ctype, 'Cl', str(NSIDE)]) +\
                                '_tot.fits')
    mw2_matrix = np.loadtxt(PATH_DICT['output_path'] + NAME_RUN + '_couplingM_w.txt')
    mwt2_matrix = np.loadtxt(PATH_DICT['output_path'] + NAME_RUN + '_couplingM_wt.txt')

    # read ell of cls from band1 as there'll always be band1
    ellcl, _ = s_d.get_ell_cl(name_pol_cov, 'band1', 'band1', return_cov = False)

    nells = len(ellcl)
    assert nells == mw2_matrix.shape[0], 'using matrix or Cl that should not be binned'

    ncombs = len(s_d.get_tracer_combinations())
    assert ncombs == (nmodes*nfreqs) * (nmodes*nfreqs+1) //2,\
             'incorrect number of auto/cross spectra'

    c_ell_all = np.zeros((nfreqs, nfreqs, nells)) + np.nan

    # populate C_ell_all
    for i in range(nfreqs):
        for j in range(i, nfreqs):

            c_ell_all[i][j] = s_d.get_ell_cl(name_pol_cov, f'band{i+1}', f'band{j+1}')[1]

            if i != j:
                c_ell_all[j][i] = s_d.get_ell_cl(name_pol_cov, f'band{i+1}', f'band{j+1}')[1]

    cov = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    cov_w = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    cov_wt = np.zeros((ncombs, nells, ncombs, nells)) + np.nan

    indices_tr = np.triu_indices(nfreqs)

    # compute prefactor of covariance
    prefactor_l = np.zeros([nells, nells]) + np.nan
    sum_m = 1/ (2 * ellcl + 1)

    for i in range(int(nells)):
        prefactor_l[:, i ] = sum_m[i] * np.ones(nells)

    # populate covariance
    for i_tr, (i_tr1,i_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        print(i_tr)
        for j_tr, (j_tr1,j_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):

            cl1 = c_ell_all[i_tr1][j_tr1]
            cl2 = c_ell_all[i_tr2][j_tr2]
            cl3 = c_ell_all[i_tr1][j_tr2]
            cl4 = c_ell_all[i_tr2][j_tr1]

            for a_ell in range(nells):
                for b_ell in range(nells):
                    # calculate mean between cl at different ells
                    c1_l12 = np.sqrt(cl1[a_ell] * cl1[b_ell])
                    c2_l12 = np.sqrt(cl2[a_ell] * cl2[b_ell])
                    c3_l12 = np.sqrt(cl3[a_ell] * cl3[b_ell])
                    c4_l12 = np.sqrt(cl4[a_ell] * cl4[b_ell])
                    # covariance between different frequency channels
                    cov[i_tr,a_ell,j_tr,b_ell] = c1_l12 * c2_l12 + c3_l12 * c4_l12

            # all factors together
            prefactor = np.multiply(cov[i_tr,:,j_tr,:], prefactor_l)
            cov_w[i_tr,:,j_tr,:] = np.multiply(mw2_matrix, prefactor)
            cov_wt[i_tr,:,j_tr,:] = np.multiply(mwt2_matrix, prefactor)

    # save to fits file
    hdu_w = fits.PrimaryHDU(cov_w)
    hdu_wt = fits.PrimaryHDU(cov_wt)

    hdu_w.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) + \
                    '_nobin_w.fits', overwrite = True)
    hdu_wt.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) +\
                    '_nobin_wt.fits', overwrite = True)


def get_effective_cov():

    '''
    gets final covariance as

    Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]

    This is still unbinned, from ell = 0 to ell = 3 * NSIDE - 1
    '''

    larr_all = np.arange(3 * NSIDE)

    # read-in precomputed covariances
    hdu_dustw = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'dust', 'Cov']) +\
                            '_nobin_w.fits')
    hdu_dustwt = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'dust', 'Cov']) +\
                            '_nobin_wt.fits')
    hdu_all = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'all', 'Cov']) + \
                            '_nobin_w.fits')

    cov_dustw = hdu_dustw[0].data
    cov_dustwt = hdu_dustwt[0].data
    cov_allw = hdu_all[0].data

    assert cov_dustw.shape[1] == len(larr_all), 'error in ell array definition'

    # Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]
    total_cov= np.add(cov_allw, np.subtract(cov_dustwt, cov_dustw))

    # save to fits file
    hdu_cov = fits.PrimaryHDU(total_cov)
    hdu_cov.writeto(PATH_DICT['output_path'] + NAME_RUN + '_nobin_fullCov.fits', overwrite = True)
