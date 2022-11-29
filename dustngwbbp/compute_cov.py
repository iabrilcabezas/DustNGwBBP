'''
compute_cov
    computes covariance matrix as expressed in the Appendix
    computes effective covariance matrix for final power spectrum
'''

import numpy as np
import sacc
from astropy.io import fits
from utils.params import NBANDS, NSIDE, LMIN, DELL, POLARIZATION_cov, PATH_DICT, NAME_RUN
from utils.SED import get_band_names
from utils.bandpowers import get_ell_arrays

band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)

lmax, _, _, _ = get_ell_arrays(LMIN, DELL, NBANDS)

def geo_mean_cl(cl, a, b):

    '''
    geometric mean of power spectrum  at two different ells (a, b)
    '''

    return np.sqrt(cl[a] * cl[b])

def compute_cov(ctype):

    '''
    
    
    '''

    # TODO: assert that cl that comes is only has 1 polarization
    name_pol_cov = 'cl_' + 2 * POLARIZATION_cov.lower()
    # READ IN CL

    s_d = sacc.Sacc.load_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cl', str(NSIDE)]) + '_tot.fits')
    mw2_matrix = np.loadtxt(PATH_DICT['output_path'] + NAME_RUN + '_couplingM_w.txt')
    mwt2_matrix = np.loadtxt(PATH_DICT['output_path'] + NAME_RUN + '_couplingM_wt.txt')

    # read ell of cls from band1 as there'll always be band1
    ellcl, _ = s_d.get_ell_cl(name_pol_cov, 'band1', 'band1', return_cov = False)
    
    nells = len(ellcl)
    assert nells == mw2_matrix.shape[0], 'using matrix or Cl that should not be binned'
   # assert nells == NBANDS), 'some binning on cl went wrong'
    ncombs = len(s_d.get_tracer_combinations())

    assert ncombs == (nmodes*nfreqs) * (nmodes*nfreqs+1) //2,\
             'incorrect number of auto/cross spectra'

    C_ell_all = np.zeros((nfreqs, nfreqs, nells)) + np.nan

    # populate C_ell_all
    for i in range(nfreqs):
        for j in range(i, nfreqs):

            C_ell_all[i][j] = s_d.get_ell_cl(name_pol_cov, f'band{i+1}', f'band{j+1}')[1]

            if i != j:
                C_ell_all[j][i] = s_d.get_ell_cl(name_pol_cov, f'band{i+1}', f'band{j+1}')[1]
    
    Cov = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    Cov_w = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    Cov_wt = np.zeros((ncombs, nells, ncombs, nells)) + np.nan

    indices_tr = np.triu_indices(nfreqs)

    prefactor_l = np.zeros([nells, nells]) + np.nan

    sum_m = 1/ (2 * ellcl + 1)

    for i in range(int(nells)):

        prefactor_l[:, i ] = sum_m[i] * np.ones(nells)

    for ii, (i1,i2) in enumerate(zip(indices_tr[0], indices_tr[1])):    
        print(ii)
        for jj, (j1,j2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        
            c1 = C_ell_all[i1][j1]
            c2 = C_ell_all[i2][j2]
            c3 = C_ell_all[i1][j2]
            c4 = C_ell_all[i2][j1]
            
            for a in range(nells):
                for b in range(nells):
                    c1_l12 = geo_mean_cl(c1, a, b)
                    c2_l12 = geo_mean_cl(c2, a, b)
                    c3_l12 = geo_mean_cl(c3, a, b)
                    c4_l12 = geo_mean_cl(c4, a, b)
                    
                    Cov[ii,a,jj,b] = c1_l12 * c2_l12 + c3_l12 * c4_l12
            
            prefactor = np.multiply(Cov[ii,:,jj,:], prefactor_l)
            Cov_w[ii,:,jj,:] = np.multiply(mw2_matrix, prefactor)
            Cov_wt[ii,:,jj,:] = np.multiply(mwt2_matrix, prefactor)


    hdu_w = fits.PrimaryHDU(Cov_w)        
    hdu_wt = fits.PrimaryHDU(Cov_wt)

    hdul_w = fits.HDUList([hdu_w])
    hdul_wt = fits.HDUList([hdu_wt])

    hdul_w.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) + '_nobin_w.fits', overwrite = True)
    hdul_wt.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) + '_nobin_wt.fits', overwrite = True)


def get_effective_cov():

    '''
    s
    '''
    larr_all = np.arange(3 * NSIDE)

    hdul_dustw = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'dust', 'Cov']) + '_nobin_w.fits')
    hdul_dustwt = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'dust', 'Cov']) + '_nobin_wt.fits')
    hdul_all = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'all', 'Cov']) + '_nobin_w.fits')
    
    Cov_dustw = hdul_dustw[0].data
    Cov_dustwt = hdul_dustwt[0].data
    Cov_allw = hdul_all[0].data

    assert Cov_dustw.shape[1] == len(larr_all), 'error in ell array definition'
    
    total_Cov= np.add(Cov_allw, np.subtract(Cov_dustwt, Cov_dustw))

    hdu_Cov = fits.PrimaryHDU(total_Cov)        
    hdul_Cov = fits.HDUList([hdu_Cov])

    hdul_Cov.writeto(PATH_DICT['output_path'] + NAME_RUN + '_nobin_fullCov.fits', overwrite = True)


