'''

'''

import numpy as np
import yaml
import sacc
from astropy.io import fits
from utils.params import Config
from utils.SED import get_band_names

config = Config(yaml.load(open('config.yaml'), yaml.FullLoader))

EXPERIMENT = config.global_param.experiment
NSIDE = config.global_param.nside
# MACHINE = config.global_param.machine
# LMIN = config.bpw_param.lmin
# DELL = config.bpw_param.dell
# NBANDS = config.bpw_param.nbands
POLARIZATION_cov = config.pol_param.pol_cov

band_names = get_band_names(EXPERIMENT)
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)

def geo_mean_cl(cl, a, b):

    '''
    geometric mean of power spectrum  at two different ells (a, b)
    '''

    return np.sqrt(cl[a] * cl[b])

def compute_cov():

    '''
    
    
    '''

    # TODO: assert that cl that comes is only has 1 polarization
    name_pol_cov = 'cl_' + 2 * POLARIZATION_cov.lower()
    # READ IN CL

    s_d = sacc.Sacc.load_fits()
    mw2_matrix = np.loadtxt()
    mwt2_matrix = np.loadtxt()

    # read ell of cls from band1 as there'll always be band1
    ellcl, _ = s_d.get_ell_cl(name_pol_cov, 'band1', 'band1')

    # bin COV MW matrix:



    nells = len(ellcl)
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

    prefactor_l = np.zeros([NSIDE * 3, NSIDE * 3]) + np.nan

    sum_m = 1/ (2 * ellcl + 1)

    for i in range(int(NSIDE*3)):

        prefactor_l[:, i ] = sum_m[i] * np.ones(NSIDE*3)

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

    hdul_w.writeto()
    hdul_wt.writeto()







