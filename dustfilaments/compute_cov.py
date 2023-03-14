'''
compute_cov

computes covariance from DustFilaments sims
'''

import numpy as np
import pymaster as nmt
import healpy as hp
from astropy.io import fits
from scipy.optimize import curve_fit
from utils_dustfil.params import DF_BASE_PATH, DF_OUTPUT_PATH, DF_END_NAME_S, DF_END_NAME_A
from utils_dustfil.params import DF_NAME_RUN, DF_NAME_SIM
from utils_dustfil.params import DF_FREQ, DF_NSIDE, DF_MASK
from utils_dustfil.params import DF_ALPHA, DF_AMP
from utils_dustfil.params import DF_LMIN, DF_LMAX, DF_DELL
from utils.sed import dl_plaw, get_band_names, get_convolved_seds
from utils.binning import rebin, cut_array
from utils.noise_calc import get_fsky
from utils.namaster import binning_bbpower, get_mask
from utils.params import POLARIZATION_cov, NAME_RUN, PATH_DICT, NSIDE, NAME_COMP
from utils.params import name_couplingmatrix_wt, name_couplingmatrix_w, COV_CORR, config
from dustngwbbp.compute_cl import import_bandpasses
from dustngwbbp.compute_couplingmatrix import compute_couplingmatrix

b_df, ell_eff_df = binning_bbpower(DF_LMIN, DF_LMAX, DF_DELL, DF_NSIDE)
nell_df = len(ell_eff_df)
mask_so_gal = get_mask(DF_NSIDE, DF_MASK, **config.mask_param.__dict__)
fsky = get_fsky()

band_names = get_band_names()

bpss       = import_bandpasses()

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncombs = len(indices_tr[0])

def model_dustfil_dust(x, Ad):

    '''
    PL spectra for dust at fixed alpha. 
    Note that we work with C_ell instead of D_ell
    '''
    dell_dust = dl_plaw(A=Ad, alpha=DF_ALPHA, ls = x)
    cell_dust = dell_dust / ( (x * (x + 1)) / (2 * np.pi) )
    return cell_dust



def compute_cell_dustfil(nseeds):

    '''
    Measures power spectrum of all DF simulations in the bandpowers of BBPower

    ** Parameters **
    nseeds: int
        number of simulations available
    '''

    # number of simulations available
    seeds = np.arange(nseeds)

    # initialize array
    cl_store_all = np.zeros((len(seeds), len(ell_eff_df))) + np.nan
    cl_store_small = np.zeros((len(seeds), len(ell_eff_df))) + np.nan

    for i, seed in enumerate(seeds):
        print(i)
        # extract sim
        filament_s = DF_BASE_PATH  + f'{seed:03}' + DF_END_NAME_S
        filament_a = DF_BASE_PATH  + f'{seed:03}' + DF_END_NAME_A
        # compute power spectrum on SO patch
        f_2_fil_s = nmt.NmtField(mask_so_gal, hp.read_map(filament_s, field = [1,2]))
        cl_22_s = nmt.compute_full_master(f_2_fil_s, f_2_fil_s, b_df)
        # extract BB spectrum
        cl_store_small[i] = cl_22_s[3] # only interested in BB mode

        f_2_fil_a = nmt.NmtField(mask_so_gal, hp.read_map(filament_a, field = [1,2]))
        cl_22_a = nmt.compute_full_master(f_2_fil_a, f_2_fil_a, b_df)
        # extract BB spectrum
        cl_store_all[i] = cl_22_a[3] # only interested in BB mode


    # save to fits file
    hdu_cl_s= fits.PrimaryHDU(cl_store_small)
    hdu_cl_s.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_small.fits', overwrite = True)

    hdu_cl_a = fits.PrimaryHDU(cl_store_all)
    hdu_cl_a.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_all.fits', overwrite = True)


def calibrate_cells():

    '''
    Scales 353 GHz covariance matrix to match dust model amplitude
    '''

    assert DF_FREQ == 353, 'SED scaling =1 for 353GHz, use it!'

    # open previously computed Cell of maps
    cl_all   = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_all.fits')[0].data
    cl_small = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_small.fits')[0].data
    # calculate <Cell> across sims
    mean_cell = np.mean(cl_all, axis  = 0)

    # START CALIBRATION

    # approximate error of <Cell> with Knox formula
    sigma_cell = np.sqrt(2 /  (( 2 * ell_eff_df + 1) * fsky * DF_DELL) ) * mean_cell
    # perform chi2 fit with power-law model
    popt, _ = curve_fit(model_dustfil_dust, ell_eff_df, mean_cell, sigma = sigma_cell )
    # extract overall normalization of power spectrum
    scaling_sims = popt[0]
    # renormalize to chosen amplitude
    cl_all = cl_all / scaling_sims * DF_AMP
    cl_small = cl_small / scaling_sims * DF_AMP

    # save to fits file
    hdu_cl_s= fits.PrimaryHDU(cl_small)
    hdu_cl_s.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_small_cal.fits', overwrite = True)

    hdu_cl_a = fits.PrimaryHDU(cl_all)
    hdu_cl_a.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_all_cal.fits', overwrite = True)

def compute_cov_fromsims(scale = 'small'):

    '''
    computes covariance with np.cov()
    '''

    assert scale == 'small', 'only small sim has true covariance'

    cl_small = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_{scale}_cal.fits')[0].data

    # compute covariance
    cov_cl_sims = np.cov(cl_small, rowvar = False)
    # extract SED scalings
    seds_dd = get_convolved_seds(band_names, bpss)[1] # dust component

    # obtain cov at all frequencies (assumes no decorrelation)
    cov_dustfil_scale_full = np.zeros((ncombs, nell_df, ncombs, nell_df )) + np.nan

    # populate covariance
    for i_tr, (i_tr1,i_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        for j_tr, (j_tr1,j_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
            # covariance between different frequency channels
            sed_scaling = seds_dd[i_tr1] * seds_dd[i_tr2] * seds_dd[j_tr1] * seds_dd[j_tr2]
            cov_dustfil_scale_full[i_tr,:,j_tr, :] = cov_cl_sims * sed_scaling
    # save to fits file
    hdu_cov = fits.PrimaryHDU(cov_dustfil_scale_full)
    hdu_cov.writeto(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_df00.fits', overwrite = True)

def compute_tildecov(scale):

    '''
    Computes \tilde{Cov}, given by eqn A20 with mean cell from sims
    '''

    cl_f0_sims = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_SIM}_{scale}_cal.fits')[0].data
    cl_f0      = np.mean(cl_f0_sims, axis  = 0) # mean across sims

    if scale == 'small':
        m2_matrix_nobin = np.loadtxt(name_couplingmatrix_w + '.txt')
    elif scale == 'all':
        m2_matrix_nobin = np.loadtxt(name_couplingmatrix_wt + '.txt')

    w2_mean = compute_couplingmatrix(**config.mask_param.__dict__)

    m2_matrix = rebin(cut_array(m2_matrix_nobin, \
                                  np.arange(3 * NSIDE), DF_LMIN, DF_LMAX), [nell_df, nell_df])

    c_ell_freqs = np.zeros((nfreqs, nfreqs, nell_df)) + np.nan
    # extract SED scalings
    seds_dd = get_convolved_seds(band_names, bpss)[1] # dust component

    # populate C_ell_all
    for i in range(nfreqs):
        for j in range(i, nfreqs):

            c_ell_freqs[i][j] = cl_f0 * seds_dd[i] * seds_dd[j]

            if i!=j:
                c_ell_freqs[j][i] = cl_f0 * seds_dd[i] * seds_dd[j]

    cov    = np.zeros((ncombs, nell_df, ncombs, nell_df)) + np.nan
    cov_wt = np.zeros((ncombs, nell_df, ncombs, nell_df)) + np.nan

    # compute prefactor of covariance
    prefactor_l = np.zeros([nell_df, nell_df]) + np.nan

    sum_m_f = 1 / (2 * ell_eff_df + 1)

    # NO COV_CORR CORRECTION?
    assert COV_CORR == 'w2', 'only implemented for w2 correction'
    if COV_CORR == 'w2':
        sum_m_f /= ( w2_mean**2 )

    for i in range(int(nell_df)):
        prefactor_l[:, i ] = sum_m_f[i] * np.ones(nell_df)

    # populate covariance
    for i_tr, (i_tr1,i_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        print(i_tr)
        for j_tr, (j_tr1,j_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):

            cl1 = c_ell_freqs[i_tr1][j_tr1]
            cl2 = c_ell_freqs[i_tr2][j_tr2]
            cl3 = c_ell_freqs[i_tr1][j_tr2]
            cl4 = c_ell_freqs[i_tr2][j_tr1]

            for a_ell in range(nell_df):
                for b_ell in range(nell_df):
                    # calculate mean between cl at different ells
                    c1_l12 = np.sqrt(cl1[a_ell] * cl1[b_ell])
                    c2_l12 = np.sqrt(cl2[a_ell] * cl2[b_ell])
                    c3_l12 = np.sqrt(cl3[a_ell] * cl3[b_ell])
                    c4_l12 = np.sqrt(cl4[a_ell] * cl4[b_ell])
                    # covariance between different frequency channels
                    cov[i_tr,a_ell,j_tr,b_ell] = c1_l12 * c2_l12 + c3_l12 * c4_l12

            # all factors together
            prefactor = np.multiply(cov[i_tr,:,j_tr,:], prefactor_l)
            cov_wt[i_tr,:,j_tr,:] = np.multiply(m2_matrix, prefactor)

    # save to fits file
    hdu_wt = fits.PrimaryHDU(cov_wt)
    hdu_wt.writeto(DF_OUTPUT_PATH + DF_NAME_RUN + f'_d00_Cov_bin_dft{scale[0]}.fits', overwrite = True)


def merge_cov(method):

    '''
    Merges covariance from large-scale modulating template (wt) and DustFil. Only dust component.

    # Blake's method is: take largest absolute value of the two
    # David: take DF sims without large scale template. add correction (Cov(tilde)) for large scales
    '''

    if method == 'blake':

        dustwt_nobin = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
                                '_nobin_wt.fits')[0].data

        dustwt_bin = np.zeros([ncombs, nell_df, ncombs, nell_df])

        for i_tr in range(ncombs):
            for j_tr in range(ncombs):
                dustwt_bin[i_tr,:, j_tr,:] = rebin(cut_array(dustwt_nobin[i_tr,:, j_tr,:], \
                                                np.arange(3 * NSIDE), DF_LMIN, DF_LMAX), [nell_df, nell_df])

        dust_dustfil = fits.open(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_df00.fits')[0].data

        merged_cov = np.where(np.abs(dust_dustfil) >= dustwt_bin, dust_dustfil, dustwt_bin)

    elif method == 'david':

        cov_s  = fits.open(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_df00.fits')[0].data
        covt_s = fits.open(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_dfts.fits')[0].data
        covt_a = fits.open(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_dfta.fits')[0].data

        merged_cov = np.add(np.subtract(covt_a, covt_s), cov_s)

    # save to fits file
    hducov = fits.PrimaryHDU(merged_cov)
    hducov.writeto(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_dfwt.fits', overwrite = True)

    return None

def compute_full_cov(type_dustcov):

    '''
    Merges dust covariance with covariance from other components
    '''

    cov_dustall = fits.open(DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_' +\
                            f'{type_dustcov}' + '.fits')[0].data

    cov_allw    = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, NAME_COMP, 'Cov']) + \
                            '_nobin_w.fits')[0].data
    # read in precomputed gaussian cov
    cov_dustw   = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
                            '_nobin_w.fits')[0].data

    cov_allw_bin    = np.zeros([ncombs, nell_df, ncombs, nell_df])
    cov_dustw_bin   = np.zeros([ncombs, nell_df, ncombs, nell_df])
    for i_tr in range(ncombs):
        for j_tr in range(ncombs):
            cov_allw_bin[i_tr,:, j_tr,:]  = rebin(cut_array(cov_allw[i_tr,:, j_tr,:], \
                                            np.arange(3 * NSIDE), DF_LMIN, DF_LMAX), [nell_df, nell_df])
            cov_dustw_bin[i_tr,:, j_tr,:] = rebin(cut_array(cov_dustw[i_tr,:, j_tr,:], \
                                            np.arange(3 * NSIDE), DF_LMIN, DF_LMAX), [nell_df, nell_df])


    # Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]
    total_cov = np.add(cov_allw_bin, np.subtract(cov_dustall, cov_dustw_bin))
    # save to fits file
    assert NAME_COMP == 'dcs', 'you are computing a full cov with nothing but dust'
    hdu_cov = fits.PrimaryHDU(total_cov)
    hdu_cov.writeto(PATH_DICT['output_path'] + '_'.join([ NAME_RUN, NAME_COMP, 'Cov_bin', type_dustcov]) + '.fits',\
                    overwrite = True)
