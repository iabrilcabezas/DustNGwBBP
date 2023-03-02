'''
compute_cov

computes covariance from DustFilaments sims
'''

import numpy as np
import pymaster as nmt
import healpy as hp
from astropy.io import fits
from scipy.optimize import curve_fit
from utils_dustfil.params import DF_BASE_PATH, DF_END_NAME, DF_OUTPUT_PATH
from utils_dustfil.params import DF_NAME_RUN
from utils_dustfil.params import DF_FREQ, DF_NSIDE, DF_MASK
from utils_dustfil.params import DF_ALPHA, DF_AMP
from utils_dustfil.params import DF_LMIN, DF_LMAX, DF_DELL
from utils.sed import dl_plaw, get_band_names, get_convolved_seds
from utils.binning import rebin, cut_array
from utils.noise_calc import get_fsky
from utils.namaster import binning_bbpower, get_mask
from utils.params import POLARIZATION_cov, NAME_RUN, PATH_DICT, NSIDE, NAME_COMP
from dustngwbbp.compute_cl import import_bandpasses

b_df, ell_eff_df = binning_bbpower(DF_LMIN, DF_LMAX, DF_DELL, DF_NSIDE)
nell_df = len(ell_eff_df)
mask_so_gal = get_mask(DF_NSIDE, DF_MASK)
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
    cl_store = np.zeros((len(seeds), len(ell_eff_df))) + np.nan

    for i, seed in enumerate(seeds):
        print(i)
        # extract sim
        filament = DF_BASE_PATH  + f'{seed:03}' + DF_END_NAME
        # compute power spectrum on SO patch
        f_2_fil = nmt.NmtField(mask_so_gal, hp.read_map(filament, field = [1,2]))
        cl_22 = nmt.compute_full_master(f_2_fil, f_2_fil, b_df)
        # extract BB spectrum
        cl_store[i] = cl_22[3] # only interested in BB mode

    # save to fits file
    hdu_cl = fits.PrimaryHDU(cl_store)
    hdu_cl.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_RUN}.fits', overwrite = True)


def calibrate_cells():

    '''
    Scales 353 GHz covariance matrix to the other experiment frequencies
    '''

    assert DF_FREQ == 353, 'SED scaling =1 for 353GHz, use it!'

    # open previously computed Cell of maps
    cl_store = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME_RUN}.fits')[0].data
    # calculate <Cell> across sims
    mean_cell = np.mean(cl_store, axis  = 0)

    # START CALIBRATION

    # approximate error of <Cell> with Knox formula
    sigma_cell = np.sqrt(2 /  (( 2 * ell_eff_df + 1) * fsky) ) * mean_cell
    # perform chi2 fit with power-law model
    popt, _ = curve_fit(model_dustfil_dust, ell_eff_df, mean_cell, sigma = sigma_cell )
    # extract overall normalization of power spectrum
    scaling_sims = popt[0]
    # renormalize to chosen amplitude
    cl_store_cal = cl_store / scaling_sims * DF_AMP
    # compute covariance
    cov_cl_sims = np.cov(cl_store_cal, rowvar = False)
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

def merge_cov():

    '''
    Merges covariance from large-scale modulating template (wt) and DustFil. Only dust component.

    # ### It will be the largest absolute value of the two
    TO-DO: CHANGE TO DAVIDS PROPOSAL FOR MERGING
    '''

    # dustwt_nobin = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
    #                         '_nobin_wt.fits')[0].data

    # dustwt_bin = np.zeros([ncombs, nell_df, ncombs, nell_df])

    # for i_tr in range(ncombs):
    #     for j_tr in range(ncombs):
    #         dustwt_bin[i_tr,:, j_tr,:] = rebin(cut_array(dustwt_nobin[i_tr,:, j_tr,:], \
    #                                         np.arange(3 * NSIDE), DF_LMIN, DF_LMAX), [nell_df, nell_df])

    # dust_dustfil = fits.open(DF_OUTPUT_PATH + DF_NAME_RUN + '_dustDFCov.fits')[0].data

    # merged_cov = np.where(np.abs(dust_dustfil) >= dustwt_bin, dust_dustfil, dustwt_bin)

    # save to fits file
   # hducov = fits.PrimaryHDU(merged_cov)
   # hducov.writeto(DF_OUTPUT_PATH + DF_NAME_RUN + '_bin_dust_dfwtCov.fits', overwrite = True)
    ### probably something like DF_OUTPUT_PATH + DF_NAME_RUN + '_d00_Cov_bin_dfwt.fits')
    print('currently invalid function')

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
    hdu_cov = fits.PrimaryHDU(total_cov)
    hdu_cov.writeto(DF_OUTPUT_PATH + DF_NAME_RUN + \
                    '_'.join( [NAME_COMP, 'Cov_bin', type_dustcov]) + '.fits',\
                    overwrite = True)
