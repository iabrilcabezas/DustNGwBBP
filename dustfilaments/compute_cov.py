'''
compute_cov

computes covariance from DustFilaments sims
'''

import numpy as np
import pymaster as nmt
import healpy as hp
from astropy.io import fits
from scipy.optimize import curve_fit
from pymaster.workspaces import compute_coupled_cell
from dustngwbbp.compute_cl import import_bandpasses
from dustngwbbp.compute_couplingmatrix import compute_couplingmatrix
from dustngwbbp.compute_cov import get_cov_fromeq
from utils.sed import dl_plaw, get_band_names, get_convolved_seds
from utils.binning import rebin, cut_array, bin_covs
from utils.noise_calc import get_fsky
from utils.namaster import binning_bbpower, get_mask
from utils.params import POLARIZATION_cov, NAME_RUN, PATH_DICT, NAME_COMP
from utils.params import name_couplingmatrix_wt, name_couplingmatrix_w, COV_CORR, config
from utils.params import DF_BASE_PATH, DF_OUTPUT_PATH, DF_END_NAME_A, DF_END_NAME_S
from utils.params import DF_NAME
from utils.params import nu0_dust as DF_FREQ
from utils.params import NSIDE, MTYPE
from utils.params import alpha_dust_BB as DF_ALPHA
from utils.params import A_dust_BB as DF_AMP
from utils.params import LMIN, DELL, NBANDS
from utils.bandpowers import get_ell_arrays

LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)
b_df, ell_eff_df = binning_bbpower(LMIN, LMAX, DELL, NSIDE)
mask_so_gal = get_mask(NSIDE, MTYPE, **config.mask_param.__dict__)
w2_mean = compute_couplingmatrix(**config.mask_param.__dict__)
fsky = get_fsky()

band_names = get_band_names()
bpss       = import_bandpasses()
# extract SED scalings
seds_dd = get_convolved_seds(band_names, bpss)[1] # dust component


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

def cell0_tofreqs(cell_array, seds):

    '''
    take Cell at f0 and extrapolate to other frequencies
    '''

    nfreq = len(seds)
    nell   = len(cell_array)

    cell_freq = np.zeros([nfreq, nfreq, nell]) + np.nan

    for i in range(nfreq):
        for j in range(nfreq):
            cell_freq[i][j] = cell_array * seds[i] * seds[j]

    for i in range(nfreq):
        for j in range(i,nfreq):
            assert np.all(np.isclose(cell_freq[i][j], cell_freq[i][j].T)), \
                'cell(f1,f2) != cell(f2,f1)'

    return cell_freq


def compute_cell_bin_dustfil(nseeds, nside):

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
        maps_s = hp.read_map(filament_s, field = [1,2])
        maps_a = hp.read_map(filament_a, field = [1,2])

        maps_s = hp.ud_grade(maps_s, nside)
        maps_a = hp.ud_grade(maps_a, nside)

        f_2_fil_s = nmt.NmtField(mask_so_gal, maps_s)
        cl_22_s = nmt.compute_full_master(f_2_fil_s, f_2_fil_s, b_df)
        # extract BB spectrum
        cl_store_small[i] = cl_22_s[3] # only interested in BB mode

        f_2_fil_a = nmt.NmtField(mask_so_gal, maps_a)
        cl_22_a = nmt.compute_full_master(f_2_fil_a, f_2_fil_a, b_df)
        # extract BB spectrum
        cl_store_all[i] = cl_22_a[3] # only interested in BB mode

    # save to fits file
    hdu_cl_s= fits.PrimaryHDU(cl_store_small)
    hdu_cl_s.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME}_bin_small.fits', overwrite = True)

    hdu_cl_a = fits.PrimaryHDU(cl_store_all)
    hdu_cl_a.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME}_bin_all.fits', overwrite = True)

def compute_cell_nobin_dustfil(nseeds, nside):

    '''
    procedure from alonso. computes cell of masked map at each ell
    '''

    # number of simulations available
    seeds = np.arange(nseeds)

    # initialize array
    cl_store_all = np.zeros((nseeds, 3 * nside)) + np.nan
    cl_store_small = np.zeros((nseeds, 3 * nside)) + np.nan

    for i, seed in enumerate(seeds):
        print(i)
        # extract sim
        filament_s = DF_BASE_PATH  + f'{seed:03}' + DF_END_NAME_S
        filament_a = DF_BASE_PATH  + f'{seed:03}' + DF_END_NAME_A
        # compute power spectrum on SO patch
        maps_s = hp.read_map(filament_s, field = [1,2])
        maps_a = hp.read_map(filament_a, field = [1,2])

        maps_s = hp.ud_grade(maps_s, nside)
        maps_a = hp.ud_grade(maps_a, nside)

        f1_s = nmt.NmtField(np.ones(hp.nside2npix(nside)), mask_so_gal * maps_s)
        f1_a = nmt.NmtField(np.ones(hp.nside2npix(nside)), mask_so_gal * maps_a)

        cl_store_small[i] = compute_coupled_cell(f1_s,f1_s)[3] / w2_mean
        cl_store_all[i]   = compute_coupled_cell(f1_a, f1_a)[3] / w2_mean

    # save to fits file
    hdu_cl_s= fits.PrimaryHDU(cl_store_small)
    hdu_cl_s.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME}_nobin_small.fits', overwrite = True)

    hdu_cl_a = fits.PrimaryHDU(cl_store_all)
    hdu_cl_a.writeto(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME}_nobin_all.fits', overwrite = True)


def calibrate_cells(bin_type):

    '''
    Scales 353 GHz covariance matrix to match dust model amplitude
    '''

    assert int(DF_FREQ) == 353, 'SED scaling =1 for 353GHz, use it!'

    # open previously computed Cell of maps
    cl_all   = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME}_{bin_type}_all.fits')[0].data
    cl_small = fits.open(DF_OUTPUT_PATH + f'cl_BB_{DF_NAME}_{bin_type}_small.fits')[0].data
    # calculate <Cell> across sims
    mean_cell = np.mean(cl_all, axis  = 0)

    # START CALIBRATION
    if bin_type == 'nobin':
        ell_array = np.arange(len(mean_cell))
        delta_ell = 1.

    elif bin_type == 'bin':
        ell_array = ell_eff_df
        delta_ell = DELL

    mask = (ell_array < 300) & (ell_array > 30)

    # approximate error of <Cell> with Knox formula
    sigma_cell = np.sqrt(2 /  (( 2 * ell_array + 1) * fsky * delta_ell) ) * mean_cell
    # perform chi2 fit with power-law model
    popt = curve_fit(model_dustfil_dust, \
                     ell_array[mask], mean_cell[mask], sigma = sigma_cell[mask] )[0]
    # extract overall normalization of power spectrum
    scaling_sims = popt[0]
    # renormalize to chosen amplitude
    cl_all = cl_all / scaling_sims * DF_AMP
    cl_small = cl_small / scaling_sims * DF_AMP

    # save to fits file
    hdu_cl_s= fits.PrimaryHDU(cl_small)
    hdu_cl_s.writeto(DF_OUTPUT_PATH + \
                     f'cl_BB_{DF_NAME}_{bin_type}_small_cal.fits', overwrite = True)

    hdu_cl_a = fits.PrimaryHDU(cl_all)
    hdu_cl_a.writeto(DF_OUTPUT_PATH + \
                     f'cl_BB_{DF_NAME}_{bin_type}_all_cal.fits', overwrite = True)

def compute_cov_fromsims(scale = 'small', bin_type = 'bin'):

    '''
    computes covariance with np.cov()
    '''

    assert scale == 'small', 'only small sim has true covariance'

    cl_small = fits.open(DF_OUTPUT_PATH +\
                         f'cl_BB_{DF_NAME}_{bin_type}_{scale}_cal.fits')[0].data

    # compute covariance
    cov_cl_sims = np.cov(cl_small, rowvar = False)

    # obtain cov at all frequencies (assumes no decorrelation)
    nel = cov_cl_sims.shape[0]
    cov_dustfil_scale_full = np.zeros((ncombs, nel, ncombs, nel )) + np.nan

    # populate covariance
    for i_tr, (i_tr1,i_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        for j_tr, (j_tr1,j_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
            # covariance between different frequency channels
            sed_scaling = seds_dd[i_tr1] * seds_dd[i_tr2] * seds_dd[j_tr1] * seds_dd[j_tr2]
            cov_dustfil_scale_full[i_tr,:,j_tr, :] = cov_cl_sims * sed_scaling
    # save to fits file
    hdu_cov = fits.PrimaryHDU(cov_dustfil_scale_full)
    hdu_cov.writeto(DF_OUTPUT_PATH + DF_NAME +\
                    f'_d00_Cov_{bin_type}_df00.fits', overwrite = True)

def compute_tildecov(scale, bin_type):

    '''
    Computes \\tilde{Cov}, given by eqn A20 with mean cell from sims
    '''
    assert bin_type == 'nobin', 'correct expression is with binning after computing cov'

    cl_f0_sims = fits.open(DF_OUTPUT_PATH + \
                           f'cl_BB_{DF_NAME}_{bin_type}_{scale}_cal.fits')[0].data
    cl_f0      = np.mean(cl_f0_sims, axis  = 0) # mean across sims

    if scale == 'small':
        mcm = np.loadtxt(name_couplingmatrix_w + '.txt')
    elif scale == 'all':
        mcm = np.loadtxt(name_couplingmatrix_wt + '.txt')

    c_ell_freqs = cell0_tofreqs(cl_f0, seds_dd)
    ell_array   = np.arange(c_ell_freqs.shape[-1]) # works for nobin, otherwise rescribed below

    if bin_type == 'bin':
        mcm = rebin(cut_array(mcm, np.arange(3 * NSIDE), LMIN, LMAX), [NBANDS, NBANDS])
        ell_array = ell_eff_df

    if bin_type == 'nobin':
        assert len(ell_array) == 3 * NSIDE

    assert COV_CORR == 'w2', 'only implemented for w2 correction'
    cov_tilde = get_cov_fromeq(ell_array, c_ell_freqs, mcm, w2_mean**2)
    hdu_wt_nobin = fits.PrimaryHDU(cov_tilde)
    hdu_wt_nobin.writeto(DF_OUTPUT_PATH + DF_NAME +\
                         f'_d00_Cov_nobin_dft{scale[0]}.fits', overwrite = True)

    # save also binned results:
    cov_tilde_bin = bin_covs(cov_tilde)
    # save to fits file
    hdu_wt = fits.PrimaryHDU(cov_tilde_bin)
    hdu_wt.writeto(DF_OUTPUT_PATH + DF_NAME + \
                   f'_d00_Cov_bin_dft{scale[0]}.fits', overwrite = True)


def merge_cov(type_cov, bin_type = 'bin'):

    '''
    Merges covariance from large-scale modulating template (wt) and DustFil. Only dust component.

    # Blake's method is: take largest absolute value of the two
    # David: take DF sims without large scale template. add correction (Cov(tilde)) for large scales
    '''

    if type_cov == 'dfwtm':

        dustwt = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
                                f'_{bin_type}_wt.fits')[0].data

        dust_df = fits.open(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_{bin_type}_df00.fits')[0].data

        merged_cov = np.where(np.abs(dust_df) >= dustwt, dust_df, dustwt)

    elif type_cov == 'dfwt':

        cov_s  = fits.open(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_{bin_type}_df00.fits')[0].data
        covt_s = fits.open(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_{bin_type}_dfts.fits')[0].data
        covt_a = fits.open(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_{bin_type}_dfta.fits')[0].data

        merged_cov = np.add(np.subtract(covt_a, covt_s), cov_s)

    else:
        print('unresolved type_cov, check naming')
        return None

    # save to fits file
    hducov = fits.PrimaryHDU(merged_cov)

    hducov.writeto(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_{bin_type}_{type_cov}.fits', overwrite = True)

    if bin_type == 'nobin':

        bin_merged_cov = bin_covs(merged_cov)
        hducov_bin = fits.PrimaryHDU(bin_merged_cov)
        hducov_bin.writeto(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_bin_{type_cov}.fits', overwrite = True)


def compute_full_cov(type_dustcov, bin_type = 'bin'):

    '''
    Merges dust covariance with covariance from other components
    '''

    cov_dustall = fits.open(DF_OUTPUT_PATH + DF_NAME + f'_d00_Cov_{bin_type}_' +\
                            f'{type_dustcov}' + '.fits')[0].data

    cov_allw    = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, NAME_COMP, 'Cov']) + \
                            f'_{bin_type}_w.fits')[0].data
    # read in precomputed gaussian cov
    cov_dustw   = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
                            f'_{bin_type}_w.fits')[0].data

    # Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]
    total_cov = np.add(cov_allw, np.subtract(cov_dustall, cov_dustw))

    # save to fits file
    assert NAME_COMP == 'dcs', 'you are computing a full cov with nothing but dust'
    hdu_cov = fits.PrimaryHDU(total_cov)
    hdu_cov.writeto(PATH_DICT['output_path'] + \
                    '_'.join([ NAME_RUN, NAME_COMP, f'Cov_{bin_type}', type_dustcov]) + '.fits',\
                    overwrite = True)

    if bin_type == 'nobin':
        bin_total_cov = bin_covs(total_cov)
        hducov_bin = fits.PrimaryHDU(bin_total_cov)
        hducov_bin.writeto(PATH_DICT['output_path'] + \
                           '_'.join([ NAME_RUN, NAME_COMP, 'Cov_bin', type_dustcov]) + '.fits',\
                    overwrite = True)
