'''
compute_cl
    computes c_ell to input into covariance calculation
    creates sacc files that contain spectrum and covariances

'''

from copy import deepcopy
import numpy as np
import sacc
from astropy.io import fits
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP
from utils.params import EXPERIMENT, NSIDE, LMIN, DELL, NBANDS, POLARIZATION_cov
import utils.noise_calc as nc
from utils.binning import rebin, cut_array
from utils.sed import get_band_names, Bpass, get_component_spectra, get_convolved_seds
from utils.bandpowers import get_ell_arrays, dell2cell_lmax

from utils_dustfil.params import DF_OUTPUT_PATH, DF_NAME_RUN

band_names = get_band_names()
LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)

ctype_dict_ncomp = {'d00': 1, 'dc0': 2, 'dcs': 3}
bincov_dict      = {'wt': 'nobin', 'w': 'nobin', 'df00': 'bin', 'dfwt': 'bin'}

def import_bandpasses():

    '''
    Returns dictionary of bandpass class for each channel in EXPERIMENT
    '''

    if EXPERIMENT == 'bicep':
        bpss = {n: Bpass(n, PATH_DICT['BK15_data'] +\
                             f'BK15_{n}_bandpass_20180920.txt') for n in band_names}
    if EXPERIMENT in ['so','cmbs4','cmbs4d']:
        bpss = {n: Bpass(n,PATH_DICT['bbpipe_path'] +\
                             f'examples/data/bandpasses/{n}.txt') for n in band_names}

    return bpss

def import_beams(ell_array):

    '''
    Returns dictionary of beams for each channel in EXPERIMENT, evaluated at each ell in ell_array

    ** Parameters **
    ell_array: np.array()
        array where to evaluate the beams
    '''

    if EXPERIMENT == 'bicep':
        beams ={band_names[i]: b for i, b in \
                    enumerate(nc.bicep_beams_exp(ell_array))}
    if EXPERIMENT in ['so','cmbs4','cmbs4d']:
        beams ={band_names[i]: b for i, b in \
                    enumerate(nc.Simons_Observatory_V3_SA_beams(ell_array))}

    return beams

def get_windows(weight):

    '''
    Returns window with binning according to weights

    ** Parameters **
    weight: 'Dl' or 'Cl'
        weight as Dell [l * (l + 1)/ 2pi] or 'Cell' [equal weights]

    ** Returns **
    windows: np.array([NBANDS, LMAX + 1])
    '''

    weight_types = ['Dl', 'Cl']
    assert weight in weight_types, 'not a type of weight!'

    windows = np.zeros([NBANDS,LMAX+1])

    cl_weights = np.ones_like(LARR_ALL)

    for b_i,(b_l0,b_lf) in enumerate(zip(LBANDS[:-1],LBANDS[1:])):

        if weight == 'Dl':
            windows[b_i,b_l0:b_lf] = (LARR_ALL * (LARR_ALL + 1)/(2*np.pi))[b_l0:b_lf]
        if weight == 'Cl':
            windows[b_i, b_l0:b_lf] = cl_weights[b_l0:b_lf]

        windows[b_i,:] /= DELL

    return windows

def add_tracers(ell_array):

    '''
    Creates sacc object and add tracers according to EXPERIMENT
    Beam of tracer is evaluated at ell = ell_array

    ** Parameters **
    ell_array: np.array()
        ell array where beams are calculated

    ** Returns **
    s_d: sacc.Sacc()
        sacc object with added tracers
    '''

    s_d = sacc.Sacc()

     # Bandpasses:
    bpss = import_bandpasses()
    # Beams
    beams = import_beams(ell_array)

    for i_band, name_b in enumerate(band_names):
        bandpass = bpss[name_b]
        beam = beams[name_b]
        s_d.add_tracer('NuMap', f'band{i_band+1}',
                        quantity='cmb_polarization',
                        spin=2,
                        nu=bandpass.nu,
                        bandpass=bandpass.bnu,
                        ell=ell_array,
                        beam=beam,
                        nu_unit='GHz',
                        map_unit='uK_CMB')

    return s_d

def add_powerspectra(s_d, bpw_freq_sig, leff, do_bin, weight = 'Cl'):

    '''
    Adds power spectra to Sacc object

    ** Parameters **
    s_d: sacc object
        object to add P(k)
    bpw_freq_sig: np.array
        power spectra that will be added
    leff: np.array
        ell each power spectra band corresponds to
    do_bin: bool
        if the power spectra is binned, provide window too
    weight: 'Cl' or 'Dl'
        binnin of window if do_bin
    '''

    if do_bin:
        windows = get_windows(weight)
        s_wins = sacc.BandpowerWindow(LARR_ALL, windows.T)

    map_names=[]
    for i_freq in range(nfreqs):
        if 'E' in POLARIZATION_cov:
            map_names.append(f'band{i_freq+1}_E')
        if 'B' in POLARIZATION_cov:
            map_names.append(f'band{i_freq+1}_B')

    for (i_tr1, i_tr2) in zip(indices_tr[0], indices_tr[1]):
        band12 = np.array([map_names[i_tr1][:-2], map_names[i_tr2][:-2]])
        pol12  = np.array([map_names[i_tr1][-1].lower(), map_names[i_tr2][-1].lower()])
        cl_type = f'cl_{pol12[0]}{pol12[1]}'

        if do_bin:
            s_d.add_ell_cl(cl_type, band12[0], band12[1], LEFF, bpw_freq_sig[i_tr1, i_tr2, :],
                         window = s_wins)
        else:
            s_d.add_ell_cl(cl_type, band12[0], band12[1], leff, bpw_freq_sig[i_tr1, i_tr2, :]) # note differences in leff


    return s_d


def add_noise(lmax, bpw_freq_sig, fsky, do_bin, weight = 'Cl'):

    '''
    Returns noise array according to experiment

    ** Parameters **
    lmax: int
        maximum ell to compute noise
    bpw_freq_sig: np.array
        provides size of noise array
    fsky:
        sky fraction of experiment (noise \propto 1/fsky)
    do_bin: bool
        if the noise should be binned
    weight: 'Cl' or 'Dl'
        binnin of window if do_bin
    '''

    if do_bin:
        windows = get_windows(weight)

    if EXPERIMENT == 'bicep':

        nell = np.zeros([nfreqs, lmax + 1])
        _, nell[:, 2:] = nc.bicep_noise(lmax)


    if EXPERIMENT == 'so':

        # N_ell
        sens=2
        knee=1
        ylf=1
        nell=np.zeros([nfreqs,lmax+1])
        _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,lmax+1,1)

    if EXPERIMENT in ['cmbs4', 'cmbs4d']:

        # N_ell
        sens=2
        knee=1
        ylf=1
        nell=np.zeros([nfreqs,lmax+1])
        _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,lmax+1,1)
        nell /= 5

    if do_bin:
        n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)
        noise_array = n_bpw
    else:
        noise_array = nell

    bpw_freq_noi=np.zeros_like(bpw_freq_sig)

    for ib in range(len(noise_array)):
        bpw_freq_noi[ib,0,ib,0,:]=noise_array[ib,:]
        if nmodes == 2:
            bpw_freq_noi[ib,1,ib,1,:]=noise_array[ib,:]

    return bpw_freq_noi


def get_bpw_freq_sig(ctype, lmax, do_bin, weight = 'Cl'):

    '''
    Computes SED of all components
    Convolves SED with instrument (bands, CMB units)
    Computes total signal in each bandpower

    ** Parameters **
    ctype: 'd00', 'dc0', 'dcs'
        type of model (dust only, dust + CMB, dust + CMB + sync)
    lmax: int
        max ell to compute SED
    do_bin: bool
        bin signal with window?
    weight: 'Cl' or 'Dl'
        binning of window
    '''

    assert POLARIZATION_cov == 'B', 'reading B components but you have specified otherwise'

    bpss = import_bandpasses()

    dl2cl = dell2cell_lmax(lmax)

    ncomp = ctype_dict_ncomp[ctype]

    dls_comp = np.zeros([ncomp,nmodes,ncomp,nmodes,lmax+1]) #[ncomp,np,ncomp,np,nl]
    dls_sed  = get_component_spectra(lmax)

    if ctype == 'd00':
        dls_comp[0,0,0,0,:] = dls_sed[1]

    if ctype == 'dc0':
        (dls_comp[1,0,1,0,:],
        dls_comp[0,0,0,0,:]) = dls_sed[1], dls_sed[3]

    if ctype == 'dcs':
        (dls_comp[1,0,1,0,:],
        dls_comp[0,0,0,0,:],
        dls_comp[2,0,2,0,:]) = dls_sed[1], dls_sed[3], dls_sed[5]

    dls_comp *= dl2cl[None, None, None, None, :]

    if do_bin:
        windows = get_windows(weight)
        bpw_comp=np.sum(dls_comp[:,:,:,:,None,:]*windows[None,None,None, None, :,:],axis=5)
        to_bpw_freq = bpw_comp
    else:
        to_bpw_freq = dls_comp

    # Convolve with bandpasses
    seds = get_convolved_seds(band_names, bpss)

    if ncomp == 1:
        seds = np.array([seds[1,:]])
    elif ncomp == 2:
        seds = seds[:2, :]

    bpw_freq_sig = np.einsum('ik,jm,iljno', seds, seds, to_bpw_freq)

    return bpw_freq_sig



def compute_cl(ctype, type_cov):

    '''
    Computes full Cl fits files for BBCompSep (signal [fiducial], noise, and total)

    ** Parameters **
    ctype: 'dc0', 'dcs'
        components used in model
    type_cov: 'w' or 'wt'
        w: gaussian covariance
        wt: adds non gaussianity with modulating template
    '''

    weight = 'Cl'
    ncomp = ctype_dict_ncomp[ctype]

    # compute power spectra
    bpw_freq_sig = get_bpw_freq_sig(ctype, LMAX, True, weight)

    # add noise
    fsky = nc.get_fsky()

    bpw_freq_tot = deepcopy(bpw_freq_sig)
    bpw_freq_noi = np.zeros_like(bpw_freq_sig)

    if ncomp > 1:
        ## add noise
        bpw_freq_noi = add_noise(LMAX, bpw_freq_sig, fsky, True, weight)
        bpw_freq_tot += bpw_freq_noi

    # correct format
    bpw_freq_sig = bpw_freq_sig.reshape([nfreqs*nmodes,nfreqs*nmodes, NBANDS])
    bpw_freq_tot = bpw_freq_tot.reshape([nfreqs*nmodes,nfreqs*nmodes, NBANDS])
    bpw_freq_noi = bpw_freq_noi.reshape([nfreqs*nmodes,nfreqs*nmodes, NBANDS])

    # Create sacc and add tracers
    print("Adding tracers")
    s_d = add_tracers(LARR_ALL)
    s_f = add_tracers(LARR_ALL)
    s_n = add_tracers(LARR_ALL)
    # Adding power spectra
    print("Adding spectra")
    s_d = add_powerspectra(s_d, bpw_freq_sig, LEFF, True, weight)
    s_f = add_powerspectra(s_f, bpw_freq_sig, LEFF, True, weight)
    s_n = add_powerspectra(s_n, bpw_freq_noi, LEFF, True, weight)

    # add covariance
    assert (ctype != 'd00'), 'adding cov to dust only?'
    ncombs = len(s_d.get_tracer_combinations())
    assert ncombs == len(indices_tr[0]), 'how many maps do you have?'

    cov_bpw = np.zeros([ncombs, NBANDS, ncombs, NBANDS])

    if type_cov == 'dfwt' or type_cov == 'df00':
        # read in cov, already binned
        cov_bpw         = fits.open(PATH_DICT['output_path'] + '_'.join([ NAME_RUN, NAME_COMP]) + '_' + \
                                '_'.join(['Cov', 'bin', type_cov]) + '.fits')[0].data

    elif type_cov == 'w' or type_cov == 'wt':
        cov_bpw_nobin   = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, NAME_COMP]) + '_' + \
                                '_'.join(['Cov', 'nobin', type_cov]) + '.fits')[0].data

        # cut and bin COV MW matrix:
        for i_tr in range(ncombs):
            for j_tr in range(ncombs):
                cov_bpw[i_tr,:, j_tr,:] = rebin(cut_array(cov_bpw_nobin[i_tr,:, j_tr,:], \
                                                np.arange(3 * NSIDE), LMIN, LMAX), [NBANDS, NBANDS])

        # save to fits file
        hducov = fits.PrimaryHDU(cov_bpw)
        hducov.writeto(PATH_DICT['output_path'] + \
                                '_'.join([NAME_RUN, NAME_COMP, 'Cov', 'bin']) + \
                                f'_{type_cov}.fits', overwrite = True)

    cov_bpw = cov_bpw.reshape([ncombs * NBANDS, ncombs * NBANDS ])

    s_d.add_covariance(cov_bpw)

    # save files
    print("Writing")
    s_d.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, type_cov]) + \
                    '_tot.fits', overwrite = True)
    s_f.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, type_cov]) +\
                     '_fid.fits', overwrite = True)
    s_n.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, type_cov]) + \
                    '_noi.fits', overwrite = True)

def compute_cl_nobin(ctype):

    '''
    Computes fiducial Cl fits file
    This Cl will go into the computation of the covariance
    No binning applied

    ** Parameters **
    ctype: 'd00', 'dc0', 'dcs'
        components used in model
    '''

    lmax = 3 * NSIDE - 1
    larr_all = np.arange(lmax + 1)

    ncomp = ctype_dict_ncomp[ctype]

    bpw_freq_sig = get_bpw_freq_sig(ctype, lmax, False)

    fsky = nc.get_fsky()

    if ncomp > 1:
        ## add noise
        bpw_freq_noi = add_noise(lmax, bpw_freq_sig, fsky, False)
        bpw_freq_sig += bpw_freq_noi

    bpw_freq_sig = bpw_freq_sig.reshape([nfreqs*nmodes,nfreqs*nmodes, lmax + 1])

    # Create sacc and add tracers
    print("Adding tracers")
    s_d = add_tracers(larr_all)
    # Adding power spectra
    print("Adding spectra")
    s_d = add_powerspectra(s_d, bpw_freq_sig, larr_all, False)

    print("Writing")
    s_d.save_fits(PATH_DICT['output_path'] +\
                    '_'.join([NAME_RUN, ctype]) + \
                    '_clnobin.fits', overwrite = True)
