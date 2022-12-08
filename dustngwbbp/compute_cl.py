'''
compute_cl
    computes c_ell to input into covariance calculation
    creates sacc files that contain spectrum and covariances

'''

from copy import deepcopy
import numpy as np
import sacc
from astropy.io import fits
from utils.params import NAME_RUN, PATH_DICT
from utils.params import EXPERIMENT, NSIDE, LMIN, DELL, NBANDS, POLARIZATION_cov
import utils.noise_calc as nc
from utils.binning import rebin, cut_array
from utils.sed import get_band_names, Bpass, get_component_spectra, get_convolved_seds
from utils.bandpowers import get_ell_arrays, dell2cell_lmax

band_names = get_band_names()
LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)

ctype_dict_ncomp = {'dust': 1, 'all': 2}

def import_bandpasses():

    '''
    Returns dictionary of bandpass class for each channel in EXPERIMENT
    '''

    if EXPERIMENT == 'bicep':
        bpss = {n: Bpass(n, PATH_DICT['BK15_data'] +\
                             f'BK15_{n}_bandpass_20180920.txt') for n in band_names}
    if EXPERIMENT == 'so':
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
    if EXPERIMENT == 'so':
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

def add_powerspectra(s_d, bpw_freq_sig, polarization, weight):

    '''
    add power spectra to sacc
    '''

    nmaps=nmodes*nfreqs
    indices_tr=np.triu_indices(nmaps)
    windows = get_windows(weight)
    s_wins = sacc.BandpowerWindow(LARR_ALL, windows.T)

    map_names=[]
    for i_freq in range(nfreqs):
        if 'E' in polarization:
            map_names.append(f'band{i_freq+1}_E')
        if 'B' in polarization:
            map_names.append(f'band{i_freq+1}_B')

    for (i_tr1, i_tr2) in zip(indices_tr[0], indices_tr[1]):
        band12 = np.array([map_names[i_tr1][:-2], map_names[i_tr2][:-2]])
        pol12  = np.array([map_names[i_tr1][-1].lower(), map_names[i_tr2][-1].lower()])
        cl_type = f'cl_{pol12[0]}{pol12[1]}'
        s_d.add_ell_cl(cl_type, band12[0], band12[1], LEFF, bpw_freq_sig[i_tr1, i_tr2, :],
                         window = s_wins)

    return s_d

def add_powerspectra_nobin(s_d, bpw_freq_sig, polarization, leff):

    '''
    Same as add_powerspectra() but without windows (Cl weighting)
    '''

    nmaps=nmodes*nfreqs
    indices_tr=np.triu_indices(nmaps)

    map_names=[]
    for i_freq in range(nfreqs):
        if 'E' in polarization:
            map_names.append(f'band{i_freq+1}_E')
        if 'B' in polarization:
            map_names.append(f'band{i_freq+1}_B')

    for (i_tr1, i_tr2) in zip(indices_tr[0], indices_tr[1]):
        band12 = np.array([map_names[i_tr1][:-2], map_names[i_tr2][:-2]])
        pol12  = np.array([map_names[i_tr1][-1].lower(), map_names[i_tr2][-1].lower()])
        cl_type = f'cl_{pol12[0]}{pol12[1]}'
        s_d.add_ell_cl(cl_type, band12[0], band12[1], leff, bpw_freq_sig[i_tr1, i_tr2, :])

    return s_d


def add_noise(weight, bpw_freq_sig, fsky):

    '''
    add noise
    '''

    windows = get_windows(weight)

    if EXPERIMENT == 'bicep':

        nell = np.zeros([nfreqs, LMAX + 1])
        _, nell[:, 2:] = nc.bicep_noise(LMAX)

        n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)

        # n_ell, n_bpw = nc.bicep_noise_fromfile(MACHINE)
        # assert np.all(np.isclose(n_ell, LEFF))
        # if weight == 'dl':
            #     n_bpw *= (n_ell * (n_ell + 1) / (2 * np.pi) )

    if EXPERIMENT == 'so':

        # N_ell
        sens=2
        knee=1
        ylf=1
        nell=np.zeros([nfreqs,LMAX+1])
        _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,LMAX+1,1)
        n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)

    bpw_freq_noi=np.zeros_like(bpw_freq_sig)

    for ib in range(len(n_bpw)):
        bpw_freq_noi[ib,0,ib,0,:]=n_bpw[ib,:]
        if nmodes == 2:
            bpw_freq_noi[ib,1,ib,1,:]=n_bpw[ib,:]


    return bpw_freq_noi

def add_noise_nobin(lmax, bpw_freq_sig, fsky):

    '''
    Same as add_noise() but without binning according to window (cl weighting, no binning)
    '''

    if EXPERIMENT == 'bicep':

        # n_ell, n_bpw = nc.bicep_noise_fromfile(MACHINE)
        # assert np.all(np.isclose(n_ell, LEFF))

        # if weight == 'dl':

        #     n_bpw *= (n_ell * (n_ell + 1) / (2 * np.pi) )
        nell = np.zeros([nfreqs, lmax + 1])
        _, nell[:, 2:] = nc.bicep_noise(lmax)

    if EXPERIMENT == 'so':

        # N_ell
        sens=2
        knee=1
        ylf=1
        nell=np.zeros([nfreqs,lmax+1])
        _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,lmax+1,1)

    bpw_freq_noi=np.zeros_like(bpw_freq_sig)

    for ib in range(len(nell)):
        bpw_freq_noi[ib,0,ib,0,:]=nell[ib,:]
        if nmodes == 2:
            bpw_freq_noi[ib,1,ib,1,:]=nell[ib,:]

    return bpw_freq_noi


def compute_cl(ctype, type_cov):

    '''
    because for cov , weight = 'cl'
    '''
    weight = 'Cl'
    weight_name = '_'.join([weight, type_cov])

    bpss = import_bandpasses()

    dl2cl = dell2cell_lmax(LMAX)

    ncomp = ctype_dict_ncomp[ctype]

    dls_comp = np.zeros([ncomp,nmodes,ncomp,nmodes,LMAX+1]) #[ncomp,np,ncomp,np,nl]
    dls_sed  = get_component_spectra(LMAX)

    if ncomp == 1:
        dls_comp[0,0,0,0,:] = dls_sed[1]

    if ncomp == 2:
        (dls_comp[1,0,1,0,:],
        dls_comp[0,0,0,0,:]) = dls_sed[1], dls_sed[3]


    dls_comp *= dl2cl[None, None, None, None, :]

    windows = get_windows(weight)

    bpw_comp=np.sum(dls_comp[:,:,:,:,None,:]*windows[None,None,None, None, :,:],axis=5)

    # Convolve with bandpasses
    seds = get_convolved_seds(band_names, bpss)

    if ncomp == 1:
        seds = np.array([seds[1,:]])

    bpw_freq_sig = np.einsum('ik,jm,iljno', seds, seds, bpw_comp)

    fsky = 1. # mask = full, nc.get_fsky()

    bpw_freq_tot = deepcopy(bpw_freq_sig)
    bpw_freq_noi = np.zeros_like(bpw_freq_sig)

    if ctype == 'all':
        ## add noise
        bpw_freq_noi = add_noise(weight, bpw_freq_sig, fsky)
        bpw_freq_tot += bpw_freq_noi

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
    s_d = add_powerspectra(s_d, bpw_freq_sig, POLARIZATION_cov, weight)
    s_f = add_powerspectra(s_f, bpw_freq_sig, POLARIZATION_cov, weight)
    s_n = add_powerspectra(s_n, bpw_freq_noi, POLARIZATION_cov, weight)

    assert (ctype == 'all'), 'adding cov to dust only?'
    ncombs = len(s_d.get_tracer_combinations())

    cov_bpw = np.zeros([ncombs, NBANDS, ncombs, NBANDS])

    nmaps=nmodes*nfreqs
    indices_tr=np.triu_indices(nmaps)

    assert len(indices_tr[0]) == ncombs, 'how many maps do you have?'

    # factor_modecount = 1./((2*leff+1)*DELL*fsky)
    # for ii, (i1, i2) in enumerate(zip(indices_tr[0], indices_tr[1])):
    #     for jj, (j1, j2) in enumerate(zip(indices_tr[0], indices_tr[1])):
    #         covar = (bpw_freq_tot[i1, j1, :]*bpw_freq_tot[i2, j2, :]+
    #                 bpw_freq_tot[i1, j2, :]*bpw_freq_tot[i2, j1, :]) * factor_modecount
    #         cov_bpw[ii, :, jj, :] = np.diag(covar)

    if type_cov == 'w':
        cov_bpw_full = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'all', 'Cov']) +\
                             '_nobin_w.fits')[0].data

    if type_cov == 'wt':
        cov_bpw_full =  fits.open(PATH_DICT['output_path'] + NAME_RUN + \
                            '_nobin_fullCov.fits')[0].data

    # cut and bin COV MW matrix:
    for i_tr in range(ncombs):
        for j_tr in range(ncombs):
            cov_bpw[i_tr,:, j_tr,:] = rebin(cut_array(cov_bpw_full[i_tr,:, j_tr,:], \
                                            np.arange(3 * NSIDE), LMIN, LMAX), [NBANDS, NBANDS])

    cov_bpw = cov_bpw.reshape([ncombs * NBANDS, ncombs * NBANDS ])
    s_d.add_covariance(cov_bpw)

    print("Writing")
    s_d.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, weight_name]) + \
                    '_tot.fits', overwrite = True)

    s_f.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, weight_name]) +\
                     '_fid.fits', overwrite = True)
    s_n.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, weight_name]) + \
                    '_noi.fits', overwrite = True)


def compute_cl_nobin(ctype):

    '''
    because for cov , weight = 'cl'
    '''
    weight = 'Cl'

    bpss = import_bandpasses()

    lmax = 3 * NSIDE - 1
    larr_all = np.arange(lmax + 1)

    dl2cl = dell2cell_lmax(lmax)
    ncomp = ctype_dict_ncomp[ctype]

    dls_comp = np.zeros([ncomp,nmodes,ncomp,nmodes,lmax+1]) #[ncomp,np,ncomp,np,nl]
    dls_sed  = get_component_spectra(lmax)

    if ncomp == 1:
        dls_comp[0,0,0,0,:] = dls_sed[1]

    if ncomp == 2:
        (dls_comp[1,0,1,0,:],
        dls_comp[0,0,0,0,:]) = dls_sed[1], dls_sed[3]

    dls_comp *= dl2cl[None, None, None, None, :]

    # Convolve with bandpasses
    seds = get_convolved_seds(band_names, bpss)

    if ncomp == 1:
        seds = np.array([seds[1,:]])

    bpw_freq_sig = np.einsum('ik,jm,iljno', seds, seds, dls_comp)

    fsky = 1. # mask = full, nc.get_fsky()

    if ctype == 'all':
        ## add noise
        bpw_freq_noi = add_noise_nobin(lmax, bpw_freq_sig, fsky)
        bpw_freq_sig += bpw_freq_noi

    bpw_freq_sig = bpw_freq_sig.reshape([nfreqs*nmodes,nfreqs*nmodes, lmax + 1])

    # Create sacc and add tracers
    print("Adding tracers")
    s_d = add_tracers(larr_all)
    # Adding power spectra
    print("Adding spectra")
    s_d = add_powerspectra_nobin(s_d, bpw_freq_sig, POLARIZATION_cov, larr_all)

    print("Writing")
    s_d.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, weight, str(NSIDE)]) + \
                    '_tot.fits', overwrite = True)
 