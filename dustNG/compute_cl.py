'''
compute cell
'''

#import sys

#sys.path.append('/global/common/software/act/python/DustNGwBBP')

import numpy as np
import sacc
import yaml
from utils.params import pathnames, Config
import utils.noise_calc as nc
from utils.SED import get_band_names, Bpass, get_component_spectra, get_convolved_seds
from utils.bandpowers import get_ell_arrays, dell2cell_lmax


config = Config(yaml.load(open('config.yaml'), yaml.FullLoader))

EXPERIMENT = config.global_param.experiment
MACHINE = config.global_param.machine
LMIN = config.bpw_param.lmin
DELL = config.bpw_param.dell
NBANDS = config.bpw_param.nbands
POLARIZATION_cov = config.pol_param.pol_cov

band_names = get_band_names(EXPERIMENT)
path_dict = dict(pathnames(MACHINE))
lmax, larr_all, lbands, leff = get_ell_arrays(LMIN, DELL, NBANDS)

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)

#sys.path.append('/global/common/software/act/python/DustNGwBBP')
# from utils.binning import cut_array, rebin

def import_bandpasses():

    '''
    Imports bandpasses
    '''

    if EXPERIMENT == 'bicep':
        bpss = {n: Bpass(n, path_dict['BK15_data'] +\
                             f'BK15_{n}_bandpass_20180920.txt') for n in band_names}
    if EXPERIMENT == 'so':
        bpss = {n: Bpass(n,path_dict['bbpipe_path'] +\
                             f'examples/data/bandpasses/{n}.txt') for n in band_names}

    return bpss

def import_beams():

    '''
    Imports beams
    '''

    if EXPERIMENT == 'bicep':
        beams ={band_names[i]: b for i, b in enumerate(nc.bicep_beams_exp(larr_all))}
    if EXPERIMENT == 'so':
        beams ={band_names[i]: b for i, b in enumerate(nc.Simons_Observatory_V3_SA_beams(larr_all))}

    return beams

def get_windows(weight):

    '''
    get windows
    '''

    weight_types = ['Dl', 'Cl']
    assert weight in weight_types, 'not a type of weight!'

    windows = np.zeros([NBANDS,lmax+1])

    cl_weights = np.ones_like(larr_all)

    for b,(l0,lf) in enumerate(zip(lbands[:-1],lbands[1:])):

        if weight == 'Dl':
            windows[b,l0:lf] = (larr_all * (larr_all + 1)/(2*np.pi))[l0:lf]
        if weight == 'Cl':
            windows[b, l0:lf] = cl_weights[l0:lf]

        windows[b,:] /= DELL

    return windows

def add_tracers():

    '''
    Creates sacc object and add tracers
    '''

    s_d = sacc.Sacc()

     # Bandpasses:
    bpss = import_bandpasses()
    # Beams
    beams = import_beams()

    for ib, n in enumerate(band_names):
        bandpass = bpss[n]
        beam = beams[n]
        s_d.add_tracer('NuMap', f'band{ib+1}',
                        quantity='cmb_polarization',
                        spin=2,
                        nu=bandpass.nu,
                        bandpass=bandpass.bnu,
                        ell=larr_all,
                        beam=beam,
                        nu_unit='GHz',
                        map_unit='uK_CMB')

    return s_d

def add_powerspectra(s_d, bpw_freq_sig, polarization):

    '''
    add power spectra to sacc
    '''

    nmaps=nmodes*nfreqs
    indices_tr=np.triu_indices(nmaps)

    map_names=[]
    for ib in range(nfreqs):
        if 'E' in polarization:
            map_names.append(f'band{ib+1}_E')
        if 'B' in polarization:
            map_names.append(f'band{ib+1}_B')

    for (i1, i2) in zip(indices_tr[0], indices_tr[1]):
        band1 = map_names[i1][:-2]
        band2 = map_names[i2][:-2]
        pol1 = map_names[i1][-1].lower()
        pol2 = map_names[i2][-1].lower()
        cl_type = f'cl_{pol1}{pol2}'
        s_d.add_ell_cl(cl_type, band1, band2, leff, bpw_freq_sig[i1, i2, :])

    return s_d


def add_noise(weight, bpw_freq_sig, fsky):

    '''
    add noise
    TODO: SO like noise for bicep
    '''

    windows = get_windows(weight)
        
    if EXPERIMENT == 'bicep':

        n_ell, n_bpw = nc.bicep_noise_fromfile(MACHINE)
        assert np.all(np.isclose(n_ell, leff))

        if weight == 'dl':
        
            n_bpw *= (n_ell * (n_ell + 1) / (2 * np.pi) )


    if EXPERIMENT == 'so':

        # N_ell
        sens=2
        knee=1
        ylf=1
        nell=np.zeros([nfreqs,lmax+1])
        _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,lmax+1,1)
        n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)
        
    bpw_freq_noi=np.zeros_like(bpw_freq_sig)
    
    for ib in range(len(n_bpw)):
        bpw_freq_noi[ib,0,ib,0,:]=n_bpw[ib,:]
        if nmodes == 2:
            bpw_freq_noi[ib,1,ib,1,:]=n_bpw[ib,:]


    return bpw_freq_noi


def compute_cl_forcov(ctype):

    '''
    because for cov , weight = 'cl'
    '''

    weight = 'Cl'

    bpss = import_bandpasses()

    dl2cl = dell2cell_lmax(lmax)

    if ctype == 'dust':
        ncomp = 1
    if ctype == 'all':
        ncomp = 2

    dls_comp = np.zeros([ncomp,nmodes,ncomp,nmodes,lmax+1]) #[ncomp,np,ncomp,np,nl]
    dls_sed  = get_component_spectra(lmax)

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
        # seds = np.zeros([1, nfreqs])
        # seds[0,:] = seds_all[1,:]
        seds = np.array([seds[1,:]])
    # if ncomp == 2:
    #     seds = seds_all

    bpw_freq_sig = np.einsum('ik,jm,iljno', seds, seds, bpw_comp)

    fsky = nc.get_fsky(MACHINE, EXPERIMENT)

    if ctype == 'all':

        ## add noise
        bpw_freq_noi = add_noise(weight, bpw_freq_sig, fsky)
        bpw_freq_sig += bpw_freq_noi

    bpw_freq_sig = bpw_freq_sig.reshape([nfreqs*nmodes,nfreqs*nmodes, NBANDS])

    # Create sacc and add tracers
    print("Adding tracers")
    s_d = add_tracers()

    # Adding power spectra
    print("Adding spectra")
    s_d = add_powerspectra(s_d, bpw_freq_sig, POLARIZATION_cov)

    return s_d
