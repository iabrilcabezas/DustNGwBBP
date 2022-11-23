'''
compute cell
'''

#import sys
import numpy as np
import sacc
from utils.params import pathnames
import utils.noise_calc as nc
from utils.SED import get_band_names, Bpass, get_component_spectra, get_convolved_seds
from utils.bandpowers import get_ell_arrays

#sys.path.append('/global/common/software/act/python/DustNGwBBP')

# from utils.binning import cut_array, rebin


# def add_noise(machine, bpw_freq, experiment, nmodes, forcov):

#     if experiment == 'bicep':

#         n_ell, n_bpw = nc.bicep_noise_fromfile(machine)
#         assert len(n_ell) == bpw_freq.shape[-1]

#         if not forcov:
#                 n_bpw *= (n_ell * (n_ell + 1) / (2 * np.pi) ) # weight, we are working with Dl

#         fsky = nc.fsky_fromnoise(machine)[0]

#         bpw_freq_noi=np.zeros_like(bpw_freq)
#         for ib,n in enumerate(n_bpw):
#             for nm in nmodes:
#                    bpw_freq_noi[ib,nm,ib,nm,:]=n_bpw[ib,:]
        
#     if experiment == 'so':
#         # N_ell
#         sens=2
#         knee=1
#         ylf=1
#         fsky=0.1
#         nell=np.zeros([bpw_freq.shape[0],lmax+1])
#         _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,lmax+1,1)
#         n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)
#         bpw_freq_noi=np.zeros_like(bpw_freq_sig)
#         for ib,n in enumerate(n_bpw):
#             bpw_freq_noi[ib,0,ib,0,:]=n_bpw[ib,:]
#             bpw_freq_noi[ib,1,ib,1,:]=n_bpw[ib,:]

def import_bandpasses(machine, experiment):

    '''
    Imports bandpasses
    '''

    band_names = get_band_names(experiment)
    path_dict = dict(pathnames(machine))

    if experiment == 'bicep':
        bpss = {n: Bpass(n, path_dict['BK15_data'] +\
                             f'BK15_{n}_bandpass_20180920.txt') for n in band_names}
    if experiment == 'so':
        bpss = {n: Bpass(n,path_dict['bbpipe_path'] +\
                             f'examples/data/bandpasses/{n}.txt') for n in band_names}

    return bpss

def import_beams(experiment, larr_all):

    '''
    Imports beams
    '''

    band_names = get_band_names(experiment)

    if experiment == 'bicep':
        beams ={band_names[i]: b for i, b in enumerate(nc.bicep_beams_exp(larr_all))}
    if experiment == 'so':
        beams ={band_names[i]: b for i, b in enumerate(nc.Simons_Observatory_V3_SA_beams(larr_all))}

    return beams

def get_windows(lmin, dell, nbands, weight):

    '''
    get windows
    '''

    weight_types = ['Dl', 'Cl']
    assert weight in weight_types, 'not a type of weight!'

    lmax, larr_all, lbands, _ = get_ell_arrays(lmin, dell, nbands)

    windows = np.zeros([nbands,lmax+1])

    cl_weights = np.ones_like(larr_all)

    for b,(l0,lf) in enumerate(zip(lbands[:-1],lbands[1:])):

        if weight == 'Dl':
            windows[b,l0:lf] = (larr_all * (larr_all + 1)/(2*np.pi))[l0:lf]
        if weight == 'Cl':
            windows[b, l0:lf] = cl_weights[l0:lf]

        windows[b,:] /= dell

    return windows

def add_tracers(machine, experiment, larr_all):

    '''
    Creates sacc object and add tracers
    '''

    s_d = sacc.Sacc()

    band_names = get_band_names(experiment)
     # Bandpasses: 
    bpss = import_bandpasses(machine, experiment)
    # Beams
    beams = import_beams(experiment, larr_all)

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

def compute_cl_forcov(machine, experiment, lmin, dell, nbands, ctype):

    '''
    
    
    '''

    nmodes = 1

    lmax, larr_all, lbands, leff = get_ell_arrays(lmin, dell, nbands)

    band_names = get_band_names(experiment)
    nfreqs = len(band_names)

    bpss = import_bandpasses(machine, experiment)

    cl2dl=larr_all*(larr_all+1)/(2*np.pi)
    dl2cl=np.zeros_like(cl2dl)
    dl2cl[1:] = 1/(cl2dl[1:])

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

    windows = get_windows(lmin, dell, nbands, weight = 'Cl')

    bpw_comp=np.sum(dls_comp[:,:,:,:,None,:]*windows[None,None,None, None, :,:],axis=5)

    # Convolve with bandpasses
    seds_all = get_convolved_seds(band_names, bpss)

    if ncomp == 1:
        seds = np.zeros([1, nfreqs])
        seds[0,:] = seds_all[1,:]
    if ncomp == 2:
        seds = seds_all

    bpw_freq_sig = np.einsum('ik,jm,iljno', seds, seds, bpw_comp)
    
    if (ctype == 'all'):
        
        if experiment == 'bicep':

            n_ell, n_bpw = nc.bicep_noise_fromfile(machine)
            assert len(n_ell) == bpw_freq_sig.shape[-1]

            fsky = nc.fsky_fromnoise(machine)[0]

            bpw_freq_noi=np.zeros_like(bpw_freq_sig)
            

        if experiment == 'so':

           # N_ell
            sens=2
            knee=1
            ylf=1
            fsky=0.1
            nell=np.zeros([nfreqs,lmax+1])
            _,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,lmax+1,1)
            n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)
            bpw_freq_noi=np.zeros_like(bpw_freq_sig)
        
        for ib,n in enumerate(n_bpw):
            bpw_freq_noi[ib,0,ib,0,:]=n_bpw[ib,:]        
        bpw_freq_sig += bpw_freq_noi

    bpw_freq_sig = bpw_freq_sig.reshape([nfreqs*nmodes,nfreqs*nmodes, nbands])

    # Create sacc and add tracers
    print("Adding tracers")
    s_d = add_tracers(machine, experiment, larr_all)

    # Adding power spectra
    print("Adding spectra")
    nmaps=nmodes*nfreqs
    indices_tr=np.triu_indices(nmaps)
    map_names=[]
    for ib, n in enumerate(band_names):
        map_names.append('band%d' % (ib+1) + '_B')

    for (i1, i2) in zip(indices_tr[0], indices_tr[1]):
        n1 = map_names[i1][:-2]
        n2 = map_names[i2][:-2]
        p1 = map_names[i1][-1].lower()
        p2 = map_names[i2][-1].lower()
        cl_type = f'cl_{p1}{p2}'
        s_d.add_ell_cl(cl_type, n1, n2, leff, bpw_freq_sig[i1, i2, :])

    return s_d