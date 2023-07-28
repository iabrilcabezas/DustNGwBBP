'''
sample_to_sacc.py

contains

sacc_from_sample: samples covariance and saves fit with g,ng cov
'''

import numpy as np
from astropy.io import fits
import sacc
from utils.params import LMIN, DELL, NBANDS, POLARIZATION_cov
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP
from utils.sed import get_band_names
from utils.bandpowers import get_ell_arrays
from dustngwbbp.compute_cl import import_bandpasses, import_beams, get_windows


band_names = get_band_names()
LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

weight = 'Cl'

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncross = len(indices_tr[0])

INPUT_CL = '/global/cfs/cdirs/act/data/iabril/BBPower/230525_sample/cells_model.fits'
COV_G_PATH = PATH_DICT['input_path'] + '_'.join([ NAME_RUN, NAME_COMP]) + \
                                '_' + '_'.join(['Cov', 'bin', 'w']) + '.fits'
COV_NG_PATH = PATH_DICT['input_path'] + '_'.join([ NAME_RUN, NAME_COMP]) + \
                                '_' + '_'.join(['Cov', 'bin', 'dfwt']) + '.fits'

counter = np.arange(len(indices_tr[0]))
dicttr = {}
for i, counter_val in enumerate(counter):
    dicttr[(indices_tr[0][i], indices_tr[1][i])] = counter_val

# Bandpasses:
bpss = import_bandpasses()
# Beams
beams = import_beams(LARR_ALL)
# windows
windows = get_windows(weight)
s_wins = sacc.BandpowerWindow(LARR_ALL, windows.T)

def sacc_from_sample(input_cl, cov_g_path, cov_ng_path, test = False):

    '''
    samples from input_cl with cov_ng and saves sacc with both covs
    '''

    # Creating Sacc files
    s_d_g = sacc.Sacc()
    s_d_ng = sacc.Sacc()

    # Adding tracers
    print("Adding tracers")
    for ib, n in enumerate(band_names):
        bandpass = bpss[n]
        beam = beams[n]
        for s in [s_d_g, s_d_ng]:
            s.add_tracer('NuMap', f'band{ib+1}',
                        quantity='cmb_polarization',
                        spin=2,
                        nu=bandpass.nu,
                        bandpass=bandpass.bnu,
                        ell=LARR_ALL,
                        beam=beam,
                        nu_unit='GHz',
                        map_unit='uK_CMB')

    # read in cov, already binned
    cov_g_bpw  = fits.open(cov_g_path)[0].data
    cov_ng_bpw  = fits.open(cov_ng_path)[0].data

    cell_model = sacc.Sacc.load_fits(input_cl)

    random_cells = np.random.default_rng().multivariate_normal(cell_model.mean,\
                                                               cov_ng_bpw.reshape([NBANDS*ncross, NBANDS*ncross]),\
                                                               size = 1)[0] #nsims

    if test:
        cell_for_sacc = np.reshape(cell_model.mean, (ncross, NBANDS))
    else:
        cell_for_sacc = np.reshape(random_cells, (ncross, NBANDS))

    # Adding power spectra
    print("Adding spectra")
    map_names=[]
    for i_freq in range(nfreqs):
        if 'E' in POLARIZATION_cov:
            map_names.append(f'band{i_freq+1}_E')
        if 'B' in POLARIZATION_cov:
            map_names.append(f'band{i_freq+1}_B')
    for ii, (i1, i2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        n1 = map_names[i1][:-2]
        n2 = map_names[i2][:-2]
        p1 = map_names[i1][-1].lower()
        p2 = map_names[i2][-1].lower()
        cl_type = f'cl_{p1}{p2}'
        for s in [s_d_g, s_d_ng]:
            s.add_ell_cl(cl_type, n1, n2, LEFF, cell_for_sacc[dicttr[(i1, i2)]], window=s_wins)
        

    cov_g_bpw = cov_g_bpw.reshape([ncross * NBANDS, ncross * NBANDS])
    s_d_g.add_covariance(cov_g_bpw)
    cov_ng_bpw = cov_ng_bpw.reshape([ncross * NBANDS, ncross * NBANDS])
    s_d_ng.add_covariance(cov_ng_bpw)

    return [s_d_g, s_d_ng]
