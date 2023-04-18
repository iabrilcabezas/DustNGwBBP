'''
add_cov2decorr
    adds gaussian and non-gaussian cov matrices to data with decorrelation
'''

from astropy.io import fits
import sacc
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP
from utils.params import POLARIZATION_cov
from utils.sed import get_band_names

band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
pols = [POLARIZATION_cov]
nmaps=nmodes*nfreqs
ncross = (nmaps * (nmaps + 1))//2

def add_cov2decorr(type_cov, weight):

    '''
    reads in cells with decorrelation, and adds covariance to saccs
    saves gaussian (w) and non gaussian (type_cov) covariance
    '''

    # load fiducial power spectrum for NAME_RUN parameters
    name_sf = PATH_DICT['output_path'] + 'cells_model.fits'

    s_fg = sacc.Sacc.load_fits(name_sf)
    s_fng = sacc.Sacc.load_fits(name_sf)

    # no. of bandpowers after ell range selection
    tr_names = sorted(list(s_fg.tracers.keys()))

    ell_b = s_fg.get_ell_cl('cl_' + 2 * pols[0].lower(), \
                                            tr_names[0], tr_names[0])[0]
    n_bpws = len(ell_b)


    # populate with gaussian and NG covariance
    cov_ng = fits.open(PATH_DICT['input_path'] + NAME_RUN + '_' + \
                            '_'.join([NAME_COMP,'Cov', 'bin']) + f'_{type_cov}.fits')[0].data
    cov_ng = cov_ng.reshape([ncross* n_bpws,ncross * n_bpws ])
    s_fng.add_covariance(cov_ng)
    # repeat:
    cov_g = fits.open(PATH_DICT['input_path'] + NAME_RUN + '_' + \
                            '_'.join([NAME_COMP,'Cov', 'bin']) + '_w.fits')[0].data
    cov_g = cov_g.reshape([ncross* n_bpws,ncross * n_bpws ])
    s_fg.add_covariance(cov_g)

    # save:
    s_fg.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, 'w']) + \
                    '_tot.fits', overwrite = True)
    s_fng.save_fits(PATH_DICT['output_path'] + '_'.join([NAME_RUN, weight, type_cov]) + \
                    '_tot.fits', overwrite = True)
