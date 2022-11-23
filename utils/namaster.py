'''
namaster
'''

import numpy as np
import healpy as hp
import pymaster as nmt
from astropy import units as u
from astropy.coordinates import SkyCoord
from utils.params import pathnames

def get_bandpowers(nside, width_ell):

    '''
    Construct linear binning scheme of bandpower

    Parameters
    ----------
    nside: int
            resolution parameter
    width_ell: int
            bandpower width

    Returns
    --------
    bin_bpw: Class
            binning of bandpower
    ell_eff: np.array()
            effective multiple associated with each bandpower
    '''

    bin_bpw = nmt.NmtBin.from_nside_linear(nside, width_ell)
    ell_eff = bin_bpw.get_effective_ells()

    return bin_bpw, ell_eff


def get_couplingmatrix(f_a, f_b, bin_bpw):

    '''
    Generates NmtWorkspace to compute and obtain the mode coupling matrix. 
    The coupling matrix depends only on the masks of the two fields to correlate

    Parameters
    ----------
    f_a, f_b: nmt.NmtField()
            fields that we correlate
    bin_bpw: Class nmt.NmtBin.from_nside_linear()
            binning of bandpowers

    Returns
    ----------
    np.array()
            mode coupling matrix

    '''

    nmtw = nmt.NmtWorkspace()
    nmtw.compute_coupling_matrix(f_a, f_b, bin_bpw)

    return nmtw.get_coupling_matrix()


def get_mask(nside, mtype, **kwargs):

    '''
    Obtain mask of type mtype

    **Parameters**
    nside: int
            resolution
    mtype: str
            type of mask. 'full' or 'bicep'
    kwargs: dict
            apo_deg: float
                    apodization scale in degrees for mtype == 'bicep'

    **Returns**

    np.array()
            mask

    '''

    npix  = hp.nside2npix(nside)

    if mtype == 'full':

        return np.ones(npix)


    if mtype == 'bicep':

        w_bicep = np.zeros(npix)

        # join three triangles to create bicep-like mask
        # coordinates in (RA, DEC) taken from Fig. 4 https://arxiv.org/pdf/2210.05684.pdf. region amplify to recover fsky bicep (0.014) after apodization
        radecs = [[[-60, -75], [60, -75], [0, -35]],
                [[60, -35],  [60, -75], [-1, -35]],
                [[-60, -35], [-60, -75], [1, -35]]]

        for radec_i in radecs:

            # 3 vertices triangle. store 3D coordinates of vertices
            vecs = np.empty((len(radec_i), 3)) + np.nan

            for j, radec_ij in enumerate(radec_i):

                # (RA, DEC)
                c_coord = SkyCoord(ra=radec_ij[0]*u.degree, dec=radec_ij[1]*u.degree)
                # (longitude, latitude) conversion
                l_gal_rad = c_coord.galactic.l.radian
                b_gal_rad = c_coord.galactic.b.radian
                # hp in (colatitude and longitude)
                vecs[j] = hp.ang2vec(np.pi/2 -  b_gal_rad, l_gal_rad)

            ipix_bicep =hp.query_polygon(nside = nside, vertices = vecs)
            w_bicep[ipix_bicep] = 1. # this is the region we want (1)

        # smooth the mask
        w_apo = nmt.mask_apodization(w_bicep, kwargs['apo_deg'], apotype="Smooth")

        return w_apo

    return None


def get_template_wmask(machine, nside, mtype, **kwargs):

    '''
    Obtains modulating template. Constructed from Appendix A instructions

    ** Parameters **
    machine: str
            Machine where code runs
    nside: int
            resolution of maps
    mtype: str
            type of mask
    kwargs: dict
            'apo_deg':      float. Apodization scale of mask in degrees
            'smooth_deg':   float. Smoothing scale of template in degrees

    ** Returns **
    template: np.array()
            modulating template
    w_mask: np.array()
            mask
    '''

    path_dict = dict(pathnames(machine))

    # map from which to create anisotropic correlated template
    p353 = hp.read_map(path_dict['planck_data']+"HFI_SkyMap_353-psb-field-IQU_2048_R3.00_full.fits")
    p353_hi = hp.ud_grade(p353, nside)
    # create large-scale modulating template by smoothing 353 map
    template_0 = hp.smoothing(p353_hi, fwhm=np.deg2rad(kwargs['smooth_deg']))

    # obtain mask
    w_mask = get_mask(nside, mtype, **kwargs)
    w2_omega_mask = np.mean(w_mask**2)

    template = template_0 * np.sqrt( w2_omega_mask / np.mean(w_mask**2 * template_0**2))

    # check
    wtilde_mask = w_mask * template
    wtilde2_omega_mask = np.mean(wtilde_mask**2)
    check = np.isclose(wtilde2_omega_mask, w2_omega_mask)
    assert check, "<w2> != <w2> ERROR"

    return (template, w_mask)
