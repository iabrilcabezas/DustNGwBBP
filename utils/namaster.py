'''
namaster
    pymaster package routines
'''

import numpy as np
import healpy as hp
import pymaster as nmt
from astropy import units as u
from astropy.coordinates import SkyCoord
from utils.params import PATH_DICT, nu0_dust
from utils.sed import fcmb

def hp_rotate(map_hp, coord):
    """Rotate healpix map between coordinate systems

    :param map_hp: A healpix map in RING ordering
    :param coord: A len(2) list of either 'G', 'C', 'E'
    Galactic, equatorial, ecliptic, eg ['G', 'C'] converts
    galactic to equatorial coordinates
    :returns: A rotated healpix map
    """
    if map_hp is None:
        return None
    if coord[0] == coord[1]:
        return map_hp
    r = hp.rotator.Rotator(coord=coord)
    new_map = r.rotate_map_pixel(map_hp)
    return new_map


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


def binning_bbpower(lmin, lmax, dell, nside):

    '''
    specifies binning scheme for namaster and also returns effective ells
    
    *Parameters*
    
    lmin: int
        minimum ell
    lmax: int
        ell max
    dell: int
        width of ell bins
    nside: int
        resolution of map
    is_Dell: bool
        weight by ell * (ell + 1) / 2pi ?
    
    * Returns * 
    b: nmt bin
    ell_eff: np.array 
    '''

    ells = np.arange(lmin, lmax, dtype = 'int32')
    weights = 1/dell * np.ones_like(ells)
    bpws = -1 + np.zeros_like(ells)

    i = 0
    while dell * (i + 1) < lmax:
        bpws[dell * i : dell * (i + 1)] = i
        i += 1

    b_df = nmt.NmtBin(nside = nside, bpws = bpws, ells = ells, weights = weights)
    ell_eff_df = b_df.get_effective_ells()

    return b_df, ell_eff_df

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
            type of mask. 'full' or 'bicep' or 'so'
    kwargs: dict
            apo_deg: float
                    apodization scale in degrees for mtype == 'bicep'

    **Returns**

    np.array()
            mask

    '''

    npix  = hp.nside2npix(nside)

    if mtype == 'full':

        # no regions masked, 'full'
        return np.ones(npix)

    if mtype == 'bicep':

        w_bicep = np.zeros(npix)

        # join three triangles to create bicep-like mask
        # coordinates in (RA, DEC) taken from Fig. 4 https://arxiv.org/pdf/2210.05684.pdf.
        # the region amplify to recover fsky bicep (0.014) after apodization
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
        w_apo = nmt.mask_apodization(w_bicep, kwargs['apo_deg'], apotype="C1")

        return w_apo

    if mtype == 'so' or mtype == 'soflat':

        # read in SO mask (Figure 8 - SO science & forecasts)
        w_so = hp.read_map(PATH_DICT['so_path'] + 'mask_apodized_david_nside512.fits')

        # RA/DEC to Galactic
        mask_so_gal = hp_rotate(w_so, coord = ['C','G'])

        w_so_hi = hp.ud_grade(mask_so_gal, nside)
    
        if mtype == 'so':

            # not smooth the mask. it says apodized already
            # w_so_apo = nmt.mask_apodization(w_so_hi, kwargs['apo_deg'], apotype="C1")
            return w_so_hi

        else:

            so_patch = w_so_hi  > 0.35
            mask_so_now = np.zeros_like(w_so_hi)
            mask_so_now[so_patch] = 1
            assert np.isclose(kwargs['apo_deg'], 5.), 'if apo_deg !=5, fsky is not 0.1'
            mask_so_apo = nmt.mask_apodization(mask_so_now, kwargs['apo_deg'] , apotype="C1") # apo_deg

            return mask_so_apo


    if mtype == 'fsky07':

        fsky07 = hp.read_map(PATH_DICT['planck_data'] + 'HFI_Mask_GalPlane-apo5_2048_R2.00.fits', field = 4) # 60% of sky (70% + apodization)
        fsky07_hi = hp.ud_grade(fsky07, nside)
        # already apodized mask

        return fsky07_hi

    return None


def get_template_wmask(nside, ttype, mtype, **kwargs):

    '''
    Obtains modulating template. Constructed from Appendix A instructions

    ** Parameters **
    machine: str
            Machine where code runs
    nside: int
            resolution of maps
    mtype: str
            type of mask
    ttype: str
            template to use
    kwargs: dict
            'apo_deg':      float. Apodization scale of mask in degrees
            'smooth_deg':   float. Smoothing scale of template in degrees

    ** Returns **
    template: np.array()
            modulating template
    w_mask: np.array()
            mask
    '''

    if ttype == 'p353':
        # map from which to create anisotropic correlated template
        templ = hp.read_map(PATH_DICT['planck_data']+"HFI_SkyMap_353-psb-field-IQU_2048_R3.00_full.fits")
    if ttype == 'd10':
        templ_rj = hp.read_map(PATH_DICT['template_path'] + 'dust_pysm3_small_scales.fits', field = 0)
        templ = templ_rj / fcmb(nu0_dust) # RJ to CMB units
    else:
        return None

    templ_hi = hp.ud_grade(templ, nside)
    # create large-scale modulating template by smoothing 353 map
    if np.isnan(kwargs['smooth_deg']):
        template_0 = templ_hi
        print('not smoothed')
    else:
        template_0 = hp.smoothing(templ_hi, fwhm=np.deg2rad(kwargs['smooth_deg']))

    # obtain mask
    w_mask = get_mask(nside, mtype, **kwargs)
    w2_omega_mask = np.mean(w_mask**2)

    # normalize template, here we see map units do not matter
    template = template_0 * np.sqrt( w2_omega_mask / np.mean(w_mask**2 * template_0**2))

    # check
    wtilde_mask = w_mask * template
    wtilde2_omega_mask = np.mean(wtilde_mask**2)
    check = np.isclose(wtilde2_omega_mask, w2_omega_mask)
    assert check, "<w2> != <w2> ERROR"

    return (template, w_mask, w2_omega_mask)
