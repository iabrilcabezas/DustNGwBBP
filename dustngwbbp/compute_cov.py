'''
compute_cov
    computes covariance matrix as expressed in the Appendix
    computes effective covariance matrix for final power spectrum
'''

import numpy as np
import sacc
from astropy.io import fits
from utils.params import NBANDS, NSIDE, LMIN, DELL, POLARIZATION_cov, COV_CORR
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP
from utils.params import name_couplingmatrix_w, name_couplingmatrix_wt
from utils.sed import get_band_names
from utils.bandpowers import get_ell_arrays
from utils.noise_calc import get_fsky
from utils.binning import bin_covs

band_names = get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)

LMAX, _, _, _ = get_ell_arrays(LMIN, DELL, NBANDS)

def get_cov_fromeq(ell_array, cell_array, mcm, corr_factor):

    '''
    computes covariance matrix from equation A20
    
    '''


    assert len(ell_array) == cell_array.shape[-1], 'ell and C_ell dimensions do not match'
    assert cell_array.shape[0] == cell_array.shape[1], 'cell frequency cross not square'
    # assert mcm.shape[0] == 3 * nside, 'mode coupling matrix is not matching required resolution'
    assert mcm.shape[0] == mcm.shape[1], 'mode coupling matrix is not square'
    assert len(ell_array) == mcm.shape[0], 'ell and M_\ell,\ellprime x-dimension do not match'
    assert isinstance(corr_factor, float), 'expecting covariance correction factor to be a float'

    nells  = len(ell_array)
    nfreq = cell_array.shape[0]
    indices_tr = np.triu_indices(nfreq)
    ncombs = len(indices_tr[0])

    prefactor_l = np.zeros([nells, nells]) + np.nan
    sum_m_f     = 1 / (2 * ell_array + 1)
    sum_m_f    /= corr_factor

    for i in range(int(nells)):
        prefactor_l[:, i ] = sum_m_f[i] * np.ones(nells)

    cov_clean  = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    cov_output = np.zeros((ncombs, nells, ncombs, nells)) + np.nan

    for i_tr, (i_tr1,i_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        print(i_tr)
        for j_tr, (j_tr1,j_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):

            cl1 = cell_array[i_tr1][j_tr1]
            cl2 = cell_array[i_tr2][j_tr2]
            cl3 = cell_array[i_tr1][j_tr2]
            cl4 = cell_array[i_tr2][j_tr1]

            for a_ell in range(nells):
                for b_ell in range(nells):
                    # calculate mean between cl at different ells
                    c1_l12 = np.sqrt(cl1[a_ell] * cl1[b_ell])
                    c2_l12 = np.sqrt(cl2[a_ell] * cl2[b_ell])
                    c3_l12 = np.sqrt(cl3[a_ell] * cl3[b_ell])
                    c4_l12 = np.sqrt(cl4[a_ell] * cl4[b_ell])
                    # covariance between different frequency channels
                    cov_clean[i_tr,a_ell,j_tr,b_ell] = c1_l12 * c2_l12 + c3_l12 * c4_l12

            assert np.all(np.isclose(cov_clean[i_tr, :, j_tr,:], cov_clean[i_tr,:,j_tr,:].T)), \
                    'cell4freq not diag'

            # all factors together
            prefactor = np.multiply(cov_clean[i_tr,:,j_tr,:], prefactor_l)
            cov_output[i_tr,:,j_tr,:]  = np.multiply(mcm, prefactor)

            assert np.all(np.isclose(cov_output[i_tr,:,j_tr,:], cov_output[i_tr,:,j_tr,:].T)) ,\
                    f'({i_tr},{j_tr}) output cov is not diagonal'

    return cov_output

def compute_cov(ctype, w2_factor):

    '''
    Compute covariance (with and without modulating template) of ctype power spectrum
    Covariances are unbinned, from ell = 0 to ell = 3 * NSIDE - 1

    This is because coupling matrix is given in that shape. and binning of covariance is only
    linear once the final covariance has been calculated (cannot bin elements first)

    Formula for the covariance is given in the Appendix

    ** Parameters **

    ctype: 'd00', 'dc0' or 'dcs'
        type of power spectrum
    '''

    assert len(POLARIZATION_cov) == 1, 'covariance formula is for spin 0 field only!'
    name_pol_cov = 'cl_' + 2 * POLARIZATION_cov.lower()

    # read in cells
    s_d = sacc.Sacc.load_fits(PATH_DICT['output_path'] + \
                              '_'.join([NAME_RUN, ctype]) + '_clnobin.fits')
    # read-in coupling matrices
    mw2_matrix = np.loadtxt(name_couplingmatrix_w + '.txt')
    mwt2_matrix = np.loadtxt(name_couplingmatrix_wt + '.txt')

    tr_names = sorted(list(s_d.tracers.keys()))

    # read ell of cls from the first band (tr_names)
    ellcl = s_d.get_ell_cl(name_pol_cov, tr_names[0], tr_names[0], return_cov = False)[0]

    nells = len(ellcl)
    assert nells == mw2_matrix.shape[0], 'using matrix or Cl that should not be binned'

    ncombs = len(s_d.get_tracer_combinations())
    assert ncombs == (nmodes*nfreqs) * (nmodes*nfreqs+1) //2,\
             'incorrect number of auto/cross spectra'

    c_ell_all = np.zeros((nfreqs, nfreqs, nells)) + np.nan

    # populate C_ell_all
    for i in range(nfreqs):
        for j in range(i, nfreqs):

            c_ell_all[i][j] = s_d.get_ell_cl(name_pol_cov, f'band{i+1}', f'band{j+1}')[1]

            if i != j:
                c_ell_all[j][i] = s_d.get_ell_cl(name_pol_cov, f'band{i+1}', f'band{j+1}')[1]

    for i in range(nfreqs):
        for j in range(i,nfreqs):
            assert np.all(np.isclose(c_ell_all[i][j] , c_ell_all[i][j].T)), \
                'cell(f1,f2) != cell(f2,f1)'


    cov = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    cov_w = np.zeros((ncombs, nells, ncombs, nells)) + np.nan
    cov_wt = np.zeros((ncombs, nells, ncombs, nells)) + np.nan

    indices_tr = np.triu_indices(nfreqs)

    # compute prefactor of covariance
    prefactor_l = np.zeros([nells, nells]) + np.nan

    sum_m_f = 1 / (2 * ellcl + 1)

    if COV_CORR == 'fsky':
        fsky = get_fsky()
        sum_m_f /= fsky

    elif COV_CORR == 'w2':
        sum_m_f /= ( w2_factor**2 )

    for i in range(int(nells)):
        prefactor_l[:, i ] = sum_m_f[i] * np.ones(nells)

    # populate covariance
    for i_tr, (i_tr1,i_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):
        print(i_tr)
        for j_tr, (j_tr1,j_tr2) in enumerate(zip(indices_tr[0], indices_tr[1])):

            cl1 = c_ell_all[i_tr1][j_tr1]
            cl2 = c_ell_all[i_tr2][j_tr2]
            cl3 = c_ell_all[i_tr1][j_tr2]
            cl4 = c_ell_all[i_tr2][j_tr1]

            for a_ell in range(nells):
                for b_ell in range(nells):
                    # calculate mean between cl at different ells
                    c1_l12 = np.sqrt(cl1[a_ell] * cl1[b_ell])
                    c2_l12 = np.sqrt(cl2[a_ell] * cl2[b_ell])
                    c3_l12 = np.sqrt(cl3[a_ell] * cl3[b_ell])
                    c4_l12 = np.sqrt(cl4[a_ell] * cl4[b_ell])
                    # covariance between different frequency channels
                    cov[i_tr,a_ell,j_tr,b_ell] = c1_l12 * c2_l12 + c3_l12 * c4_l12

            assert np.all(np.isclose(cov[i_tr, :, j_tr,:], cov[i_tr,:,j_tr,:].T)), \
                    'cell4freq not diag'

            # all factors together
            prefactor = np.multiply(cov[i_tr,:,j_tr,:], prefactor_l)
            cov_w[i_tr,:,j_tr,:] = np.multiply(mw2_matrix, prefactor)
            cov_wt[i_tr,:,j_tr,:] = np.multiply(mwt2_matrix, prefactor)

            assert np.all(np.isclose(cov_w[i_tr,:,j_tr,:], cov_w[i_tr,:,j_tr,:].T)) ,\
                    f'({i_tr},{j_tr}) w cov is not diagonal'
            assert np.all(np.isclose(cov_wt[i_tr,:,j_tr,:], cov_wt[i_tr,:,j_tr,:].T)) ,\
                    f'({i_tr},{j_tr}) wt cov is not diagonal'


    # save to fits file
    hdu_w = fits.PrimaryHDU(cov_w)
    hdu_w.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) + \
                    '_nobin_w.fits', overwrite = True)
    # save binned:
    cov_w_bin = bin_covs(cov_w)
    hdu_w_bin = fits.PrimaryHDU(cov_w_bin)
    hdu_w_bin.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) + \
                    '_bin_w.fits', overwrite = True)

    if ctype == 'd00':

        # save binned:
        cov_wt_bin = bin_covs(cov_wt)
        hdu_wt = fits.PrimaryHDU(cov_wt)
        hdu_wt_bin = fits.PrimaryHDU(cov_wt_bin)

        hdu_wt.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) +\
                    '_nobin_wt.fits', overwrite = True)

        hdu_wt_bin.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, ctype, 'Cov']) +\
                    '_bin_wt.fits', overwrite = True)

def get_effective_cov(bin_type):

    '''
    gets final covariance as

    Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]

    This is still unbinned, from ell = 0 to ell = 3 * NSIDE - 1
    '''

    larr_all = np.arange(3 * NSIDE)

    # read-in precomputed covariances
    hdu_dustw = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
                            f'_{bin_type}_w.fits')
    hdu_dustwt = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
                            f'_{bin_type}_wt.fits')
    hdu_all = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, NAME_COMP, 'Cov']) + \
                            f'_{bin_type}_w.fits')

    cov_dustw = hdu_dustw[0].data
    cov_dustwt = hdu_dustwt[0].data
    cov_allw = hdu_all[0].data

    if bin_type == 'nobin':
        assert cov_dustw.shape[1] == len(larr_all), 'error in ell array definition'

    # Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]
    total_cov= np.add(cov_allw, np.subtract(cov_dustwt, cov_dustw))

    # save to fits file
    hdu_cov = fits.PrimaryHDU(total_cov)
    hdu_cov.writeto(PATH_DICT['output_path'] + '_'.join([NAME_RUN, NAME_COMP, 'Cov']) + \
                    f'_{bin_type}_wt.fits', overwrite = True)

def get_crazy_cov():#covtype):

    '''

    Adds crazy covariance of dust to full Covariance

    ** Parameters **
    covtype: str, 'uniform' or 'offset'
    '''
    # N_offset = 30
    # corr_uniform = 0.1
    # corr_offset = 0.5
    # # read in precomputed gaussian cov
    # hdu_dustw = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, 'd00', 'Cov']) +\
    #                         '_nobin_w.fits')
    # hdu_all = fits.open(PATH_DICT['output_path'] + '_'.join([NAME_RUN, NAME_COMP, 'Cov']) + \
    #                         '_nobin_w.fits')

    # cov_dustw = hdu_dustw[0].data
    # cov_allw = hdu_all[0].data

    # cross, Nell = cov_dustw.shape[:2]

    # cov_dustNG = np.zeros_like(cov_dustw)

    # for x in range(cross):
    #     for y in range(cross):
    #         np.fill_diagonal(cov_dustNG[x,:,y,:], np.diag(cov_dustw[x,:,y,:]))

    #         if covtype == 'uniform':

    #             for i in range(Nell - 1):
    #                 for j in range(i+1, Nell):
    #                     cov_dustNG[x,i, y,j] = corr_uniform * np.sqrt( cov_dustNG[x,i,y,i] * cov_dustNG[x,j,y,j] )
    #                     cov_dustNG[x,j,y,i]  = cov_dustNG[x,i,y,j]

    #         if covtype == 'offset':

    #             for n in range(N_offset):
    #                 for i,j in zip(range(Nell - 1), range(n+1, Nell)):
    #                     cov_dustNG[x,i,y,j] = corr_offset * np.sqrt( cov_dustNG[x,i,y,i] * cov_dustNG[x,j,y,j] )
    #                     cov_dustNG[x,j,y,i]  = cov_dustNG[x,i,y,j]

    # # Cov = Cov(all, gaussian) + [ Cov(dust, non gaussian) - Cov(dust, gaussian) ]
    # total_cov= np.add(cov_allw, np.subtract(cov_dustNG, cov_dustw))
    # # save to fits file
    # hdu_cov = fits.PrimaryHDU(total_cov)
    # hdu_cov.writeto(PATH_DICT['output_path'] + NAME_RUN + '_nobin_fullCov.fits', overwrite = True)
    print('currently invalid function')
