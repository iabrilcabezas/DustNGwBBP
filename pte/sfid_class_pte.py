'''
sfid_class_pte
    module to compute PTE values from gaussian and non-gaussian covariances simulations
'''

import numpy as np
from astropy.io import fits
import sacc
from scipy import stats

from utils.params import NAME_CELLS, PATH_DICT, NAME_RUN, NAME_COUPLINGM, NAME_COMP
from utils.params import NSIDE, POLARIZATION_cov

from utils.binning import rebin, cut_array
from utils.sed import get_band_names

WEIGHT = 'Cl'
band_names = get_band_names()

class SfClass():

    def __init__(self, bands,lmin_bbp, lmax_bbp):

        name_sf = PATH_DICT['output_path'] + '_'.join([NAME_RUN, WEIGHT, 'wt']) + '_fid.fits'

        self.s_f = sacc.Sacc.load_fits(name_sf)
        self.s_fg = sacc.Sacc.load_fits(name_sf)
        self.s_fng = sacc.Sacc.load_fits(name_sf)

        cov_ng_full = fits.open(PATH_DICT['output_path'] + NAME_RUN +'_nobin_fullCov.fits')[0].data
        cov_g_full  = fits.open(PATH_DICT['output_path'] + \
                            '_'.join([NAME_CELLS, NAME_COUPLINGM, NAME_COMP, 'Cov']) +\
                            '_nobin_w.fits')[0].data

        for obj in [self.s_f, self.s_fg, self.s_fng]:
            obj.remove_selection(ell__gt=lmax_bbp)
            obj.remove_selection(ell__lt=lmin_bbp)
            obj.keep_selection('cl_bb')

        tr_names = sorted(list(self.s_f.tracers.keys()))

        if bands == 'all':
            self.tr_names = sorted(list(self.s_f.tracers.keys()))
        else:
            self.tr_names = bands

        self.pols = [POLARIZATION_cov]
        self.nfreqs = len(self.tr_names)
        nfreqs = len(tr_names)
        self.npol = len(self.pols)
        self.nmaps = self.nfreqs * self.npol
        nmaps = nfreqs * self.npol
        self.index_ut = np.triu_indices(self.nmaps)
        self.ncross = (self.nmaps * (self.nmaps + 1)) // 2
        ncross = (nmaps * (nmaps + 1))//2
        self.pol_order = dict(zip(self.pols, range(self.npol)))

        self.ell_b, _ = self.s_f.get_ell_cl('cl_' + 2 * self.pols[0].lower(), \
                                            self.tr_names[0], self.tr_names[0])
        self.n_bpws = len(self.ell_b)

        cov_g = np.zeros([ncross, self.n_bpws, ncross, self.n_bpws])
        cov_ng =  np.zeros([ncross, self.n_bpws, ncross, self.n_bpws])
        # cut and bin COV MW matrix:
        for i_tr in range(ncross):
            for j_tr in range(ncross):

                cov_g[i_tr,:, j_tr,:] = rebin(cut_array(cov_g_full[i_tr,:, j_tr,:], \
                                                    np.arange(3 * NSIDE), lmin_bbp, lmax_bbp), \
                                                    [self.n_bpws, self.n_bpws])
                cov_ng[i_tr,:, j_tr,:] = rebin(cut_array(cov_ng_full[i_tr,:, j_tr,:], \
                                                    np.arange(3 * NSIDE), lmin_bbp, lmax_bbp), \
                                                    [self.n_bpws,self.n_bpws])

        cov_g = cov_g.reshape([ncross * self.n_bpws, ncross * self.n_bpws ])
        cov_ng = cov_ng.reshape([ncross* self.n_bpws,ncross * self.n_bpws ])

        self.s_fg.add_covariance(cov_g)
        self.s_fng.add_covariance(cov_ng)

        self.dof = self.get_dof()

    def _freq_pol_iterator(self):
        icl = -1
        for b1 in range(self.nfreqs):
            for p1 in range(self.npol):
                m1 = p1 + self.npol * b1
                for b2 in range(b1, self.nfreqs):
                    if b1 == b2:
                        p2_r = range(p1, self.npol)
                    else:
                        p2_r = range(self.npol)
                    for p2 in p2_r:
                        m2 = p2 + self.npol * b2
                        icl += 1
                        yield b1, b2, p1, p2, m1, m2, icl

    def get_dof(self):

        dof = self.n_bpws * self.ncross

        return dof

    def name_band2trac(self):

        values = list(self.s_f.tracers.keys())

        return dict(zip(band_names,values))

    def get_cellandcov(self):

        v2d_fid = np.zeros([self.n_bpws, self.ncross])
        cv2d_g = np.zeros([self.n_bpws, self.ncross, self.n_bpws, self.ncross])
        cv2d_ng = np.zeros([self.n_bpws, self.ncross, self.n_bpws, self.ncross])

        # Parse into the right ordering
        itr1 = self._freq_pol_iterator()

        for b1, b2, p1, p2, m1, m2, ind_vec in itr1:
            t1 = self.tr_names[b1]
            t2 = self.tr_names[b2]
            pol1 = self.pols[p1].lower()
            pol2 = self.pols[p2].lower()
            cl_typ = f'cl_{pol1}{pol2}'
            ind_a = self.s_f.indices(cl_typ, (t1, t2))
            if len(ind_a) != self.n_bpws:
                raise ValueError("All power spectra need to be "
                                 "sampled at the same ells")
            _, v2d_fid[:, ind_vec] = self.s_f.get_ell_cl(cl_typ, t1, t2)

            itr2 = self._freq_pol_iterator()
            for b1b, b2b, p1b, p2b, m1b, m2b, ind_vecb in itr2:
                t1b = self.tr_names[b1b]
                t2b = self.tr_names[b2b]
                pol1b = self.pols[p1b].lower()
                pol2b = self.pols[p2b].lower()
                cl_typb = f'cl_{pol1b}{pol2b}'
                ind_b = self.s_f.indices(cl_typb, (t1b, t2b))
                cv2d_g[:, ind_vec, :, ind_vecb] = self.s_fg.covariance.covmat[ind_a][:, ind_b]
                cv2d_ng[:, ind_vec, :, ind_vecb] = self.s_fng.covariance.covmat[ind_a][:, ind_b]

        cell_array = v2d_fid.flatten('F')

        bbcovar_g = cv2d_g.reshape([self.n_bpws * self.ncross,
                                     self.n_bpws * self.ncross])
        bbcovar_ng = cv2d_ng.reshape([self.n_bpws * self.ncross,
                                     self.n_bpws * self.ncross])
        invcov_g = np.linalg.solve(bbcovar_g,
                                      np.identity(len(bbcovar_g)))
        invcov_ng = np.linalg.solve(bbcovar_ng,
                                      np.identity(len(bbcovar_ng)))

        return cell_array, bbcovar_g, bbcovar_ng, invcov_g, invcov_ng

def calc_chi2(cell_random, cell_array, invcov, dof_chi2):

    deltax = cell_random - cell_array
    chi2 = np.linalg.multi_dot([deltax, invcov, deltax])
    pval = stats.chi2.sf(chi2, df = dof_chi2)

    return (chi2, pval)

def ptechi2_gvsng(nsims, Sf):

    cell_array, _, bbcovar_ng, invcov_g, invcov_ng = Sf.get_cellandcov()
    dof_chi2 = Sf.dof - 1

    chi_ng_array = np.zeros(nsims) + np.nan
    p_ng_array = np.zeros(nsims) + np.nan
    chi_g_array = np.zeros(nsims) + np.nan
    p_g_array = np.zeros(nsims) + np.nan

    random_cells = np.random.default_rng().multivariate_normal(cell_array, bbcovar_ng, size = nsims)

    for i in range(nsims):
        chi_g_array[i], p_g_array[i] = calc_chi2(random_cells[i], cell_array, invcov_g, dof_chi2)
        chi_ng_array[i], p_ng_array[i] = calc_chi2(random_cells[i], cell_array, invcov_ng, dof_chi2)

    return chi_g_array, chi_ng_array, p_g_array, p_ng_array
