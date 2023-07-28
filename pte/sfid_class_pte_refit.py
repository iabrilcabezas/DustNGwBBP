'''
sfid_class_pte
    module to compute PTE values from gaussian and non-gaussian covariances simulations
'''

import numpy as np
import sacc
from utils.params import POLARIZATION_cov
from utils.sed import get_band_names

WEIGHT = 'Cl'
band_names = get_band_names()

class SfClass_refit():

    '''
    For a given NAME_RUN, reads in fiducial Cell value and precomputed G and NG covariances
    Selects bands, pol and ell range specified by user
    Computes Cell and (inv)covariances

    Methods
    -------
    get_dof()
        Return degrees of freedom of [flattened] power spectrum

    name_band2trac()
        Returns dictionary that matches tracer names (band1, band2, ...)
        to experiment band names ('UHF1', 'UHF2', ...)

    get_cellandcov()
        According to user-specified bands, ell range and polarization, extracts
        - flattened fiducial Cell [Cl = (Clfreq1l1, Clfreq1l2, … Clfreq1lN, Clfreq2, … ClfreqN)]
        - gaussian covariance
        - non-gaussian covariances
        and also returns the inverse of the last two
    '''

    def __init__(self, nosim, bands,lmin_bbp, lmax_bbp):

        # load fiducial power spectrum for NAME_RUN parameters
        name_sf_model = f'/global/cfs/cdirs/act/data/iabril/BBPower/230725/sims/{nosim}/w/cells_model.fits'
        name_sf_random = f'/global/cfs/cdirs/act/data/iabril/BBPower/230725/sims/{nosim}/w/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_tot.fits'
        
        self.s_fg = sacc.Sacc.load_fits(name_sf_random)
        self.s_mg = sacc.Sacc.load_fits(name_sf_model)
        # select polarization and ell range
        for obj in [self.s_fg, self.s_mg]:
            obj.remove_selection(ell__gt=lmax_bbp)
            obj.remove_selection(ell__lt=lmin_bbp)
            obj.keep_selection('cl_bb')

        # name of bands [full]
        tr_names = sorted(list(self.s_fg.tracers.keys()))

        # name of bands [user selection]
        if bands == 'all':
            self.tr_names = sorted(list(self.s_fg.tracers.keys()))

        else:
            self.tr_names = sorted(list(bands))

        self.pols = [POLARIZATION_cov]   # polarization
        self.nfreqs = len(self.tr_names) # number of bands [user selection]
        nfreqs = len(tr_names)           # total no. of bands
        self.npol = len(self.pols)       # no. of polarizations
        self.nmaps = self.nfreqs * self.npol    # total no. of combinations [user]
        nmaps = nfreqs * self.npol              # total no. of combinations [full]
        self.index_ut = np.triu_indices(self.nmaps) # index covariance matrix
        self.ncross = (self.nmaps * (self.nmaps + 1)) // 2  # no. combinations [user]
        ncross = (nmaps * (nmaps + 1))//2                   # no. combinations [full]
        self.pol_order = dict(zip(self.pols, range(self.npol))) # arange polarization

        # no. of bandpowers after ell range selection
        self.ell_b = self.s_fg.get_ell_cl('cl_' + 2 * self.pols[0].lower(), \
                                            self.tr_names[0], self.tr_names[0])[0]
        self.n_bpws = len(self.ell_b)

        # degrees of freedom (no. ells * no. maps)
        self.dof = self.get_dof()

    def _freq_pol_iterator(self):
        '''From https://github.com/simonsobs/BBPower/blob/master/bbpower/compsep.py'''
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

        '''
        Return degrees of freedom of [flattened] power spectrum:
        Cl = (Clfreq1l1, Clfreq1l2, … Clfreq1lN, Clfreq2, … ClfreqN)
        '''
        # no. ell bandpowers * no. of maps
        dof = self.n_bpws * self.ncross

        return dof

    def name_band2trac(self):

        '''
        Returns dictionary that matches tracer names (band1, band2, ...)
        to experiment band names ('UHF1', 'UHF2', ...)
        '''

        values = list(self.s_fg.tracers.keys())

        return dict(zip(band_names,values))

    def get_cellandcov(self):

        '''
        According to user-specified bands, ell range and polarization, extracts
        - flattened fiducial Cell [Cl = (Clfreq1l1, Clfreq1l2, … Clfreq1lN, Clfreq2, … ClfreqN)]
        - gaussian covariance
        - non-gaussian covariances
        and also returns the inverse of the last two

        ** Returns **
        cell_array: np.array
            flattened fiducial Cell
        bbcovar_g: np.array
            gaussian covariance
        bbcovar_ng: np.array
            non-gaussian covariance
        invcov_g: np.array
            inverse of bbcovar_g
        invcov_ng: np.array
            inverse of bbcovar_ng
        '''

        v2d_fg = np.zeros([self.n_bpws, self.ncross])
        v2d_mg = np.zeros([self.n_bpws, self.ncross])
        cv2d_g = np.zeros([self.n_bpws, self.ncross, self.n_bpws, self.ncross])
        
        # From https://github.com/simonsobs/BBPower/blob/master/bbpower/compsep.py
        itr1 = self._freq_pol_iterator() # Parse into the right ordering
        for b1, b2, p1, p2, m1, m2, ind_vec in itr1:
            t1 = self.tr_names[b1]
            t2 = self.tr_names[b2]
            pol1 = self.pols[p1].lower()
            pol2 = self.pols[p2].lower()
            cl_typ = f'cl_{pol1}{pol2}'
            ind_a = self.s_fg.indices(cl_typ, (t1, t2))
            if len(ind_a) != self.n_bpws:
                raise ValueError("All power spectra need to be "
                                 "sampled at the same ells")
            v2d_fg[:, ind_vec] = self.s_fg.get_ell_cl(cl_typ, t1, t2)[1]
            v2d_mg[:, ind_vec] = self.s_mg.get_ell_cl(cl_typ, t1, t2)[1]

            itr2 = self._freq_pol_iterator()
            for b1b, b2b, p1b, p2b, m1b, m2b, ind_vecb in itr2:
                t1b = self.tr_names[b1b]
                t2b = self.tr_names[b2b]
                pol1b = self.pols[p1b].lower()
                pol2b = self.pols[p2b].lower()
                cl_typb = f'cl_{pol1b}{pol2b}'
                ind_b = self.s_fg.indices(cl_typb, (t1b, t2b))
                cv2d_g[:, ind_vec, :, ind_vecb] = self.s_fg.covariance.covmat[ind_a][:, ind_b]

        # flatten all arrays
        cell_array_fg = v2d_fg.flatten() # bug 'F')
        cell_array_mg = v2d_mg.flatten() # bug 'F')

        bbcovar_g = cv2d_g.reshape([self.n_bpws * self.ncross,
                                     self.n_bpws * self.ncross])
        invcov_g = np.linalg.solve(bbcovar_g,
                                      np.identity(len(bbcovar_g)))
        #print(cell_array_fg)
        #print(cell_array_mg)
        #print(invcov_g)
        return cell_array_fg, cell_array_mg, bbcovar_g, invcov_g

def calc_chi2(cell_random, cell_array, invcov):

    '''
    Calculates chi2 value:
        chisq = (Cl-<Cl>)Cov^{-1} (Cl-<Cl>)
    and pte (probability to exceed):
        probability of obtaining a chi2 value larger than that measured data chi2

    ** Parameters **
    cell_random: np.array
        array containing all randomly sampled cells
    cell_array: np.array
        fiducial cell, mean used to randomly sample cell_random
    invcov: np.array
        inverse covariance
    dof_chi2: int
        number of degrees of freedom of data for chi2, that is, dof - nparams

    ** Returns **
    chi2: float
        chisquare value
    pval: float
        PTE value
    '''

    deltax = cell_random - cell_array
    chi2 = np.linalg.multi_dot([deltax, invcov, deltax])

    return chi2

def get_chi2(sfclassobj):

    '''
    Computes chi2 for Cell
    calculated with a gaussian covariance

    ** Parameters **
    Sf: SfClass
        populated SfClass object

    ** Returns **
    chi_ng_array: np.array(nsims)
        chi2 values for each sim calculated with NG covariance
    '''

    # obtain cell array and covariances
    cell_array_fg, cell_array_mg, _, invcov_g = sfclassobj.get_cellandcov()
    dof_chi2 = sfclassobj.dof

    chi_ng = calc_chi2(cell_array_fg, cell_array_mg, invcov_g)

    return chi_ng, dof_chi2