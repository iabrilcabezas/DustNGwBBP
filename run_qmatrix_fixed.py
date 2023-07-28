'''
better version than run_qmatrix.py because it contains no bug. Q should be recalculated for each sim
'''

import numpy as np
import sacc
import yaml
from utils.sed import get_convolved_seds_varpar, get_band_names
from utils.bandpowers import get_ell_arrays
from dustngwbbp.compute_cl import get_windows, import_bandpasses
import utils.noise_calc as nc
from utils.params import LMIN, DELL, NBANDS, POLARIZATION_cov

band_names = get_band_names()
LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

nfreqs = len(band_names)
nmodes = len(POLARIZATION_cov)
nmaps = nfreqs * nmodes
ncross = (nmaps * (nmaps + 1)) // 2

indices_tr = np.triu_indices(nfreqs)
ncombs = len(indices_tr[0])
assert ncross == ncombs

NCOMP = 3
TYPE_COV = 'dfwt'
WEIGHT = 'Cl'

# for matrices index afterwards
counter = np.arange(len(indices_tr[0]))
dicttr = {}
for i, counter_val in enumerate(counter):
    dicttr[(indices_tr[0][i], indices_tr[1][i])] = counter_val

fsky = nc.get_fsky()
bpss = import_bandpasses()

sens=2
knee=1
ylf=1

nell=np.zeros([nfreqs,LMAX+1])
_,nell[:,2:],_=nc.Simons_Observatory_V3_SA_noise(sens,knee,ylf,fsky,LMAX+1,1, atm_noise = True)

windows = get_windows(WEIGHT)
n_bpw=np.sum(nell[:,None,:]*windows[None,:,:],axis=2)
noise_array = n_bpw

invnoise_ell = np.zeros((NBANDS, nfreqs, nfreqs)) + np.nan

for i in range(NBANDS):
    invnoise_ell[i,:,:]  = np.linalg.inv( np.diag( noise_array[:,i] ) )  # nfreq x nfreq


def Q_from_S(S):
    
    A_ell = np.zeros(( NBANDS, NCOMP, NCOMP)) + np.nan
    B_ell = np.zeros(( NBANDS, NCOMP, nfreqs)) + np.nan
    for i in range(NBANDS):
        A_ell[i,:,:] = np.linalg.inv(np.matmul( S, np.matmul(invnoise_ell[i], S.transpose()) ))
        B_ell[i,:,:] = np.matmul( S, invnoise_ell[i] )
    Q_ell = np.zeros(( NBANDS, nfreqs )) + np.nan
    for i in range(NBANDS):
        for j in range(nfreqs):
            Q_ell[i,j] = sum(A_ell[i,0,:] * B_ell[i,:,j])

    return Q_ell

def clean_cov(Q_ell, cov):
    cov_clean = np.zeros((NBANDS, NBANDS))

    for l1 in range(NBANDS): 
        for l2 in range(NBANDS):
            suma = 0
            for n1 in range(nfreqs):
                for n2 in range(nfreqs):
                    cros1 = (n1,n2)
                    cros1_val = dicttr[tuple(sorted(cros1))]
                    for n3 in range(nfreqs):
                        for n4 in range(nfreqs):
                            cros2 = (n3,n4)
                            cros2_val = dicttr[tuple(sorted(cros2))]
                            qval = (Q_ell[l1, n1] * Q_ell[l1, n2] * Q_ell[l2, n3] * Q_ell[l2, n4]) 
                            suma  += qval * cov[ cros1_val, l1, cros2_val, l2]
            cov_clean[l1,l2] = suma
        
    return cov_clean

def get_clean_cells(Q, nosim, covtype, cell_name ):
    
    cell_sampled = sacc.Sacc.load_fits(path_sims + f'{nosim}/{covtype}/{cell_name}').mean
    cell_sampled = cell_sampled.reshape((ncross, NBANDS))

    # reshape into cell(nu, nuprime) shape
    cell_nunup_ell = np.zeros((NBANDS, nfreqs, nfreqs))
    for i in range(NBANDS):
        # https://stackoverflow.com/questions/17527693/transform-the-upper-lower-triangular-part-of-a-symmetric-matrix-2d-array-into
        X = np.zeros((nfreqs, nfreqs))
        X[np.triu_indices(X.shape[0], k = 0)] = cell_sampled[:,i]
        X = X + X.T - np.diag(np.diag(X))
        cell_nunup_ell[i] = X
        
    # clean:
    clean_cell = np.zeros(NBANDS) + np.nan
    for i in range(NBANDS):
        clean_cell[i] = np.matmul( (Q[i,:]).transpose(), np.matmul(cell_nunup_ell[i], Q[i,:]) )
    
    return clean_cell

def get_clean_chi2(nosim, covtype):
    
    with open(path_sims + f'{nosim}/{covtype}/bestfitcell.yml') as f:
        data = yaml.load(f, yaml.SafeLoader)

    beta_d = float(data['BBCompSep']['fg_model']['component_1']['sed_parameters']['beta_d'][-1][0])
    temp_d = float(data['BBCompSep']['fg_model']['component_1']['sed_parameters']['temp_d'][-1][0])
    beta_s = float(data['BBCompSep']['fg_model']['component_2']['sed_parameters']['beta_s'][-1][0])

    S = get_convolved_seds_varpar(band_names, bpss, beta_d, temp_d, beta_s)
    Q = Q_from_S(S)

    Cov = sacc.Sacc.load_fits(path_sims + f'{nosim}/{covtype}/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_{covtype}_tot.fits').covariance.covmat
    Cov  = Cov.reshape([ncross, NBANDS, ncross, NBANDS])

    cov_cleancmb = clean_cov(Q, Cov)
    invcov_cleancmb = np.linalg.solve(cov_cleancmb, np.identity(len(cov_cleancmb)))

    clean_sample_cell = get_clean_cells(Q, nosim, covtype, f'so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_{covtype}_tot.fits')
    clean_best_fit = get_clean_cells(Q, nosim, covtype, 'cells_model.fits')

    deltax = clean_sample_cell - clean_best_fit
    chi2 = np.linalg.multi_dot([deltax, invcov_cleancmb, deltax])

    return chi2

path_sims = '/global/cfs/cdirs/act/data/iabril/BBPower/230725/sims/'
nsims = int(1e4)
i0 = 0

chi2_array = np.zeros((nsims, 2))

for i_ns, ns in enumerate(range(i0, nsims)):
    if ns%500 == 0:
        print(ns)
    for i_cov, covs in enumerate(['w', 'dfwt']):
        chi2_array[i_ns, i_cov] = get_clean_chi2(ns, covs)
        
np.savetxt(path_sims + f'/../chi2_array_cleaned_{i0}_{nsims}_v2.txt', chi2_array)
