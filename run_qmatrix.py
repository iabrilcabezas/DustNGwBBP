'''
this function is wrong because Q need to be regenerated for each sim
'''
import numpy as np
from pte.qmatrix import clean_cells

save_path = '/global/cfs/cdirs/act/data/iabril/BBPower/230525_sample/'
covs_name = 'so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_'

Q_ell = np.loadtxt(save_path + 'Qmatrix_clean.txt')
cov_g = np.loadtxt(save_path + 'Cov_g_clean.txt')
invcov_g = np.linalg.solve(cov_g, np.identity(len(cov_g)))
cov_ng = np.loadtxt(save_path + 'Cov_ng_clean.txt')
invcov_ng = np.linalg.solve(cov_ng, np.identity(len(cov_ng)))

nsims = int(1e5)

for i in range(6400, int(1e4)): #, nsims):

    if i%500 ==0:
        print(i)

    covs_path = save_path + f'sims/{i}/'
    sample_cell_path = covs_path + 'w/' + covs_name + 'w_tot.fits'
    
    clean_cells(Q_ell, sample_cell_path, covs_path + 'clean_sample_cell.txt')

# same for cells_model:
for i in range(6400, int(1e4)): #, nsims):

    if i%500 ==0:
        print(i)

    for wt in ['w', 'dfwt']:
        cells_path = save_path + f'sims/{i}/{wt}/'
        
        clean_cells(Q_ell, cells_path + 'cells_model.fits', cells_path + 'clean_cellbestfit.txt')

# compute chi2
invcovs = [invcov_g, invcov_ng]

chi2_array = np.zeros((2, nsims))

for i in range(nsims):
    
    path = save_path + f'sims/{i}/'

    for j, wt in enumerate(['w', 'dfwt']):
        
        cdata  =  np.loadtxt(path + 'clean_sample_cell.txt')
        cmodel =  np.loadtxt(path + wt + '/clean_cellbestfit.txt')
        
        deltax = cdata - cmodel
        chi2_array[j, i] = np.linalg.multi_dot([deltax, invcovs[j], deltax])

np.savetxt(save_path + f'chi2_array_g_{nsims}.txt', chi2_array[0])
np.savetxt(save_path + f'chi2_array_ng_{nsims}.txt', chi2_array[1])
