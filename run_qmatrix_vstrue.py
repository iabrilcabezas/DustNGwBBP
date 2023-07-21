
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

# sample cells already cleaned -- great

# need to clean true value:
clean_cells(Q_ell, save_path + 'cells_model.fits', save_path + 'clean_celltrue.txt')
print('cleaned true')
# compute chi2
invcovs = [invcov_g, invcov_ng]

chi2_array = np.zeros((2, nsims))

# load the best fit
cmodel =  np.loadtxt(save_path + '/clean_celltrue.txt')

for i in range(nsims):

    if (i %500) ==0:
        print(i)
    path = save_path + f'sims/{i}/'

    for j, wt in enumerate(['w', 'dfwt']):

        cdata  =  np.loadtxt(path + 'clean_sample_cell.txt')

        deltax = cdata - cmodel
        chi2_array[j, i] = np.linalg.multi_dot([deltax, invcovs[j], deltax])

np.savetxt(save_path + f'chi2_array_true_g_{nsims}.txt', chi2_array[0])
np.savetxt(save_path + f'chi2_array_true_ng_{nsims}.txt', chi2_array[1])
