'''
first step
'''
from utils.params import config, path_dict, MACHINE, NSIDE, EXPERIMENT
from dustNG.compute_cl import compute_cl
from dustNG.compute_couplingmatrix import compute_couplingmatrix
from dustNG.compute_cov import compute_cov, get_effective_cov

compute_couplingmatrix(**config.mask_param.__dict__)

s_dust = compute_cl('dust', 'Cl', False)
s_all = compute_cl('all', 'Cl', False)

compute_cov('dust')
compute_cov('all')

get_effective_cov()

s_all_Dl = compute_cl('all', 'Dl', True)
