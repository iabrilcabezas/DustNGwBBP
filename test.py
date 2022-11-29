'''
first step
'''
from utils.params import config
from dustngwbbp.compute_cl import compute_cl_nobin, compute_cl
from dustngwbbp.compute_couplingmatrix import compute_couplingmatrix
from dustngwbbp.compute_cov import compute_cov, get_effective_cov

compute_couplingmatrix(**config.mask_param.__dict__)

compute_cl_nobin('dust')
compute_cl_nobin('all')

compute_cov('dust')
compute_cov('all')

get_effective_cov()

compute_cl('all', 'w')
compute_cl('all','wt')
