'''
first step
'''
from utils.params import config
from utils.params import NAME_COMP, EXPERIMENT
from utils.bbpw_script import write_config_yml_script
from dustngwbbp.compute_cl import compute_cl_nobin, compute_cl
from dustngwbbp.compute_couplingmatrix import compute_couplingmatrix
from dustngwbbp.compute_cov import compute_cov, get_effective_cov

compute_couplingmatrix(**config.mask_param.__dict__)

compute_cl_nobin('d00')
compute_cl_nobin(NAME_COMP)

compute_cov('d00')
compute_cov(NAME_COMP)

get_effective_cov()

compute_cl(NAME_COMP,'w')
compute_cl(NAME_COMP,'wt')
                                                                                                                                                       
nside_bp = 256
niters = int(1e4)
dict_bb = {'lmin': 30, 'lmax':300}

if EXPERIMENT == 'so':
    bands_iter = ['all', ['LF1', 'LF2'], ['MF1', 'MF2'], ['UHF1', 'UHF2']]
else:
    bands_iter = ['all']

for biter in bands_iter:
    dict_bb['bands'] = biter
    for template in ['w', 'wt']:
        for mom in [True, False]:
            write_config_yml_script(template, nside_bbpw = nside_bp, dict_compsep= dict_bb,
                                    moments = mom, n_iters=niters)
            if NAME_COMP == 'dcs':
                write_config_yml_script(template, nside_bbpw = nside_bp, dict_compsep= dict_bb, 
                                        cross = True, moments = mom, n_iters=niters)
