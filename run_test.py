'''
first step
'''
from utils.params import config
from utils.params import NAME_COMP
from utils_bbpw.bbpw_script import write_config_yml_script
from utils_bbpw.params import dict_params_bbpw, dict_ells_bbpw, niterbands
from dustngwbbp.compute_cl import compute_cl_nobin, compute_cl
from dustngwbbp.compute_couplingmatrix import compute_couplingmatrix
from dustngwbbp.compute_cov import compute_cov, get_effective_cov # , get_crazy_cov

w2_mean = compute_couplingmatrix(**config.mask_param.__dict__)

compute_cl_nobin('d00')
compute_cl_nobin(NAME_COMP)

compute_cov('d00', w2_mean)
compute_cov(NAME_COMP, w2_mean)

get_effective_cov()

# get_crazy_cov('offset')

compute_cl(NAME_COMP,'w')
compute_cl(NAME_COMP,'wt')

for biter in niterbands:
    dict_ells_bbpw['bands'] = biter
    for template in ['w', 'wt']:
        write_config_yml_script(template, params_bbpw = dict_params_bbpw,
                                    dict_compsep= dict_ells_bbpw,
                                    cros = False, mmt = False) # dc0 or dcs, simple
        if NAME_COMP == 'dcs':
            for moment in [True, False]:
                write_config_yml_script(template, params_bbpw = dict_params_bbpw,
                                    dict_compsep= dict_ells_bbpw,
                                    cros = True, mmt = moment)
