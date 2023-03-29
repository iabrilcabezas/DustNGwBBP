'''
generates all outputs needed to run BBCompsep
'''

from utils.params import config
from utils.params import NAME_COMP, DF_NSIMS, NSIDE
from utils_bbpw.bbpw_script import write_config_yml_script
from utils_bbpw.params import dict_params_bbpw, dict_ells_bbpw, niterbands
from dustngwbbp.compute_cl import compute_cl_nobin, compute_cl
from dustngwbbp.compute_couplingmatrix import compute_couplingmatrix
from dustngwbbp.compute_cov import compute_cov, get_effective_cov
from dustfilaments.compute_cov import calibrate_cells, compute_full_cov, merge_cov
from dustfilaments.compute_cov import compute_cov_fromsims, compute_tildecov
from dustfilaments.compute_cov import compute_cell_nobin_dustfil

w2_mean = compute_couplingmatrix(**config.mask_param.__dict__)

compute_cl_nobin('d00')
compute_cl_nobin(NAME_COMP)

compute_cov('d00', w2_mean)
compute_cov(NAME_COMP, w2_mean)

get_effective_cov('nobin')
get_effective_cov('bin')

compute_cl(NAME_COMP,'w')
compute_cl(NAME_COMP,'wt')

## dustfilaments:
compute_cell_nobin_dustfil(DF_NSIMS, NSIDE)
calibrate_cells('nobin')

compute_cov_fromsims('small', 'nobin')

compute_tildecov('small', 'nobin')
compute_tildecov('all', 'nobin')

for TYPE_DCOV in ['df00', 'dfwt', 'dfwtm']:

    if TYPE_DCOV != 'df00':
        merge_cov(TYPE_DCOV, 'nobin')

    compute_full_cov(TYPE_DCOV, 'nobin')

    compute_cl('dcs',TYPE_DCOV)

# write config scripts

for TYPE_DCOV in ['df00', 'dfwt', 'dfwtm']:

    for biter in niterbands:
        dict_ells_bbpw['bands'] = biter
        # dict_ells_bbpw['bands'] = 'all'

        for mom, decs in zip([False, False, True],[False, True, False]):
            write_config_yml_script(TYPE_DCOV, params_bbpw = dict_params_bbpw,
                                    dict_compsep= dict_ells_bbpw,
                                    cros = mom, mmt = mom, dec = decs)

for biter in niterbands:
    dict_ells_bbpw['bands'] = biter
    #dict_ells_bbpw['bands'] = 'all'

    for template in ['w', 'wt']:
        for decs in [True, False]:
            write_config_yml_script(template, params_bbpw = dict_params_bbpw,
                                    dict_compsep= dict_ells_bbpw,
                                    cros = False, mmt = False, dec = decs) # dc0 or dcs, simple
        if NAME_COMP == 'dcs':
            write_config_yml_script(template, params_bbpw = dict_params_bbpw,
                                    dict_compsep= dict_ells_bbpw,
                                    cros = True, mmt = True, dec = False)
