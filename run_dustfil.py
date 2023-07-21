'''
runs computation of covariance from sims, and merges to exisiting cov from other components
'''

#from dustfilaments.compute_cov import compute_cell_bin_dustfil
from dustfilaments.compute_cov import calibrate_cells, compute_full_cov, merge_cov
from dustfilaments.compute_cov import compute_cov_fromsims, compute_tildecov
from dustfilaments.compute_cov import compute_cell_nobin_dustfil
from utils.params import DF_NSIMS, NSIDE
from dustngwbbp.compute_cl import compute_cl
from utils_bbpw.bbpw_script import write_config_yml_script
from utils_bbpw.params import dict_ells_bbpw, dict_params_bbpw

REGENERATE = False # = False only adds DF covariance to the output scheme

if REGENERATE:

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

        dict_ells_bbpw['bands'] = 'all'

        for mom, decs in zip([False, False, True],[False, True, False]):
            write_config_yml_script(TYPE_DCOV, params_bbpw = dict_params_bbpw,
                                    dict_compsep= dict_ells_bbpw,
                                    cros = mom, mmt = mom, dec = decs)

else:

    # dust doesnt change:
    TYPE_DCOV = 'dfwt'
    compute_full_cov(TYPE_DCOV, 'nobin')
    compute_cl('dcs',TYPE_DCOV)

    dict_ells_bbpw['bands'] = 'all'

    mom, decs = False, False 
    write_config_yml_script(TYPE_DCOV, params_bbpw = dict_params_bbpw,
                            dict_compsep= dict_ells_bbpw,
                            cros = mom, mmt = mom, dec = decs)
