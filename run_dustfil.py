'''
runs computation of covariance from sims, and merges to exisiting cov from other components
'''

from dustfilaments.compute_cov import compute_cell_dustfil
from dustfilaments.compute_cov import calibrate_cells, compute_full_cov, merge_cov
from dustfilaments.compute_cov import compute_cov_fromsims, compute_tildecov
from utils_dustfil.params import DF_NSIMS, TYPE_DCOV
from dustngwbbp.compute_cl import compute_cl
from utils_bbpw.bbpw_script import write_config_yml_script
from utils_bbpw.params import dict_ells_bbpw, dict_params_bbpw

# compute_cell_dustfil(DF_NSIMS)
calibrate_cells()

compute_cov_fromsims('small')
compute_tildecov('small')
compute_tildecov('all')

merge_cov('david')

compute_full_cov(TYPE_DCOV)

compute_cl('dcs',TYPE_DCOV)

dict_ells_bbpw['bands'] = 'all'

write_config_yml_script(TYPE_DCOV, params_bbpw = dict_params_bbpw,
                            dict_compsep= dict_ells_bbpw,
                            cros = False, mmt = False, dec = False) # dc0 or dcs, simple

write_config_yml_script(TYPE_DCOV, params_bbpw = dict_params_bbpw,
                            dict_compsep= dict_ells_bbpw,
                            cros = False, mmt = False, dec = True) # dc0 or dcs, simple

write_config_yml_script(TYPE_DCOV, params_bbpw = dict_params_bbpw,
                            dict_compsep= dict_ells_bbpw,
                            cros = True, mmt = True, dec = False) # dc0 or dcs, simple
