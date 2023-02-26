'''
runs computation of covariance from sims, and merges to exisiting cov from other components
'''

from dustfilaments.compute_cov import compute_cell_dustfil
from dustfilaments.compute_cov import calibrate_cells, merge_cov, compute_full_cov
from utils_dustfil.params import DF_NSIMS
from dustngwbbp.compute_cl import compute_cl
from utils_bbpw.bbpw_script import write_config_yml_script
from utils_bbpw.params import dict_ells_bbpw, dict_params_bbpw

compute_cell_dustfil(DF_NSIMS)
calibrate_cells()
merge_cov()
compute_full_cov()

compute_cl('dcs','df')

dict_ells_bbpw['bands'] = 'all'
write_config_yml_script('df', params_bbpw = dict_params_bbpw,
                            dict_compsep= dict_ells_bbpw,
                            cros = False, mmt = False) # dc0 or dcs, simple