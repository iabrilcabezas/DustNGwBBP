'''
first step
'''
import yaml

from utils.params import Config, pathnames
from dustNG.compute_cl import compute_cl_forcov
from dustNG.compute_couplingmatrix import compute_couplingmatrix

config = Config(yaml.load(open('config.yaml'), yaml.FullLoader))

MACHINE = config.global_param.machine
NSIDE   = config.global_param.nside
EXPERIMENT = config.global_param.experiment

path_dict = dict(pathnames(MACHINE))
## TODO: do i need the last part??? i changed dell_nmt place
args_dict = {**config.mask_param.__dict__, **config.bpw_param.__dict__}

compute_couplingmatrix(**args_dict)

s_dust = compute_cl_forcov('dust')

print("Writing")
s_dust.save_fits(path_dict['output_path'] + '_'.join([str(NSIDE), EXPERIMENT]) +\
                 "_cls_dust_bb.fits", overwrite=True)


