'''
script to sample from cov and generate sacc
'''
import os
from sample_pte.sample_to_sacc import sacc_from_sample
from utils.params import PATH_DICT, NAME_RUN, NAME_COMP

INPUT_CL = '/global/cfs/cdirs/act/data/iabril/BBPower/230525_sample/cells_model.fits'
COV_G_PATH = PATH_DICT['input_path'] + '_'.join([ NAME_RUN, NAME_COMP]) + \
                                '_' + '_'.join(['Cov', 'bin', 'w']) + '.fits'
COV_NG_PATH = PATH_DICT['input_path'] + '_'.join([ NAME_RUN, NAME_COMP]) + \
                                '_' + '_'.join(['Cov', 'bin', 'dfwt']) + '.fits'

TYPE_COV = 'dfwt'
WEIGHT = 'Cl'

for i in range(int(5e4), int(1e5)):
    print(i)

    newpath = PATH_DICT['output_path'] + f'sims/{i}/'

    if not os.path.exists(newpath):
        os.makedirs(newpath)
        os.makedirs(newpath + '/w/')
        os.makedirs(newpath + f'/{TYPE_COV}/')
        os.makedirs(newpath + '/w/LF/')
        os.makedirs(newpath + '/w/MF/')
        os.makedirs(newpath + '/w/UHF/')

    s_d_g, s_d_ng = sacc_from_sample(INPUT_CL, COV_G_PATH, COV_NG_PATH, test = False)

    s_d_g.save_fits(PATH_DICT['output_path'] + f'sims/{i}/w/' + \
                   '_'.join([NAME_RUN, WEIGHT, 'w']) + \
                   '_tot.fits', overwrite = True)
    s_d_ng.save_fits(PATH_DICT['output_path'] + f'sims/{i}/{TYPE_COV}/' + \
                   '_'.join([NAME_RUN, WEIGHT, TYPE_COV]) + \
                   '_tot.fits', overwrite = True)
