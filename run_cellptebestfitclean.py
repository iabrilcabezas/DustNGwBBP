'''
workflow to call best fits and write .yml files corresponding to those
'''

from pte.clean2cmb import dictmodel_fromparams

path = '/global/cfs/cdirs/act/data/iabril/BBPower/230725/sims/' # 03_pte/'
w_types = ['w', 'dfwt']

for i in range(int(5e4), int(1e5)):

    if i%500 == 0:
        print(i)

    for wt in w_types:
        dictmodel_fromparams(path + f'{i}/' + wt + '/' )
