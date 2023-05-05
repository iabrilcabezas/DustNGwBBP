from pte.clean2cmb import dictmodel_fromparams

path = '/global/cfs/cdirs/act/data/iabril/BBPower/230503_pte/'
w_types = ['w', 'dfwt']

for i in range(1000):

    for wt in w_types:
        dictmodel_fromparams(path + f'{i}/' + wt + '/' )
