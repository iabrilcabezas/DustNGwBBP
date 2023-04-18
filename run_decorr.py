'''
prepares decorr data to be analyzed with BBCompSep
'''

from decorr.add_cov2decorr import add_cov2decorr

TYPE_COV = 'dfwt'
WEIGHT = 'Cl'

add_cov2decorr(TYPE_COV,WEIGHT)
