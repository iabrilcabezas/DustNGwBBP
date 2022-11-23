'''
params
'''

def get_namerun(kwargs_dict):

    '''
    Returns string that contains all arguments used in the code run

    ** Parameters **
    kwargs_dict: dict
            'width_l'
            'apo_deg'
            'smooth_deg'

    ** Returns **
    str
    '''

    name_run = '_'.join(map(str, kwargs_dict.values()))

    return name_run

def pathnames(machine):

    '''
    Returns path to data depending on machine

    ** Parameters **
    machine: str
            machine where code runs ('cori' or 'perl')

    **Returns**
    dict
        stores path to data
        'planck_path': Planck data
    '''

    path_dict = {}

    bicep_BK15 = '/global/cfs/cdirs/act/data/iabril/BICEP/'
    bbpipe_path = '/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower/'

    if machine == 'cori':

    #   figure_path = '/global/cscratch1/sd/iabril/plots/'
        planck_path = '/global/cscratch1/sd/iabril/PlanckData/'
        path_dict['planck_data'] = planck_path

    if machine == 'perl':

    #   figure_path = '/pscratch/sd/i/iabril/plots/'
        planck_path = '/pscratch/sd/i/iabril/data/PlanckData/'
        path_dict['planck_data'] = planck_path

    path_dict['BK15_data']   = bicep_BK15
    path_dict['bbpipe_path'] = bbpipe_path

    return path_dict
