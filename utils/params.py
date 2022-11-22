def params2str(NSIDE, experiment, mtype):

    return '_'.join([str(NSIDE), experiment, mtype])

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

    if machine == 'cori':

    #   figure_path = '/global/cscratch1/sd/iabril/plots/'
        planck_path = '/global/cscratch1/sd/iabril/PlanckData/'

    if machine == 'perl':

    #   figure_path = '/pscratch/sd/i/iabril/plots/'
        planck_path = '/pscratch/sd/i/iabril/data/PlanckData/'

    path_dict['planck_data'] = planck_path

    return path_dict
