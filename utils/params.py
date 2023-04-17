'''
params
    config file read in and definition of global variables

'''
import yaml

def get_namerun(kwargs_dict):

    '''
    Converts to string the input dictionary

    ** Parameters **
    kwargs_dict: dict
        dictionary to be converted to string

    ** Returns **
    str
    '''

    return '_'.join(map(str, kwargs_dict.values()))

def get_namecouplingm(nside, ttype, mtype, **kwargs_mask):

    '''
    Returns string to name results of compute_couplingmatrix()

    ** Parameters **
    nside: int
        resolution of maps
    dell_nmt: int
        width of ell bins in namaster bandpower bins
    mtype: str
        type of mask
    ttype: str
        type of template
    apo_deg: float or NaN
        apodization scale of mask

    * Returns *
    str
    '''

    values = [nside, kwargs_mask['dell_nmt'], ttype, mtype, \
              kwargs_mask['apo_deg'], kwargs_mask['smooth_deg']]


    return '_'.join(map(str, values))

def get_namecells(experiment, nside, polarization, ctype):

    '''
    Returns string name for compute_cl_nobin()

    * Parameters *
    experiment: 'so', 'bicep', 'cmbs4','cmbs4d', 'lbrd'
        experiment to determine bandpasses, noise levels, sky coverage (fsky)
    nside: int
        resolution of maps
    polarization: 'e', 'b' or 'eb'
        polarization channels to use
    ctype: 'd00','dc0','dcs'
        components (d)ust+(c)mb+(s)ynchrotron

    * Returns *
    str
    '''

    values = [experiment, nside, polarization, ctype]

    return '_'.join(map(str, values))

def get_namecomp(comp_dict):

    '''
    Returns name string to specify components
    '''

    assert comp_dict['dust'] is True, 'no dust in your model!'

    if comp_dict['sync']:
        return 'dcs'

    return 'dc0'


def pathnames():

    '''
    Returns path to data depending on machine

    ** Parameters **
    machine: str
            machine where code runs ('cori' or 'perl')

    ** Returns **
    dict: stores path to data
        'planck_path': Planck data
        'camb_cmb_lens_nobb': CMB power spectrum, including lensing, without bb
        'output_path': location of file output storage
        'BK15_data': BK15 data (bandpasses, noise)
        'bbpipe_path': path to bbpower code. contains SO bandpasses and noise in examples folder
        'so_path': path where SO mask is stored
    '''

    dict_path = {}

    bicep_bk15 = config.path_param.bicep_bk15
    bbpipe_path = config.path_param.bbpipe
    output_path = config.path_param.output
    camb_cmb_lens_nobb_path = config.path_param.camb_cmb_lens_nobb
    so_path = config.path_param.so
    template_path = config.path_param.template

    if MACHINE == 'cori':

        dict_path['planck_data'] = config.path_param.planck_path_cori

    if MACHINE == 'perl':

        dict_path['planck_data'] = config.path_param.planck_path_perl

    dict_path['BK15_data']   = bicep_bk15
    dict_path['bbpipe_path'] = bbpipe_path
    dict_path['output_path'] = output_path
    dict_path['camb_cmb_lens_nobb'] = camb_cmb_lens_nobb_path
    dict_path['so_path'] = so_path
    dict_path['template_path'] = template_path

    return dict_path

class GlobalConfig:

    '''
    2nd level class info on config file
    contains global config parameters

    machine: 'cori' or 'perl'
        machine where the code runs
    experiment: 'so' or 'bicep'
        experiment we simulate
    nside: 2**N
        resolution of maps

    '''

    def __init__(self, param):
        self.experiment = param['experiment']
        self.nside      = param['nside']
        self.cov_corr   = param['cov_correction']
        self.template   = param['template']

class BpwConfig:

    '''
    2nd level class info on config file
    contains bandpower parameters

    lmin: int
        min ell
    nbands: int
        number of ell bands
    dell: int
        width of ell bands

    '''
    def __init__(self, param):
        self.lmin = param['lmin']
        self.nbands = param['nbands']
        self.dell = param['dell']

class PolConfig:

    '''
    2nd level class info of config file
    contains polarization parameters

    pol_cov: 'E', 'B' or 'EB'
        polarization channels to use

    '''
    def __init__(self, param):
        self.pol_cov = param['pol_cov']

class PathConfig:

    '''
    2nd level class of config file
    contains path info

    output: str
        path to location of file output storage
    bicep_bk15: str
        path to BK15 data (bandpasses, noise)
    bbpipe: str
        path to bbpower code
    planck_path_cori: str
        path to Planck data on cori machine
    planck_path_perl: str
        path to Planck data on perl machine
    camb_cmb_lens_nobb: str
        CMB power spectrum, including lensing, without bb
    so: str
        path to SO mask
    '''

    def __init__(self, param):
        self.output = param['output_path']
        self.bicep_bk15 = param['bicep_BK15']
        self.bbpipe  = param['bbpipe_path']
        self.planck_path_cori = param['planck_path_cori']
        self.planck_path_perl = param['planck_path_perl']
        self.camb_cmb_lens_nobb = param['camb_cmb_lens_nobb']
        self.so = param['so_path']
        self.template = param['template_path']

class MaskConfig:

    '''
    2nd level class of config file
    contains info on mask_param
        mask_type: str
            type of mask
        apo_deg: float
            apodization scale of mask [degrees]
        smooth_deg: float
            smoothing scale of template [degrees]
        dell_nmt: int
            width of ell bins in namaster
    '''

    def __init__(self, param):
        self.mask_type  = param['mask_type']
        self.apo_deg    = param['apo_deg']
        self.smooth_deg = param['smooth_deg']
        self.dell_nmt   = param['dell_nmt']


class CosmoConfig:

    '''
    2nd level class of config file
    contains info on cosmo_params
        CMB: CMB parametrization
        dust: dust parametrization
        model: model parameters
    '''

    def __init__(self, param):
        self.CMB = param['CMB']
        self.dust = param['dust']
        self.sync = param['sync']
        self.model = param['model']
        # self.nparams = param['nparams']

class BandConfig:

    '''
    2nd level class of config file
    contains info on band_names
    so: list
        so band names
    bicep: list
        bicep band names
    '''

    def __init__(self, param):
        self.so = param['SO']
        self.bicep = param['BICEP']

class CompConfig:

    '''
    2nd level class of config file
    contains info on model components used
    dust: bool
        contains dust
    sync: bool
        contains sync
    '''
    def __init__(self, param):
        self.dust = param['dust']
        self.sync = param['sync']
        # self.cross = param['cross']

class ExtConfig:

    '''
    2nd level class of config file
    contains info on external parameters
    machine: str, 'perl' or 'cori'
        machine where code is running
    '''

    def __init__(self, param):
        self.machine = param['machine']

class DustFilaments:

    '''
    2nd level class of config file
    Contains parameters for DustFilaments sims
    '''

    def __init__(self, param):
        self.nsims = param['nsims']
        self.paths = param['paths']

class Config:

    '''
    1st level class for config file

    global_param: global parameters
    bpw_param: bandpower parameters
    pol_param: polarization parameters
    path_param: path to files
    mask_param: mask parameters
    cosmo_param: parameters on CMB and dust P(k)
    band_names: name of bands for different experiments
    components: model components (dust, sync) used
    external: external parameters
    '''

    def __init__(self, param):
        self.global_param = GlobalConfig(param['global'])
        self.bpw_param    = BpwConfig(param['bandpowers'])
        self.pol_param    = PolConfig(param['polarization'])
        self.path_param   = PathConfig(param['paths'])
        self.mask_param   = MaskConfig(param['mask'])
        self.cosmo_param  = CosmoConfig(param['cosmology'])
        self.band_names   = BandConfig(param['band_names'])
        self.components   = CompConfig(param['components'])
        self.external     = ExtConfig(param['external'])
        self.dustfil      = DustFilaments(param['dust_filaments'])


with open('config.yaml', 'r', encoding = 'utf-8') as config_file:
    config = Config(yaml.load(config_file, yaml.FullLoader))


MACHINE          = config.external.machine
NSIDE            = config.global_param.nside
EXPERIMENT       = config.global_param.experiment
COV_CORR         = config.global_param.cov_corr
TTYPE            = config.global_param.template
MTYPE            = config.mask_param.mask_type
DELL_NMT         = config.mask_param.dell_nmt
LMIN             = config.bpw_param.lmin
DELL             = config.bpw_param.dell
NBANDS           = config.bpw_param.nbands
POLARIZATION_cov = config.pol_param.pol_cov

band_names_config = config.band_names

cosmo_params    = config.cosmo_param
cmb_params      = cosmo_params.CMB
dust_params     = cosmo_params.dust
sync_params     = cosmo_params.sync
model_params    = cosmo_params.model

# NPARAMS = cosmo_params.nparams
A_dust_BB = dust_params['A_dust_BB']
EB_dust = dust_params['EB_dust']
alpha_dust_EE = dust_params['alpha_dust_EE']
alpha_dust_BB = dust_params['alpha_dust_BB']
beta_dust = dust_params['beta_dust']
temp_dust = dust_params['temp_dust']
nu0_dust = dust_params['nu0_dust']
Alens = cmb_params['Alens']
lnorm_PL = model_params['lnorm_PL']

A_sync_BB = sync_params['A_sync_BB']
EB_sync =  sync_params['A_sync_BB']
alpha_sync_EE = sync_params['alpha_sync_BB']
alpha_sync_BB =  sync_params['alpha_sync_BB']
beta_sync =  sync_params['beta_sync']
nu0_sync = sync_params['nu0_sync']

DF_NSIMS = config.dustfil.nsims
# construct full path
DF_PATHS = config.dustfil.paths

DF_BASE_PATH  = DF_PATHS['base'] + DF_PATHS['sim'] + DF_PATHS['base_name']
DF_END_NAME_S = '_Small' + DF_PATHS['mid_name'] + f'{int(nu0_dust)}' + DF_PATHS['end_name']
DF_END_NAME_A = '_All'   + DF_PATHS['mid_name'] + f'{int(nu0_dust)}' + DF_PATHS['end_name']
DF_OUTPUT_PATH = DF_PATHS['base'] + DF_PATHS['sim'] + DF_PATHS['output']

PATH_DICT = dict(pathnames())

namerun_dict = {**config.global_param.__dict__, **config.bpw_param.__dict__, \
                **config.mask_param.__dict__, **config.pol_param.__dict__}

NAME_COUPLINGM = get_namecouplingm(NSIDE, TTYPE, MTYPE, **config.mask_param.__dict__)
NAME_COMP = get_namecomp(config.components.__dict__)
NAME_CELLS = get_namecells(EXPERIMENT, NSIDE, POLARIZATION_cov, NAME_COMP)
NAME_RUN  = get_namerun(namerun_dict) + '_' + NAME_COMP

name_couplingmatrix_w = PATH_DICT['output_path'] + NAME_COUPLINGM + '_couplingM_w' # .txt'
name_couplingmatrix_wt = PATH_DICT['output_path'] + NAME_COUPLINGM + '_couplingM_wt' #.txt'

DF_NAME = NAME_RUN + f'{int(nu0_dust)}_{abs(alpha_dust_BB)}_{A_dust_BB}'

print(NAME_RUN)
print(DF_OUTPUT_PATH)

# print(PATH_DICT)
# print(NAME_COMP)
# print(NAME_COUPLINGM)
# print(NAME_CELLS)
