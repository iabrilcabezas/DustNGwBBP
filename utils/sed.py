'''
sed
    computes component spectra and SEDs
'''

import numpy as np
from scipy.constants import physical_constants
from utils.params import PATH_DICT, EXPERIMENT, band_names_config
from utils.params import A_dust_BB, EB_dust, alpha_dust_BB, alpha_dust_EE
from utils.params import beta_dust, temp_dust, nu0_dust
from utils.params import Alens, lnorm_PL
from utils.params import A_sync_BB, EB_sync, alpha_sync_BB, alpha_sync_EE
from utils.params import beta_sync, nu0_sync

h_JHz = physical_constants['Planck constant'][0]
k_JK  = physical_constants['Boltzmann constant'][0]
T_CMB =  2.72548 # fFixser 2009

def get_band_names():

    '''
    Define name of frequency bands according to experiment

    ** Parameters**
    experiment: str
                'so' or 'bicep' or 'cmbs4', 'cmbs4d'

    **Returns**
    list of band names
    '''
    if EXPERIMENT in ['so', 'cmbs4', 'cmbs4d', 'lbrd']:
        band_names = band_names_config.so
    if EXPERIMENT == 'bicep':
        band_names = band_names_config.bicep

    return band_names

### from https://github.com/simonsobs/BBPower/blob/master/examples/utils.py

#CMB spectrum
def fcmb(nu):

    '''
    spectral energy density in CMB units

    ** parameters **
    nu: float
        frequency in GHz
    '''

    x = h_JHz * 1e9/(k_JK*T_CMB) *nu
    return np.exp(x)*(x/(np.expm1(x)))**2

#All spectra
def comp_sed(nu,nu0,beta,temp,typ):

    '''
    SED for CMB and dust components
    '''

    if typ == 'cmb':
        return fcmb(nu)
    if typ == 'dust':
        x_to=0.04799244662211351*nu/temp
        x_from=0.04799244662211351*nu0/temp
        return (nu/nu0)**(1+beta)*np.expm1(x_from)/np.expm1(x_to)*fcmb(nu0)
    elif typ == 'sync':
        return (nu/nu0)**beta*fcmb(nu0)
    return None

#Component power spectra
def dl_plaw(A,alpha,ls):
    '''power law of index alpha'''
    return A*((ls+0.001)/lnorm_PL)**alpha

def read_camb(fname, lmax):
    ''''reads camb data for CMB P(k)'''
    larr_all = np.arange(lmax+1)
    l,dtt,dee,dbb,dte = np.loadtxt(fname,unpack=True)
    l = l.astype(int)
    msk = l <= lmax
    l = l[msk]
    dltt = np.zeros(len(larr_all))
    dltt[l] = dtt[msk]
    dlee = np.zeros(len(larr_all))
    dlee[l] = dee[msk]
    dlbb = np.zeros(len(larr_all))
    dlbb[l] = dbb[msk]
    dlte = np.zeros(len(larr_all))
    dlte[l] = dte[msk]
    return dltt,dlee,dlbb,dlte


#Bandpasses
class Bpass:
    '''bandpass class object'''
    def __init__(self,name,fname):
        self.name = name
        self.nu,self.bnu = np.loadtxt(fname, usecols = (0,1), unpack= True)
        self.dnu = np.zeros_like(self.nu)
        self.dnu[1:] = np.diff(self.nu)
        self.dnu[0] = self.dnu[1]
        # CMB units
        norm = np.sum(self.dnu*self.bnu*self.nu**2*fcmb(self.nu))
        self.bnu /= norm

    def convolve_sed(self,f):
        '''convolves SED to get normalization'''
        sed = np.sum(self.dnu*self.bnu*self.nu**2*f(self.nu))
        return sed

def get_component_spectra(lmax):
    '''gets component spectra of CMB and dust, E and B polarization'''

    larr_all = np.arange(lmax+1)

    dls_dust_ee=dl_plaw(A_dust_BB*EB_dust,alpha_dust_EE,larr_all)
    dls_dust_bb=dl_plaw(A_dust_BB,alpha_dust_BB,larr_all)

    dls_sync_ee=dl_plaw(A_sync_BB*EB_sync,alpha_sync_EE,larr_all)
    dls_sync_bb=dl_plaw(A_sync_BB,alpha_sync_BB,larr_all)

    _,dls_cmb_ee,dls_cmb_bb,_=read_camb( PATH_DICT['camb_cmb_lens_nobb'], lmax)

    return (dls_dust_ee, dls_dust_bb,
            dls_cmb_ee, Alens*dls_cmb_bb,
            dls_sync_ee, dls_sync_bb)

def get_convolved_seds(names, bpss):
    '''convolves SED with bandpasses'''
    nfreqs = len(names)
    seds = np.zeros([3,nfreqs])
    for ib, n in enumerate(names):
        b = bpss[n]
        seds[0,ib] = b.convolve_sed(lambda nu : comp_sed(nu,None,None,None,'cmb'))
        seds[1,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_dust,beta_dust,temp_dust,'dust'))
        seds[2,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_sync,beta_sync,None,'sync'))
    return seds

def get_convolved_seds_varpar(names, bpss, beta_d, temp_d, beta_s):

    '''
    parameters SED read-in, instead of default sky model
    '''

    nfreqs = len(names)
    seds = np.zeros([3,nfreqs])
    for ib, n in enumerate(names):
        b = bpss[n]
        seds[0,ib] = b.convolve_sed(lambda nu : comp_sed(nu,None,None,None,'cmb'))
        seds[1,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_dust,beta_d,temp_d,'dust'))
        seds[2,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_sync,beta_s,None,'sync'))
    return seds
