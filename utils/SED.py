'''
SED
'''
import numpy as np
from utils.params import path_dict, band_names_config, A_dust_BB, EB_dust, alpha_dust_BB, alpha_dust_EE, beta_dust, temp_dust, nu0_dust, Alens, lnorm_PL

def get_band_names(experiment):

    '''
    Define name of frequency bands according to experiment

    ** Parameters**
    experiment: str
                'so' or 'bicep'

    **Returns**
    list of band names
    '''
    if experiment == 'so':
        band_names = band_names_config.so
    if experiment == 'bicep':
        band_names = band_names_config.bicep

    return band_names

#CMB spectrum
def fcmb(nu):
    x = 0.017608676067552197*nu
    ex = np.exp(x)
    return ex*(x/(ex-1))**2

#All spectra
def comp_sed(nu,nu0,beta,temp,typ):
    if typ == 'cmb':
        return fcmb(nu)
    elif typ == 'dust':
        x_to=0.04799244662211351*nu/temp
        x_from=0.04799244662211351*nu0/temp
        return (nu/nu0)**(1+beta)*np.expm1(x_from)/np.expm1(x_to)*fcmb(nu0)
    return None

#Component power spectra
def dl_plaw(A,alpha,ls):
    return A*((ls+0.001)/lnorm_PL)**alpha

def read_camb(fname, lmax):
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
class Bpass(object):
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
        sed = np.sum(self.dnu*self.bnu*self.nu**2*f(self.nu))
        return sed

def get_component_spectra(lmax):
    
    larr_all = np.arange(lmax+1)
    
    dls_dust_ee=dl_plaw(A_dust_BB*EB_dust,alpha_dust_EE,larr_all)
    dls_dust_bb=dl_plaw(A_dust_BB,alpha_dust_BB,larr_all)
    
    _,dls_cmb_ee,dls_cmb_bb,_=read_camb( path_dict['camb_cmb_lens_nobb'], lmax)
    
    return (dls_dust_ee, dls_dust_bb,
            dls_cmb_ee, Alens*dls_cmb_bb)

def get_convolved_seds(names, bpss):
    nfreqs = len(names)
    seds = np.zeros([2,nfreqs])
    for ib, n in enumerate(names):
        b = bpss[n]
        seds[0,ib] = b.convolve_sed(lambda nu : comp_sed(nu,None,None,None,'cmb'))
        seds[1,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_dust,beta_dust,temp_dust,'dust'))
    return seds
