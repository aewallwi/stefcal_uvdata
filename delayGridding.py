import numpy as np
import numpy.fft as fft
import scipy.signal as signal
from astropy.cosmology import Planck15 as pl15
pi=np.pi
c=299792458. #speed of light
f21=1.42040575177e9 #21 cm frequency
DH=c/(pl15.H0.value*1e3)#hubble distance Mpc
DHcm=DH*3.08567758e24#hubble distance in cm
kb=1.38e-23 #Boltzmann constant
jy=1e-26 #Jy to SI
mH=1.6726231e-24 #mass of hydrogen atom in kg. 
littleh=pl15.H0.value/100.


'''
X(f) where k_perp = 2*pi*u/X
'''
def X(f):
    z=f21/f-1
    return pl15.comoving_distance(z).value
def Y(f):
    z=f21/f-1
    return DH*(1+z)*(1+z)*de(z)/f21
def de(z):
    return 1./np.sqrt(pl15.Om0*(1+z)**3.+1.-pl15.Om0)
def f2z(f):
    return f21/f-1.
#convert from surface temperature to surface brightness in jy
def i2t(f,i):
    return i*(c/f)**2/(2*kb)*jy
def u2kperp(u,z):
    return u*2.*pi/pl15.comoving_distance(z).value

def eta2kpara(eta,z):
    return (2.*pi*f21)/(DH*(1.+z)**2.*de(z))*eta

#method to compute delay-transform of visibilities. 
def delayTransformVisibilities(visList,df):
    '''
     ARGS:
     visList: nfreqxnvis complex matrix of visibilities
     df: float indicating channel width
     RETURNS:
     ndelayxnvis grid of delay-transformed visibilities
    '''
    nf=visList.shape[0]
    nv=visList.shape[1]
    windowGrid=np.zeros_like(visList)
    for vNum in range(nv):
        windowGrid[:,vNum]=signal.blackmanharris(nf)
        windowGrid[:,vNum]/=np.sqrt(np.mean(windowGrid[:,vNum]**2.))
    return df*fft.fftshift(fft.fft(fft.fftshift(visList*windowGrid,axes=[0]),axis=0),axes=[0])



#method to grid delay-transformed visibilities
def gridDelay(uAIn,tAIn,ftVis,nCells,uMax,visFlags=None):
    '''
    ARGS:
    uAin, nvis vector of u vales for individual delay-transformed visibilities
    tAin, ndelay vector of delay values for individual delay-transformed visibilities
    ftVis, ndelay x nvis matrix of delay-transformed visibilities
    nCells, number of bins in u-Axis. 
    uMax, maximum u for gridded delay-transformed visibilities
    visFlags, nvis vector of booleans indicating delay-transformed visibilities that you don't want to be included in gridding. 
    RETURNS:
    ncell vector of u values, ndelay/2 (averages negative/positive delays) axis of tau values, ndelayxncell gridded delay-transformed visibilities 
    '''
    dGrid=np.zeros((ftVis.shape[0]/2,nCells))
    counts=np.zeros(nCells).astype(int)
    nu=ftVis.shape[1]
    nt=ftVis.shape[0]
    if(visFlags is None):
        visFlags=np.empty(nu,dtype=bool)
        visFlags[:]=False
    for uNum in range(nu):
         binNum=int(np.round(uAIn[uNum]/(uMax/nCells)))
         if(binNum<nCells and not(np.any(np.isnan(ftVis[:,uNum]))) and not(np.any(np.isinf(ftVis[:,uNum]))) and not(visFlags[uNum])):
             counts[binNum]+=1
             dGrid[:,binNum]+=np.abs(ftVis[nt/2:,uNum])**2.
             dGrid[:,binNum]+=np.abs(np.flipud(ftVis[:nt/2,uNum]))**2.
    for bNum in range(nCells):
        if(counts[bNum]>0):
            dGrid[:,bNum]/=counts[bNum]
    uAxis=(np.arange(0,nCells)+.5)*uMax/nCells
    tAxis=tAIn[nt/2:]
    return uAxis,tAxis,dGrid


#method to convert delay-transformed and squared gridded visibilities to power spectrum
def delaySq2Ps(vs,f,sigma,band):
    '''
    ARGS:
    vs: arbitrary dimension grid of delay-transformed and squared visibilities
    f: center frequency of delay-transform measurement
    sigma: standard deviation of central lobe of primary beam (0.45*lambda/D)
    band: noise-equivalent bandwidth of delay-transform
    '''
    return littleh**3.*1e6*vs*X(f)**2.*Y(f)/(band*pi*sigma**2.)*i2t(f,1)**2.

#just use this function to do everything
def delayTransformAndGrid(visList,bList,f0,df,nCells,uMax,dAnt,flags):
    '''
    ARGS:
    visList: nfreqxnvis list of complex visibilities
    bList: nvis list of baseline lengths
    f0: center frequency (Hz)
    df: channel width (Hz)
    nCells: number of grid cells you want to have in the k_perp direction
    uMax: maximal u value that you want to grid. 
    dAnt: aperture diameters
    Returns:
    nCell kperp values (hMpc^-1), nfreq/2 kpara values (hMpc^-1), nfreq/2xnCell grid of cosmological power spectrum (mK^2 h^{-3}Mpc^3)
    '''
    dtV=delayTransformVisibilities(visList,df)
    uAxis,tAxis,ftV2Grid=gridDelay(bList*f0/c,fft.fftshift(fft.fftfreq(visList.shape[0],df)),dtV,nCells,uMax,visFlags=flags)
    kperpAxis=u2kperp(uAxis,f2z(f0))/littleh
    kparaAxis=eta2kpara(tAxis,f2z(f0))/littleh
    return kperpAxis,kparaAxis,delaySq2Ps(ftV2Grid,f0,0.45*(c/f0)/dAnt,len(tAxis)*df)
    
    
