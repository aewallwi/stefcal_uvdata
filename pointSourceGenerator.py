import numpy as np
PI=np.pi
KAPPA=3200#6c dn/ds params hard coded
BETA=2.51
GAMMA=1.75
S0=.88
#smin is the minimum flux
#use 6C source counts
#draw alpha from distribution with mean 0.5 and std 0.5
def power_law(alpha,xmin,xmax,x):
    return (1-alpha)/(xmax**(-alpha+1)-xmin**(-alpha+1))*x**-alpha
#draw power law distribution between xmin and xmax with P(x) ~ x^-alpha
#with size samples
def draw_power_law(xmin,xmax,alpha,size=1):
    y=np.random.rand(size)
    return ((xmax**(-alpha+1)-xmin**(-alpha+1))*y+xmin**(-alpha+1))**(1/(-alpha+1))
def drawSources(area,smin,smax):

    if(smin>S0):
        coeff=KAPPA*S0**BETA*(smin**(1.-BETA)-smax**(1.-BETA))/(BETA-1.)
        n1=int(np.random.poisson(int(area*coeff)))
        n2=0
        fluxes=draw_power_law(smin,smax,BETA,size=n1)
    if(smin<S0):
        coeff1=KAPPA*S0**GAMMA*(S0**(1-GAMMA)-smin**(1-GAMMA))/(1-GAMMA)
        coeff2=KAPPA*S0/(BETA-1.)
        n1=int(np.random.poisson(int(area*coeff1)))
        n2=int(np.random.poisson(int(area*coeff2)))
        fluxes=np.hstack([draw_power_law(smin,S0,GAMMA,size=n1),draw_power_law(S0,smax,BETA,size=n2)])
    spectral_inds=np.random.normal(loc=-0.5,scale=0.5,size=n1+n2)
    return fluxes,spectral_inds
#lAngleList=np.array([0])
def drawRandomSources(smin,smax=1e3,area=2*PI):    
    flux,spIndices=drawSources(area,smin,smax)#draw sources greater than 1 Jy from 2 PI sr
    nsrcs=len(flux)
    print 'number of sources=%d'%(nsrcs)
    phiList=np.random.rand(int(nsrcs))*PI*2.
    thetaList=np.arccos(np.random.rand(int(nsrcs))-1)
    lAngleList=np.sin(thetaList)*np.cos(phiList)
    mAngleList=np.sin(thetaList)*np.sin(phiList)
    select=flux<=smax
    return np.array([lAngleList[select],mAngleList[select],flux[select],spIndices[select]]).T

#    computeVisList(lAngles,mAngles,uList,vList,fluxList,output)
