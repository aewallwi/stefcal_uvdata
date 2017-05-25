#create a simulation of many point sources for MWA
#drawn from power law distribution
import pointSourceGenerator as psg
import numpy as np
PI=np.pi

def sc_2_j2000(theta,phi):
    """
    convert from spherical radians to dir string
    """
    deg_ra=np.degrees(theta)
    hr_ra=np.floor(deg_ra/15.)
    min_ra=np.floor((deg_ra/15.-hr_ra)*60.)
    sec_ra=((deg_ra/15.-hr_ra)*60.-min_ra)*60.

    deg_dec_full=np.degrees(theta-PI/2.)
    deg_sign=np.sign(deg_dec_full)
    deg_dec=np.floor(np.abs(deg_dec_full))
    min_dec=np.floor((np.abs(deg_dec_full)-deg_dec)*60.)
    sec_dec=(60.*(np.abs(deg_dec_full)-deg_dec)-min_dec)*60.
    deg_dec*=deg_sign
    if deg_sign*1>0:
        return "J2000 %02dh%02dm%05.2fs %02dd%02dm%05.2fs"%(hr_ra,min_ra,sec_ra,
                                                            deg_dec,min_dec,sec_dec)
    else:
        return "J2000 %02dh%02dm%05.2fs %03dd%02dm%05.2fs"%(hr_ra,min_ra,sec_ra,
                                                            deg_dec,min_dec,sec_dec)
phase_center="J2000 10h00m00.8s 34d04m43.0s"
#draw random sources from a sphere. 
SMIN=.1#minimum flux (Jy)
SCUT=1.0#minimum calibration catalog flux (Jy)
psList=psg.drawRandomSources(smin=SMIN,area=4*PI)
#add components
cl.done()
for snum in range(psList.shape[0]):
    src=psList[snum]
    dirstr=sc_2_j2000(src[0],src[1])
    print("(theta=%06.4f,phi=%06.4f) = %s"%(src[0],src[1],dirstr))
    cl.addcomponent(dir=dirstr,
                    fluxunit="Jy",
                    freq="150MHz",
                    shape="point",
                    index=src[-1],
                    spectrumtype="spectral index",
                    flux=src[-2])
cl.rename("point_source_cat_%e2.cl"%SMIN)
cl.done()
psList=psList[psList[:,-2]>SCUT,:]
for snum in range(psList.shape[0]):
    src=psList[snum]
    dirstr=sc_2_j2000(src[0],src[1])
    print("(theta=%06.4f,phi=%06.4f) = %s"%(src[0],src[1],dirstr))
    cl.addcomponent(dir=dirstr,
                    fluxunit="Jy",
                    freq="150MHz",
                    shape="point",
                    index=src[-1],
                    spectrumtype="spectral index",
                    flux=src[-2])
cl.rename("point_source_cat_%e2.cl"%SCUT)
cl.done()
default=("simobserve")
project="../data/mwa_sim_%e2Jy_data.py"%SMIN
complist="point_source_cat_%e2.cl"%SMIN
compwidth="100kHz," for m in range(200)




