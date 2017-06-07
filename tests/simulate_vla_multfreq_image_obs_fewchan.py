#************************************************************
#simulation code to create a single point source at zenith
#with spectral index of 0. 
import numpy as np

phase_center="J2000 10h00m00.08s 34d04m43.0s"
cl.done()


default('simobserve')
integration='2s'
direction=phase_center
obsmode="int"
antennalist="vla.cfg"
totaltime="2s"
thermalnoise="tsys-atm"
antennalist='vla.a.cfg'

cl.done()
cl.addcomponent(dir=phase_center,
                flux=100.0,
                fluxunit="Jy",
                freq="100MHz",
                shape="point",
                index=0.0)
cl.rename("point_100.cl")
cl.done()

complist="point_100.cl"
project="point_source_sim_mfs_20chan"
incenter="150MHz"
inwdith="100kHz"
skymodel="dummy_multifreq_20chan.fits"
inbright="0Jy/pixel"
indirection=phase_center
incell="0.1 arcsec"
simobserve()
    



