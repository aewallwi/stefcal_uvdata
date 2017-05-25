#create a single zenith point source.
#remove directories
#exec('rm -rf point_source_sim')
#exec('rm -rf point.cl')
import numpy as np

phase_center="J2000 10h00m00.08s 34d04m43.0s"
df=0.01
nf=20
freqs = 1.4+np.arange(-nf/2,nf/2)*df
cl.done()

#cl.rename("point_mfs.cl")
#cl.done()

default('simobserve')
integration='2s'
project="point_source_sim_mfs"
complist="../data/point.cl"
compcenter=";".join(["%.2fGHz"%freqs[m] for m in range(nf)])
compwidth=";".join(["%.2fGHz"%df for m in range(nf)])
direction=phase_center
obsmode="int"
antennalist="vla.cfg"
totaltime="60s"
thermalnoise="tsys-atm"
antennalist='vla.a.cfg'
simobserve()
tb.open('point_source_sim_nfs/point_source_sim_nfs.vla.a.ms/SPECTRAL_WINDOW')
print tb.getcol('CHAN_FREQ')


