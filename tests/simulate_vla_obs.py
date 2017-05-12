#create a single zenith point source.
#remove directories
#exec('rm -rf point_source_sim')
#exec('rm -rf point.cl')


phase_center="J2000 10h00m00.08s 34d04m43.0s"

cl.done()
cl.addcomponent(dir=phase_center,
                flux=1.0,
                fluxunit="Jy",
                freq="1.4GHz",
                shape="point")
cl.rename("point.cl")
cl.done()

default('simobserve')
integration='2s'
project="point_source_sim"
complist="point.cl"
compwidth="400MHz"
direction=phase_center
obsmode="int"
antennalist="vla.cfg"
totaltime="60s"
thermalnoise="tsys-atm"
antennalist='vla.a.cfg'
simobserve()

