import numpy as np
import pointSourceGenerator as psg
PI=np.pi
SMIN=0.1
SCUT=1.0

psList=psg.drawRandomSources(smin=SMIN,area=4*PI)

deg_ra=np.degrees(psList[1])
deg_dec=np.degrees(psList[0]-PI/2.)
nps=len(deg_ra)
output=np.vstack([deg_ra,deg_dec,psList[2],psList[3],np.zeros(nps),np.zeros(nps),np.zeros(nps)]).T
np.savetxt('point_source_list_SMIN_%.2f.txt'%SMIN,output)

output=output[output[:,2]>SCUT,:]
np.savetxt('point_source_list_SMIN_%.2f.txt'%SCUT,output)
