import numpy as np

from numba import jit
            
    
@jit(nopython=True)
def calibrate_2_basic_scalar(dataMatrix,modelMatrix,gains,gbuffer,eps,weights,refant=0,n_phase_iter=5):
    '''
    2-basic stefcal routine from Salvini et al.
    Args:
    data: array of data, complex N(N-1)/2 visibilities 
    model: array of complex N(N-1)/2 model visibilities
    gains: initial guesses for antenna gains, vector
    weights: weights for each visibility vector
    epsilon: convergence condition
    gbuffer: buffer to store last gain values
    '''
    _eps=1.
    niter=0
    gnew=0.
    niter=0
#    while(not(_eps<eps and niter<n_phase_iter)):
    while(_eps>eps or niter < n_phase_iter):
        #zmat=np.einsum('ij,jk',np.diag(gains.conj()),modelMatrix)
        zmat=np.dot(np.diag(np.conj(gains)),modelMatrix)
        for j in range(len(gains)):
            gbuffer[j]=gains[j]
            zmat_weighted=np.dot(np.diag(np.conj(gains)),modelMatrix*weights)
            #gnew=np.dot(np.conj(zmat[:,j]),dataMatrix[:,j])/np.dot(np.conj(zmat[:,j]),zmat[:,j])
            gnew=np.dot(np.conj(zmat_weighted[:,j]),dataMatrix[:,j])/np.dot(np.conj(zmat_weighted[:,j]),zmat[:,j])
            #gnew=np.dot(np.conj(zmat_weighted[:,j]),dataMatrix[:,j]*weights[:,j])/np.dot(np.conj(zmat_weighted[:,j]),zmat_weighted[:,j])
            #update gains
            #print 'gnew='+str(gnew)
            #print 'gbuffer[j]='+str(gbuffer[j])
            #print gbuffer.shape
            if(niter<n_phase_iter):
                gnew=gnew*np.abs(gbuffer[j])/np.abs(gnew)
            gains[j]=gbuffer[j]/2.+gnew/2.
            #zmat=np.einsum('ij,jk',np.diag(gains.conj()),modelMatrix)
            #zmat=np.dot(np.diag(np.conj(gains)),modelMatrix)
        niter+=1
        #print _eps
        for j in range(len(gains)):
            gains[j]*=np.conj(gains[refant])/np.abs(gains[refant])
        
        _eps=np.abs(1.-gains/gbuffer).max()
    #print 'converged in %d iterations'%(niter)
    return gains
            
    
