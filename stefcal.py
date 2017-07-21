import numpy as np
import copy
import stefcal_utils as utils
DEBUG=False
#from numba import jit

'''
This method no longer does not do any flagging and can be vastly streamlined in light of this.
'''
def stefcal_scaler(data_matrix,model_matrix,weights_matrix,flag_matrix,
                    refant=0,n_phase_iter=5,n_cycles=1,min_bl_per_ant=2,
                    eps=1e-10,min_ant_times=1,perterb=0.):
    '''
    2-basic stefcal algorithm described in Salvini et al. 
    Args:
        data_matrix: Ntime x Nant x Nant complex numpy array of measured visibilities 
                    where the (i,j) entry is V_{ij}^{meas}= g_i^* g_j V_{ij}^{true}
        model_matrix: Ntime x Nant x Nant complex numpy array of model visibilities 
                     where the (i,j) entry is V_{ij}^{true}
        weights_matrix: Ntime x Nant x Nant numpy array of real weights. (i,j) entry
                       indicates how V_{ij}^{meas} will be weighted to solve for
                       g_i in each iterative step. 
        flag_matrix: Ntime x Nant x Nant numpy array of data flags. The (i,j) entry
                    indicates whether the (i,j) visibility is to be excluded from 
                    calibration. Symmetric matrix is enforced. 
        refant: integer representing the reference antenna 
        n_cycles: integer specifies the number of cycles of stefcal that will be run.
        min_bl_per_ant: minimum effective number of baselines per antenna. Computed
                     from visibility weights. Antennas with less than the minimum will 
                     be flagged. 
        eps: float specifying the maximal change in gain parameters 
             before termination.
        min_ant_times: integer minimum number of times an antenna must not be flagged
                       for it to remain unflagged. 
        trim_neff: specify whether to try trimming baselines responsible for small
                   min_bl_per_ant. 
   Returns: 
        ant_flags, NtimesxNant np array of boolean flags indicating antennas flagged over all times
        flag_matrix, Ntime x Nant x Nant np array of boolean flags indicating 
                   baselines flagged. 
        niter, n_cycles np array of integers. Denotes number of iterations 
               that were evaluated in each cycle. 
        gains, nAnt numpy array of complex numbers with estimated gains
    '''
    #start with type checks
    assert weights_matrix.dtype==np.float64
    assert data_matrix.dtype==np.complex128 or data_matrix.dtype==np.complex64
    assert model_matrix.dtype==np.complex128 or data_matrix.dtype==np.complex64
    assert flag_matrix.dtype==np.bool
    assert weights_matrix.shape==data_matrix.shape
    assert model_matrix.shape==weights_matrix.shape
    assert flag_matrix.shape==model_matrix.shape
    
    nAnt=data_matrix.shape[1]
    nTimes=data_matrix.shape[0]
    #start by flagging antennas with insufficient baselines
    weights_matrix=copy.copy(weights_matrix)
    #ant_flags=np.empty(data_matrix.shape[:2],dtype=bool);ant_flags[:]=False
    #for nt in range(nTimes):
    #   ant_flags[nt][utils.compute_neff(weights_matrix[nt])<min_bl_per_ant]=True
    #for nt in range(nTimes):
    #    print weights_matrix[nt][ant_flags[nt],:]

    weights_matrix[flag_matrix]=0.
    #if trim_neff:
    #    for nt in range(data_matrix.shape[0]):
    #        nf,ant_flags[nt],weights_matrix[nt],flag_matrix[nt],_=flag_neff(weights_matrix[nt],
    #                                                                        flag_matrix[nt],
    #                                                                        min_bl_per_ant)
            #print('nflagged='+str(nf))
            #print('ant_flags='+str(ant_flags[nt]))
            #print('new_weights='+str(weights_matrix[nt]))
            #print('any new flags?'+str(np.any(flag_matrix[nt])))
            #print('new_flags='+str(flag_matrix[nt]))
    #print ant_flags
    ant_flags_combined=np.empty(nAnt,dtype=bool);ant_flags_combined[:]=False
    #for m in range(nAnt):
    #    ant_flags_combined[m]=len(ant_flags[np.invert(ant_flags[:,m])])<min_ant_times

    #print('ant_flags_combined='+str(ant_flags_combined))
    
    #print('ant_flags_combined='+str(ant_flags_combined))
    antNumbers=np.arange(nAnt).astype(int)[np.invert(ant_flags_combined)]#only calibrate the antennas
    #that exceed the minimal number of baselines on the minimum required number of times
    #now run stefcal
    _eps=1.
    niter=np.zeros(n_cycles)
    gnew=0

    nAntG=len(antNumbers)
    gains=np.ones(nAnt,dtype=complex)
    gainsG=np.ones(nAntG,dtype=complex)
    gainsG_temp=np.ones_like(gainsG)
    data_matrixG=np.zeros((nTimes,nAntG,nAntG),dtype=complex)
    model_matrixG=np.zeros_like(data_matrixG)
    flag_matrixG=np.empty(data_matrixG.shape,dtype=bool)
    weights_matrixG=np.zeros((nTimes,nAntG,nAntG),dtype=float)
    #make sure that numerator and denominator are same length as gainsG
    gNumerator=np.zeros_like(gainsG)
    gDenominator=np.zeros_like(gainsG)
    gnew=np.zeros_like(gainsG)
    #initialize matrices
    for n in range(nAntG):
        for m in range(n):
            data_matrixG[:,m,n]=data_matrix[:,antNumbers[m],antNumbers[n]]
            data_matrixG[:,n,m]=data_matrix[:,antNumbers[n],antNumbers[m]]
            model_matrixG[:,m,n]=model_matrix[:,antNumbers[m],antNumbers[n]]
            model_matrixG[:,n,m]=model_matrix[:,antNumbers[n],antNumbers[m]]
            weights_matrixG[:,m,n]=weights_matrix[:,antNumbers[m],antNumbers[n]]
            weights_matrixG[:,n,m]=weights_matrix[:,antNumbers[n],antNumbers[m]]
            flag_matrixG[:,n,m]=flag_matrix[:,antNumbers[n],antNumbers[m]]
    #run stefcal cycles
    #if DEBUG:
    #    print('data_matrix='+str(data_matrixG[0,1,:]))
    #    print('model_matrix='+str(model_matrixG[0,1,:]))
    #    print('weights_matrix='+str(weights_matrixG[0,1,:]))
    for cycle in range(n_cycles):
        #if DEBUG:
        #    print('cycle='+str(cycle))
            
        perterbation=np.random.randn(len(gainsG))*perterb
        gainsG_temp[:]=1.+perterbation
        gainsG[:]=1.+perterbation

        while _eps>eps or niter[cycle]<=n_phase_iter:
            #print('_eps='+str(_eps))
            gNumerator[:]=0j
            gDenominator[:]=0j
            gainsG_temp[:]=gainsG[:]
            for nt in range(nTimes):
                #print('data='+str(data_matrixG[nt][0]))
                #print('model='+str(model_matrixG[nt][0]))
                
                zmat=np.diag(np.conj(gainsG)).dot(model_matrixG[nt])
                if DEBUG:
                    print('gainsG.shape='+str(gainsG.shape))
                    print('model_matrixG[nt].shape='+str(model_matrix[nt].shape))
                    print('weights_matrixG[nt].shape='+str(weights_matrixG[nt].shape))
                zmat_w=np.diag(np.conj(gainsG)).dot(model_matrixG[nt]*weights_matrixG[nt])
                gNumerator+=np.sum(np.conj(zmat_w)*data_matrixG[nt],axis=0)
                gDenominator+=np.sum(np.conj(zmat_w)*zmat,axis=0)
            #if DEBUG:
            gainsG=gNumerator/gDenominator
            if(niter[cycle]<n_phase_iter):
                gainsG=gainsG*np.abs(gainsG_temp)/np.abs(gainsG)
            gainsG=gainsG_temp/2.+gainsG/2. 
            gainsG*=np.conj(gainsG[refant])/np.abs(gainsG[refant])
            _eps_angle=np.abs(np.angle(gainsG)-np.angle(gainsG_temp)).max()
            _eps_abs=np.abs(gainsG-gainsG_temp).max()
            if(_eps_abs>_eps_angle):
                _eps=_eps_abs
            else:
                _eps=_eps_angle
            #if DEBUG:
            #    print('_eps='+str(_eps))
            #    if np.isnan(_eps):
            #        print('gDenominator='+str(gDenominator))
            #        print('len(gDenominator)='+str(len(gDenominator)))
            #        print('gNumerator='+str(gNumerator))
            #        print('len(gNumerator)='+str(len(gNumerator)))

            niter[cycle]+=1.
        for n in range(nAntG):
            gains[antNumbers[n]]*=gainsG[n]
            for m in range(n):
                for nt in range(nTimes):
                    data_matrixG[nt,m,n]/=(gainsG[n]*np.conj(gainsG[m]))
                    data_matrixG[nt,n,m]/=(gainsG[m]*np.conj(gainsG[n]))
    return niter,gains
                
    
        
    
    
