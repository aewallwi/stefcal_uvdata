import numpy as np
import copy
DEBUG=True
#from numba import jit



def compute_neff(weights_matrix):
    '''
    computes the number of effective baselines per antenna for a weighting matrix wMat
    Args:
        weights_matrix, Nant x Nant matrix of floats
    Returns:
        Nant numpy array of floats giving the effective number of baselines per antenna. 
    '''
    
    nEff=np.zeros(weights_matrix.shape[0])
    for antNum in range(weights_matrix.shape[0]):
        nEff[antNum]=np.abs(weights_matrix[antNum,:]).sum()/np.abs(weights_matrix[antNum,:]).max()
    return nEff

def flag_neff(weights_matrix,flagsMatrix=None,threshold=2):
    '''
    given a weighting matrix, flag baselines until the effective
    of antennas for each antenna is greater than some threshold.
    Args: 
        weights_matrix, Nant x Nant numpy array of floats giving weights of 
        each vis in stefcal. 
        flagsMatrix, Nant x Nant numpy array of booleans giving pre-existing
                     flags
    '''
    assert weights_matrix.dtype==float
    assert weights_matrix.shape[0]==weights_matrix.shape[1]
    
    weights_matrix_c=copy.deepcopy(weights_matrix)
    nAnt= len(weights_matrix)
    nVis=nAnt*(nAnt-1)/2
    
    if(flagsMatrix is None):
        flagsMatrix_c=np.empty(weights_matrix_c.shape,dtype=bool);flagsMatrix_c[:]=False
    else:
        assert flagsMatrix.dtype==bool
        assert flagsMatrix.shape==weights_matrix.shape
        flagsMatrix_c=copy.copy(flagsMatrix)
        weights_matrix_c[flagsMatrix_c]=0.
    #flagsMatrixL is a list of visibility flags
    
    flagsMatrixL=np.empty(nVis,dtype=bool);flagsMatrixL[:]=False
    
    nEff=computeNeff(weights_matrix_c)
    nFlag=0
    nvis=0

    while(np.any(nEff<threshold)):
        for i in range(weights_matrix_c.shape[0]):
            if(nEff[i]<threshold):
                maxInd=np.where(weights_matrix_c[i,:]==weights_matrix_c[i,:].max())[0][0]
                weights_matrix_c[i,maxInd]=0.
                weights_matrix_c[maxInd,i]=0.
                flagsMatrix_c[maxInd,i]=True
                flagsMatrix_c[i,maxInd]=True
                nFlag+=1        
        nEff=computeNeff(weights_matrix_c)
    for i in range(len(vFlags)):
        for j in range(i):
            flagsMatrixL[nvis]=flagsMatrix_c[i,j]
            nvis+=1
    antFlag=np.empty(weights_matrix.shape[0],dtype=bool)
    antFlag[:]=False
    for i in range(weights_matrix.shape[0]):
        antFlag[i]=np.all(weights_matrix_c[i,:]==0)
    return nFlag,antFlag,weights_matrix_c,flagsMatrix_c,flagsMatrixL
                

    
def stefcal_scaler(data_matrix,model_matrix,weights_matrix,flag_matrix,
                    refant=0,n_phase_iter=5,n_cycles=1,min_bl_per_ant=2,
                    eps=1e-10,min_ant_times=1,trim_neff=False,perterb=0.):
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
        ant_flags, Nant np array of boolean flags indicating antennas flagged
        flag_matrix, Ntime x Nant x Nant np array of boolean flags indicating 
                   baselines flagged. 
        niter, n_cycles np array of integers. Denotes number of iterations 
               that were evaluated in each cycle. 
        gains, nAnt numpy array of complex numbers with estimated gains
    '''
    #start with type checks
    assert weights_matrix.dtype==np.float64
    assert data_matrix.dtype==np.complex64
    assert model_matrix.dtype==np.complex64
    assert flag_matrix.dtype==np.bool
    assert weights_matrix.shape==data_matrix.shape
    assert model_matrix.shape==weights_matrix.shape
    assert flag_matrix.shape==model_matrix.shape
    
    nAnt=data_matrix.shape[1]
    nTimes=data_matrix.shape[0]
    #start by flagging antennas with insufficient baselines
    weights_matrix=copy.copy(weights_matrix)
    weights_matrix[flag_matrix]=0.
    ant_flags=np.empty(data_matrix.shape[:2],dtype=bool);ant_flags[:]=False


    for nt in range(nTimes):
        ant_flags[nt][compute_neff(weights_matrix)<min_bl_per_ant]=True
    
    if trim_neff:
        for nt in range(data_matrix.shape[0]):
            _,ant_flags[nt],weights_matrix[nt],flag_matrix[nt],_=flag_neff(weights_matrix[nt],
                                                                         weights_matrix[nt],
                                                                         min_bl_per_ant)
    ant_flags_combined=np.empty(nAnt,dtype=bool);ant_flags_combined[:]=False
    for m in range(nAnt):
        ant_flags_combined[m]=len(ant_flags[np.invert(ant_flags[:,m])])>min_ant_times
    antNumbers=np.arange(nAnt).astype(int)[np.invert(ant_flags_combined)]#the numbers of antennas not flagged
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
    weights_matrixG=np.zeros((nTimes,nAntG,nAntG),dtype=float)
    gNumerator=np.zeros_like(gains)
    gDenominator=np.zeros_like(gains)
    gnew=np.zeros_like(gains)
    #initialize matrices
    for n in range(nAntG):
        for m in range(n):
            data_matrixG[:,m,n]=data_matrix[:,antNumbers[m],antNumbers[n]]
            data_matrixG[:,n,m]=data_matrix[:,antNumbers[n],antNumbers[m]]
            model_matrixG[:,m,n]=model_matrix[:,antNumbers[m],antNumbers[n]]
            model_matrixG[:,n,m]=model_matrix[:,antNumbers[n],antNumbers[m]]
            weights_matrixG[:,m,n]=weights_matrix[:,antNumbers[m],antNumbers[n]]
            weights_matrixG[:,n,m]=weights_matrix[:,antNumbers[n],antNumbers[m]]
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
            gNumerator[:]=0.
            gDenominator[:]=0.
            gainsG_temp[:]=gainsG[:]
            for nt in range(nTimes):
                zmat=np.diag(np.conj(gainsG)).dot(model_matrixG[nt])
                zmat_w=np.diag(np.conj(gainsG)).dot(model_matrixG[nt]*weights_matrixG[nt])
                gNumerator+=np.sum(np.conj(zmat_w)*data_matrixG[nt],axis=0)
                gDenominator+=np.sum(np.conj(zmat_w)*zmat,axis=0)
            #if DEBUG:
            gainsG=gNumerator/gDenominator
            if(niter[cycle]<n_phase_iter):
                gainsG=gainsG*np.abs(gainsG_temp)/np.abs(gainsG)
            gainsG=gainsG_temp/2.+gainsG/2. 
            gainsG*=np.conj(gainsG[refant])/np.abs(gainsG[refant])
            _eps=np.abs(gainsG-gainsG_temp).max()
            if DEBUG:
                print('_eps='+str(_eps))
                if np.isnan(_eps):
                    print('gDenominator='+str(gDenominator))
                    print('len(gDenominator)='+str(len(gDenominator)))
                    print('gNumerator='+str(gNumerator))
                    print('len(gNumerator)='+str(len(gNumerator)))

            niter[cycle]+=1.

        for nt in range(nTimes):
            for m in range(nAntG):
                gains[antNumbers[m]]*=gainsG[m]
                for n in range(m):
                    data_matrixG[nt,m,n]/=(gainsG[n]*np.conj(gainsG[m]))
                    data_matrixG[nt,n,m]/=(gainsG[m]*np.conj(gainsG[n]))
    return ant_flags,ant_flags_combined,niter,gains
                
    
        
    
    
