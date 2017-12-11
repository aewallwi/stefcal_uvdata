
import numpy as np
import copy
from pyuvdata import UVCal
#************************************************************
#Generate a corrected uvdata set. 
#************************************************************
def correct_vis(uvdata,uvcal,applyGains=False):
    """
    Generate a corrected uvdata set 
    Args: uvdata, data set to apply gains to
          uvcal, calibration solution to apply
    Returns:
          uvcal object identical to uvdata except
          data_array has had gains divided out
          and calibration flags have been applied to 
          data
    """
    corrected_vis=copy.deepcopy(uvdata)
    #merge flags
    #print('data array before application='+str(corrected_vis.data_array))
    #print('gains='+str(uvcal.gain_array))
    #print('data type='+str(corrected_vis.data_array.dtype))
    #print('ant max='+str(uvdata.antenna_numbers.max()))
    #print('ant nums='+str(uvdata.antenna_numbers))
    unique_ants=np.unique(np.vstack([uvdata.ant_1_array,uvdata.ant_2_array]))
    ant_dict={}
    ant_dict_inv={}
    for antind,antnum in enumerate(unique_ants):
        ant_dict[antind]=antnum
        ant_dict_inv[antnum]=antind
    for antind in ant_dict.keys():
        ant=ant_dict[antind]
        for tnum in range(corrected_vis.Ntimes):
            selection=np.logical_and(corrected_vis.ant_2_array==ant,
                                     corrected_vis.time_array==uvcal.time_array[tnum])
            selection_conj=np.logical_and(corrected_vis.ant_1_array==ant,
                                          corrected_vis.time_array==uvcal.time_array[tnum])
            for pol in range(corrected_vis.Npols):
                for chan in range(corrected_vis.Nfreqs):
                    if applyGains:
                        corrected_vis.data_array[selection_conj,0,chan,pol]*=np.conj(uvcal.gain_array[antind,0,chan,tnum,pol])
                        corrected_vis.data_array[selection,0,chan,pol]*=uvcal.gain_array[antind,0,chan,tnum,pol]
                    else:
                        corrected_vis.data_array[selection_conj,0,chan,pol]/=np.conj(uvcal.gain_array[antind,0,chan,tnum,pol])
                        corrected_vis.data_array[selection,0,chan,pol]/=uvcal.gain_array[antind,0,chan,tnum,pol]
                    added_flags=[uvcal.flag_array[antind,0,chan,tnum,pol] for m in range(len(selection[selection]))]
                    corrected_vis.flag_array[selection,0,chan,pol]=np.logical_or(corrected_vis.flag_array[selection,0,chan,pol],added_flags)
    #print('data array after application='+str(corrected_vis.data_array))
    return corrected_vis

#************************************************************
#generate weights from baseline lengths
#************************************************************
def generate_gaussian_weights(sigma_w,uvmodel,modelweights=False,regularizer=1e-6,trim_neff=True):
    """
    Generate gaussian weights for visibilities with 
    w=exp(-|uvw|^2/(2 sigma_w^2)
    Args:
    sigma_w, standard deviation weighting parameter (meters) (float)
    uvmodel, uvdata object containing data to generate weights from
    modelweights, if True, multiply each weight by amplitude of visibility squared (optimal weighting for diagonal thermal noise). 
    """
    bl_lengths=np.linalg.norm(uvmodel.uvw_array,axis=1)
    avg_amp=np.exp(-bl_lengths**2./(2.*sigma_w**2.))
    avg_amp[avg_amp<regularizer]=regularizer
    bl_flags=np.empty(uvmodel.data_array.shape,dtype=bool)
    bl_flags[:]=False
    u_times=np.unique(uvmodel.time_array)
    ant_flags=np.empty((uvmodel.Ntimes,uvmodel.Nants_data),dtype=bool)
    ant_flags[:]=False
    avg_flags=np.empty(uvmodel.Nblts,dtype=bool)
    avg_flags[:]=False
    if trim_neff:
        for nt,time in enumerate(u_times):
            selection=uvmodel.time_array==time
            n_flag,ant_flags[nt],avg_amp[selection],avg_flags[selection]=flag_neff(avg_amp[selection],
                                                                                   bl_flags[selection,0,0,0],
                                                                                   mode='blt_list',
                                                                                   ant1List=uvmodel.ant_1_array[selection],
                                                                                   ant2List=uvmodel.ant_2_array[selection])
    #make sure to flag all channels and pols for a flagged antenna.
    for pol in range(uvmodel.Npols):
        for chan in range(uvmodel.Nfreqs):
            for spw in range(uvmodel.Nspws):
                bl_flags[:,spw,chan,pol]=avg_flags
                


        
            
    if modelweights:
        wbase=np.abs(uvmodel.data_array)**2./(np.abs(uvmodel.data_array).max())**2.
    else:
        wbase=np.ones(uvmodel.data_array.shape)
    for pol in range(uvmodel.Npols):
        for spw in range(uvmodel.Nspws):
            for chan in range(uvmodel.Nfreqs):
                wbase[:,spw,chan,pol]*=avg_amp
    return n_flag,ant_flags,bl_flags,wbase



def flag_ants(weight_array,threshhold=2.,mode='matrix',ant1List=None,ant2List=None):
    '''
    Flags antennas with neff smaller than a threshold. blt_list mode assumes a single time.
    Args: 
         weight_array: either an NantxNant matrix of antenna cross correlations with flagged antenna weights set to zero or a N_blt array of cross correlations
         thresshold: the minimum effective number of baselines per antenna
         mode: specify whether the weights are in NantxNant matrix form or Nant*(Nant-1)/2 list of baselines. 
    '''
    assert weight_array.dtype==np.float64
    assert mode in ['matrix','blt_list']
    ant_flags=np.empty(matrix.shape[0],dtype=bool);ant_flags[:]=False
    ant_flags=compute_neff(weight_array,mode=mode,
                           ant1List=ant1List,ant2List=ant2List)<threshhold
    return ant_flags


def compute_neff(weights_array,mode='matrix',ant1List=None,ant2List=None):
    '''
    computes the number of effective baselines per antenna for a weighting matrix wMat
    Args:
        weights_array, Nant x Nant matrix of floats
    Returns:
        Nant numpy array of floats giving the effective number of baselines per antenna. 
    '''
    assert mode in ['matrix','blt_list']
    if mode == 'matrix':
        nEff=np.zeros(weights_array.shape[0])
        for antnum in range(weights_array.shape[0]):
            if np.any(weights_array[antnum,:]!=0.):
                nEff[antnum]=np.abs(weights_array[antnum,:]).sum()/np.abs(weights_array[antnum,:]).max()
            else:
                nEff[antnum]=0.
        return nEff
    elif mode == 'blt_list':
        assert not(ant1List is None) and not(ant2List is None)
        assert len(ant1List)==len(ant2List) 
        u_ant=np.unique(np.vstack([ant1List,ant2List]))
        nEff=np.zeros(len(u_ant))
        n_ant=len(u_ant)
        ant_index={}
        for index,antnum in enumerate(u_ant):
            ant_index[antnum]=index
        for antnum in u_ant:
            selection=np.logical_xor(ant1List==antnum,ant2List==antnum)
            maxval=np.max(np.abs(weights_array[selection]))
            nEff[ant_index[antnum]]=np.sum(np.abs(weights_array[selection]))/maxval
        return nEff
            

def flag_neff(weights_array,flags_array=None,threshold=2.,mode='matrix',ant1List=None,ant2List=None):
    '''
    given a weighting matrix, flag baselines until the effective
    of antennas for each antenna is greater than some threshold.
    Args: 
        weights_array, Nant x Nant numpy array of floats giving weights of 
        each vis in stefcal or Nvis vector of floats (blt_list mode). 
        flags_array, Nant x Nant numpy array of booleans giving pre-existing
                     flags or Nvis list of pre-existing flags (blt_list mode)
        mode, string specififying whether weights/flags are in blt_list format or
              matrix format. 
        antList1, Nvis vector of integers for antenna 1 of array (bltList mode only)
        antList2, Nvis vector of integers in antenna 2 of array (bltList mode only)
    '''
    weights_array_c=copy.deepcopy(weights_array)
    if(flags_array is None):
        flags_array_c=np.empty(weights_array_c.shape,dtype=bool);flags_array_c[:]=False
    else:
        assert flags_array.dtype==bool
        #print('flags_array.shape='+str(flags_array.shape))
        #print('weights_array.shape='+str(weights_array.shape))
        assert flags_array.shape==weights_array.shape
        flags_array_c=copy.copy(flags_array)
        weights_array_c[flags_array_c]=0.
    assert mode in ['matrix', 'blt_list']
    
    if mode=='matrix':
        assert weights_array.dtype==float
        assert weights_array.shape[0]==weights_array.shape[1]
    
        nAnt= len(weights_array)
        nVis=nAnt*(nAnt-1)/2
    
        #flags_array_l is a list of visibility flags
        flags_array_l=np.empty(nVis,dtype=bool);flags_array_l[:]=False
        
        nEff=compute_neff(weights_array_c,mode)
        nFlag=0
        nvis=0

        while(np.any(nEff<threshold)):
            for i in range(weights_array_c.shape[0]):
                if(nEff[i]<threshold):
                    maxInd=np.where(weights_array_c[i,:]==weights_array_c[i,:].max())[0][0]
                    weights_array_c[i,maxInd]=0.
                    weights_array_c[maxInd,i]=0.
                    flags_array_c[maxInd,i]=True
                    flags_array_c[i,maxInd]=True
                    nFlag+=1        
                    nEff=compute_neff(weights_array_c,mode)
        for i in range(nAnt):
            for j in range(i):
                flags_array_l[nvis]=flags_array_c[i,j]
                nvis+=1
        antFlag=np.empty(weights_array.shape[0],dtype=bool)
        antFlag[:]=False
        for i in range(weights_array.shape[0]):
            antFlag[i]=np.all(weights_array_c[i,:]==0)
        return nFlag,antFlag,weights_array_c,flags_array_c,flags_array_l
    elif mode=='blt_list':
        assert not(ant1List is None) and not(ant2List is None)
        assert len(ant1List)==len(ant2List); assert len(ant1List)==len(weights_array)
        assert len(weights_array)==len(flags_array_c)
        weights_array_c=copy.deepcopy(weights_array)
        u_ants=np.unique(np.vstack([ant1List,ant2List]))
        ant_index={}
        for ind,antnum in enumerate(u_ants):
            ant_index[antnum]=ind

        nAnt=len(u_ants)
        nEff=compute_neff(weights_array_c,mode,ant1List,ant2List)
        nFlag=0

        while(np.any(nEff<threshold)):
            for antnum in u_ants:
                if nEff[ant_index[antnum]]<threshold:
                    selection=np.logical_xor(ant1List==antnum,
                                             ant2List==antnum)
                    temp=weights_array_c[selection]
                    temp_f=flags_array_c[selection]
                    maxind=np.where(temp==temp.max())[0][0]
                    temp[maxind]=0.
                    temp_f[maxind]=True
                    weights_array_c[selection]=temp
                    flags_array_c[selection]=temp_f
                    nFlag+=1
                    nEff=compute_neff(weights_array_c,mode,ant1List,ant2List)
        #print('neff.min='+str(nEff.min()))
        #print('nflag='+str(nFlag))
            #print('nEff='+str(nEff))
        antFlag=np.empty(nAnt,dtype=bool)
        antFlag[:]=False
        for antnum in u_ants:
            selection=np.logical_xor(ant1List==antnum,
                                     ant2List==antnum)
            antFlag[ant_index[antnum]]=np.all(weights_array[selection]==0.)
        return nFlag,antFlag,weights_array_c,flags_array_c
        
    
