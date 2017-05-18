DEBUG=True
from pyuvdata import UVData
import numpy as np
from pyuvdata import UVCal
import pyuvdata.parameter as uvp
from cal_flag_weights import CalFlagWeights
import stefcal
import pickle
from stefcal_meta import StefcalMeta
import uuid


#************************************************************
#Need to store state IDs for different times
#************************************************************

class StefcalUVData():
    """
    Defines a class for performing stefcal on uvdata sets.
    Attributes:
    model_vis: uvdata of model
    measured_vis: uvdata of measurements
    cal_flag_weights: CalFlagWeights object storing flags and weights
    cal_solution: uvcal object
    meta_params: StefCalMeta object containing information on calibration
                 run.
    _noise_tavg: uv_parameter for noise (uncalib) on each frequency channel. 
    _noise_favg: uv_parameter for noise (uncalib) at each time
    id: unique uuid corresponding to calibration run. 
    """
    
    def __init__(self,refant=0,n_phase_iter=5,
                 n_cycles=1,min_bl_per_ant=2,eps=1e-10,
                 min_ant_times=1,trim_neff=False,spw=0,
                 t_avg=1):
        self.model_vis=UVData()
        self.measured_vis=UVData()
        self.uvcal=UVCal()
        self.meta_params=StefcalMeta()

        
        self.meta_params.refant=refant
        self.meta_params.n_phase_iter=n_phase_iter
        self.meta_params.n_cycles=n_cycles
        self.meta_params.id=str(uuid.uuid4())
        self.meta_params.n_phase_iter=n_phase_iter
        self.meta_params.min_bl_per_ant=min_bl_per_ant
        self.meta_params.eps=eps
        self.meta_params.spw=spw
        self.meta_params.t_avg=t_avg
        self.meta_params.trim_neff=trim_neff
        self.meta_params.min_ant_times=min_ant_times

        self.cal_flag_weights=CalFlagWeights(self.meta_params.id)

        
    def _compare_model_data(self,model,data,compare_flags=False,compare_weights=False,compare_phase=False):
        """
        compare properties of uvdata sets and calweights. 
        Args:
             model, CalFlagWeight or UVData object to be compared to 
             data, CalFlagWeight. 
             compare_flags: directly compare the flags
             compare_weights: directly compare nsample_array 
             compare_phase: if two data sets, make sure they are phased in the same way
        """
        checkpass=True
        if model.Nfreqs!=data.Nfreqs:
            checkpass=False
            print("Nfreqs not the same")
        if model.Npols!=data.Npols:
            checkpass=False
            print("Njones not the same")
        if model.Ntimes!=data.Ntimes:
            checkpass=False
            print("Ntimes not the same")
        if model.Nspws!=data.Nspws:
            checkpass=False
            print("Nspws not the same")
        if np.any(model.freq_array!=data.freq_array):
            checkpass=False
            print("freq_range not the same")            
        if np.any(model.time_array!=data.time_array):
            checkpass=False
            print("time_range not the same")
        if model.telescope_name!=data.telescope_name:
            checkpass=False
            print("telescope_name not the same")
        if model.Nants_data!=data.Nants_data:
            checkpass=False
            print("Nants_data not the same")
        if model.Nants_telescope!=data.Nants_telescope:
            checkpass=False
            print("Nants_telescope not the same")
        if model.antenna_names!=data.antenna_names:
            checkpass=False
            print("antenna_names not the same")
        if np.any(model.antenna_numbers!=data.antenna_numbers):
            checkpass=False
            print("antenna_numbers not the same")
        if model.channel_width!=data.channel_width:
            checkpass=False
            print("channel_width not the same")
        if np.any(model.polarization_array!=data.polarization_array):
            checkpass=False
            print("polarization_array not the same")
        if np.any(model.time_array!=data.time_array):
            checkpass=False
            print("time_array not the same")
        if model.integration_time!=data.integration_time:
            checkpass=False
            print("integration times not the same")
        #if model.x_orientation!=data.x_orientation:
        #    checkpass=False
        #    print("integration times not the same")
        if compare_flags:
            if np.any(model.flag_array!=data.flag_array):
                checkpass=False
                print("flag_array not the same")
        if compare_weights:
            if np.any(model.nsample_array!=data.nsample_array):
                checkpass=False
                print("nsample_array not the same")
                
        if (model.phase_type!=data.phase_type or
            model.phase_center_ra!=data.phase_center_ra or
            model.phase_center_dec!=data.phase_center_dec):
            checkpass = False
            print("phase centers not the same")
        if np.any(model.antenna_positions!=data.antenna_positions):
            checkpass=False
            print("antenna positions not the same") 
                
            
        return checkpass
                            
    def _check_consistency(self):
        """
        check consistency between data and model uvdata objects
        also check whether the calweights object has the correct
        dimension
        """
        if not(self._compare_model_data(self.model_vis,self.measured_vis)):
               raise ValueError("model_vis not consistent with measured_vis")
        
    #if(self.compare_properties(self.model_vis,self.flag_weights)):
            #raise ValueError("flag weights non consistent with visibilities")
            
    def _load_vis_ms(self,msname):
        '''
        read in ms visibilities
        args: msName, path to .ms directory
        '''
        #need to throw an error here if the model column does not exist
        self.model_vis.read_ms(msname,columnname='MODEL')
        self.measured_vis.read_ms(msname,columnname='DATA')
        
    def _load_fhd(self,dataname):
        '''
        read in fhd visibility files
        args: dataname, name of the fhd file
        '''
        self.model_vis.read_fhd(fhdname,use_model=True)
        self.measured_vis.read_fhd(fhdname)
        
    def _load_uvfits(self,dataname, modelname):
        '''
        read in uvfits files 
        args: dataname, name of uvfits data file
              modelname, name of uvfits model file
        '''
        self.model_vis.read_uvfits(modelname)
        self.measured_vis.read_uvfits(dataname)
        #import_uvfits appears to use 0-indexed uvfits
        self.model_vis.ant_1_array+=1
        self.model_vis.ant_2_array+=1
        self.measured_vis.ant_1_array+=1
        self.measured_vis.ant_2_array+=1
        self.model_vis.antenna_numbers+=1
        self.measured_vis.antenna_numbers+=1
        if DEBUG:
            print('antenna_numbers='+str(self.model_vis.antenna_numbers))
        
    def _load_miriad(self,dataname,modelname):
        '''
        read in miriad files
        args: dataname, name of miriad data file
              modelname, name of miriad model file
        '''
        self.model_vis.read_mirad(modelname)
        self.measured_vis.read_miraid(dataname)
        
    def _compute_noise(self,mode='dTIMEmTIMEBL'):#,minmax_percentile=.05):
        """
        compute thermal noise levels for each baseline and frequency
        using adjacent time differencing and adjacent channel differencing.
        Args:
             mode='FREQ', or 'TIME' determines over which axis to diff and compute
             averages.
        """
        #assert mode in ['TIME','FREQ','dTIMEmBL','dFREQmBL']
        assert mode == 'dTIMEmTIMEBL'
        """
        if mode=='FREQ':
            #take diff in freq and average over freq
            self.meta_params.noise_favg=np.zeros((self.measured_vis.Nblts,self.measured_vis.Npols))
            for pol in range(self.measured_vis.Npols):
                for blt in range(self.measured_vis.Nblts):
                    data_select=self.measured_vis.data_array[blt,:,self.spw,pol].squeeze()
                    flag_select=self.cal_flag_weights.flag_array[blt,:,self.spw,pol].squeeze()
                    #only take diffs for adjacent channels,
                    mlist=np.array([])
                    for chan in range(1,self.measured_vis.Nfreq,2):
                        if not(flag_select[chan]) and not(flag_select[chan-1]):
                            mlist=np.append(mlist,np.abs(data_select[chan]-data_select[chan-1])**2.)
                   # minmax=np.percentile(mlist,[minmax_percentile,1-minmax_percentile])
                   #self.noise_favg[blt,pol]=np.mean(mlist[np.logical_and(mlist>=minmax[0],mlist<=minmax[1])])/2.
                   self.meta_params.noise_favg[blt,pol]=np.median(mlist)/(2.*np.log(2.))
        elif mode=='TIME':
            #take diff in time and average over time
            self.meta_params.noise_tavg=np.zeros((self.measured_vis.NBls,self.measured_vis.Nfreqs,self.Njones))
            for pol in range(self.measured_vis.Npols):
                for blnum in range(self.measured_vis.Nbls):
                    bl_selection=self.baseline_array==self.baseline_array[blnum]
                    ndiffs=0
                    for chan in range(self.measured_vis.Nfreqs):
                        data_select=self.measured_vis.data_array[bl_selection,self.spw,chan,pol]
                        flag_select=self.cal_flag_weights.flag_array[bl_selection,self.spw,chan,pol]
                        #only take diffs for adjacent channels,
                        mlist=np.array([])
                        for tnum in range(1,self.measured_vis.Ntimes,2):
                            if not(flag_select[tnum]) and not(flag_select[tnum-1]):
                                mlist=np.append(mlist, np.abs(data_select[tnum]-data_select[tnum-1])**2.)
                        #minmax=np.percentile(mlist,[minmax_percentile,1-minmax_percentile])
                        #self.noise_tavg[blnum,chan,pol]=np.mean(mlist[np.logical_and(mlist>=minmax[0],mlist<=minmax[1])])/2.
                        self.meta_params.noise_tavg[blnum,chan,pol]=np.median(mlist)/(2.*np.log(2.))
        """
        if mode=='dTIMEmTIMEBL':
            #take diff in time and average in time,baselines,and pols; good for measurement sets
            #with a small number of time samples.
            self.meta_params.noise_tblavg=np.zeros((self.measured_vis.Nfreqs))
            for chan in range(self.measured_vis.Nfreqs):
                mlist=([])
                for pol in range(self.measured_vis.Npols):
                    for blnum in range(self.measured_vis.Nbls):
                        bl_selection=self.measured_vis.baseline_array==self.measured_vis.baseline_array[blnum]
                        data_select=self.measured_vis.data_array[bl_selection,self.meta_params.spw,chan,pol]
                        flag_select=self.measured_vis.flag_array[bl_selection,self.meta_params.spw,chan,pol]
                        for tnum in range(1,self.measured_vis.Ntimes,2):
                            if not(flag_select[tnum]) and not(flag_select[tnum-1]):
                                mlist=np.append(mlist,np.abs(data_select[tnum]-data_select[tnum-1])**2.)
                #minmax=np.percentile(mlist,[minmax_percentile,1-minmax_percentile])
                #self.noise_tblavg[chan]=np.mean(mlist[np.logical_and(mlist>=minmax[0],mlist<=minmax[1])])/2.
                self.meta_params.noise_tblavg[chan]=np.median(mlist)/(2.*np.log(2.))   
                        
    def _compute_chiSQ(self,mode='dTIMEmTIMEBL'):
        """
        Compute a chi-square/DoF for each antenna gain solution by evaluationg
        \chi^2_i = \sum_{j,t} |v_{ij}-y_{i-j}g_i^*g_j |^2/\sigma_{ij}^2
        where \sigma_{ij}^2 at each frequency and baseline is estimated 
        """
        #assert mode in ['TIME','dFREQmBL']
        self._compute_noise(self,mode)
        #compute chi-squared values for each antenna, pol, and time.
        for pol in range(self.measured_vis.Npols):
            for chan in range(self.measured_vis.Nfreqs):
                for antnum in self.measured_vis.antenna_numbers:
                    for tnum in self.measured_vis.antenna_numbers:
                        selection=np.logical_and(self.measured_vis.ant_1_array==\
                                                 antnum1,
                                                 self.measured_vis.time_array==\
                                                 self.uvcal.time_array[tnum])
                        flag_select=self.measured_vis.flag_array[selection,self.meta_params.spw,
                                                                 chan,pol]
                        data_select=self.measured_vis.data_array[selection,self.meta_params.spw,
                                                                 chan,pol][np.invert(flag_select)]
                        model_select=self.measured_vis.data_array[selection,self.meta_params.spw,
                                                                  chan,pol][np.invert(flag_select)]
                        weight_select=self.cal_flag_weights.weights_array[selection,self.meta_params.spw,
                                                            chan,pol]
                        weight_select/=weight_select.sum()
                        ant2_select=self.ant_2_array[selection][np.invert(flag_select)]
                        gain_select=self.uvcal.gain_array[ant2_select,
                                                          self.meta_params.spw,
                                                          chan,
                                                          tnum,
                                                          pol]
                        this_gain=self.uvcal.gain-array[antnum,
                                                        self.meta_params.spw,
                                                        chan,
                                                        tnum,
                                                        pol]
                        self.meta_params.chi_square_per_ant[antnum,chan,tnum,pol]=\
                        np.sum(np.abs((data_select-this_gain*gain_select*model_select)*weight_select)**2.)/self.meta_params.noise_tblavg[chan]
                        neff=np.sum(weights_select)/np.max(weights_select)
                        self.meta_params.dof_per_ant[antnum,chan,tnum,pol]=neff-1
    
    def _read_files(self,data,mode,flag_weights_fromdata,flagweightsfile=None,model=None):
        '''
        read in all files including data,model,calibration weights
        args: 
           data, data file (fhd,uvfits,measurement,miriad set etc...)
           mode, specify type of data/model files
           model, model file (uvfits, miriad)
        '''
        assert mode in ['MS','UVFITS','MIRIAD','FHD']
        if mode=='MS':
            self._load_ms(data)
            self.meta_params.model_file=data
        elif mode=='UVFITS':
            assert model
            self._load_uvfits(data,model)
            self.meta_params.model_file=model
        elif mode=='MIRIAD':
            assert model
            self.__load_miriad(data,model)
            self.meta_params.model_file=model
        elif mode=='FHD':
            self._load_fhd(data)
        if(flag_weights_fromdata):
            self.cal_flag_weights.from_file(self.measured_vis,mode='UVDATA')
            self.meta_params.flag_weights_file=data
        else:
            assert flagweightsfile
            self.cal_flag_weights.from_file(flagweightsfile)
            self.meta_params.flag_weights_file=flagweightsfile
        self.meta_params.Niterations=np.zeros((self.meta_params.n_cycles,
                                               self.model_vis.Nfreqs,
                                               self.model_vis.Npols),dtype=int)
        self.meta_params.data_file=data
        self.meta_params.Nfreqs=self.model_vis.Nfreqs
        self.meta_params.Nants_data=self.model_vis.Nants_data
        self.meta_params.Njones=self.model_vis.Npols
        self.meta_params.Ntimes=self.model_vis.Ntimes
        self.meta_params.Nants_telescope=self.model_vis.Nants_telescope
        self.meta_params.chi_square_per_ant=np.zeros((self.measured_vis.Nants_data,
                                                      self.measured_vis.Nfreqs,
                                                      self.measured_vis.Ntimes,
                                                      self.measured_vis.Npols))
        self.meta_params.dof_per_ant=np.zeros_like(self.meta_params.chi_square_per_ant)
        self.meta_params.noise_tblavg=np.zeros(self.measured_vis.Nfreqs)
        self.meta_params.Ntime_steps=int(np.ceil(self.measured_vis.Ntimes/self.meta_params.t_avg))
        self._check_consistency()
        #compute noise matrices
        self._compute_noise()
        #initializes uvcal object properties
        self.uvcal.Njones=self.measured_vis.Npols
        self.uvcal.Nfreqs=self.measured_vis.Nfreqs
        self.uvcal.Ntimes=self.measured_vis.Ntimes
        self.uvcal.Nspws=self.measured_vis.Nspws
        self.uvcal.time_range=(self.measured_vis.time_array.min(),self.measured_vis.time_array.max())
        self.uvcal.telescope_name=self.measured_vis.telescope_name
        self.uvcal.Nants_data=self.measured_vis.Nants_data
        self.uvcal.Nants_telescope=self.measured_vis.Nants_telescope
        self.uvcal.ant_array=np.unique(self.measured_vis.ant_1_array)
        self.uvcal.antenna_names=self.measured_vis.antenna_names
        self.uvcal.antenna_numbers=self.measured_vis.antenna_numbers
        self.uvcal.freq_array=self.measured_vis.freq_array
        self.uvcal.channel_width=self.measured_vis.channel_width
        self.uvcal.jones_array=self.measured_vis.polarization_array
        self.uvcal.time_array=np.unique(self.measured_vis.time_array)
        self.uvcal.integration_time=self.measured_vis.integration_time
        self.uvcal.x_orientation='east' #stefcal is going to treat x as east.
        self.uvcal.cal_type='gain'
        self.uvcal.git_origin_cal='calibrated with stefcal_uvdata version %s with run id %s'\
                                   %(self.meta_params.stefcal_version_str,self.meta_params.id)
        #initialize uvcal gain and quality arrays
        self.uvcal.gain_array=9e99*np.ones((self.meta_params.Nants_data,1,self.model_vis.Nfreqs,self.model_vis.Ntimes,self.model_vis.Npols),dtype=complex)
        self.uvcal.flag_array=np.empty((self.meta_params.Nants_data,1,self.model_vis.Nfreqs,self.model_vis.Ntimes,self.model_vis.Npols),dtype=bool)

        
        

        
    def from_ms(self,msfile,flag_weights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from a measurement set
        args:
             msfile, name of measurement set cirectory
             flag_weights_fromdata, true if you want to determine 
                                   the flagweights object from the measurement set
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(msfile,mode='MS',
                         flag_weights_fromdata=flag_weights_fromdata,
                         flagweightsfile=flagweightsfile)
    def from_miriad(self,miriaddata,miriadmodel,flag_weights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from miriad files
        args:
             miriaddata, name of miriad data file
             miriadmodel, name of mirad model file
             flag_weights_fromdata, true if you want to determine 
                                   the flagweights object from the miriad data
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(miriaddata,mode='MIRIAD',
                         flag_weights_fromdata=flag_weights_fromdata,
                         flagweightsfile=flagweightsfile,
                         model=miriadmodel)
    def from_fhd(self,fhdfile,flag_weights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from a fhd
        args:
             fhdfile, name of fhd file
             flag_weights_fromdata, true if you want to determine 
                                   the flagweights object from the measurement set
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(fhdfile,mode='FHD',
                         flag_weights_fromdata=flag_weights_fromdata,
                         flagweightsfile=flagweightsfile)
        
    def from_uvfits(self,uvfitsdata,uvfitsmodel,flag_weights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from miriad files
        args:
             uvfitsdata, name of uvfits data file
             uvfitsmodel, name of uvfits model file
             flag_weights_fromdata, true if you want to determine 
                                   the flagweights object from the miriad data
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(uvfitsdata,mode='UVFITS',
                         flag_weights_fromdata=flag_weights_fromdata,
                         flagweightsfile=flagweightsfile,
                         model=uvfitsmodel)

    def _blt_list_2_matrix(self,blt_list,t_steps,hermit=True):
        """
        convert a baseline-time ordered list into a 
        Ntimes x NAnt x NAnt matrix 
        Args:
            blt_list, baseline-time ordered numpy array 
            t_steps, integer time-steps defining Ntimes axis. 
            hermit, boolean, set true if output[a,b]=conj(output[b,a])
        """
        timeList=np.unique(self.measured_vis.time_array)
        antList=np.unique(self.measured_vis.antenna_numbers)
        nT=len(t_steps)
        nAnt=len(antList)
        if DEBUG:
            print('number of antennas ='+str(nAnt))
        output=np.empty((nT,nAnt,nAnt),dtype=blt_list.dtype)
        for t_step in range(nT):
            selection=self.measured_vis.time_array==timeList[t_step]
            bl_list=blt_list[selection]
            for antNum2 in range(nAnt):
                ant2=antList[antNum2]
                for antNum1 in range(antNum2):
                    ant1=antList[antNum1]
                    selection=np.logical_and(self.measured_vis.time_array\
                                             ==timeList[t_step],
                                             np.logical_and(self.measured_vis.ant_1_array==ant1,
                                                            self.measured_vis.ant_2_array==ant2))
                    output[t_step,antNum1,antNum2]=blt_list[selection]
                    if hermit:
                        output[t_step,antNum1,antNum2]=np.conj(blt_list[selection])
                    else:
                        output[t_step,antNum1,antNum2]=blt_list[selection]
        return output

    #Function to reinitialize from files and id
    
    
    def save_state(self,output_root,clobber=False):
        """
        saves cal_flag_weights and meta information to .pkl files. 
        allows for reinitialization of stefcal at any future time
        """
        self.meta_params.to_file(output_root+'_meta',clobber)
        self.cal_flag_weights.to_file(output_root+'_cfws',clobber)
        #save calibration solution
        if clobber:
            pickle.dump(self.uvcal,open(output_root+'_'+self.meta_params.id+'_cal',"wb"))
        

    def load_state(self,input_root):
        """
        load state from meta and cal flag weights file. 
        """
        self.meta_params.from_file(input_root+'_meta')
        self.cal_flag_weights.from_file(input_root+'_cfws',mode='CFWS')
        self.uvcal=pickle.load(open(input_root+'_'+self.meta_params.id+'_cal',"rb"))
                    
                               
    def set_refant(self,refant):
        """
        set refant, int
        """
        self.meta_params.refant=refant

    def set_n_phase_iter(self,n_phase_iter):
        """
        set n_phase_iter, int
        """
        self.meta_params.self_n_phase_iter=n_phase_iter
        

    def set_n_cycles(self,n_cycles):
        """
        set n_cycles, int
        """
        self.meta_params.n_cycles=n_cycles
        self.meta_params.Nitrations=np.zeros((self.meta_params.n_cycles,
                                              self.meta_params.Nfreqs,
                                              self.meta_params.Njones),dtype=int)
    def set_min_bl_per_ant(self,min_bl_per_ant):
        """
        set min_bl_per_ant, int
        """
        self.min_bl_per_ant=min_bl_per_ant
    def set_eps(self,eps):
        """
        set eps, float
        """
        self.eps=eps
        
    def set_min_ant_times(self,min_ant_times):
        """
        set min_ant_times, int
        """
        self.min_ant_times=min_ant_times
    def set_trim_neff(self,trim_neff):
        """
        set trim_neff, bool
        """
        self.trim_neff=trim_neff

    def set_params(self,param_dict):
        """
        set stefcal parameters 
        Args:
            param_dict: dictionary of params that can include any one of
            trim_neff: boolean, determine whether you want stefcal to flag baselines 
                         that are causing small neff for each antenna and slow convergence.
            min_ant_times: int, minimum number of time-steps within a given calibration interval 
                           for an antenna to not be flagged for that interval
            eps: float, minimum values of 
            min_bl_per_ant: int, minimum neff per antenna required to not flag it from aa cal solution. 
            n_cycles: int, number of cycles to run stefcal for (a cycle is a full set of stefcal iterations 
                      to bring maximum difference between successive iterations below eps.
            n_phase_iter: int, number of iterations to execute per cycle in which on the phase is solved for
            refant: int, the reference antenna
        """
        input_keys=param_dict.keys()
        if 'trim_neff' in input_keys:
            self.set_trim_neff(param_dict['trim_neff'])
        if 'min_ant_times' in input_keys:
            self.set_min_ant_times(param_dict['min_ant_times'])
        if 'eps' in input_keys:
            self.set_eps(param_dict['eps'])
        if 'min_bl_per_ant' in input_keys:
            self.set_min_bl_per_ant(param_dict['min_bl_per_ant'])
        if 'n_cycles' in input_keys:
            self.set_n_cycles(param_dict['n_cycles'])
        if 'n_phase_iter' in input_keys:
            self.set_n_phase_iter(param_dict['n_phase_iter'])
        if 'refant' in input_keys:
            self.set_refant(param_dict['refant'])
           
        
    def stefcalibrate(self,parallelized=False):
        '''
        Run stefcal
        Args: 
            parallilized, choose this if you want calibration to be parallelized across frequency channels
        '''
        #loop through polarization
        for pol in range(self.measured_vis.Npols):
            #for each channel and time-averaging step, 
            for chan in range(self.measured_vis.Nfreqs):
                for tstep in range(self.meta_params.Ntime_steps):
                    t_steps=range(tstep*self.meta_params.t_avg,np.min([(tstep+1)*self.meta_params.t_avg,
                                                                       self.measured_vis.Ntimes-1]))
                    data_mat=self._blt_list_2_matrix(self.measured_vis.data_array[:,self.meta_params.spw,chan,pol].squeeze(),t_steps,hermit=True)
                    model_mat=self._blt_list_2_matrix(self.model_vis.data_array[:,self.meta_params.spw,chan,pol].squeeze(),t_steps,hermit=True)
                    weights_mat=self._blt_list_2_matrix(self.cal_flag_weights.weights_array[:,self.meta_params.spw,chan,pol].squeeze(),t_steps,hermit=True)
                    flags_mat=self._blt_list_2_matrix(self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,pol].squeeze(),t_steps,hermit=False)
                    ant_flags,flag_matrix,niter,gains=stefcal.stefcal_scaler(data_mat,model_mat,
                                                                             weights_mat,
                                                                             flags_mat,
                                                                             self.meta_params.refant,
                                                                             self.meta_params.n_phase_iter,
                                                                             self.meta_params.n_cycles,
                                                                             self.meta_params.min_bl_per_ant,
                                                                             self.meta_params.eps,
                                                                             self.meta_params.min_ant_times,
                                                                             self.meta_params.trim_neff)
                    if DEBUG:
                        print('gain_array.shape='+str(self.uvcal.gain_array.shape))
                        print('spw='+str(self.meta_params.spw))
                        print('chan='+str(chan))
                        print('t_steps='+str(t_steps))
                        print('pol='+str(pol))
                        print('gains.shape='+str(gains.shape))
                        print('ant_flags.shape='+str(ant_flags.shape))
                        print('flag_array.shape='+str(self.uvcal.flag_array.shape))
                        print('cfw.flag_array.shape='+str(self.cal_flag_weights.flag_array.shape))
                        print('flag_matrix.shape='+str(flag_matrix.shape))
                    for ts in t_steps:
                        self.uvcal.gain_array[:,self.meta_params.spw,chan,ts,pol]=gains
                        self.uvcal.flag_array[:,self.meta_params.spw,chan,ts,pol]=ant_flags
                        self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,ts,pol]=np.logical_or(flag_matrix,self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,ts,pol])
                        #Need to translate back into blt list!
                        self.meta_params.Niterations[ts,chan,pol,:]=niter
            #compute chi-squares
            self._compute_chiSQ()
            #get data, model, weights, flags
            self.quality_array[:,self.meta_params.spw,chan,t_steps,pol]=_compute_chiSQ
