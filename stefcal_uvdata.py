DEBUG=False
from pyuvdata import UVData
import copy
import numpy as np
from pyuvdata import UVCal
import pyuvdata.parameter as uvp
from cal_flag_weights import CalFlagWeights
import stefcal
import pickle
from stefcal_meta import StefcalMeta
import uuid
import stefcal_utils as utils


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
                 min_ant_times=1,spw=0,
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
        #self.meta_params.trim_neff=trim_neff
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
        #****************************************
        #There is a bug in uvfits.py that is breaking
        #antenna names!
        #****************************************
        #if model.antenna_names!=data.antenna_names:
        #    checkpass=False
        #    print("antenna_names not the same")
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


    def _load_vis_ms(self,msname):
        '''
        read in ms visibilities
        args: msName, path to .ms directory
        '''
        #need to throw an error here if the model column does not exist
        self.model_vis.read_ms(msname,columnname='MODEL')
        self.measured_vis.read_ms(msname,columnname='DATA')

    def _load_ms(self,dataname):
        '''
        read in visiblity files for ms files
        args: dataname, name of measurement set
        '''
        self.model_vis.read_ms(dataname,data_column='MODEL')
        self.measured_vis.read_ms(dataname,data_column='DATA')
    def _load_fhd(self,dataname,modelname):
        '''
        read in fhd visibility files
        args: dataname, name of the fhd file
        '''
        self.measured_vis.read_fhd(dataname)
        self.model_vis.read_fhd(modelname,use_model=True)
        
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
        assert mode in ['dTIMEmTIMEBL','dFREQmTIMEBL']
        if mode=='dTIMEmTIMEBL' and self.measured_vis.Ntimes>2:
            #take diff in time and average in time,baselines,and pols; good for measurement sets
            #with a small number of time samples.
            self.meta_params.noise_tblavg=np.zeros((self.measured_vis.Nfreqs))
            for chan in range(self.measured_vis.Nfreqs):
                mlist=np.array([])
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
                if self.meta_params.noise_tblavg[chan]==0.:
                    self.meta_params.noise_tblavg[chan]=9e99
        elif mode=='dFREQmTIMEBL' or self.measured_vis.Ntimes<2:
            #also able to take differences in frequencies
            self.meta_params.noise_tblavg=np.zeros((self.measured_vis.Nfreqs))
            for chan in range(1,self.measured_vis.Nfreqs):
                mlist=np.array([])
                for pol in range(self.measured_vis.Npols):
                    for blnum in range(self.measured_vis.Nbls):
                        bl_selection=self.measured_vis.baseline_array==self.measured_vis.baseline_array[blnum]
                        data_select=self.measured_vis.data_array[bl_selection,self.meta_params.spw,chan,pol]
                        flag_select=self.measured_vis.flag_array[bl_selection,self.meta_params.spw,chan,pol]
                        data_select_pc=self.measured_vis.data_array[bl_selection,self.meta_params.spw,chan-1,pol]
                        flag_select_pc=self.measured_vis.flag_array[bl_selection,self.meta_params.spw,chan-1,pol]
                        
                        for tnum in range(self.measured_vis.Ntimes):
                            if not(flag_select[tnum]) and not (flag_select_pc[tnum]):
                                mlist=np.append(mlist,np.abs(data_select[tnum]-data_select_pc[tnum])**2.)
                self.meta_params.noise_tblavg[chan]=np.median(mlist)/(2.*np.log(2.))
                if self.meta_params.noise_tblavg[chan]==0.:
                    self.meta_params.noise_tblavg[chan]=9e99
            self.meta_params.noise_tblavg[0]=self.meta_params.noise_tblavg[1]#set zero bin to 1 bin if using freq avg.
                        
    def _compute_chiSQ(self,mode='dTIMEmTIMEBL'):
        """
        Compute a chi-square/DoF for each antenna gain solution by evaluationg
        \chi^2_i = \sum_{j,t} |v_{ij}-y_{i-j}g_i^*g_j |^2/\sigma_{ij}^2
        where \sigma_{ij}^2 at each frequency and baseline is estimated 
        """
        ant_dict={}
        ant_dict_inv={}
        u_ants=np.unique(np.vstack([self.measured_vis.ant_1_array,self.measured_vis.ant_2_array]))
        for antind in range(len(u_ants)):
            ant_dict[antind]=u_ants[antind]
            ant_dict_inv[u_ants[antind]]=antind
        #assert mode in ['TIME','dFREQmBL']
        self._compute_noise(mode=mode)
        #compute chi-squared values for each antenna, pol, and time.
        for tnum in range(self.measured_vis.Ntimes):
            selection=self.measured_vis.time_array==self.uvcal.time_array[tnum]
            for pol in range(self.measured_vis.Npols):
                for chan in range(self.measured_vis.Nfreqs):
                    flag_select=self.cal_flag_weights.flag_array[selection,self.meta_params.spw,
                                                                 chan,pol]
                    ant1_select=np.array([ant_dict_inv[m] for m in self.measured_vis.ant_1_array[selection]]).astype(int)
                    ant2_select=np.array([ant_dict_inv[m] for m in self.measured_vis.ant_2_array[selection]]).astype(int)
                    #ant1_select=ant_dict_inv[self.measured_vis.ant_1_array[selection]]
                    #ant2_select=ant_dict_inv[self.measured_vis.ant_2_array[selection]]
                    autos=ant1_select==ant2_select
                    flag_select=np.logical_or(flag_select,autos)
                    ant1_select=ant1_select[np.invert(flag_select)]
                    ant2_select=ant2_select[np.invert(flag_select)]
                    data_select=self.measured_vis.data_array[selection,self.meta_params.spw,
                                                             chan,pol][np.invert(flag_select)]
                    model_select=self.measured_vis.data_array[selection,self.meta_params.spw,
                                                              chan,pol][np.invert(flag_select)]
                    weight_select=self.cal_flag_weights.weights_array[selection,self.meta_params.spw,
                                                                      chan,pol][np.invert(flag_select)]
                    gain1_select=self.uvcal.gain_array[ant1_select,
                                                       self.meta_params.spw,
                                                       chan,
                                                       tnum,
                                                       pol]
                    gain2_select=self.uvcal.gain_array[ant2_select,
                                                       self.meta_params.spw,
                                                       chan,
                                                       tnum,
                                                       pol]
                    self.meta_params.chi_square[chan,tnum,pol]=np.sum(np.abs((data_select-np.conj(gain1_select)*gain2_select*model_select)*weight_select/weight_select.sum())**2.)
                    self.meta_params.chi_square[chan,tnum,pol]/=self.meta_params.noise_tblavg[chan]
                    if len(weight_select)>0:
                        neff=np.sum(weight_select)/np.max(weight_select)
                    else:
                        neff=0.
                    self.meta_params.dof[chan,tnum,pol]=neff-1
                    for antnum in u_ants:
                        selection_ant=np.logical_or(np.logical_xor(self.measured_vis.ant_1_array==antnum,
                                                                   self.measured_vis.ant_2_array==antnum),
                                                     self.measured_vis.time_array==\
                                                     self.uvcal.time_array[tnum])
                        flag_ant_select=self.cal_flag_weights.flag_array[selection_ant,self.meta_params.spw,
                                                                         chan,pol]
                        data_ant_select=self.measured_vis.data_array[selection_ant,self.meta_params.spw,
                                                                 chan,pol][np.invert(flag_ant_select)]
                        model_ant_select=self.measured_vis.data_array[selection_ant,self.meta_params.spw,
                                                                  chan,pol][np.invert(flag_ant_select)]
                        weight_ant_select=self.cal_flag_weights.weights_array[selection_ant,self.meta_params.spw,
                                                            chan,pol][np.invert(flag_ant_select)]
                        ant2_ant_select=np.array([ant_dict_inv[ant2] for ant2 in self.measured_vis.ant_2_array[selection_ant][np.invert(flag_ant_select)]]).astype(int)
                        ant1_ant_select=np.array([ant_dict_inv[ant1] for ant1 in self.measured_vis.ant_1_array[selection_ant][np.invert(flag_ant_select)]]).astype(int)
                        #ant2_ant_select=ant_dict[self.measured_vis.ant_2_array[selection_ant][np.invert(flag_ant_select)]]
                        #ant1_ant_select=ant_dict[self.measured_vis.ant_1_array[selection_ant][np.invert(flag_ant_select)]]
                        gain2_ant_select=self.uvcal.gain_array[ant2_ant_select,
                                                               self.meta_params.spw,
                                                               chan,
                                                               tnum,
                                                               pol]
                        gain1_ant_select=self.uvcal.gain_array[ant1_ant_select,
                                                               self.meta_params.spw,
                                                               chan,
                                                               tnum,
                                                               pol]

                        self.meta_params.chi_square_per_ant[ant_dict_inv[antnum],chan,tnum,pol]=\
                        np.sum(np.abs((data_ant_select-np.conj(gain1_ant_select)*gain2_ant_select*model_ant_select)*weight_ant_select/weight_ant_select.sum())**2.)
                        self.meta_params.chi_square_per_ant[ant_dict_inv[antnum],chan,tnum,pol]/=self.meta_params.noise_tblavg[chan]
                        if len(weight_ant_select)>0:
                            neff_ant=np.sum(weight_ant_select)/np.max(weight_ant_select)
                        else:
                            neff_ant=0.
                        self.meta_params.dof_per_ant[ant_dict_inv[antnum],chan,tnum,pol]=neff_ant-1
        self.uvcal.quality_array=self.meta_params.chi_square_per_ant/self.meta_params.dof_per_ant
    def _read_files(self,data,mode,flag_weights_fromdata,flagweightsfile=None,model=None,selection={}):
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
            self._load_fhd(data,model)
        #perform selection operations:
        if not 'antenna_nums' in selection.keys():
            selection['antenna_nums']=None
        if not 'antenna_names' in selection.keys():
            selection['antenna_names']=None
        if not 'ant_pairs_nums' in selection.keys():
            selection['ant_pairs_nums']=None
        if not 'frequencies' in selection.keys():
            selection['frequencies']=None
        if not 'freq_chans' in selection.keys():
            selection['freq_chans']=None
        if not 'times' in selection.keys():
            selection['times']=None
        if not 'polarizations' in selection.keys():
            selection['polarizations']=None
        if not 'blt_inds' in selection.keys():
            selection['blt_inds']=None
        self.model_vis.select(antenna_nums=selection['antenna_nums'],
                              ant_pairs_nums=selection['ant_pairs_nums'],
                              frequencies=selection['frequencies'],
                              freq_chans=selection['freq_chans'],
                              times=selection['times'],
                              polarizations=selection['polarizations'],
                              blt_inds=selection['blt_inds'])
        self.measured_vis.select(antenna_nums=selection['antenna_nums'],
                                 ant_pairs_nums=selection['ant_pairs_nums'],
                                 frequencies=selection['frequencies'],
                                 freq_chans=selection['freq_chans'],
                                 times=selection['times'],
                                 polarizations=selection['polarizations'],
                                 blt_inds=selection['blt_inds'])
        
        if(flag_weights_fromdata):
            self.cal_flag_weights.from_file(self.measured_vis,mode='UVDATA')
            self.meta_params.flag_weights_file=data
        else:
            assert flagweightsfile
            self.cal_flag_weights.from_file(flagweightsfile)
            self.meta_params.flag_weights_file=flagweightsfile
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
        self.meta_params.chi_square=np.zeros((self.measured_vis.Nfreqs,
                                              self.measured_vis.Ntimes,
                                              self.measured_vis.Npols))
        self.meta_params.dof=np.zeros_like(self.meta_params.chi_square)
        self.meta_params.dof_per_ant=np.zeros_like(self.meta_params.chi_square_per_ant)
        self.meta_params.noise_tblavg=np.zeros(self.measured_vis.Nfreqs)
        self.meta_params.Ntime_steps=int(np.ceil(self.measured_vis.Ntimes/self.meta_params.t_avg))
        self.meta_params.Niterations=np.zeros((self.meta_params.Ntime_steps,
                                               self.meta_params.n_cycles,
                                               self.model_vis.Nfreqs,
                                               self.model_vis.Npols),dtype=int)

        self.cal_flag_weights.ant_flag_array=np.empty((self.measured_vis.Nants_data,
                                                       self.measured_vis.Nspws,
                                                       self.measured_vis.Ntimes,
                                                       self.measured_vis.Npols),dtype=bool)
        self.cal_flag_weights.ant_flag_array[:]=False
        
        self._check_consistency()
        #compute noise matrices
        self._compute_noise()
        #initializes uvcal object properties
        #self.model_vis.data_array=self.model_vis.data_array.astype('complex128')
        #self.measured_vis.data_array=self.measured_vis.data_array.astype('complex128')
        self.uvcal=self.uvcal_from_data()

        
    
    def uvcal_from_data(self):
        """
        Generate an empty uvcal object from visibility parameters
        """
        uvc=UVCal()
        uvc.Njones=self.measured_vis.Npols
        uvc.Nfreqs=self.measured_vis.Nfreqs
        uvc.Ntimes=self.measured_vis.Ntimes
        uvc.Nspws=self.measured_vis.Nspws
        uvc.time_range=(self.measured_vis.time_array.min(),self.measured_vis.time_array.max())
        uvc.telescope_name=self.measured_vis.telescope_name
        uvc.Nants_data=self.measured_vis.Nants_data
        uvc.Nants_telescope=self.measured_vis.Nants_telescope
        uvc.ant_array=np.unique(self.measured_vis.ant_1_array)
        uvc.antenna_names=self.measured_vis.antenna_names
        uvc.antenna_numbers=self.measured_vis.antenna_numbers
        uvc.freq_array=self.measured_vis.freq_array
        uvc.channel_width=self.measured_vis.channel_width
        uvc.jones_array=self.measured_vis.polarization_array
        uvc.time_array=np.unique(self.measured_vis.time_array)
        uvc.integration_time=self.measured_vis.integration_time
        uvc.x_orientation='east'#always
        uvc.cal_type='gain'
        uvc.quality_array=np.zeros((self.measured_vis.Nants_data,
                                    1,
                                    self.measured_vis.Nfreqs,
                                    self.measured_vis.Ntimes,
                                    self.measured_vis.Npols))
        uvc.git_origin_cal='calibrated with stefcal_uvdata version %s with run id %s'\
                            %(self.meta_params.stefcal_version_str,self.meta_params.id)
        uvc.gain_array=np.ones((self.meta_params.Nants_data,1,self.model_vis.Nfreqs,self.model_vis.Ntimes,self.model_vis.Npols),dtype=complex)
        uvc.flag_array=np.empty((self.meta_params.Nants_data,1,self.model_vis.Nfreqs,self.model_vis.Ntimes,self.model_vis.Npols),dtype=bool)
        uvc.flag_array[:]=False
        return uvc
        

        
    def from_ms(self,msfile,flag_weights_fromdata,
                flagweightsfile=None,selection={}):
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
                         flagweightsfile=flagweightsfile,selection=selection)
    def from_miriad(self,miriaddata,miriadmodel,
                    flag_weights_fromdata,flagweightsfile=None,
                    select_data=False,selection={}):
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
                         model=miriadmodel,
                         selection=selection)
    def from_fhd(self,fhddata,fhdmodel,
                 flag_weights_fromdata,
                 flagweightsfile=None,
                 selection={}):
        """
        initialize stefcal from a fhd
        args:
             fhdfile, name of fhd file
             flag_weights_fromdata, true if you want to determine 
                                   the flagweights object from the measurement set
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(fhddata,mode='FHD',
                         flag_weights_fromdata=flag_weights_fromdata,
                         flagweightsfile=flagweightsfile,
                         model=fhdmodel,selection=selection)
        
    def from_uvfits(self,uvfitsdata,uvfitsmodel,
                    flag_weights_fromdata,
                    flagweightsfile=None,
                    selection={}):
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
                         model=uvfitsmodel,selection=selection)

    def _matrix_2_blt_list(self,blt_matrix,ant_dict,ant_flags):
        """
        convert a matrix of Ntimes x Nant x Nant 
        to a basline-time ordered list
        Args: 
            blt_matrix, NtimesxNantxNant matrix
            ant_dict, mapping from antenna indices to antenna numbers.
        Returns: 
            blt_list, baseline-time ordered list
        """
        iscomplex=(blt_matrix.dtype=='complex64' or blt_matrix.dtype=='complex128')
        timeList=np.unique(self.measured_vis.time_array)
        ant_dict_c=copy.copy(ant_dic)
        nAnt=len(ant_dict)
        ant_flags_all=np.all(ant_flags,axis=0)
        ant_dict={}
        for ind,key in enumerate(ant_dict_c.keys()):
            ant_dict[ind]=ant_dict_c[key]
        nAnt=len(nAnt)
        output=np.empty(self.measured_vis.Nblts,dtype=blt_matrix.dtype)
        for t_step in range(self.measured_vis.Ntimes):
            t_selection=self.measured_vis.time_array==timeList[t_step]
            for antNum2 in range(nAnt):
                ant2=ant_dict[antNum2]
                for antNum1 in range(antNum2):
                    ant1=ant_dict[antNum1]
                    a_selection=np.logical_and(self.measured_vis.ant_1_array==ant1,
                                               self.measured_vis.ant_2_array==ant2)
                    if not(np.any(a_selection)):
                        a_selection=np.logical_and(self.measured_vis.ant_1_array==ant2,
                                                   self.measured_vis.ant_2_array==ant1)
                        conjugate=True
                    else:
                        conjugate=False
                            
                    selection=np.logical_and(a_selection,t_selection)
                    output[selection]=blt_matrix[t_step,antNum1,antNum2]
                    if conjugate:
                        output[selection]=np.conj(output[selection])
        return output
                    
        
    def _blt_list_2_matrix(self,blt_list,t_steps,ant_dict,ant_flags,hermit=True):
        """
        convert a baseline-time ordered list into a 
        Ntimes x NAnt x NAnt matrix 
        Args:
            blt_list, baseline-time ordered numpy array 
            t_steps, integer time-steps defining Ntimes axis. 
            ant_dict, dictionary from antenna indices to antenna numbers. 
            hermit, boolean, set true if output[a,b]=conj(output[b,a])
        """
        iscomplex=(blt_list.dtype=='complex64' or blt_list.dtype=='complex128')
        timeList=np.unique(self.measured_vis.time_array)
        ant_dict_c=copy.copy(ant_dict)
        nT=len(t_steps)
        nAnt=len(ant_dict)
        ant_flags_all=np.all(ant_flags,axis=0)
        #filter bad antennas
        #rebuild index
        ant_dict={}
        for ind,key in enumerate(ant_dict_c.keys()):
            ant_dict[ind]=ant_dict_c[key]
        nAnt=len(ant_dict)
        output=np.empty((nT,nAnt,nAnt),dtype=blt_list.dtype)
        if str(output.dtype)=='bool':
            output[:]=True
        if str(output.dtype)=='complex64' or str(output.dtype)=='float64':
            output[:]=0.
        if output.dtype==np.int:
            output[:]=0
        for ts,t_step in enumerate(t_steps):
            t_selection=self.measured_vis.time_array==timeList[t_step]
            for antNum2 in range(nAnt):
                ant2=ant_dict[antNum2]
                for antNum1 in range(antNum2):
                    ant1=ant_dict[antNum1]
                    a_selection=np.logical_and(self.measured_vis.ant_1_array==ant1,
                                               self.measured_vis.ant_2_array==ant2)
                    if not(np.any(a_selection)):
                        a_selection=np.logical_and(self.measured_vis.ant_1_array==ant2,
                                                   self.measured_vis.ant_2_array==ant1)
                        conjugate=True
                    else:
                        conjugate=False
                    selection=np.logical_and(a_selection,t_selection)
                    output[ts,antNum1,antNum2]=blt_list[selection][0]
                    if hermit:
                        output[ts,antNum2,antNum1]=np.conj(blt_list[selection])[0]
                    else:
                        output[ts,antNum2,antNum1]=blt_list[selection][0]
                    if conjugate and iscomplex:
                        output[ts,antNum2,antNum1]=np.conj(output[ts,antNum2,antNum1])
                        output[ts,antNum1,antNum2]=np.conj(output[ts,antNum1,antNum2])
                                                           
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

                
    def set_weights(self,weights):
        """
        set weights array
        """
        assert(self.cal_flag_weights.weights_array.shape==\
               weights.shape)
        assert(weights.dtype==np.float64)
        self.cal_flag_weights.weights_array=weights

    def set_flags(self,flags):
        """
        set flags array
        """
        assert(self.cal_flag_weights.flag_array.shape==\
               flags.shape)
        assert flags.dtype==bool
        self.cal_flag_weights.flag_array=flags
    def set_ant_flags(self,ant_flags):
        assert(ant_flags.shape==self.cal_flag_weights.ant_flag_array.shape)
        assert ant_flags.dtype==bol
        self.cal_flag_weights.ant_flag_array=ant_flags
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
        self.meta_params.n_phase_iter=n_phase_iter
        

    def set_n_cycles(self,n_cycles):
        """
        set n_cycles, int
        """
        self.meta_params.n_cycles=n_cycles
        self.meta_params.Niterations=np.zeros((self.meta_params.Ntimes,
                                               self.meta_params.n_cycles,
                                               self.meta_params.Nfreqs,
                                               self.meta_params.Njones),dtype=int)
    def set_min_bl_per_ant(self,min_bl_per_ant):
        """
        set min_bl_per_ant, int
        """
        self.meta_params.min_bl_per_ant=min_bl_per_ant
    def set_eps(self,eps):
        """
        set eps, float
        """
        self.meta_params.eps=eps

    def set_tavg(self,t_avg):
        """
        set number of time samples to average over.
        """
        self.meta_params.t_avg=t_avg
        self.meta_params.Ntime_steps=int(np.ceil(self.measured_vis.Ntimes/self.meta_params.t_avg))
        self.meta_params.Niterations=np.zeros((self.meta_params.Ntime_steps,
                                               self.meta_params.n_cycles,
                                               self.model_vis.Nfreqs,
                                               self.model_vis.Npols),dtype=int)
    def set_min_ant_times(self,min_ant_times):
        """
        set min_ant_times, int
        """
        self.meta_params.min_ant_times=min_ant_times
    #def set_trim_neff(self,trim_neff):
    #    """
    #    set trim_neff, bool
    #    """
    #    self.meta_params.trim_neff=trim_neff

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
            flags: bool, Nblts x Nspws x Nfreqs x Npols
        """
        input_keys=param_dict.keys()
        #if 'trim_neff' in input_keys:
        #    self.set_trim_neff(param_dict['trim_neff'])
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
        if 't_avg' in input_keys:
            self.set_tavg(param_dict['t_avg'])
        if 'weights' in input_keys:
            self.set_weights(param_dict['weights'])
        if 'flags' in input_keys:
            self.set_flags(param_dict['flags'])
    def stefcalibrate(self,perterb=0,parallelized=False):
        '''
        Run stefcal
        Args: 
            parallilized, choose this if you want calibration to be parallelized across frequency channels
        '''
        #antinds=self.measured_vis.antenna_numbers
        #if antinds.max()==len(antinds):
        #    indsub=1
        #else:
        #    indsub=0
        u_ant=np.unique(np.vstack([self.measured_vis.ant_1_array,self.measured_vis.ant_2_array]))
        ant_dict={}
        for antind,antnum in enumerate(u_ant):
            ant_dict[antind]=antnum
        #loop through polarization
        for pol in range(self.measured_vis.Npols):
            #for each channel and time-averaging step, 
            for chan in range(self.measured_vis.Nfreqs):
                full_flags=np.empty((self.measured_vis.Ntimes,self.measured_vis.Nants_data,self.measured_vis.Nants_data),dtype=bool)
                for tstep in range(self.meta_params.Ntime_steps):
                    #print('calibrating pol=%d,chan=%d,tstep=%d'%(pol,chan,tstep))
                    t_steps=range(tstep*self.meta_params.t_avg,np.min([(tstep+1)*self.meta_params.t_avg,
                                                                       self.measured_vis.Ntimes]))
                    data_mat=self._blt_list_2_matrix(self.measured_vis.data_array[:,self.meta_params.spw,chan,pol].squeeze(),
                                                     t_steps,ant_dict,ant_flags=self.cal_flag_weights.ant_flag_array[:,self.meta_params.spw,t_steps,pol])
                    model_mat=self._blt_list_2_matrix(self.model_vis.data_array[:,self.meta_params.spw,chan,pol].squeeze(),
                                                      t_steps,ant_dict,ant_flags=self.cal_flag_weights.ant_flag_array[:,self.meta_params.spw,t_steps,pol])
                    weights_mat=self._blt_list_2_matrix(self.cal_flag_weights.weights_array[:,self.meta_params.spw,chan,pol].squeeze(),
                                                        t_steps,ant_dict,ant_flags=self.cal_flag_weights.ant_flag_array[:,self.meta_params.spw,t_steps,pol])
                    flags_mat=self._blt_list_2_matrix(self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,pol].squeeze(),
                                                      t_steps,ant_dict,ant_flags=self.cal_flag_weights.ant_flag_array[:,self.meta_params.spw,t_steps,pol])
                    #if DEBUG:
                    #    print('t_steps='+str(t_steps))
                    #    print('weights_mat.shape='+str(weights_mat.shape))
                    #    print('weights_mat='+str(weights_mat[0,0,:]))
                    #    print('data_mat='+str(data_mat[0,0,:]))
                    #    print('flags_mat='+str(flags_mat[0,0,:]))
                    #    print('model_mat='+str(model_mat[0,0,:]))
                    #print('calibrating pol=%d,chan=%d,tstep=%d'%(pol,chan,tstep))
                    niter,gains=stefcal.stefcal_scaler(data_mat,model_mat,
                                                       weights_mat,
                                                       flags_mat,
                                                       refant=self.meta_params.refant,
                                                       n_phase_iter=self.meta_params.n_phase_iter,
                                                       n_cycles=self.meta_params.n_cycles,
                                                       min_bl_per_ant=self.meta_params.min_bl_per_ant,
                                                       eps=self.meta_params.eps,
                                                       min_ant_times=self.meta_params.min_ant_times,
                                                       #trim_neff=self.meta_params.trim_neff,
                                                       perterb=perterb)
                    #if DEBUG:
                        #print('gains.shape='+str(gains.shape))
                    print('success calibrating pol=%d,chan=%d,tstep=%d'%(pol,chan,tstep))
                    self.meta_params.Niterations[tstep,:,chan,pol]=niter
                    for tsn,ts in enumerate(t_steps):
                        aflag=np.invert(self.cal_flag_weights.ant_flag_array[:,self.meta_params.spw,pol]).flatten()
                        self.uvcal.gain_array[aflag,self.meta_params.spw,chan,ts,pol]=gains
                        #self.uvcal.flag_array[:,self.meta_params.spw,chan,ts,pol]=ant_flags[tsn]
                        #Need to translate back into blt list!
                        #if DEBUG:
                            #print('shape niterations='+str(self.meta_params.Niterations.shape))
                        #full_flags[ts,:,:]=flag_matrix[tsn]
                    #print('any flags in flag_matrix?'+str(np.any(full_flags)))
                #self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,pol]=np.logical_or(self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,pol],
                 #                                                                           self._matrix_2_blt_list(full_flags))
                #print('any flags?='+str(np.any(self.cal_flag_weights.flag_array[:,self.meta_params.spw,chan,pol])))


            
            #compute chi-squares
        self._compute_chiSQ()
        #get data, model, weights, flags
