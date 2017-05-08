from pyuvdata import UVData
from pyuvdata import UVCal
from calFlagWeights import CalFlagWeights
import stefcal
import uuid

class StefCalMeta():
    """
    Defines a container class for stefcal meta parameters 
    contains stefcal parameters along with filepaths
    to data, model, and flag_weights files. 
    """
    def __init__(self,id,refant=0,n_phase_iter=5,
                 n_cycles=1,min_bl_per_ant=2,eps=1e-10,
                 min_ant_times=1,trim_neff=False):
        self.refant=refant
        self.n_phase_iter=n_phase_iter
        self.min_bl_per_ant=min_bl_per_ant
        self.eps=eps
        self.min_ant_times=min_ant_times
        self.trim_neff=trim_neff
        self.id=id
        self.data_file=""
        self.flag_weights_file=""
        self.model_file=""
        


class StefCalUVData():
    """
    Defines a class for performing stefcal on uvdata sets.
    Attributes:
    model_vis: uvdata of model
    measured_vis: uvdata of measurements
    cal_flag_weights: CalFlagWeights object storing flags and weights
    cal_solution: uvcal object
    meta_params: StefCalMeta object containing information on calibration
                 run
    id: unique uuid corresponding to calibration run. 
    """
    def __init__(self,refant=0,n_phase_iter=5,
                 n_cycles=1,min_bl_per_ant=2,eps=1e-10,
                 min_ant_times=1,trim_neff=False):
        self.id=uuid.uuid4()
        self.model_vis=UVData()
        self.measured_vis=UVData()
        self.cal_flag_weights=CalFlagWeights(self.id)
        self.cal_solution=UVCal()
        self.meta_params=StefCalMeta(refant,
                                     n_phase_iter,
                                     n_cycles,
                                     min_bl_per_ant,
                                     eps,
                                     min_ant_times,
                                     trim_neff,
                                     self.id)
        
    def _compare_properties(self,property1,property2,compare_flags=False,compare_weights=False,compare_phase=False):
        """
        compare properties of uvdata sets and calweights. 
        Args:
             property1, CalFlagWeight or UVData object to be compared to 
             property2, CalFlagWeight. 
             compare_flags: directly compare the flags
             compare_weights: directly compare nsample_array 
             compare_phase: if two data sets, make sure they are phased in the same way
        """
        checkpass=True
        if property1.Nfreqs!=property2.Nfreqs:
            checkpass=False
            print("Nfreqs not the same")
        if property1.Njones!=property2.Njones:
            checkpass=False
            print("Njones not the same")
        if property1.Ntimes!=property2.Ntimes:
            checkpass=False
            print("Ntimes not the same")
        if property1.Nspws!=property2.Nspws:
            checkpass=False
            print("Nspws not the same")
        if property1.freq_range!=property2.freq_range:
            checkpass=False
            print("freq_range not the same")            
        if property1.time_range!=property2.time_range:
            checkpass=False
            print("time_range not the same")
        if property1.telescope_name!=property2.telescope_name:
            checkpass=False
            print("telescope_name not the same")
        if property1.Nants_data!=property2.Nants_data:
            checkpass=False
            print("Nants_data not the same")
        if property1.Nants_telescope!=property2.Nants_telescope:
            checkpass=False
            print("Nants_telescope not the same")
        if property1.antenna_names!=property2.antenna_names:
            checkpass=False
            print("antenna_names not the same")
        if property1.antenna_numbers!=property2.antenna_numbers:
            checkpass=False
            print("antenna_numbers not the same")
        if property1.channel_width!=property2.channel_width:
            checkpass=False
            print("channel_width not the same")
        if property1.jones_array!=property2.jones_array:
            checkpass=False
            print("jones_array not the same")
        if property1.time_array!=property2.time_array:
            checkpass=False
            print("time_array not the same")
        if property1.integration_time!=property2.integration_time:
            checkpass=False
            print("integration times not the same")
        if property1.x_orientation!=property2.x_orientation:
            checkpass=False
            print("integration times not the same")
        if compare_flags:
            if property1.flag_array!=property2.flag_array:
                checkpass=False
                print("flag_array not the same")
        if compare_weights:
            if property1.nsample_array!=property2.nsample_array:
                checkpass=False
                print("nsample_array not the same")
        if compare_phase:
            if (property1.phase_type!=property2.phase_type or
                property1.phase_center_ra!=property2.phase_center_ra or
                property1.phase_center_dec!=property2.phase_center_dec):
                checkpass = False
                print("phase centers not the same")
        if compare_ant_positions:
            if property1.antenna_positions!=property2.antenna_positions:
                checkpass=False
                print("antenna positions not the same") 
                
            
        return checkpass
                            
        
    def _check_consistency(self):
        """
        check consistency between data and model uvdata objects
        also check whether the calweights object has the correct
        dimension
        """
        if not(self._compare_properties(self.model_vis,self.measured_vis,
                                        compare_flags=True,
                                        compare_weights=True,
                                        compare-phase=True,
                                        compare_ant_positions=True)):
            raise ValueError("model_vis not consistent with measured_vis")
        
        if(self.compare_properties(self.model_vis,self.flag_weights)):
            raise ValueError("flag weights non consistent with visibilities")
            
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
    def _load_miriad(self,dataname,modelname):
        '''
        read in miriad files
        args: dataname, name of miriad data file
              modelname, name of miriad model file
        '''
        self.model_vis.read_mirad(modelname)
        self.measured_vis.read_miraid(dataname)
        
    def _read_files(self,data,mode,flagweights_fromdata,flagweightsfile=None,model=None):
        '''
        read in all files including data,model,calibration weights
        args: 
           data, data file (fhd,uvfits,measurement,miriad set etc...)
           mode, specify type of data/model files
           model, model file (uvfits, miriad)
        '''
        assert mode in ['MS','UVFITS','MIRIAD','FHD']
        if mode=='MS':
            self.__load_ms(data)
            self.meta_params.model_file=data
        elif mode=='UVFITS':
            assert model
            self.__load_uvfits(data,model)
            self.meta_params.model_file=model
        elif mode=='MIRIAD':
            assert model
            self.__load_miriad(data,model)
            self.meta_params.model_file=model
        elif mode=='FHD':
            self.__load_fhd(data)
        if(flagweights_fromdata):
            self.cal_flag_weights.from_data(self.model_data)
        else:
            assert flagweightsfile
            self.cal_flag_weights.from_file(flagweightsfile)
        self.meta_params.data_file=data
        self._check_consistency()
    def from_ms(self,msfile,flagweights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from a measurement set
        args:
             msfile, name of measurement set cirectory
             flagweights_fromdata, true if you want to determine 
                                   the flagweights object from the measurement set
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(msfile,mode='MS',
                         flagweights_fromdata=flagweights_fromdata,
                         flagweightsfile=flagweightsfile)
    def from_miriad(self,miriaddata,miriadmodel,flagweights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from miriad files
        args:
             miriaddata, name of miriad data file
             miriadmodel, name of mirad model file
             flagweights_fromdata, true if you want to determine 
                                   the flagweights object from the miriad data
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(miriaddata,mode='MIRIAD',
                         flagweights_fromdata=flagweights_fromdata,
                         flagweightsfile=flagweightsfile,
                         model=miriadmodel)
    def from_fhd(self,fhdfile,flagweights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from a fhd
        args:
             fhdfile, name of fhd file
             flagweights_fromdata, true if you want to determine 
                                   the flagweights object from the measurement set
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(fhdfile,mode='FHD',
                         flagweights_fromdata=flagweights_fromdata,
                         flagweightsfile=flagweightsfile)
        
    def from_uvfits(self,uvfitsdata,uvfitsmodel,flagweights_fromdata,flagweightsfile=None):
        """
        initialize stefcal from miriad files
        args:
             uvfitsdata, name of uvfits data file
             uvfitsmodel, name of uvfits model file
             flagweights_fromdata, true if you want to determine 
                                   the flagweights object from the miriad data
             flagweightsfile, initialize flagweights from external file. 
        """
        self._read_files(uvfitsdata,mode='UVFITS',
                         flagweights_fromdata=flagweights_fromdata,
                         flagweightsfile=flagweightsfile,
                         model=uvfitsmodel)
        
    
    def stefcalibrate(self,parallelized=False):
        '''
        Run stefcal
        Args: 
            parallilized, choose this if you want calibration to be parallelized across frequency channels
        '''
        
