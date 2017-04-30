from pyuvdata import UVData
from pyuvdata import UVCal





class CalFlagWeights():
    """
    Defines a class for storing weights and flags in calibration
    Properties: 
    """
    
    def __init__():
class Stefcal():
    """
    Defines a class for performing stefcal on uvdata sets.
    Attributes:
    _model_vis: uvdata of model
    _measured_vis: uvdata of measurements
    _cal_flag_weights: CalFlagWeights object storing flags and weights
    _cal_solution: uvcal object
    """
    def __init__():
        self._model_vis=UVData()
        self._measured_vis=UVData()
        self._cal_flag_weights=CalFlagWeights()
        self._cal_solution=UVCal()
    def _check_consistency(self):
        """
        check consistency between data and model uvdata objects
        also check whether the calweights object has the correct
        dimension
        """
        
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
        elif mode=='UVFITS':
            assert model
            self.__load_uvfits(data,model)
        elif mode=='MIRIAD':
            assert model
            self.__load_miriad(data,model)
        elif mode=='FHD':
            self.__load_fhd(data)
        if(flagweights_fromdata):
            self.cal_flag_weights.fromdata(self.model_data)
        else:
            assert flagweightsfile
            self.cal_flag_weights.read_file(flagweightsfile)
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
        
        
        
        
        
        
