import numpy as np
import pyuvdata.parameter as uvp
import pyuvdata.version as stfversion
#Need to change this to a version.py on this repository!
import copy
import pickle
import os

class CalFlagWeights():
    """
    Defines a class for storing weights and flags in calibration
    Properties: 
    """
    
    def __init__(self,id):
        self.id=id #unique identifier. 
        radian_tol = 10 * 2 * np.pi * 1e-3 / (60.0 * 60.0 * 360.0)
        desc=('Array of boolean flags of data that was flagged for calibration'
             ' but not necessarily elsewhere in data reduction')
        self._flag_array=uvp.UVParameter('flag_array',description=desc,
                                         form=('Nblts','Nspws','Nfreqs','Npols'),
                                         expected_type=np.bool)
        desc=('Array of weights applied to data only in calibration.'
              'shape: (Nblts,Nspws,Nfreqs,Npols)')
        self._weights_array=uvp.UVParameter('weights_array',description=desc,
                                           form=('Nblts','Nspws','Nfreqs','Npols'),
                                            expected_type=np.float)        
        desc='unique string identifier associating calflagweights with cal file and meta file'
        self._id=uvp.UVParameter('id',description=desc,
                                 expected_type=str)
            
    def from_file(self,data,mode='CFWS'):
        '''
        initialize flagweights object from calFlagWeights file
        args: name of data file to read in
        '''
        assert mode in ['UVDATA','CFWS']
        if mode=='CFWS':
            data=pickle.load(open(data,"rb"))
            self.weights_array=copy.copy(data.weights_array)
        else:
            self.weights_array=copy.copy(data.nsample_array).astype(np.float64)
        self.flag_array=copy.copy(data.flag_array)
        del(data)
        


    def to_file(self,datafile,clobber=False):
        '''
        write calweights object to a file using pickle
        '''
        if clobber or not(os.path.exists(datafile)):
            pickle.dump(self,open(datafile,"wb"))
        else:
            print("%s already exists. Use clobber=True to overwrite"%datafile)
            

            
