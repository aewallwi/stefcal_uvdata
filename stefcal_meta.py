import pyuvdata.version as stfversion
import pickle
import numpy as np
import copy
class StefcalMeta():
    """
    Defines a container class for stefcal meta parameters 
    contains stefcal parameters, chi-squares, noise
    calculations along with filepaths
    to data, model, and flag_weights files. 
    """
    def __init__()#self,id,refant=0,n_phase_iter=5,
                 #n_cycles=1,min_bl_per_ant=2,eps=1e-10,
                 #min_ant_times=1,trim_neff=False,spw=0,
                 #t_avg=1):


        desc='Reference antenna'
        self._refant=uvp.UVParameter('refant',description=desc,
                                     expected_type=int)
        
        desc='Number of iterations in which only the phase is fitted.'
        self._n_phase_iter=uvp.UVParameter('n_phase_iter',description=desc,
                                           expected_type=int)
        
        desc='minimum baselines per antenna for it to remain unflagged'
        self._min_bl_per_ant=uvp.UVParameter('min_bl_per_ant',description=desc,
                                             expected_type=int)
        
        desc='stopping criterion for stefcal'
        self._eps=uvp.UVParameter('eps',description=desc,
                                  expected_type=str)
        
        desc='minimum unflagged times or an antenna to not flag it'
        self._min_ant_times=uvp.UVParameter('min_ant_times',description=desc,
                                            expected_type=int)
        
        desc='flag whether neff trimming was employed'
        self._trim_neff=uvp.UVParameter('trim_neff',description=desc,
                                        expected_type=bool)
        
        desc='unique identifier linking meta-data to calibration solutions and'
             'flag weights file'
        self._id=uvp.UVParameter('id',description=desc,
                                 expected_type=str)
        
        desc='file of flags and weights'
        self._flag_weights_file=uvp.UVParameter('flag_weights_file',description=desc,
                                                expected_type=str)

        desc='file of model visibilities',
        self._model_file=uvp.UVParameter('model_file',
                                         description=desc,
                                         expected_type=str)

        desc='file for measured visibilities'
        self._measurement_file=uvp.UVParameter('measurement_file',
                                               description=desc,
                                               expected_type=str)

        self._spw=uvp.UVParameter('spw',description='SPW calibrated',
                                  expected_type=int)

        desc='Number of time steps to average calibration over'
        self._t_avg=uvp.UVParameter('t_avg',
                                    description=desc,
                                    expected_type=int)

        self._n_ycles=uvp.UVParameter('n_cycles',
                                      description='Number of cycles of stefcal performed',
                                      expected_type=int)
        
        self._Niterations=uvp.UVParameter('iterations',
                                         description='Number of iterations of stefcal performed',
                                         form=('Ncycles','Nfreqs','Njones')
                                         expected_type=int)
        
        self._Nfreqs = uvp.UVParameter('Nfreqs',
                                       description='Number of frequency channels',
                                       expected_type=int)
        
                                        expected_type=int)
        self._Ntimes = uvp.UVParameter('Ntimes',
                                       description='Number of times',
                                       expected_type=int)

        self._Nfreqs = uvp.UVParameter('Nfreqs',
                                       description='Number of frequency channels',
                                       expected_type=int)
        desc = ('Number of antennas with data present (i.e. number of unique '
                'entries in ant_array). May be smaller ' +
                'than the number of antennas in the telescope')
        self._Nants_data = uvp.UVParameter('Nants_data', description=desc,
                                           expected_type=int)

        desc = ('Number of antennas in the array. May be larger ' +
                'than the number of antennas with data')
        self._Nants_telescope = uvp.UVParameter('Nants_telescope',
                                                description=desc,
                                                expected_type=int)

        
        desc="Chi-Squares for each antenna gain solution."
        
        self._chi_square_per_ant=uvp.UVParameter('chi_square_per_ant',description=desc,
                                                 form=('Nants_data','Nfreqs',
                                                       'Ntimes','Njones'),
                                                 expected_type=float,required=False)
        
        desc='keeps track of degrees of freedom per gain solution.'
        self._dof_per_ant=uvp.UVParameter('dof_per_ant',description=desc,
                                          form=('Nants_data','Nfreqs',
                                                'Ntimes','Njones'),
                                          expected_type=float)
                    
        
        #self._noise_tavg=uvp.UVParameter('noise_tavg',
        #                                 description='noise levels in uncalibrated'
        #                                 'visibilities computed by taking differences'
        #                                 'in frequency and median over'
        #                                 ' all times',
        #                                 form=('Nbls','Nfreqs','Njones'),
        #                                 expected_type=np.float,required=False)
        
        #self._noise_favg=uvp.UVParameter('noise_favg',
        #                                 description='noise levels in uncalibrated'
        #                                 'visibilities computed by taking differences'
        #                                 'in time and median over'
        #                                 ' all frequency',
        #                                 form=('Nblts','Njones'),
        #                                 expected_type=np.float,required=False)

        

        self._noise_tblavg=uvp.UVParameter('noise_tblavg',
                                           description='noise levels in uncalibrated'
                                           'visibilities computed by taking differences'
                                           'in time and and taking median over frequency'
                                           'baseline, polarization, and time',
                                          form=('Nfreqs'),
                                          expected_type=np.float,required=False)
        self.stefcal_version_str+=(' Git origin: ' + stfversion.git.origin +
                                   '. Git hash: ' + stfversion.get_hash +
                                   '. Git branch: ' + stfversion.git_branch+
                                   '. Git description: ' + stfversion.get_description)

        def to_file(self,datafile,clobber=False):
            """
            write meta object to a file using pickle
            """
            if clobber or not(os.path.exists(datafile)):
                pickle.dump(self,open(datafile,"wb"))
            else:
                print("%s already exists. Use clobber=True to overwrite"%datafile)
                
        def from_file(self,datafile):
            """
            read meta object from pickled file
            """
            data=pickle.load(open(datafile,"rb"))
            self.refant=copy.copy(data.refant)
            self.n_phase_iter=copy.copy(data.n_phase_iter)
            self.min_bl_per_ant=copy.copy(data.min_bl_per_ant)
            self.eps=copy.copy(data.eps)
            self.min_ant_times=copy.copy(data.min_ant_times)
            self.trim_neff=copy.copy(data.trim_neff)
            self.id=copy.copy(data.id)
            self.flag_weights_file=copy.copy(data.flag_weights_file)
            self.model_file=copy.copy(data.model_file)
            self.measurement_file=copy.copy(data.measurement_file)
            self.model_file=copy.copy(data.model_file)
            self.spw=copy.copy(data.spw)
            self.t_avg=copy.copy(data.t_avg)
            self.n_cycles=copy.copy(data.n_cycles)
            self.Niterations=copy.copy(data.Niterations)
            self.Nfreqs=copy.copy(data.Nfreqs)
            self.Ntimes=copy.copy(data.Ntimes)
            self.Nfreqs=copy.copy(data.Nfreqs)
            self.Nants_data=copy.copy(data.Nants_data)
            self.Nants_telescope=copy.copy(data.Nants_telescope)
            self.chi_square_per_ant=copy.copy(data.chi_square_per_ant)
            self.dof_per_ant=copy.copy(data.dof_per_ant)
            self.noise_tblavg=copy.copy(data.noise_tblavg)
            self.stefcal_version_str=copy.copy(data.stefcal_version_str)
            del(data)
