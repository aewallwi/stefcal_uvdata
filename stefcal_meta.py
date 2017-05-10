
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

        self._refant=uvp.UVParameter('refant',
                                     description='Reference antenna number',
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
        
        self._chi_squares=uvp.UVParameter('chi_squares',description=desc,
                                          form=('Nants_data','Nfreqs',
                                                'Ntimes','Njones'),
                                          expected_type=float,required=False)
        
        self._noise_tavg=uvp.UVParameter('noise_tavg',
                                         description='noise levels in uncalibrated'
                                         'visibilities computed by taking differences'
                                         'in frequency and median over'
                                         ' all times',
                                         form=('Nbls','Nfreqs','Njones'),
                                         expected_type=np.float,required=False)
        
        self._noise_favg=uvp.UVParameter('noise_favg',
                                         description='noise levels in uncalibrated'
                                         'visibilities computed by taking differences'
                                         'in time and median over'
                                         ' all frequency',
                                         form=('Nblts','Njones'),
                                         expected_type=np.float,required=False)

        

        self._noise_tblavg=uvp.UVParameter('noise_tblavg',
                                           description='noise levels in uncalibrated'
                                           'visibilities computed by taking differences'
                                           'in time and and taking median over frequency'
                                           'baseline, polarization, and time',
                                          form=('Nfreqs'),
                                          expected_type=np.float,required=False)
        
