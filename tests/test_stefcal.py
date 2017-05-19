"""Tests for stefcal uvdata"""
import nose.tools as nt
import numpy as np
import os
import stefcal_uvdata as suv
from pyuvdata import UVCal
from stefcal_uvdata import StefcalUVData
from data import DATA_PATH

def test_read_write():
    """
    Test stefcal on noiseless data set with unity gains
    and a perfect sky model. 
    """
    stefcal_perfect=StefcalUVData()
    measured_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    model_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    #generate_flag_weights file from data
    stefcal_perfect.from_uvfits(measured_vis,model_vis,
                                flag_weights_fromdata=True)
    #set parameters
    stefcal_perfect.set_params({'trim_neff':False,
                                'min_ant_times':1,
                                'eps':1e-10,
                                'min_bl_per_ant':2,
                                'n_cycles':1,
                                'n_phase_iter':5,
                                'refant':0})
    #test writing out flag_weights_file
    print('writing')
    stefcal_perfect.save_state(DATA_PATH+'/test/perfect_cal',clobber=True)
    print('writing successful')
    print('reading')
    stefcal_perfect.load_state(DATA_PATH+'/test/perfect_cal')
    print('reading succesful')

def test_matrix_reordering():
    """
    Test whether method that converts from matrix to baseline-times 
    reverses correctly
    """
    stefcal_test=StefcalUVData()
    measured_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    model_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    stefcal_test.from_uvfits(measured_vis,model_vis,
                             flag_weights_fromdata=True)
    time_array_test=stefcal_test._matrix_2_blt_list(stefcal_test._blt_list_2_matrix(stefcal_test.measured_vis.time_array,range(stefcal_test.measured_vis.Ntimes)))
    nt.assert_true(np.all(time_array_test==stefcal_test.measured_vis.time_array))


def test_perfect_calibration_random_gains():
    """
    test stefcal on noiseless data set with gains that have
    reflections and perfect sky model. Difference between 
    true and estimated gains should be small.
    """
    eps=1e-10
    stefcal_perfect=StefcalUVData()
    measured_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    model_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    stefcal_perfect.from_uvfits(measured_vis,model_vis,
                                flag_weights_fromdata=True)
    #set params
    stefcal_perfect.set_params({'trim_neff':False,
                                'min_ant_times':1,
                                'eps':eps,
                                'min_bl_per_ant':2,
                                'n_cycles':1,
                                'n_phase_iter':0,
                                'refant':0})
    gshape=stefcal_perfect.uvcal.gain_array.shape
    gain_array=2.+np.random.randn(*gshape)+1j*np.random.randn(*gshape)
    gains_random=stefcal_perfect.uvcal_from_data()
    gains_random.gain_array=1./gain_array
    stefcal_perfect.measured_vis=suv.correct_vis(stefcal_perfect.measured_vis,gains_random)
    #calibrate
    stefcal_perfect.stefcalibrate(perterb=1e-2)
    #compare calibration solution with gain array
    nt.assert_true(np.max(np.abs(stefcal_perfect.uvcal.gain_array-gain_array))<=10.*eps)
    
    
    
def test_perfect_calibration_unity_gains():
    """
    Test stefcal on noiseless data set with unity gains
    and a perfect sky model. 
    """
    eps=1e-10
    stefcal_perfect=StefcalUVData()
    measured_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    model_vis=os.path.join(DATA_PATH,'point_source_sim_model.uvfits')
    #generate_flag_weights file from data
    stefcal_perfect.from_uvfits(measured_vis,model_vis,
                                flag_weights_fromdata=True)
    #set parameters
    stefcal_perfect.set_params({'trim_neff':False,
                                'min_ant_times':1,
                                'eps':epse,
                                'min_bl_per_ant':2,
                                'n_cycles':1,
                                'n_phase_iter':0,
                                'refant':0})    
    stefcal_perfect.stefcalibrate(perterb=1e-2)
    #verify that solutions are no greater than eps
    

    
    #test that the number of degress of freedom
    #per antenna is equal to the number of antennas - 1
    #Chi square should be infinite    
#test_matrix_reordering()
#test_perfect_calibration()
test_perfect_calibration_random_gains()    
    
