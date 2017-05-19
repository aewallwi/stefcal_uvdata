"""Tests for stefcal uvdata"""
import nose.tools as nt
import os
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


def test_perfect_calibration():
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
                                'n_phase_iter':0,
                                'refant':0})
    
    stefcal_perfect.stefcalibrate(perterb=1e-2)
    print('stef-calibrated with %s iterations'%stefcal_perfect.meta_params.Niterations)
    #test that the number of degress of freedom
    #per antenna is equal to the number of antennas - 1

    #Chi square should be infinite
    

test_perfect_calibration()
    
    
