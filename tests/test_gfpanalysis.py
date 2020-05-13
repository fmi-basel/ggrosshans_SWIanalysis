#test the analysis functions
import pytest
import os
import sys
import pandas as pd
import numpy as np
sys.path.append('../')
import processing.procswi as processing_SWI

#import example lethargus data
DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'example_files/')

def lethargus_testcase():
    '''Utility to load the test case. Could be extended to load
    not only one but several test cases.
    '''
    return pd.read_csv(DATA_PATH + "Lethargus_pYPH5_EV.csv")
  
@pytest.fixture
def gfp_testcase():
    '''Utility to load the test case. Could be extended to load
    not only one but several test cases.
    '''
    return pd.read_csv(DATA_PATH + "Kymograph_Quantification.csv")


lethargus = lethargus_testcase()
intmolts, molts = processing_SWI.lethargus_analysis(lethargus)

dev_length = 231

gfpdata_original = pd.read_csv(DATA_PATH + "Kymograph_Quantification.csv")
gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")
gfp_adjusted = processing_SWI.adjust_gfp(gfpdata, intmolts, dev_length)

#set up f_clean for test_run_hilbert(). This can run correctly if test_interpolation runs correctly
f_clean = processing_SWI.interpolate_gfp(gfpdata, intmolts, dev_length)
f_clean_df = pd.DataFrame(f_clean).T
f_clean_df_clean = f_clean_df.loc[:,:] #normally, in the processing, valid worms are selected
#tests
def test_interpolation():
    f_clean = processing_SWI.interpolate_gfp(gfpdata, intmolts, dev_length)
    assert f_clean[3][20] == 192.0

def test_setup_valid_worms():
    valid_worms = processing_SWI.set_up_valid_worms(gfpdata)
    assert len(valid_worms) == 16

def test_calculate_durations():
    larval_stage_dur, molt_dur, intermolt_dur = processing_SWI.calculate_durations(lethargus, 
                                                                            intmolts, molts)
    assert np.sum(molt_dur["variable"].unique() == ['M1', 'M2', 'M3', 'M4']) == 4




def test_run_hilbert():
    lowcut = 1/14 #default in procswi
    highcut = 1/5 #default in procswi
    fs = 6
    my_phase, PeriodoverTime, phase_melt = processing_SWI.run_hilbert(f_clean_df_clean, 
                                                        lethargus, lowcut, 
                                                        highcut, 
                                                        fs, 
                                                        dev_length)
    assert np.max(my_phase) <= np.pi and np.min(my_phase) >= -np.pi
    assert len(my_phase) 
