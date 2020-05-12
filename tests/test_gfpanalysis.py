#test the analysis functions
import pytest
import os
import sys
import pandas as pd
import numpy as np
sys.path.append('../')
import processing.procswi as processing_SWI

#import example lethargus data
data_path = "C:/Users/hausyann/source/repos/ggrosshans_SWIanalysis/example_files/"

lethargus = pd.read_csv(data_path + "Lethargus_pYPH5_EV.csv")
intmolts, molts = processing_SWI.lethargus_analysis(lethargus)

dev_length = 231

gfpdata_original = pd.read_csv(data_path + "Kymograph_Quantification.csv")
gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")
gfp_adjusted = processing_SWI.adjust_gfp(gfpdata, intmolts, dev_length)

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

#set up f_clean for test_run_hilbert(). This can run correctly if test_interpolation runs correctly
f_clean = processing_SWI.interpolate_gfp(gfpdata, intmolts, dev_length)
f_clean_df = pd.DataFrame(f_clean).T
f_clean_df_clean = f_clean_df.loc[:,:] #normally, in the processing, valid worms are selected


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
