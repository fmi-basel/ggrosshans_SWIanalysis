#tests for SWI
import pytest
import os
import sys
import pandas as pd
sys.path.append('../')
import processing.procswi as processing_SWI

#import example lethargus data
DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'example_files/')
@pytest.fixture
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


def test_lethargus_intmolts(lethargus_testcase):
    intmolts, molts = processing_SWI.lethargus_analysis(lethargus_testcase)
    assert intmolts[1,1,12] == 128.0
    assert molts[1,1,12] == 137.0
    assert len(intmolts[1,0,:]) == 16

def test_gfp_adjusting(lethargus_testcase, gfp_testcase):
    gfpdata_original = gfp_testcase
    gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")

    intmolts, molts = processing_SWI.lethargus_analysis(lethargus_testcase)
    dev_length= 231
    gfp_adjusted = processing_SWI.adjust_gfp(gfpdata, intmolts, dev_length)
    assert gfp_adjusted[1][50] == 62.0
    assert len(gfpdata.columns) == 16


