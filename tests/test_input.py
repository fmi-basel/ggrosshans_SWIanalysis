#tests for SWI
import pytest
import os
import sys
import pandas as pd
sys.path.append('../')

import processing.procswi as processing_SWI


#import example lethargus data
data_path = "C:/Users/hausyann/source/repos/ggrosshans_SWIanalysis/example_files/"

#test load_gfp:
def load_gfp(uploaded_file_GFP):
    uploaded_file_GFP = st.sidebar.file_uploader("Choose file for GFP data", type="csv")
    if uploaded_file_GFP is not None:
        gfpdata_original = pd.read_csv(uploaded_file_GFP)
        gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")
        return gfpdata


def test_lethargus_intmolts():
    lethargus = pd.read_csv(data_path + "Lethargus_pYPH5_EV.csv")
    intmolts, molts = processing_SWI.lethargus_analysis(lethargus)
    assert intmolts[1,1,12] == 128.0

def test_lethargus_molts():
    lethargus = pd.read_csv(data_path + "Lethargus_pYPH5_EV.csv")

    intmolts, molts = processing_SWI.lethargus_analysis(lethargus)
    assert molts[1,1,12] == 137.0


def test_lethargus_length():
    lethargus = pd.read_csv(data_path + "Lethargus_pYPH5_EV.csv")

    intmolts, molts = processing_SWI.lethargus_analysis(lethargus)
    assert len(intmolts[1,0,:]) == 16

def test_gfp_adjusting():
    gfpdata = pd.read_csv(data_path + "Kymograph_Quantification.csv")
    lethargus = pd.read_csv(data_path + "Lethargus_pYPH5_EV.csv")
    intmolts, molts = processing_SWI.lethargus_analysis(lethargus)
    dev_length = 231

    gfp_adjusted = adjust_gfp(gfpdata, intmolts, dev_length)
    assert 

