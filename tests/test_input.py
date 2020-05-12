#tests for SWI
sys.path.append('../')
import pytest
import os
import sys
import pandas as pd
import processing.procswi as processing_SWI

#import example lethargus data
data_path = "C:/Users/hausyann/source/repos/ggrosshans_SWIanalysis/example_files/"

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
    gfpdata_original = pd.read_csv(data_path + "Kymograph_Quantification.csv")
    gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")
    lethargus = pd.read_csv(data_path + "Lethargus_pYPH5_EV.csv")
    intmolts, molts = processing_SWI.lethargus_analysis(lethargus)
    dev_length= 231
    gfp_adjusted = processing_SWI.adjust_gfp(gfpdata, intmolts, dev_length)
    assert gfp_adjusted[1][50] == 62.0
