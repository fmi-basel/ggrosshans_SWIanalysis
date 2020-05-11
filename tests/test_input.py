#tests for SWI
import pytest
import os
import sys
import pandas as pd
sys.path.append('../')
import processing.procswi as processing_SWI


#import example lethargus data
lethargus_data = pd.read_csv("../Lethargus_pYPH5_EV.csv")
gfp_data = pd.read_csv("../Kymograph_Quantification.csv")

def test_lethargus_intmolts(lethargus_data):
    intmolts, molts = processing_SWI.lethargus_analysis(lethargus_data)
    assert intmolts[1,1,12] == 128.0

def test_lethargus_molts(lethargus_data):
    intmolts, molts = processing_SWI.lethargus_analysis(lethargus_data)
    assert molts[1,1,12] == 137.0
"""
def test_easy():
    a = 5
    b = 3
    assert a + b == 8


def func(x):
    return x + 1

def test_answer():
    assert func(3) == 4
"""
