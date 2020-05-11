#tests for SWI
import pytest
import os
import sys
import pandas as pd
from processing.procswi import *


#import example lethargus data
lethargus_data = pd.read_csv("../example_files/Lethargus_pYPH5_EV.csv")



class Test_class:
    def test_lethargus_input(lethargus_data):
        


