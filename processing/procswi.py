import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import scipy.stats as sc
from scipy.signal import butter, filtfilt
from scipy.signal import hilbert
from uncertainties import ufloat
import time
import io
import os
import streamlit as st
from PIL import Image


#Lethargus analysis including progress bar
@st.cache(suppress_st_warning=True)
def lethargus_analysis(leth):

    # Add a placeholder
    latest_iteration = st.empty()
    bar = st.progress(0)
    intmolts = np.zeros((2,5,len(leth.columns)-1))
    molts = np.zeros((2,4,len(leth.columns)-1))

    for i in np.arange(1,len(leth.columns)):
        tempintmolts = np.zeros((10,2))
        toggle = 0
        count = -1
        for t in np.arange(1,len(leth)):

            if leth.iloc[t,i] == 1 and toggle == 0:
                toggle = 1
                count = count+1
                tempintmolts[count,0] = t
        
        
            if leth.iloc[t,i] == 0 and toggle == 1:
                toggle = 0
                tempintmolts[count,1] = t-1
                
        tempintmolts[count,1] = t
        lintermolts = pd.DataFrame(np.diff(tempintmolts, axis=1)+1)
        
        sorted_lintermolts = pd.DataFrame(lintermolts.sort_values([0], ascending = False))

        if tempintmolts.shape[0] >4:
            intmolts[0,0:5,i-1] = np.sort(np.take(tempintmolts[:,0],sorted_lintermolts.index.values[0:5]))
            intmolts[1,0:5,i-1] = np.sort(np.take(tempintmolts[:,1],sorted_lintermolts.index.values[0:5]))

        else:
            intmolts[0,0:5,i-1] = np.nan
            intmolts[1,0:5,i-1] = np.nan  
        
        
        for ii in np.arange(0,4):
            molts[0, ii, i-1] = intmolts[1,ii,i-1]+1
            molts[1, ii, i-1] = intmolts[0,ii+1,i-1]-1

        # Update the progress bar with each iteration.
        latest_iteration.text(f'Worm number {i}')
        bar.progress(int(np.round(100*i/(len(leth.columns)-1),0)))
        time.sleep(0.1)
    return intmolts, molts


#show image in app
@st.cache()
def show_worm_image():
    image = Image.open(os.path.dirname(__file__) + "/SWI_chambers_OP50_202020215.png")
    return image


#load raw data
@st.cache()
def load_lethargus(uploaded_file_leth):
    uploaded_file_leth = st.sidebar.file_uploader("Choose file for lethargus data", type="csv")
    if uploaded_file_leth is not None:
        leth = pd.read_csv(uploaded_file_leth)
        return leth

@st.cache()
def load_gfp(uploaded_file_GFP):
    uploaded_file_GFP = st.sidebar.file_uploader("Choose file for GFP data", type="csv")
    if uploaded_file_GFP is not None:
        gfpdata_original = pd.read_csv(uploaded_file_GFP)
        gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")
        return gfpdata



#adjust gfpdata
@st.cache()
def adjust_gfp(gfpdata, intmolts, dev_length):
    gfp_adj = []
    for i in np.arange(0,len(gfpdata.columns)): 
        gfp_adj.append(gfpdata.iloc[(int(intmolts[0,0,i])):(int(intmolts[0,0,i])+dev_length),i])
    gfp_adj = np.asarray(gfp_adj)
    return gfp_adj



#interpolate gfp data
@st.cache()
def interpolate_gfp(gfpdata, intmolts, dev_length):
    for i in np.arange(0,len(gfpdata.columns)):
        f = gfpdata.interpolate(method = 'linear', axis =0, limit = 60, limit_direction = 'backward')

    f_clean = []
    for i in np.arange(0, len(gfpdata.columns)):
        f_clean.append(f.iloc[int(intmolts[0,0,i]):int(intmolts[0,0,i]+dev_length),i].values)
    return f_clean

#select valid worms
@st.cache(allow_output_mutation=True)
def set_up_valid_worms(gfpdata):
    valid_worms_1 = np.repeat(1,len(gfpdata.columns))
    return valid_worms_1

#larval stage durations
@st.cache()
def calculate_durations(leth_clean, intmolts_clean, molts_clean):

    L1_int_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        L1_int_wt.append((intmolts_clean[1,0,i] - intmolts_clean[0,0,i])/6)

    L2_int_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        L2_int_wt.append((intmolts_clean[1,1,i] - intmolts_clean[0,1,i])/6)

    L3_int_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        L3_int_wt.append((intmolts_clean[1,2,i] - intmolts_clean[0,2,i])/6)

    L4_int_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        L4_int_wt.append((intmolts_clean[1,3,i] - intmolts_clean[0,3,i])/6)


    #molt durations
    M1_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        M1_wt.append((intmolts_clean[0,1,i] - intmolts_clean[1,0,i])/6)

    M2_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        M2_wt.append((intmolts_clean[0,2,i] - intmolts_clean[1,1,i])/6)

    M3_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        M3_wt.append((intmolts_clean[0,3,i] - intmolts_clean[1,2,i])/6)

    M4_wt = []
    for i in np.arange(0,len(leth_clean.columns)-1):  
        M4_wt.append((intmolts_clean[0,4,i] - intmolts_clean[1,3,i])/6)

    #larval stage durations
    L1_dur_wt = []
    L2_dur_wt = []
    L3_dur_wt = []
    L4_dur_wt = []

    for i in np.arange(0,len(leth_clean.columns)-1):
        L1_dur_wt.append((intmolts_clean[0,1,i]- intmolts_clean[0,0,i])/6)
        L2_dur_wt.append((intmolts_clean[0,2,i] - intmolts_clean[0,1,i])/6)
        L3_dur_wt.append((intmolts_clean[0,3,i] - intmolts_clean[0,2,i])/6)
        L4_dur_wt.append((intmolts_clean[0,4,i] - intmolts_clean[0,3,i])/6)

    larval_stage_dur = pd.DataFrame([L1_dur_wt, L2_dur_wt, L3_dur_wt, L4_dur_wt], index=["L1", "L2", "L3", "L4"]).T.melt()
    molt_dur = pd.DataFrame([M1_wt, M2_wt, M3_wt, M4_wt], index = ["M1", "M2", "M3", "M4"]).T.melt()
    intermolt_dur = pd.DataFrame([L1_int_wt, L2_int_wt, L3_int_wt, L4_int_wt], index = ["IM1", "IM2", "IM3", "IM4"]).T.melt()

    return larval_stage_dur, molt_dur, intermolt_dur





#Scale the data to each individual larval stage according to mean length of larval stage
@st.cache()
def scale_data(intmolts_clean, molts_clean, gfpdata_clean):
    mean_L1 = int(np.round(np.mean(intmolts_clean[0,1,:] - intmolts_clean[0,0,:])))
    mean_L2 = int(np.round(np.mean(intmolts_clean[0,2,:] - intmolts_clean[0,1,:])))
    mean_L3 = int(np.round(np.mean(intmolts_clean[0,3,:] - intmolts_clean[0,2,:])))
    mean_L4 = int(np.round(np.mean(intmolts_clean[0,4,:] - intmolts_clean[0,3,:])))
        
    scaled_L1 = np.zeros((mean_L1, len(gfpdata_clean.columns)))
    scaled_L2 = np.zeros((mean_L2, len(gfpdata_clean.columns)))
    scaled_L3 = np.zeros((mean_L3, len(gfpdata_clean.columns)))
    scaled_L4 = np.zeros((mean_L4, len(gfpdata_clean.columns)))

    for i in np.arange(0,len(gfpdata_clean.columns)):
        scaled_L1[:,i] = np.interp(np.arange(0,mean_L1), np.arange(0,(intmolts_clean[0,1,i] - intmolts_clean[0,0,i])), gfpdata_clean.iloc[int(intmolts_clean[0,0,i]):int(intmolts_clean[0,1,i]),i])
        scaled_L2[:,i] = np.interp(np.arange(0,mean_L2), np.arange(0,(intmolts_clean[0,2,i] - intmolts_clean[0,1,i])), gfpdata_clean.iloc[int(intmolts_clean[0,1,i]):int(intmolts_clean[0,2,i]),i])
        scaled_L3[:,i] = np.interp(np.arange(0,mean_L3), np.arange(0,(intmolts_clean[0,3,i] - intmolts_clean[0,2,i])), gfpdata_clean.iloc[int(intmolts_clean[0,2,i]):int(intmolts_clean[0,3,i]),i])
        scaled_L4[:,i] = np.interp(np.arange(0,mean_L4), np.arange(0,(intmolts_clean[0,4,i] - intmolts_clean[0,3,i])), gfpdata_clean.iloc[int(intmolts_clean[0,3,i]):int(intmolts_clean[0,4,i]),i])
    
    scaled_all = pd.DataFrame({"scaled_L1": scaled_L1, "scaled_L2": scaled_L2, "scaled_L3": scaled_L3, "scaled_L4": scaled_L4})

    new_molt = molts_clean[:,:,:]-intmolts_clean[0,0:4,:]

    new_molt_mean = np.zeros((4,2))

    for i in np.arange(0,4):
        new_molt_mean[i,0] = np.mean(new_molt[0,i,:])
        new_molt_mean[i,1] = np.mean(new_molt[1,i,:])

    new_molt_std = np.zeros((4,2))

    for i in np.arange(0,4):
        new_molt_std[i,0] = np.std(new_molt[0,i,:])
        new_molt_std[i,1] = np.std(new_molt[1,i,:])

    return scaled_all, new_molt_mean, new_molt_std



#Hilbert transform
def butter_bandpass(lowcut, highcut, fs, order=1):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=1):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data, padtype='constant')
    return y

@st.cache()
def run_hilbert(f_clean_df_clean, leth_clean, lowcut, highcut, fs, dev_length):
    data = f_clean_df_clean

    y = []
    for i in np.arange(0,len(f_clean_df_clean.columns)):
        y.append(butter_bandpass_filter((data.iloc[:,i]-np.mean(data.iloc[:,i])), lowcut, highcut, fs, order=1))

    analytic_signal = []
    amplitude_envelope = []
    instantaneous_phase = []
    instantaneous_frequency = []
    inst_phase_deriv = []
    PeriodoverTime = []
    my_phase = []
    for i in np.arange(0,len(f_clean_df_clean.columns)): 
        analytic_signal.append(hilbert(y[i]))
        amplitude_envelope.append(np.abs(analytic_signal[i]))
        instantaneous_phase.append(np.unwrap(np.angle(analytic_signal[i])))
        instantaneous_frequency.append((np.diff(instantaneous_phase[i]) / (2.0*np.pi) * fs))
        inst_phase_deriv.append(np.diff(instantaneous_phase[i]))
        PeriodoverTime.append((2*np.pi)/inst_phase_deriv[i])    
        my_phase.append(np.angle(analytic_signal[i]))

    worm_names = leth_clean.columns

    phase = pd.DataFrame(my_phase).T
    phase.columns = worm_names[1:]
    phase_melt = phase.melt()
    phase_melt.columns = ["Worm", "phase"]
    phase_melt["Timepoint"] = np.tile(np.arange(1,dev_length+1), len(my_phase))
    
    return my_phase, PeriodoverTime, phase_melt


@st.cache()
def phase_molt_entry_exit(gfpdata_clean, intmolts_clean, my_phase):
    molt_entry_phase = []
    molt_exit_phase = []

    for i in np.arange(0,len(gfpdata_clean.columns)):
        for n in np.arange(0,4):
        
            entry_tp = int(intmolts_clean[1,n,i] - intmolts_clean[0,0,i])
            exit_tp = int(intmolts_clean[0,n+1,i] - intmolts_clean[0,0,i])
            molt_entry_phase.append(my_phase[i][entry_tp])
            molt_exit_phase.append(my_phase[i][exit_tp])
    molt_entr_ph_L1 = []
    molt_entr_ph_L2 = []
    molt_entr_ph_L3 = []
    molt_entr_ph_L4 = []

    molt_exit_ph_L1 = []
    molt_exit_ph_L2 = []
    molt_exit_ph_L3 = []
    molt_exit_ph_L4 = []

    for i in np.arange(0,len(gfpdata_clean.columns)):
        molt_entr_ph_L1.append(molt_entry_phase[4*i]) 
        molt_entr_ph_L2.append(molt_entry_phase[4*i+1])
        molt_entr_ph_L3.append(molt_entry_phase[4*i+2])
        molt_entr_ph_L4.append(molt_entry_phase[4*i+3])
        
        molt_exit_ph_L1.append(molt_exit_phase[4*i])
        molt_exit_ph_L2.append(molt_exit_phase[4*i+1])
        molt_exit_ph_L3.append(molt_exit_phase[4*i+2])
        molt_exit_ph_L4.append(molt_exit_phase[4*i+3])

    return molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4

    #correct (switch) phase 
@st.cache()
def correct_phase(molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4):
    corr_molt_entr_ph_L1 = [] #the following code switched the phase in case some data points run over 2pi
    for i in np.arange(0,len(molt_entr_ph_L1)):
        if np.sign(np.median(molt_entr_ph_L1)) != np.sign(molt_entr_ph_L1[i]):
            corr_molt_entr_ph_L1.append(molt_entr_ph_L1[i]+2*np.pi) 
        else:
            corr_molt_entr_ph_L1.append(molt_entr_ph_L1[i])

    corr_molt_exit_ph_L1 = []
    for i in np.arange(0,len(molt_exit_ph_L1)):
        if np.sign(np.median(molt_exit_ph_L1)) != np.sign(molt_exit_ph_L1[i]):
            corr_molt_exit_ph_L1.append(molt_exit_ph_L1[i]+2*np.pi)
        else:
            corr_molt_exit_ph_L1.append(molt_exit_ph_L1[i])

    corr_molt_entr_ph_L2 = []
    for i in np.arange(0,len(molt_entr_ph_L2)):
        if np.sign(np.median(molt_entr_ph_L2)) != np.sign(molt_entr_ph_L2[i]):
            corr_molt_entr_ph_L2.append(molt_entr_ph_L2[i]+2*np.pi)
        else:
            corr_molt_entr_ph_L2.append(molt_entr_ph_L2[i])
    
    corr_molt_exit_ph_L2 = []
    for i in np.arange(0,len(molt_exit_ph_L2)):
        if np.sign(np.median(molt_exit_ph_L2)) != np.sign(molt_exit_ph_L2[i]):
            corr_molt_exit_ph_L2.append(molt_exit_ph_L2[i]+2*np.pi)
        else:
            corr_molt_exit_ph_L2.append(molt_exit_ph_L2[i])

    corr_molt_entr_ph_L3 = [] #the following code switched the phase in case some data points run over 2pi
    for i in np.arange(0,len(molt_entr_ph_L3)):
        if np.sign(np.median(molt_entr_ph_L3)) != np.sign(molt_entr_ph_L3[i]):
            corr_molt_entr_ph_L3.append(molt_entr_ph_L3[i]+2*np.pi) 
        else:
            corr_molt_entr_ph_L3.append(molt_entr_ph_L3[i])

    corr_molt_exit_ph_L3 = []
    for i in np.arange(0,len(molt_exit_ph_L3)):
        if np.sign(np.median(molt_exit_ph_L3)) != np.sign(molt_exit_ph_L3[i]):
            corr_molt_exit_ph_L3.append(molt_exit_ph_L3[i]+2*np.pi)
        else:
            corr_molt_exit_ph_L3.append(molt_exit_ph_L3[i])

    corr_molt_entr_ph_L4 = []
    for i in np.arange(0,len(molt_entr_ph_L4)):
        if np.sign(np.median(molt_entr_ph_L4)) != np.sign(molt_entr_ph_L4[i]):
            corr_molt_entr_ph_L4.append(molt_entr_ph_L4[i]+2*np.pi)
        else:
            corr_molt_entr_ph_L4.append(molt_entr_ph_L4[i])
    
    corr_molt_exit_ph_L4 = []
    for i in np.arange(0,len(molt_exit_ph_L4)):
        if np.sign(np.median(molt_exit_ph_L4)) != np.sign(molt_exit_ph_L4[i]):
            corr_molt_exit_ph_L4.append(molt_exit_ph_L4[i]+2*np.pi)
        else:
            corr_molt_exit_ph_L4.append(molt_exit_ph_L4[i])

    return corr_molt_entr_ph_L1, corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4

@st.cache()
def sort_lethargus(leth_clean, intmolts_clean):
    leth_corr = []
    end = []
    for i in np.arange(1,len(leth_clean.columns)):
        end.append(len(leth_clean)-int(intmolts_clean[0,0,i-1]))
    max_len = np.min(end)

    for i in np.arange(1,len(leth_clean.columns)):
        leth_corr.append(leth_clean.iloc[int(intmolts_clean[0,0,i-1]):(int(intmolts_clean[0,0,i-1])+max_len),i])

    leth_new = np.asarray(leth_corr)

    return leth_new[np.argsort(intmolts_clean[1,0,:]-intmolts_clean[0,0,:]),:]

@st.cache()
def period_per_LS(PeriodoverTime, gfpdata_clean, intmolts_clean):
    period_L2 = []
    period_L3 = []
    period_L4 = []

    for i in np.arange(0,len(gfpdata_clean.columns)):
            entry_tp_2 = int(intmolts_clean[0,1,i] - intmolts_clean[0,0,i])
            exit_tp_2 = int(intmolts_clean[0,2,i] - intmolts_clean[0,0,i])
            entry_tp_3 = int(intmolts_clean[0,2,i] - intmolts_clean[0,0,i])
            exit_tp_3 = int(intmolts_clean[0,3,i] - intmolts_clean[0,0,i])
            entry_tp_4 = int(intmolts_clean[0,3,i] - intmolts_clean[0,0,i])
            exit_tp_4 = int(intmolts_clean[0,4,i] - intmolts_clean[0,0,i])
            period_L2.append(np.nanmean(PeriodoverTime[i][entry_tp_2:exit_tp_2])/6)
            period_L3.append(np.nanmean(PeriodoverTime[i][entry_tp_3:exit_tp_3])/6)
            period_L4.append(np.nanmean(PeriodoverTime[i][entry_tp_4:exit_tp_4])/6)

    return period_L2, period_L3, period_L4

@st.cache(allow_output_mutation=True)
def error_prop(L2_dur_wt, L3_dur_wt, L4_dur_wt, L2_int_wt, L3_int_wt, L4_int_wt, period_L2, period_L3, period_L4):
    MEAN_periods_and_LS = [[np.mean(L2_dur_wt), np.mean(L3_dur_wt), np.mean(L4_dur_wt)],[np.mean(period_L2), np.mean(period_L3),  np.mean(period_L4)], [np.mean(L2_int_wt), np.mean(L3_int_wt), np.mean(L4_int_wt)]]
    STD_periods_and_LS = [[np.std(L2_dur_wt), np.std(L3_dur_wt),np.std(L4_dur_wt)], [np.std(period_L2), np.std(period_L3),  np.std(period_L4)], [np.std(L2_int_wt), np.std(L3_int_wt), np.std(L4_int_wt)]]
    prop_err_exit = []
    prop_err_entry = []
    for i in np.arange(0,3):
        LS = ufloat(MEAN_periods_and_LS[0][i], STD_periods_and_LS[0][i])
        per = ufloat(MEAN_periods_and_LS[1][i], STD_periods_and_LS[1][i])
        IM = ufloat(MEAN_periods_and_LS[2][i], STD_periods_and_LS[2][i])
        prop_err_exit.append((2*np.pi)/per*LS)
        prop_err_entry.append((2*np.pi)/per*IM)

    return prop_err_entry, prop_err_exit

@st.cache()
def use_phases(molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4):
    molt_ph_entry = pd.DataFrame([molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
    molt_ph_entry["Molt"] = "entry"
    
    molt_ph_exit = pd.DataFrame([molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
    molt_ph_exit["Molt"] = "exit"

    molt_phase_both = molt_ph_entry.append(molt_ph_exit)
    return molt_phase_both

@st.cache()
def use_corrected_phases(corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4):
    molt_ph_entry = pd.DataFrame([corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
    molt_ph_entry["Molt"] = "entry"
    
    molt_ph_exit = pd.DataFrame([corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
    molt_ph_exit["Molt"] = "exit"

    molt_phase_both_corr = molt_ph_entry.append(molt_ph_exit)
    return molt_phase_both_corr

@st.cache()
def combine_stds(molt_ph_entry, prop_err_entry, molt_ph_exit, prop_err_exit):
    std_phases_entry = molt_ph_entry[["variable", "value"]].groupby("variable").std().iloc[1:3,:]
    error_prop_std_entry = []
    for i in np.arange(0,2):
        error_prop_std_entry.append(prop_err_entry[i].std_dev)

    std_phases_entry["error_prop_std"] = error_prop_std_entry
    std_phases_entry["ratio sd_obs/sd_exp"] = std_phases_entry["value"]/std_phases_entry["error_prop_std"]
    std_phases_entry["entry_or_exit"] = "entry"

    std_phases_exit = molt_ph_exit[["variable", "value"]].groupby("variable").std().iloc[1:3,:]
    error_prop_std_exit = []
    for i in np.arange(0,2):
        error_prop_std_exit.append(prop_err_exit[i].std_dev)

    std_phases_exit["error_prop_std"] = error_prop_std_exit
    std_phases_exit["ratio sd_obs/sd_exp"] = std_phases_exit["value"]/std_phases_exit["error_prop_std"]
    std_phases_exit["entry_or_exit"] = "exit"

    std_phases_all = std_phases_exit.append(std_phases_entry)
    
    std_phases_all["observation"] = std_phases_all.index
    return std_phases_all











