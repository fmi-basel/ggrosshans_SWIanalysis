# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 14:30:58 2018

@author: hausyann

This is a GUI for single worm imaging data.

"""

import tkinter as tk
from tkinter.ttk import Button
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import imp
import scipy.stats as sc
from scipy.signal import butter, filtfilt
from scipy.signal import hilbert
from uncertainties import ufloat






LARGE_FONT = ("Verdana", 12)
dt = 0.16666666 # sampling interval in hours


#dev_length = 231 #in time points, has to be adjusted for individual experiments


class SWI_analysis(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)
        container.pack(side="top", fill = "both", expand = True)
        
        
        self.frames = {}
        
        frame = StartPage(container, self)
        
        self.frames[StartPage] =    frame
        
        frame.grid(row = 0, column = 0, sticky = "nsew")
        
        self.show_frame(StartPage)
        
        
    def show_frame(self, cont):
        
        frame = self.frames[cont]
        frame.tkraise()
        
        
#Define functions here:

def leth_ana():

    root1 = tk.Tk()
    root1.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
    global leth
    leth = pd.DataFrame(pd.read_csv(root1.filename))

    global intmolts
    intmolts = np.zeros((2,5,len(leth.columns)-1))
    global molts
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

    root1.geometry("100x100")
    messagebox.showinfo("Lethargus analysis status", "Lethargus analysis done")
    root1.destroy()

def gfp_ana(self):

    root1 = tk.Tk()
    root1.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
    leth = pd.DataFrame(pd.read_csv(root1.filename))

    root2 = tk.Tk()
    root2.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))

    gfpdata = pd.read_csv(root2.filename)
    #gfpdata = gfpdata.pivot("Frame","Position","Mean (Statistic Features)")        
    gfpdata = gfpdata.pivot("Frame","Position","Intensity_BGsub")        

    #add user input for dev_length input
    root3 = tk.Tk()
    root3.title("Input developmental length")
    canvas1 = tk.Canvas(root3, width = 400, height = 300)
    canvas1.pack()
    entry1 = tk.Entry(root3) 
    canvas1.create_window(200, 140, window=entry1)

    dev_length = entry1.get()
    

    for i in np.arange(0,len(gfpdata.columns)):
        f = gfpdata.interpolate(method = 'linear', axis =0, limit = 60, limit_direction = 'backward')
  
    f_clean = []
    for i in np.arange(0, len(gfpdata.columns)):
        f_clean.append(f.iloc[int(intmolts[0,0,i]):int(intmolts[0,0,i]+dev_length),i].values)
    

        
    #Scale the data to each individual larval stage according to mean length of larval stage
    mean_L1 = int(np.round(np.mean(intmolts[0,1,:] - intmolts[0,0,:])))
    mean_L2 = int(np.round(np.mean(intmolts[0,2,:] - intmolts[0,1,:])))
    mean_L3 = int(np.round(np.mean(intmolts[0,3,:] - intmolts[0,2,:])))
    mean_L4 = int(np.round(np.mean(intmolts[0,4,:] - intmolts[0,3,:])))
        
    scaled_L1 = np.zeros((mean_L1, len(gfpdata.columns)))
    scaled_L2 = np.zeros((mean_L2, len(gfpdata.columns)))
    scaled_L3 = np.zeros((mean_L3, len(gfpdata.columns)))
    scaled_L4 = np.zeros((mean_L4, len(gfpdata.columns)))

    for i in np.arange(0,len(gfpdata.columns)):
        scaled_L1[:,i] = np.interp(np.arange(0,mean_L1), np.arange(0,(intmolts[0,1,i] - intmolts[0,0,i])), gfpdata.iloc[int(intmolts[0,0,i]):int(intmolts[0,1,i]),i])
        scaled_L2[:,i] = np.interp(np.arange(0,mean_L2), np.arange(0,(intmolts[0,2,i] - intmolts[0,1,i])), gfpdata.iloc[int(intmolts[0,1,i]):int(intmolts[0,2,i]),i])
        scaled_L3[:,i] = np.interp(np.arange(0,mean_L3), np.arange(0,(intmolts[0,3,i] - intmolts[0,2,i])), gfpdata.iloc[int(intmolts[0,2,i]):int(intmolts[0,3,i]),i])
        scaled_L4[:,i] = np.interp(np.arange(0,mean_L4), np.arange(0,(intmolts[0,4,i] - intmolts[0,3,i])), gfpdata.iloc[int(intmolts[0,3,i]):int(intmolts[0,4,i]),i])


    new_molt = molts[:,:,:]-intmolts[0,0:4,:]

    new_molt_mean = np.zeros((4,2))

    for i in np.arange(0,4):
        new_molt_mean[i,0] = np.mean(new_molt[0,i,:])
        new_molt_mean[i,1] = np.mean(new_molt[1,i,:])

    new_molt_std = np.zeros((4,2))

    for i in np.arange(0,4):
        new_molt_std[i,0] = np.std(new_molt[0,i,:])
        new_molt_std[i,1] = np.std(new_molt[1,i,:])


    order = 1
    def butter_bandpass(lowcut, highcut, fs, order=order):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a


    def butter_bandpass_filter(data, lowcut, highcut, fs, order=order):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = filtfilt(b, a, data, padtype='constant')
        return y


    # Sample rate and desired cutoff frequencies (in 1/h).
    fs = 6
    lowcut = 1/14 #14hour period
    highcut = 1/5 #5hour period
    
    data = f_clean

    y = []
    for i in np.arange(0,len(f_clean)):
        y.append(butter_bandpass_filter((data[i]-np.mean(data[i])), lowcut, highcut, fs, order=order))

    analytic_signal = []
    amplitude_envelope = []
    instantaneous_phase = []
    instantaneous_frequency = []
    inst_phase_deriv = []
    PeriodoverTime = []
    my_phase = []
    for i in np.arange(0,len(f_clean)): 
        analytic_signal.append(hilbert(y[i]))
        amplitude_envelope.append(np.abs(analytic_signal[i]))
        instantaneous_phase.append(np.unwrap(np.angle(analytic_signal[i])))
        instantaneous_frequency.append((np.diff(instantaneous_phase[i]) / (2.0*np.pi) * fs))
        inst_phase_deriv.append(np.diff(instantaneous_phase[i]))
        PeriodoverTime.append((2*np.pi)/inst_phase_deriv[i])    
        my_phase.append(np.angle(analytic_signal[i]))


    
    molt_entry_phase = []
    molt_exit_phase = []

    for i in np.arange(0,len(gfpdata.columns)):
        for n in np.arange(0,4):
        
            entry_tp = int(intmolts[1,n,i] - intmolts[0,0,i])
            exit_tp = int(intmolts[0,n+1,i] - intmolts[0,0,i])
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

    for i in np.arange(0,len(gfpdata.columns)):
        molt_entr_ph_L1.append(molt_entry_phase[4*i]) #select M1 molt entry phase, +np.pi because wavelets start at -pi (=0 degree) over 0 (=180degree) to pi (=360 degree)
        molt_entr_ph_L2.append(molt_entry_phase[4*i+1])
        molt_entr_ph_L3.append(molt_entry_phase[4*i+2])
        molt_entr_ph_L4.append(molt_entry_phase[4*i+3])
        
        molt_exit_ph_L1.append(molt_exit_phase[4*i])
        molt_exit_ph_L2.append(molt_exit_phase[4*i+1])
        molt_exit_ph_L3.append(molt_exit_phase[4*i+2])
        molt_exit_ph_L4.append(molt_exit_phase[4*i+3])

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

    #Larval stage durations
    L1_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L1_int_wt.append((intmolts[1,0,i] - intmolts[0,0,i])/6)

    L2_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L2_int_wt.append((intmolts[1,1,i] - intmolts[0,1,i])/6)



    L3_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L3_int_wt.append((intmolts[1,2,i] - intmolts[0,2,i])/6)

    L4_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L4_int_wt.append((intmolts[1,3,i] - intmolts[0,3,i])/6)

    intermolt_dur = [L1_int_wt, L2_int_wt, L3_int_wt, L4_int_wt]

#molt durations

    M1_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M1_wt.append((intmolts[0,1,i] - intmolts[1,0,i])/6)

    M2_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M2_wt.append((intmolts[0,2,i] - intmolts[1,1,i])/6)

    M3_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M3_wt.append((intmolts[0,3,i] - intmolts[1,2,i])/6)

    M4_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M4_wt.append((intmolts[0,4,i] - intmolts[1,3,i])/6)

    molt_dur = [M1_wt, M2_wt, M3_wt, M4_wt]


    #larval stage durations
    L1_dur_wt = []
    L2_dur_wt = []
    L3_dur_wt = []
    L4_dur_wt = []

    for i in np.arange(0,len(leth.columns)-1):
        L1_dur_wt.append((intmolts[0,1,i]- intmolts[0,0,i])/6)
        L2_dur_wt.append((intmolts[0,2,i] - intmolts[0,1,i])/6)
        L3_dur_wt.append((intmolts[0,3,i] - intmolts[0,2,i])/6)
        L4_dur_wt.append((intmolts[0,4,i] - intmolts[0,3,i])/6)

    larval_stage_dur = [L1_dur_wt, L2_dur_wt, L3_dur_wt, L4_dur_wt]
    
    leth_corr = []
    end = []
    for i in np.arange(1,len(leth.columns)):
        end.append(len(leth)-int(intmolts[0,0,i-1]))
    max_len = np.min(end)
    
    for i in np.arange(1,len(leth.columns)):
        leth_corr.append(leth.iloc[int(intmolts[0,0,i-1]):(int(intmolts[0,0,i-1])+max_len),i])
    
    leth_new = np.asarray(leth_corr)
    
    leth_new_sorted = leth_new[np.argsort(intmolts[1,0,:]-intmolts[0,0,:]),:]
    

    period_L2 = []
    period_L3 = []
    period_L4 = []

    sem_L2 = []
    sem_L3 = []
    sem_L4 = []

    for i in np.arange(0,len(gfpdata.columns)):
            entry_tp_2 = int(intmolts[0,1,i] - intmolts[0,0,i])
            exit_tp_2 = int(intmolts[0,2,i] - intmolts[0,0,i])
            entry_tp_3 = int(intmolts[0,2,i] - intmolts[0,0,i])
            exit_tp_3 = int(intmolts[0,3,i] - intmolts[0,0,i])
            entry_tp_4 = int(intmolts[0,3,i] - intmolts[0,0,i])
            exit_tp_4 = int(intmolts[0,4,i] - intmolts[0,0,i])
            period_L2.append(np.nanmean(PeriodoverTime[i][entry_tp_2:exit_tp_2])/6)
            period_L3.append(np.nanmean(PeriodoverTime[i][entry_tp_3:exit_tp_3])/6)
            period_L4.append(np.nanmean(PeriodoverTime[i][entry_tp_4:exit_tp_4])/6)
            #sem_L1_1.append(sc.sem(periods_all[i][entry_tp_1_1:exit_tp_1_1])/6)
            sem_L2.append(sc.sem(PeriodoverTime[i][entry_tp_2:exit_tp_2])/6)
            sem_L3.append(sc.sem(PeriodoverTime[i][entry_tp_3:exit_tp_3])/6)
            sem_L4.append(sc.sem(PeriodoverTime[i][entry_tp_4:exit_tp_4])/6)

    #Error propagation of the phase calling
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
 



    gfp_adj = []
    for i in np.arange(0,len(gfpdata.columns)): 
        gfp_adj.append(gfpdata.iloc[(int(intmolts[0,0,i])):(int(intmolts[0,0,i])+dev_length),i])
    gfp_adj = np.asarray(gfp_adj)
    
    plt.style.use("classic")

    f = plt.figure(figsize=(15,5), dpi=150)
    
    a1 = f.add_subplot(231)
    linewidth=0.3
    
    
    for i in np.arange(0,len(gfpdata.columns)): 
        gfp_data_adjusted = gfpdata.iloc[(int(intmolts[0,0,i])):(int(intmolts[0,0,i])+dev_length),i]
        a1.plot(np.arange(0,dev_length)/6, gfp_data_adjusted, color="black", linewidth=linewidth, alpha = 0.3)

        for n in np.arange(0,4):
            molt_tp = np.arange((molts[0,n,i])-int(intmolts[0,0,i]),(molts[1,n,i]+1-int(intmolts[0,0,i])))
            gfp_dMolt = gfpdata.iloc[np.arange((molts[0,n,i]),(molts[1,n,i]+1)),i]
            a1.plot((molt_tp/6), gfp_dMolt, color = "red", linewidth = linewidth, alpha = 0.3)        
        a1.set_title("GFP intensities, rel to hatch", fontsize=12)
        a1.set_xlabel("Time after hatch (h)", size=15)
        a1.set_ylim(np.min(np.min(gfpdata)), np.max(np.max(gfpdata)))
        a1.set_xlim(0, dev_length/6)
        a1.set_ylabel("GFP intensities (a.u.)", size=15)        
        a1.set_facecolor("None")
        a1.tick_params(axis='both', which='major', labelsize=15)
        a1.spines['right'].set_visible(False)
        a1.spines['top'].set_visible(False)
        a1.yaxis.set_ticks_position('left')
        a1.xaxis.set_ticks_position('bottom')
    a1.plot(np.arange(0,dev_length)/6, np.nanmean(gfp_adj, axis=0))
    a1.fill_between(np.arange(0,dev_length)/6, np.nanmean(gfp_adj, axis=0)-np.nanstd(gfp_adj, axis=0), np.nanmean(gfp_adj, axis=0)+np.nanstd(gfp_adj, axis=0), color="royalblue", alpha=0.4)
    

    
    cmap = cm.magma
    norm = Normalize(vmin=-3.2, vmax=3.2)
    
    #for scatter
    a2 = f.add_subplot(232)
    size=0.3
    for i in np.arange(0,len(gfpdata.columns)):
        a2.scatter(np.arange(0,dev_length)/6, np.repeat(i, dev_length), color = cmap(norm(my_phase[i])), s=size)
        for n in np.arange(0,4):
            molt_tp = np.arange((molts[0,n,i])-int(intmolts[0,0,i]),(molts[1,n,i]+1-int(intmolts[0,0,i])))
            gfp_dMolt = gfpdata.iloc[np.arange((molts[0,n,i]),(molts[1,n,i]+1)),i]
            plt.scatter((molt_tp/6), np.repeat(i,len(molt_tp)), color = cmap(norm(my_phase[i][int(molt_tp[0]):int(molt_tp[-1]+1)])), s=size+3)        
    a2.set_xlim(0,dev_length/6)
    #a2.set_ylim(np.min(np.min(gfpdata)), 1500)
    a2.set_xlim(0, dev_length/6)
    a2.set_title("Phases of GFP intensities", fontsize=12)
    a2.set_xlabel("Time after hatch (h)", size=15)
    a2.set_ylabel("GFP intensity (a.u.)", size=15)
    a2.tick_params(axis='both', which='major', labelsize=15)
    a2.set_facecolor("None")
    a2.spines['right'].set_visible(False)
    a2.spines['top'].set_visible(False)
    a2.yaxis.set_ticks_position('left')
    a2.xaxis.set_ticks_position('bottom')
    

    a3 = f.add_subplot(233)
    size=30
    molt_ph_entry = [corr_molt_entr_ph_L1,molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4]
    box3 = a3.boxplot(molt_ph_entry, patch_artist=True, showfliers = True, showcaps=False, boxprops=dict(facecolor="grey"))

    molt_ph_exit = [corr_molt_exit_ph_L1,molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4]
    box3_1 = a3.boxplot(molt_ph_exit, patch_artist=True, showfliers = True, showcaps=False, boxprops=dict(facecolor="skyblue"))

    a3.set_title("Molt entry and exit phase", fontsize=12)
    a3.set_xlabel("Molts", size=15)
    a3.set_ylabel("Phase (rad)", size=15)
    a3.set_ylim(-np.pi,4)
    a3.set_xticklabels([1,2,3,4])
    a3.set_xticklabels(["L1", "L2","L3", "L4"], minor=True)
    a3.tick_params(axis='both', which='major', labelsize=15)
    
    colors = ['grey', 'grey', 'grey', 'grey']
    for whisker in box3['whiskers']:
        whisker.set(color="black")
    for patch, color in zip(box3['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_alpha(1)
    colors_3_1 = ['blue', 'blue', 'blue', 'blue']
    for whisker in box3_1['whiskers']:
        whisker.set(color="black")
    for patch, color in zip(box3_1['boxes'], colors_3_1):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_alpha(1)
    a3.set_facecolor("None")
    a3.spines['right'].set_visible(False)
    a3.spines['top'].set_visible(False)
    a3.yaxis.set_ticks_position('left')
    a3.xaxis.set_ticks_position('bottom')
    a3.legend([box3["boxes"][0], box3_1["boxes"][0]], ['Molt entry phase', 'Molt exit phase'], loc="best", prop={'size': 10}, frameon=False)  
    

    plt.style.use("bmh")
    plt.style.use("classic")

    a1_1 = f.add_subplot(235)
    for i in np.arange(0,2):
        a1_1.bar(i+1.85,np.std(molt_ph_entry[i+1])/prop_err_entry[i].std_dev, color="grey", width=0.3, label="entry")
        a1_1.bar(i+2.15,np.std(molt_ph_exit[i+1])/prop_err_exit[i].std_dev, color="blue", width=0.3, label="exit")
    a1_1.set_facecolor("None")
    a1_1.spines['right'].set_visible(False)
    a1_1.spines['top'].set_visible(False)
    a1_1.yaxis.set_ticks_position('left')
    a1_1.xaxis.set_ticks_position('bottom')
    a1_1.set_ylabel("sd obs / sd exp")
    a1_1.set_xlabel("larval stage")
    a1_1.hlines(1,1.5,3.5)
    a1_1.legend(loc=2, fontsize="xx-small", frameon=False)  
    

    periods_and_LS = [L2_dur_wt, period_L2, L3_dur_wt, period_L3, L4_dur_wt, period_L4]
                       
    a4 = f.add_subplot(234)
    box4 = a4.boxplot(periods_and_LS, patch_artist=True, showfliers = True, showcaps=False)
    a4.set_ylim(5,14)
    #a4.set_xlim(0, dev_length/6)
    a4.set_title("period length for larval stages", fontsize=12)
    a4.set_ylabel("GFP period (h)", size=15)
    a4.set_xlabel("larval stage", size=15)
    a4.set_xticklabels(['L2 dur', 'L2 period' 'L3 dur', 'L3 period', 'L4 dur', 'L4 period'], minor=True, rotation=90, fontsize=10)
    a4.set_xticks( [1,2,3,4,5,6] )
    a4.tick_params(axis='both', which='major', labelsize=15)
    colors = ['grey', 'grey', 'grey', 'grey', 'grey', 'grey']
    for whisker in box4['whiskers']:
        whisker.set(color="black")
    for patch, color in zip(box4['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_alpha(0.5)
    a4.spines['right'].set_visible(False)
    a4.spines['top'].set_visible(False)
    a4.yaxis.set_ticks_position('left')
    a4.xaxis.set_ticks_position('bottom')
    a4.set_facecolor("None")
    
    
    f.tight_layout()
        
    canvas = FigureCanvasTkAgg(f, self)
    canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, self)
    toolbar.update()
    canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


    root2.destroy()
 

def stage_dur(self):

    root3 = tk.Tk()
    root3.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
    leth = pd.DataFrame(pd.read_csv(root3.filename))

    #intermolt durations

    L1_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L1_int_wt.append((intmolts[1,0,i] - intmolts[0,0,i])/6)

    L2_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L2_int_wt.append((intmolts[1,1,i] - intmolts[0,1,i])/6)



    L3_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L3_int_wt.append((intmolts[1,2,i] - intmolts[0,2,i])/6)

    L4_int_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        L4_int_wt.append((intmolts[1,3,i] - intmolts[0,3,i])/6)

    intermolt_dur = [L1_int_wt, L2_int_wt, L3_int_wt, L4_int_wt]

#molt durations

    M1_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M1_wt.append((intmolts[0,1,i] - intmolts[1,0,i])/6)

    M2_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M2_wt.append((intmolts[0,2,i] - intmolts[1,1,i])/6)

    M3_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M3_wt.append((intmolts[0,3,i] - intmolts[1,2,i])/6)

    M4_wt = []
    for i in np.arange(0,len(leth.columns)-1):  
        M4_wt.append((intmolts[0,4,i] - intmolts[1,3,i])/6)

    molt_dur = [M1_wt, M2_wt, M3_wt, M4_wt]


    #larval stage durations
    L1_dur_wt = []
    L2_dur_wt = []
    L3_dur_wt = []
    L4_dur_wt = []
    
    for i in np.arange(0,len(leth.columns)-1):
        L1_dur_wt.append((intmolts[0,1,i]- intmolts[0,0,i])/6)
        L2_dur_wt.append((intmolts[0,2,i] - intmolts[0,1,i])/6)
        L3_dur_wt.append((intmolts[0,3,i] - intmolts[0,2,i])/6)
        L4_dur_wt.append((intmolts[0,4,i] - intmolts[0,3,i])/6)

    larval_stage_dur = [L1_dur_wt, L2_dur_wt, L3_dur_wt, L4_dur_wt]
 
    leth_corr = []
    end = []
    for i in np.arange(1,len(leth.columns)):
        end.append(len(leth)-int(intmolts[0,0,i-1]))
    max_len = np.min(end)
    
    for i in np.arange(1,len(leth.columns)):
        leth_corr.append(leth.iloc[int(intmolts[0,0,i-1]):(int(intmolts[0,0,i-1])+max_len),i])
    
    leth_new = np.asarray(leth_corr)
    for i in (np.arange(0,4)):
        print("Intermolt Median = " + str(np.median(intermolt_dur[i])))
        print("Molt Median = " + str(np.median(molt_dur[i])))
        print("Larval stage Median = " + str(np.median(larval_stage_dur[i])))
        
    leth_new_sorted = leth_new[np.argsort(intmolts[1,0,:]-intmolts[0,0,:]),:]
    
    plt.style.use("classic")
    f = Figure(figsize=(9,7), dpi=100)


    a1 = f.add_subplot(221)
    
    box1= a1.boxplot(intermolt_dur, patch_artist=True, showfliers = True, showcaps=False)
    a1.set_title("Intermolt duration", size=20)
    a1.set_ylabel("Duration [h]", size=15)
    a1.set_xlabel("Intermolts", size=15)
    a1.set_ylim(0,15)
    a1.set_facecolor("None")
    a1.tick_params(axis='both', which='major', labelsize=14)
    colors = ['red', 'red', 'red', 'red']
    for patch, color in zip(box1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

    a1.spines['right'].set_visible(False)
    a1.spines['top'].set_visible(False)
    a1.yaxis.set_ticks_position('left')
    a1.xaxis.set_ticks_position('bottom')
    
    a2 = f.add_subplot(222)
    box2 = a2.boxplot(molt_dur, patch_artist=True, showfliers = True, showcaps=False)
    a2.set_title("Molt duration", size=20)
    a2.set_ylabel("Duration [h]", size=15)
    a2.set_xlabel("Molts", size=15)
    a2.set_ylim(0,5)
    a2.set_facecolor("None")
    a2.tick_params(axis='both', which='major', labelsize=14)
    colors = ['blue', 'blue', 'blue', 'blue']
    for patch, color in zip(box2['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
    a2.spines['right'].set_visible(False)
    a2.spines['top'].set_visible(False)
    a2.yaxis.set_ticks_position('left')
    a2.xaxis.set_ticks_position('bottom')

    a3 = f.add_subplot(223)
    box3 = a3.boxplot(larval_stage_dur, patch_artist=True, showfliers = True, showcaps=False)
    a3.set_title("Larval stage duration", size=20)
    a3.set_ylabel("Duration [h]", size=15)
    a3.set_xlabel("Larval stage", size=15)
    a3.set_ylim(0,15)
    a3.set_facecolor("None")
    a3.tick_params(axis='both', which='major', labelsize=14)
    colors = ['green', 'green', 'green', 'green']
    for patch, color in zip(box3['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
    a3.spines['right'].set_visible(False)
    a3.spines['top'].set_visible(False)
    a3.yaxis.set_ticks_position('left')
    a3.xaxis.set_ticks_position('bottom')


    a4 = f.add_subplot(224)
    heatplot = a4.imshow(leth_new_sorted, cmap='plasma_r', interpolation="None", aspect=6)
    a4.set_title("Larval development heatmap", size=20)
    a4.set_ylabel("Worms", size=15)
    a4.set_xlabel("Time points", size=15)
    a4.set_facecolor("None")
    a4.tick_params(axis='both', which='major', labelsize=14)    
    f.tight_layout()
    
    canvas = FigureCanvasTkAgg(f, self)
    #canvas.show()
    canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, self)
    toolbar.update()
    canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


    root3.destroy()


      
class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text = "Single worm imaging analysis", font=LARGE_FONT)
        
        label.pack(pady=10, padx=10)
        
        
        button1 = tk.Button(self, text = "Analyze lethargus data", command = lambda:leth_ana())
        button1.pack()
        
        button2 = tk.Button(self, text = "Analyze GFP", command = lambda:gfp_ana(self))
        button2.pack()
        
        button3 = tk.Button(self, text = "Analyze stage durations", command = lambda:stage_dur(self))
        button3.pack()


app = SWI_analysis()
app.mainloop()
