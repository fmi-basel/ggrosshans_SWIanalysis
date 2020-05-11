import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import imp
import scipy.stats as sc
from scipy.signal import butter, filtfilt
from scipy.signal import hilbert
from uncertainties import ufloat
import streamlit as st
import time
from bokeh.plotting import figure, output_file, show
import matplotlib
import io
import os
from PIL import Image


#plotting functions
#@st.cache()
def bokeh_plot_gfp_raw(gfpdata, gfp_adj, molts, intmolts, dev_length):
    p = figure(
    title='Raw GFP data with molts in red',
    x_axis_label='Time of larval development (h)',
    y_axis_label='GFP intensities (a.u.)')
    for i in np.arange(0,len(gfpdata.columns)): 
        p.line(np.arange(0,dev_length)/6, gfp_adj[i], legend_label='Single worm GFP trace', line_width=2, line_color="black", line_alpha = 0.4)
        for n in np.arange(0,4):
            molt_tp = np.arange((molts[0,n,i])-int(intmolts[0,0,i]),(molts[1,n,i]+1-int(intmolts[0,0,i])))
            gfp_dMolt = gfpdata.iloc[np.arange((molts[0,n,i]),(molts[1,n,i]+1)),i]
            p.line((molt_tp/6), gfp_dMolt, line_color = "red", legend_label = "molt", line_width = 2, alpha = 0.3)        
            p.xaxis.major_label_text_font_size = "15pt"
            p.xaxis.axis_label_text_font_size = "15pt"
            p.yaxis.major_label_text_font_size = "15pt"
            p.yaxis.axis_label_text_font_size = "15pt"
    return p

#plot individual GFP intensities in sidebar to select only valid worms
def plot_ind_worm_sidebar(gfpdata, f_clean, dev_length, molts, intmolts, worm_to_check, linewidth_sidebar):
    fig, ax = plt.subplots(1,1)
    ax.plot(np.arange(0,dev_length)/6, f_clean[worm_to_check], color="black", linewidth=linewidth_sidebar, alpha = 0.8)

    for n in np.arange(0,4):
        molt_tp = np.arange((molts[0,n,worm_to_check])-int(intmolts[0,0,worm_to_check]),(molts[1,n,worm_to_check]+1-int(intmolts[0,0,worm_to_check])))
        gfp_dMolt = f_clean[worm_to_check][int(molts[0,n,worm_to_check])-int(intmolts[0,0,worm_to_check]):int(molts[1,n,worm_to_check]+1)-int(intmolts[0,0,worm_to_check])]
        ax.plot((molt_tp/6), gfp_dMolt, color = "red", linewidth = linewidth_sidebar+1, alpha = 1)        
    ax.axvline(dev_length/6, 0, np.max(np.max(f_clean[worm_to_check])), color="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    st.sidebar.pyplot()


    #plot f_clean (interpolated data) to investigate molts and developmental length in more detail

def plot_gfp_and_molts(gfpdata, dev_length, f_clean_df_clean, alpha_mean, alpha_gfp, alpha_molt, alpha_std, molts_clean, intmolts_clean, labelsize, linewidth):
    f = plt.figure(figsize=(8,4), dpi=150)

    a1 = f.add_subplot(111)
    linewidth=0.3
    a1.plot(np.arange(0,dev_length)/6, f_clean_df_clean.mean(axis=1), color = "blue", alpha = alpha_mean/100)
    a1.fill_between(np.arange(0,dev_length)/6, f_clean_df_clean.mean(axis=1)-f_clean_df_clean.std(axis=1), f_clean_df_clean.mean(axis=1)+f_clean_df_clean.std(axis=1), color="royalblue", alpha=(alpha_std/100))
    
    for i in np.arange(0,len(f_clean_df_clean.columns)): 
        a1.plot(np.arange(0,dev_length)/6, f_clean_df_clean.iloc[:,i], color="black", linewidth=linewidth, alpha = alpha_gfp/100)
        for n in np.arange(0,4):
            molt_tp = np.arange((molts_clean[0,n,i])-int(intmolts_clean[0,0,i]),(molts_clean[1,n,i]+1-int(intmolts_clean[0,0,i])))
            gfp_dMolt = f_clean_df_clean.iloc[:,i][int(molts_clean[0,n,i])-int(intmolts_clean[0,0,i]):int(molts_clean[1,n,i]+1)-int(intmolts_clean[0,0,i])]
            a1.plot((molt_tp/6), gfp_dMolt, color = "red", linewidth = linewidth+1, alpha = alpha_molt/100)        
            #print(molt_tp)
        a1.axvline(dev_length/6, 0, np.max(np.max(f_clean_df_clean.iloc[:,i])), color="black")
        a1.set_title("GFP intensities, interpolated and relative to hatch", fontsize=10)
        a1.set_xlabel("Time of larval development (h)", size=labelsize)
        a1.set_ylim(np.min(np.min(f_clean_df_clean)), np.max(np.max(f_clean_df_clean)))
        a1.set_xlim(0, len(gfpdata)/6)
        a1.set_ylabel("GFP intensities (a.u.)", size=labelsize)        
        a1.set_facecolor("None")
        a1.tick_params(axis='both', which='major', labelsize=labelsize)
        a1.spines['right'].set_visible(False)
        a1.spines['top'].set_visible(False)
        a1.yaxis.set_ticks_position('left')
        a1.xaxis.set_ticks_position('bottom')
    plt.tight_layout()

    if st.checkbox("click for saving figure 1"):
        save_fig2 = str(st.text_input("location/filename to save the figure:", ""))
        plt.savefig(save_fig2)
        st.pyplot()
    else:
        st.pyplot()


def plot_durations(molt_dur, intermolt_dur, larval_stage_dur, y_lim_low_molt, y_lim_high_molt, y_lim_low_intermolt, y_lim_high_intermolt, y_lim_low_larvalstage, y_lim_high_larvalstage):
    fig, ax = plt.subplots(1,3, figsize=(21,7), dpi=150)
    width = 0.8
    labelsize = 30
    sns.boxplot(x = "variable", y = "value", data = molt_dur, width = width, palette = "Blues", ax = ax[0])
    sns.boxplot(x = "variable", y = "value", data = intermolt_dur, width = width, palette = "Greens", ax = ax[1])
    sns.boxplot(x = "variable", y = "value", data = larval_stage_dur, width = width, palette = "Greys", ax = ax[2])
    ax[0].set_ylim(y_lim_low_molt, y_lim_high_molt)
    ax[1].set_ylim(y_lim_low_intermolt, y_lim_high_intermolt)
    ax[2].set_ylim(y_lim_low_larvalstage, y_lim_high_larvalstage)
    ax[0].set_ylabel("duration (h)", size = labelsize)
    ax[1].set_ylabel("")
    ax[2].set_ylabel("")

    for axis in ax:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(labelsize=labelsize)
        axis.set_xlabel("")
    fig.tight_layout(w_pad = 5)

    if st.checkbox("click for saving figure 3"):
        save_fig3 = str(st.text_input("location/filename to save the figure:", ""))
        plt.savefig(save_fig3)
        st.pyplot()
    else:
        st.pyplot()
     

#plot phases and phase at molt entry / exit
def plot_phases(molt_phases_normal, molt_phases_corr, cmap, norm, ylim_phase_min, ylim_phase_max, gfpdata_clean, dev_length, my_phase, molts_clean, intmolts_clean, corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4, molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4):
    f = plt.figure(figsize=(8,3), dpi=150)

    a2 = f.add_subplot(121)
    size=0.3
    for i in np.arange(0,len(gfpdata_clean.columns)):
        sc = a2.scatter(np.arange(0,dev_length)/6, np.repeat(i, dev_length), color = cmap(norm(my_phase[i])), cmap = "magma", s=size)
        for n in np.arange(0,4):
            molt_tp = np.arange((molts_clean[0,n,i])-int(intmolts_clean[0,0,i]),(molts_clean[1,n,i]+1-int(intmolts_clean[0,0,i])))
            gfp_dMolt = gfpdata_clean.iloc[np.arange((molts_clean[0,n,i]),(molts_clean[1,n,i]+1)),i]
            plt.scatter((molt_tp/6), np.repeat(i,len(molt_tp)), color = cmap(norm(my_phase[i][int(molt_tp[0]):int(molt_tp[-1]+1)])), s=size+3)        
    a2.set_xlim(0,dev_length/6)
    #a2.set_ylim(np.min(np.min(gfpdata_clean)), 1500)
    a2.set_xlim(0, dev_length/6)
    a2.set_xlabel("Time after hatch (h)", size=10)
    a2.set_ylabel("GFP intensity (a.u.)", size=10)
    a2.tick_params(axis='both', which='major', labelsize=10)
    a2.set_facecolor("None")
    a2.spines['right'].set_visible(False)
    a2.spines['top'].set_visible(False)
    a2.yaxis.set_ticks_position('left')
    a2.xaxis.set_ticks_position('bottom')


    if st.checkbox("flip phases (in case they are at boundary)", value = False):
        
        #plot
        ax = f.add_subplot(122)
        sns.boxplot(x = "variable", y = "value", data = molt_phases_corr,
        hue="Molt", palette = "Blues", ax = ax)
        ax.set_ylabel("Phase (rad)")
        ax.set_xlabel("molts_clean")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False)
        ax.set_ylim(ylim_phase_min, ylim_phase_max)
        plt.tight_layout()

        st.pyplot()

    else:
        molt_phases_normal
        
        #plot
        ax = f.add_subplot(122)
        sns.boxplot(x = "variable", y = "value", data = molt_phases_normal,
        hue="Molt", palette = "Blues", ax = ax)
        ax.set_ylabel("Phase (rad)")
        ax.set_xlabel("molts_clean")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False)
        ax.set_ylim(ylim_phase_min, ylim_phase_max)
        plt.tight_layout()
        if st.checkbox("click for saving figure 4"):
            save_fig4 = str(st.text_input("location/filename to save the figure:", ""))
            plt.savefig(save_fig4)
            st.pyplot()
        else:
            st.pyplot()


def plot_period_and_error_prop(larval_stage_dur, period_L2, period_L3, period_L4, std_phases_all):
    f = plt.figure(figsize=(8,3), dpi=150)
    y_lim_low_LS_and_PER = st.sidebar.number_input("y-axis lower limit of larval stage and period", 0,100, 0)
    y_lim_high_LS_and_PER = st.sidebar.number_input("y-axis upper limit of larval stage and period", 0,100, 17)
    
    periods_and_LS = pd.DataFrame([larval_stage_dur.loc[larval_stage_dur["variable"]=="L2",:]["value"].values, period_L2, larval_stage_dur.loc[larval_stage_dur["variable"]=="L3",:]["value"].values, period_L3, larval_stage_dur.loc[larval_stage_dur["variable"]=="L4",:]["value"].values, period_L4], index = ["L2_dur", "L2_period", "L3_dur", "L3_period", "L4_dur", "L4_period"]).T.melt()
    a1 = f.add_subplot(121)
    sns.boxplot(x = "variable", y = "value", data = periods_and_LS,
    palette = "Greys", ax = a1)
    a1.set_ylim(y_lim_low_LS_and_PER, y_lim_high_LS_and_PER)
    a1.set_ylabel("Durations (h)")
    a1.spines["top"].set_visible(False)
    a1.spines["right"].set_visible(False)
    
    a1_1 = f.add_subplot(122)
    sns.barplot(x = "observation", y = "ratio sd_obs/sd_exp", data = std_phases_all, 
    color="observation", hue = "entry_or_exit", palette = "Blues_r", ax = a1_1)
    a1_1.legend(loc=2, fontsize="small", frameon=False)  
    a1_1.spines["top"].set_visible(False)
    a1_1.spines["right"].set_visible(False)
    a1_1.set_ylim(0,1.1)
    plt.tight_layout()

    if st.checkbox("click for saving figure 5"):
        save_fig5 = str(st.text_input("location/filename to save the figure:", ""))
        plt.savefig(save_fig5)
        st.pyplot()
    else:
        st.pyplot()











