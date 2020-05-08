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
from procswi import * #these are the custom functions for the single worm imaging analysis
from plotting_SWI import *

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

linewidth = 15
linewidth_sidebar = 2
labelsize = 15


@st.cache()
def show_worm_image():
    image = Image.open(os.path.dirname(__file__) + "/SWI_chambers_OP50_202020215.png")
    return image

image = show_worm_image()
st.image(image, caption='', use_column_width=True)

st.title(":microscope: :snake: :microscope: :snake: SWI Analyzer :snake: :microscope: :snake: :microscope:")


if st.sidebar.checkbox("Show useful information on how to work with this tool"):
    str_1 = "This tool is designed to analyze single worm imaging data. It relies on lethargus data in the format of a csv file where the columns indicate the individual worms and the rows represent time points. "
    str_2 = "The lethargus data is usually generated manually by checking whether the worm pumps (=1) or not (=0) at each time point. "
    str_3 = "The GFP data should be in the tidy format, where Timepoints are in the first column (Frame), worms in the second (Position) and GFP intensity in the third column (Intensity_BGsub). "
    str_4 = "The tool is structured in an interactive way through the sidebar. The idea is to go through the steps from top to bottom by first loading in the raw data and then performing the lethargus analysis. "
    str_5 = "After the lethargus analysis the raw GFP data is plotted with the lethargus times in red. It is then possible to clean the data from worms that did not survive the assay or escape during the imaging. "
    str_6 = "All further calculations and plots will be done using the cleaned data. Finally, it is possible to save the figures as pdf or png and the results in csv format to a desired location."
    st.info(str_1 + str_2 + str_3 + str_4 + str_5 + str_6)


#Load the data
st.sidebar.title("Import the data")

@st.cache()
def load_lethargus(uploaded_file_leth):
    return pd.read_csv(uploaded_file_leth)

uploaded_file_leth = st.sidebar.file_uploader("Choose file for lethargus data", type="csv")
if uploaded_file_leth is not None:
    leth = load_lethargus(uploaded_file_leth)

uploaded_file_GFP = st.sidebar.file_uploader("Choose file for GFP data", type="csv")
if uploaded_file_GFP is not None:
    gfpdata_original = pd.read_csv(uploaded_file_GFP)
    gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")

#leth = pd.read_csv("G:/user/Yannick.Hauser/Resource Analysis Paper/Results_Figures/Data/SWI/pYPH70_20180309/Lethargus_pYPH70_EV_clean.csv")
#gfpdata_original = pd.read_csv("G:/user/Yannick.Hauser/Resource Analysis Paper/Results_Figures/Data/SWI/pYPH70_20180309/Kymograph_Quantification_20180309_BGsub_easy_clean.csv")
#gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")

st.sidebar.title("Display raw data")

if st.sidebar.checkbox("show lethargus data"):
    st.header("Raw data input")
    st.subheader("Lethargus data:")
    st.dataframe(leth)

if st.sidebar.checkbox("show GFP data"):
    st.header("Raw data input")
    st.subheader("GFP data:")
    st.dataframe(gfpdata)

#run the lethargus analysis

st.sidebar.title("Lethargus analysis")

if st.sidebar.checkbox("start lethargus analysis"):
    intmolts, molts = lethargus_analysis(leth)
    st.success("Lethargus analysis successful!")


#choose developmental length
    st.sidebar.subheader("choose developmental length (in time points) for GFP analysis")
    dev_length = st.sidebar.number_input("", 0, 360 - int(np.max(intmolts[0,0,:])), 231)

st.sidebar.title("Raw data plotting")
if st.sidebar.checkbox("plot GFP data with molts"):

    st.subheader("Raw GFP data with molts in red")
    st.write("This will be a nice description for the plot")

    #adjusted gfp data to start from hatch for every worm individually
    gfp_adj = adjust_gfp(gfpdata, intmolts, dev_length)

    #plot GFP raw data
    #figure_raw_gfp = bokeh_plot_gfp_raw(gfpdata, gfp_adj, molts, intmolts, dev_length)
    #st.bokeh_chart(figure_raw_gfp, use_container_width=True)

    #######run GFP analysis

    # calculations (can be improved with dataframes etc)
    # interpolate gfp data  
    f_clean = interpolate_gfp(gfpdata, intmolts, dev_length)

    # select valid worms
    valid_worms = set_up_valid_worms(gfpdata)

    if st.sidebar.checkbox("filter out bad worms"):
        worm_to_check = st.sidebar.number_input("worm number", 1 , step = 1) - 1
        
        #plot individual worm in sidebar
        evaluated = plot_ind_worm_sidebar(gfpdata, f_clean, dev_length, molts, intmolts, worm_to_check, linewidth_sidebar)
        evaluated = st.sidebar.checkbox("valid worm", value = True)
        valid_worms[worm_to_check] = evaluated
        valid_worms_df = pd.DataFrame({"worm": gfpdata.columns, "valid": valid_worms})
    
    #filter gfpdata and lethargus for valid samples
    leth_clean = leth.iloc[:,1:].loc[:,valid_worms==True]
    leth_clean.insert(loc=0, column='Timepoint', value=leth["Timepoint"])

    f_clean_df = pd.DataFrame(f_clean).T
    f_clean_df_clean = f_clean_df.loc[:,valid_worms==True]
    intmolts_clean = intmolts[:,:,valid_worms==True]
    molts_clean = molts[:,:,valid_worms==True]
    gfpdata_clean = gfpdata.loc[:,valid_worms == True]
    
    #plot f_clean (interpolated data) to investigate molts and developmental length in more detail
    if st.sidebar.checkbox("Cleaned data: highlight molt times in comparison to developmental length"):
        st.subheader("Cleaned GFP data with molts")
        
        def color_zero_red(val):
    
            color = 'red' if val == 0 else 'black'
            return 'color: %s' % color

        st.dataframe(valid_worms_df.T.style.applymap(color_zero_red))

        alpha_gfp = st.sidebar.slider("transparency of single worm GFP", 0,100,40)
        alpha_molt = st.sidebar.slider("transparency of molts", 0,100,40)
        alpha_mean = st.sidebar.slider("transparency of mean", 0,100,90)
        alpha_std = st.sidebar.slider("transparency of standard deviation", 0,100,20)
        #plot
        plot_gfp_and_molts(gfpdata, dev_length, f_clean_df_clean, alpha_mean, alpha_gfp, alpha_molt, alpha_std, molts_clean, intmolts_clean, labelsize, linewidth)

    
    st.sidebar.title("Plot developmental durations")

    #Larval stage durations
    larval_stage_dur, molt_dur, intermolt_dur = calculate_durations(leth_clean, intmolts_clean, molts_clean)


    if st.sidebar.checkbox("plot molt, intermolt and larval stage duration"):
        st.subheader("Molt, Intermolt and Larval stage durations in hours")
        y_lim_low_molt = st.sidebar.number_input("molt y-axis lower limit", 0,100, 0)
        y_lim_high_molt = st.sidebar.number_input("molt y-axis upper limit", 0,100, 5)
        y_lim_low_intermolt = st.sidebar.number_input("intermolt y-axis lower limit", 0,100, 0)
        y_lim_high_intermolt = st.sidebar.number_input("intermolt y-axis upper limit", 0,100, 17)
        y_lim_low_larvalstage = st.sidebar.number_input("larval stage y-axis lower limit", 0,100, 0)
        y_lim_high_larvalstage = st.sidebar.number_input("larval stage y-axis upper limit", 0,100, 17)
        #plot
        plot_durations(molt_dur, intermolt_dur, larval_stage_dur, y_lim_low_molt, y_lim_high_molt, y_lim_low_intermolt, y_lim_high_intermolt, y_lim_low_larvalstage, y_lim_high_larvalstage)



        
    #Scale the data to each individual larval stage according to mean length of larval stage (not working yet)
    #scaled_all, new_molt_mean, new_molt_std = scale_data(intmolts_clean, molts_clean, gfpdata_clean)
    
    #hilbert analysis
    st.sidebar.title("Hilbert analysis for phases at molt entry and exit")

    if st.sidebar.checkbox("plot phase from Hilbert at molt entry and exit"):

        # Sample rate and desired cutoff frequencies (in 1/h).
        st.sidebar.subheader("parameters: only for advanced users!")
        if st.sidebar.checkbox("set parameters for BW filter and hilbert transform"):
            st.sidebar.warning("only for advanced users")
            fs = st.sidebar.number_input("sample frequency (tp/h)", 1,12, 6) #default = 6
            lowcut = 1/st.sidebar.number_input("upper limit period (h)", 5,20, 14) #default = 14
            highcut = 1/st.sidebar.number_input("lower limit period (h)", 1,10, 5) #default = 5
        
        else:
            fs = 6 #default:6 
            lowcut = 1/14 #default 1/14 --> 14 hour period
            highcut = 1/5 #default 1/5 --> 5 hour period

        my_phase, PeriodoverTime, phase_melt = run_hilbert(f_clean_df_clean, leth_clean, lowcut, highcut, fs, dev_length)
 
        #obtain phase at molt entry and exit
        molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4 = phase_molt_entry_exit(gfpdata_clean, intmolts_clean, my_phase)
        
        #correct (switch) phase 
        corr_molt_entr_ph_L1, corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4 = correct_phase(molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4)
        
        #sort lethargus (optional, not used for plotting)
        leth_new_sorted = sort_lethargus(leth_clean, intmolts_clean)

        #calculate period per larval stage
        period_L2, period_L3, period_L4 = period_per_LS(PeriodoverTime, gfpdata_clean, intmolts_clean)
        
        #Error propagation of the phase calling
        prop_err_entry, prop_err_exit = error_prop(larval_stage_dur.loc[larval_stage_dur["variable"]=="L2",:]["value"].values, larval_stage_dur.loc[larval_stage_dur["variable"]=="L3",:]["value"].values, larval_stage_dur.loc[larval_stage_dur["variable"]=="L4",:]["value"].values, intermolt_dur.loc[intermolt_dur["variable"]=="IM2",:]["value"].values, intermolt_dur.loc[intermolt_dur["variable"]=="IM3",:]["value"].values, intermolt_dur.loc[intermolt_dur["variable"]=="IM4",:]["value"].values, period_L2, period_L3, period_L4)
        
        #plot 2
        st.subheader("Phases at molt entry and exit")

        molt_phases_normal = use_phases(molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4)
        molt_phases_corr = use_corrected_phases(corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4)

        st.sidebar.subheader("adjust y-axis limits")

        ylim_phase_min = st.sidebar.number_input("phase ylim minimum", -4, 4, -4)
        ylim_phase_max = st.sidebar.number_input("phase ylim minimum", -4, 4, 4)
        cmap = cm.magma
        norm = Normalize(vmin=-3.2, vmax=3.2)
        
        
        #plot scatter
        plot_phases(molt_phases_normal, molt_phases_corr, cmap, norm, ylim_phase_min, ylim_phase_max, gfpdata_clean, dev_length, my_phase, molts_clean, intmolts_clean, corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4, molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4)
        #period and larval stage duration and error propagated plots

        #prepare error propagated data
        st.subheader("Larval stage vs osc period and error propagation")

        if st.checkbox("use flipped phases for error propagation", value = False):
            std_phases_all = combine_stds(molt_phases_corr.loc[molt_phases_corr["Molt"]=="entry",:], prop_err_entry, molt_phases_corr.loc[molt_phases_corr["Molt"]=="exit",:], prop_err_exit)
        else:
            std_phases_all = combine_stds(molt_phases_normal.loc[molt_phases_normal["Molt"]=="entry",:], prop_err_entry, molt_phases_normal.loc[molt_phases_normal["Molt"]=="exit",:], prop_err_exit)

        #plot
        plot_period_and_error_prop(larval_stage_dur, period_L2, period_L3, period_L4, std_phases_all)

    st.sidebar.title("export results")
    if st.sidebar.checkbox("Save results in following directory"):
        
        params_hilbert = pd.DataFrame({"sample frequency (tp)": fs, 
                                        "lower period limit (h)": 1/highcut,
                                        "upper period limit (h)": 1/lowcut}, index = ["value"]).T

        
        save_dir_data = st.sidebar.text_input("add location", "")
        my_phase.to_csv(save_dir_data + "phase.csv") 
        larval_stage_dur.to_csv(save_dir_data + "Larval_stage_durations.csv") 
        molt_dur.to_csv(save_dir_data + "Molt_durations.csv")
        intermolt_dur.to_csv(save_dir_data + "Intermolt_durations.csv")
        params_hilbert.to_csv(save_dir_data + "parameters_used_for_hilbert_analysis.csv")
        valid_worms_df.to_csv(save_dir_data + "valid_worms.csv")
        #write molting time points per worm
        #write phases at molt entry / exit per worm
        #write error prop data
