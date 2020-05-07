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
from bokeh.plotting import figure
from procswi import * #these are the custom functions for the single worm imaging analysis

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

    #plot
    p = figure(
    title='Raw GFP data with molts in red',
    x_axis_label='Time of larval development (h)',
    y_axis_label='GFP intensities (a.u.)')
    for i in np.arange(0,len(gfpdata.columns)): 
        p.line(np.arange(0,dev_length)/6, gfp_adj[i], legend='Single worm GFP trace', line_width=2, line_color="black", line_alpha = 0.4)
        for n in np.arange(0,4):
            molt_tp = np.arange((molts[0,n,i])-int(intmolts[0,0,i]),(molts[1,n,i]+1-int(intmolts[0,0,i])))
            gfp_dMolt = gfpdata.iloc[np.arange((molts[0,n,i]),(molts[1,n,i]+1)),i]
            p.line((molt_tp/6), gfp_dMolt, line_color = "red", legend = "molt", line_width = 2, alpha = 0.3)        
            p.xaxis.major_label_text_font_size = "15pt"
            p.xaxis.axis_label_text_font_size = "15pt"
            p.yaxis.major_label_text_font_size = "15pt"
            p.yaxis.axis_label_text_font_size = "15pt"
    st.bokeh_chart(p, use_container_width=True)

    
    #######run GFP analysis

    # calculations (can be improved with dataframes etc)
    # interpolate gfp data  
    f_clean = interpolate_gfp(gfpdata, intmolts, dev_length)

    # select valid worms
    valid_worms = set_up_valid_worms(gfpdata)
#up to here with functions
    if st.sidebar.checkbox("filter out bad worms"):
        worm_to_check = st.sidebar.number_input("worm number", 1 , step = 1) - 1
        
        #plot individual worm in sidebar
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
    
    #plot f_clean (interpolated data) to investigate molts and developmental length in 
    #more detail
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
        st.sidebar.subheader("adjust y-axis limits")

        ylim_phase_min = st.sidebar.number_input("phase ylim minimum", -4, 4, -4)
        ylim_phase_max = st.sidebar.number_input("phase ylim minimum", -4, 4, 4)
        cmap = cm.magma
        norm = Normalize(vmin=-3.2, vmax=3.2)

        #for scatter
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


        if st.checkbox("flip phases (in case they are at boundary)"):
            molt_phases = use_corrected_phases(molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4, molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4)
            
            #plot
            ax = f.add_subplot(122)
            sns.boxplot(x = "variable", y = "value", data = molt_phases,
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
            molt_phases = use_phases(corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4, corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4)
            
            #plot
            ax = f.add_subplot(122)
            sns.boxplot(x = "variable", y = "value", data = molt_phases,
            hue="Molt", palette = "Blues", ax = ax)
            ax.set_ylabel("Phase (rad)")
            ax.set_xlabel("molts_clean")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.legend(frameon=False)
            ax.set_ylim(ylim_phase_min, ylim_phase_max)
            plt.tight_layout()

            st.pyplot()

        #period and larval stage duration and error propagated plots

        #prepare error propagated data 
        std_phases_all = combine_stds(molt_phases.loc[molt_phases["Molt"]=="entry",:], prop_err_entry, molt_phases.loc[molt_phases["Molt"]=="exit",:], prop_err_exit)

        #plot

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

        st.pyplot()


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
