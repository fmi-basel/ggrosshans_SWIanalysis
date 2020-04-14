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


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

linewidth = 15
linewidth_sidebar = 2
labelsize = 15
#functions:

#lethargus analysis
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

leth = pd.read_csv("G:/user/Yannick.Hauser/Resource Analysis Paper/Results_Figures/Data/SWI/pYPH70_20180309/Lethargus_pYPH70_EV_clean.csv")
gfpdata_original = pd.read_csv("G:/user/Yannick.Hauser/Resource Analysis Paper/Results_Figures/Data/SWI/pYPH70_20180309/Kymograph_Quantification_20180309_BGsub_easy_clean.csv")
gfpdata = gfpdata_original.pivot("Frame","Position", "Intensity_BGsub")

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


    st.sidebar.subheader("choose developmental length (in time points) for GFP analysis")
    dev_length = st.sidebar.number_input("", 0, 360 - int(np.max(intmolts[0,0,:])), 231)

st.sidebar.title("Raw data plotting")
if st.sidebar.checkbox("plot GFP data with molts"):

    st.subheader("Raw GFP data with molts in red")
    st.write("This will be a nice description for the plot")

    #adjusted gfp data
    gfp_adj = []
    for i in np.arange(0,len(gfpdata.columns)): 
        gfp_adj.append(gfpdata.iloc[(int(intmolts[0,0,i])):(int(intmolts[0,0,i])+dev_length),i])
    gfp_adj = np.asarray(gfp_adj)

    #plot
    p = figure(
    title='Raw GFP data with molts in red',
    x_axis_label='Time of larval development (h)',
    y_axis_label='GFP intensities (a.u.)')
    for i in np.arange(0,len(gfpdata.columns)): 
        gfp_data_adjusted = gfpdata.iloc[(int(intmolts[0,0,i])):(int(intmolts[0,0,i])+dev_length),i]
        p.line(np.arange(0,dev_length)/6, gfp_data_adjusted, legend='Single worm GFP trace', line_width=2, line_color="black", line_alpha = 0.4)
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
    @st.cache()
    def interpolate_gfp(gfpdata):
        for i in np.arange(0,len(gfpdata.columns)):
            f = gfpdata.interpolate(method = 'linear', axis =0, limit = 60, limit_direction = 'backward')

        f_clean = []
        for i in np.arange(0, len(gfpdata.columns)):
            f_clean.append(f.iloc[int(intmolts[0,0,i]):int(intmolts[0,0,i]+dev_length),i].values)
        return f_clean
    
    f_clean = interpolate_gfp(gfpdata)

    # select valid worms
    @st.cache(allow_output_mutation=True)
    def set_up_valid_worms():
        valid_worms_1 = np.repeat(1,len(gfpdata.columns))
        return valid_worms_1

    valid_worms = set_up_valid_worms()

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
    #@st.cache()
    #def calculate_durations(leth_clean, intmolts_clean, molts_clean):

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

    #return larval_stage_dur, molt_dur, intermolt_dur
    
    #larval_stage_dur, molt_dur, intermolt_dur = calculate_durations(leth_clean, intmolts_clean, molts_clean)
    

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
     


        
    #Scale the data to each individual larval stage according to mean length of larval stage
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


    new_molt = molts_clean[:,:,:]-intmolts_clean[0,0:4,:]

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

        data = f_clean_df_clean

        y = []
        for i in np.arange(0,len(f_clean_df_clean.columns)):
            y.append(butter_bandpass_filter((data.iloc[:,i]-np.mean(data.iloc[:,i])), lowcut, highcut, fs, order=order))

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
            molt_entr_ph_L1.append(molt_entry_phase[4*i]) #select M1 molt entry phase, +np.pi because wavelets start at -pi (=0 degree) over 0 (=180degree) to pi (=360 degree)
            molt_entr_ph_L2.append(molt_entry_phase[4*i+1])
            molt_entr_ph_L3.append(molt_entry_phase[4*i+2])
            molt_entr_ph_L4.append(molt_entry_phase[4*i+3])
            
            molt_exit_ph_L1.append(molt_exit_phase[4*i])
            molt_exit_ph_L2.append(molt_exit_phase[4*i+1])
            molt_exit_ph_L3.append(molt_exit_phase[4*i+2])
            molt_exit_ph_L4.append(molt_exit_phase[4*i+3])



    #correct (switch) phase 

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
        
        
        leth_corr = []
        end = []
        for i in np.arange(1,len(leth_clean.columns)):
            end.append(len(leth_clean)-int(intmolts_clean[0,0,i-1]))
        max_len = np.min(end)

        for i in np.arange(1,len(leth_clean.columns)):
            leth_corr.append(leth_clean.iloc[int(intmolts_clean[0,0,i-1]):(int(intmolts_clean[0,0,i-1])+max_len),i])

        leth_new = np.asarray(leth_corr)

        leth_new_sorted = leth_new[np.argsort(intmolts_clean[1,0,:]-intmolts_clean[0,0,:]),:]


        period_L2 = []
        period_L3 = []
        period_L4 = []

        sem_L2 = []
        sem_L3 = []
        sem_L4 = []

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

            molt_ph_entry = pd.DataFrame([corr_molt_entr_ph_L1,corr_molt_entr_ph_L2, corr_molt_entr_ph_L3, corr_molt_entr_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
            molt_ph_entry["Molt"] = "entry"
            
            molt_ph_exit = pd.DataFrame([corr_molt_exit_ph_L1, corr_molt_exit_ph_L2, corr_molt_exit_ph_L3, corr_molt_exit_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
            molt_ph_exit["Molt"] = "exit"

            molt_phases = molt_ph_entry.append(molt_ph_exit)
            
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

            molt_ph_entry = pd.DataFrame([molt_entr_ph_L1, molt_entr_ph_L2, molt_entr_ph_L3, molt_entr_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
            molt_ph_entry["Molt"] = "entry"
            
            molt_ph_exit = pd.DataFrame([molt_exit_ph_L1, molt_exit_ph_L2, molt_exit_ph_L3, molt_exit_ph_L4], index = ["M1", "M2", "M3", "M4"]).T.melt()
            molt_ph_exit["Molt"] = "exit"

            molt_phases = molt_ph_entry.append(molt_ph_exit)
            
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
        
        #plot

        f = plt.figure(figsize=(8,3), dpi=150)
        y_lim_low_LS_and_PER = st.sidebar.number_input("y-axis lower limit of larval stage and period", 0,100, 0)
        y_lim_high_LS_and_PER = st.sidebar.number_input("y-axis upper limit of larval stage and period", 0,100, 17)
        
        periods_and_LS = pd.DataFrame([L2_dur_wt, period_L2, L3_dur_wt, period_L3, L4_dur_wt, period_L4], index = ["L2_dur", "L2_period", "L3_dur", "L3_period", "L4_dur", "L4_period"]).T.melt()
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
        phase.to_csv(save_dir_data + "phase.csv") 
        larval_stage_dur.to_csv(save_dir_data + "Larval_stage_durations.csv") 
        molt_dur.to_csv(save_dir_data + "Molt_durations.csv")
        intermolt_dur.to_csv(save_dir_data + "Intermolt_durations.csv")
        params_hilbert.to_csv(save_dir_data + "parameters_used_for_hilbert_analysis.csv")
        valid_worms_df.to_csv(save_dir_data + "valid_worms.csv")
        #write molting time points per worm
        #write phases at molt entry / exit per worm
        #write error prop data






