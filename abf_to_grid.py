import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import pyabf
from pathlib import Path
from scipy.signal import find_peaks 
from scipy.io import loadmat
from asymmetry import abf_to_data as ad
from asymmetry import utils


# a function to use stim times, slice the cell data, find the max, 
# and plot the max value on the corresponding coordinate taken from 
# pattern_order dataframe column 3 on a grid

def plot_heatmap(data_path, protocol_path, recording_type='continuous', stim_times=None, expected_stims=100, cut_trace=False, pulse_slice=[0,100], ISI=20040,
                stim_threshold=100, prop='peak', clip_spikes=True, clip_value=30, Fs=2e4, palette='rocket',protocol_order=''):
    
    # load data
    # if file extension is mat then load mat file
    # if file extension is abf then load abf file
    if data_path.suffix == '.mat':
        data = loadmat(data_path)
        cell = data['data'][0][0]
        stim = data['data'][0][1]
    elif data_path.suffix == '.abf':
        data = ad.abf_to_data(data_path, baseline_subtraction=True,
                    signal_scaling=1, sampling_freq=2e4, filter_type='', filter_cutoff=1000,
                    data_order="channelwise", plot_data=True)

    pattern_order = pd.read_csv(protocol_path, sep=' ', header=None)
    # if recording_type == 'continuous':
    cell = data[0][0].reshape(-1)
    stim = data[0][1].reshape(-1)
    print("Cell data shape: ", cell.shape)

    if cut_trace:
        expected_stims = pulse_slice[1] - pulse_slice[0]

    if stim_times is None:
        _, props = find_peaks(stim, height=stim_threshold, width=10, distance=ISI-1000)
        print("Number of stimulations: ", len(props['left_ips']))
        stim_times = props['left_ips']
        print(stim_times/20000)
        if (not cut_trace) & (len(stim_times) != expected_stims):
            print("Number of stimulations do not match the expected number.")
            if np.abs(len(stim_times) - expected_stims) > 10:
                print('Difference is too large. Breaking.')
                # plot stim time diff
                plt.figure()
                plt.plot(np.diff(stim_times))
                return stim_times, None, None, None, None, None, None, None
            else:
                print("@@ Creating artificial stim times from first pulse and IPI")
                #  creating artificial stim times from first pulse and IPI")
                first_stim = stim_times[0]
                stim_times = np.arange(first_stim, first_stim + expected_stims*ISI, ISI)


    # if cut_trace is true then slice the cell data with pulse_slice
    if cut_trace:
        stim_times = stim_times[pulse_slice[0]:pulse_slice[1]]

    grid_size = pattern_order.iloc[0, 2]
    peakdict = []

    for i, stim_time in enumerate(stim_times):
        if protocol_order == 'random':
            spot_loc = pattern_order.iloc[i, 3] - 1
        else:
            spot_loc = i 
        
        x, y = spot_loc // grid_size, spot_loc % grid_size

        cell_slice = cell[int(stim_time):int(stim_time+ISI)]
        if cell_slice.shape[0] < ISI:
            print("Cell slice is smaller than ISI. Skipping this stimulation.")
            continue
        baseline = np.mean(cell[int(stim_time-100):int(stim_time)])
        cell_slice = cell_slice - baseline
        peak = np.max(cell_slice)
        auc  = np.trapz(cell_slice[:1000]) / Fs
        
        peakdict.append({'x':x, 'y':y, 'peak':peak, 'auc':auc})
        # print(i, x, y, peak)
    
    # convert peakdict to df
    peak_df = pd.DataFrame(peakdict)
    # histogram of response values in new figure
    fig1, [ax1_1,ax1_2] = plt.subplots(1,2)
    sns.histplot(peak_df['peak'], bins=20, kde=True, ax=ax1_1)
    sns.histplot(peak_df['auc'], bins=20, kde=True, ax=ax1_2)
    ax1_1.set_title('Distribution of peak values')
    ax1_1.set_xlabel('Peak value (mV)')
    ax1_2.set_title('Distribution of AUC values')
    ax1_2.set_xlabel('AUC (mV*ms)')
    
    # if clip spikes is true than any peak value above 30 should be clipped to 30
    if clip_spikes:
        print("Clipping spikes above 30 mV")
        peak_df['peak'] = peak_df['peak'].apply(lambda k: clip_value if k > 30 else k)
        peak_df['auc'] = peak_df['auc'].apply(lambda k: 3 if k > 3 else k)
    
    df_pivot_peak = peak_df.pivot(index='x', columns='y', values='peak')
    df_pivot_auc  = peak_df.pivot(index='x', columns='y', values='auc')
    
    fig2_1, ax2_1 = plt.subplots(1,1, figsize=(6,5), layout='tight')
    sns.heatmap(df_pivot_peak, ax=ax2_1, cmap=palette, annot=False, vmin=np.min(df_pivot_peak.values), vmax=np.max(df_pivot_peak.values))
    ax2_1.set_title('Peak Heatmap')
    # remove x and y labels, ticks and grid
    ax2_1.set_xticks([])
    ax2_1.set_yticks([])
    ax2_1.set_xticklabels([])
    ax2_1.set_yticklabels([])
    # axs labels
    ax2_1.set_xlabel('X', fontsize=18)
    ax2_1.set_ylabel('Y', fontsize=18)
    # add colorbar label
    cbar = ax2_1.collections[0].colorbar
    unit = 'mV'
    cbar.set_label(f'Peak Response ({unit})')
    print("Peak Heatmap", np.min(df_pivot_peak.values), np.max(df_pivot_peak.values))

    fig2_2, ax2_2 = plt.subplots(1,1, figsize=(6,5), layout='tight')
    sns.heatmap(df_pivot_auc,  ax=ax2_2, cmap=palette, annot=False, vmin=np.min(df_pivot_auc.values), vmax=np.max(df_pivot_auc.values))
    ax2_2.set_title('AUC Heatmap')
    # remove x and y labels, ticks and grid
    ax2_2.set_xticks([])
    ax2_2.set_yticks([])
    ax2_2.set_xticklabels([])
    ax2_2.set_yticklabels([])
    # axs labels
    ax2_2.set_xlabel('X', fontsize=18)
    ax2_2.set_ylabel('Y', fontsize=18)
    # add colorbar label
    cbar = ax2_2.collections[0].colorbar
    unit = 'mV*ms'
    cbar.set_label(f'AuC Response ({unit})')
    print("AUC Heatmap", np.min(df_pivot_auc.values), np.max(df_pivot_auc.values))

    return stim_times, peak_df, df_pivot_peak, df_pivot_auc, fig1, [fig2_1, fig2_2], cell, stim
