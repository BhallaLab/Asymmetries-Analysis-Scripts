import sys
import os
import datetime
import importlib
import pathlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.signal import filter_design
from scipy.signal import butter, bessel, decimate, sosfiltfilt
from scipy.signal import find_peaks, peak_widths
from scipy import stats

from asymmetry import utils
from asymmetry import pattern_index

def create_edge_colormap():
    # get only the middle row
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap
   
    nodes = [0, 0.28, 0.52, 0.78, 1.0]
    colors = np.array([[109,11,212], [53,42,212], [28,83,202], [0,163,178], [92,228,112]])/256

    edge = LinearSegmentedColormap.from_list("edge", list(zip(nodes, colors)))
    edge_r = edge.reversed()

    if not 'edge' in mpl.colormaps:
        mpl.colormaps.register(cmap=edge)
        mpl.colormaps.register(cmap=edge_r)

create_edge_colormap()


def simplify_axes(axes, splines_to_keep=['bottom','left'], axis_offset=10, remove_ticks=True, xtick_locs=[], xtick_labels=[], ytick_locs=[], ytick_labels=[]):
    '''simplify axis properties to remove clutter like ticks, ticklabels, spines, etc.'''
    # check if ax is a list of axes or a numpy array
    if not isinstance(axes, (list, np.ndarray)):
        axes = [axes]

    for ax in axes:        
        # remove spines
        for side in ['top', 'right', 'left', 'bottom']:
            if side not in splines_to_keep: # remove top and right
                ax.spines[side].set_visible(False)
                ax.set_xticks([]) if side=='bottom' else ''
                ax.set_yticks([]) if side=='left' else ''
                # set xticklabels and yticklabels to empty list
                ax.set_xticklabels([]) if side=='bottom' else ''
                ax.set_yticklabels([]) if side=='left' else ''
            else:    # keep ticks on all splines
                ax.spines[side].set_linewidth(0.5)
                ax.spines[side].set_position(('outward', axis_offset))
                ax.set_xticks(xtick_locs, labels=xtick_labels)
                ax.set_yticks(ytick_locs, labels=ytick_labels)
                # set xlim and ylim
                ax.set_xlim([min(xtick_locs), max(xtick_locs)])
                ax.set_ylim([min(ytick_locs), max(ytick_locs)])

        #remove title       
        ax.set_title('')

        # remove legend box and location top right
        ax.legend(frameon=False, loc='upper right')

    return axes


def split_axes(ax, which_axes=['left', 'bottom'], offset=10):
    _adjust_spines(ax, which_axes, offset_distance=offset)


def _adjust_spines(ax, visible_spines, offset_distance=10):
    '''adjust spines to be inside or outside the plot'''
    
    # check if ax is a list of axes or a numpy array
    if not isinstance(ax, (list, np.ndarray)):
        axs = [ax]

    for axx in axs:
        axx.label_outer(remove_inner_ticks=True)

        for loc, spine in axx.spines.items():
            print(loc, spine)
            if loc in visible_spines:
                spine.set_position(('outward', offset_distance))  # outward by 10 points
            else:
                spine.set_visible(False)


def add_floating_scalebar(ax, scalebar_origin=[0,0], xlength=1.0, ylength=1.0, labelx='', labely='', unitx='', unity='', fontsize=12, color='black', linewidth=2, pad=0.01, show_labels=False):
    """Simplifies a matplotlib axes object and adds a floating scalebar.
    Args:
        ax: matplotlib axes object
        x: x position of the scalebar in data coordinates
        y: y position of the scalebar in data coordinates
        xl: length of the x-axis of scalebar in data coordinates
        yl: length of the y-axis of scalebar in data coordinates
        labelx: label of the a-axis of scalebar
        labely: label of the y-axis of scalebar
        unitx: units of the x-axis of scalebar
        unity: units of the y-axis of scalebar
        fontsize: fontsize of the label
        color: color of the scalebar
        linewidth: linewidth of the scalebar
        pad: padding between the scalebar and the label
    Returns:
        None
    Example:
        add_floating_scalebar(ax, 1, 1, 0.1, 0.2, '0.1', '0.2', 's', 'mV')

    """
    x,y = scalebar_origin
    xl = xlength
    yl = ylength

    # draw a line for x axis
    ax.plot([x, x+xl], [y, y]   , color=color, linewidth=linewidth)

    # draw a line for y axis
    ax.plot([x, x],    [y, y+yl], color=color, linewidth=linewidth)
    
    if show_labels:    
        # write x axis label
        ax.text(x+xl/2, y-2*pad, labelx+' '+unitx, fontsize=fontsize, horizontalalignment='center', verticalalignment='top')
        # write y axis label
        ax.text(x-pad, y+yl/2, labely+' '+unity, fontsize=fontsize, horizontalalignment='right', verticalalignment='center', rotation=90) # alignment of the rotated text as a block
    

def plot_abf_data(dataDict, label=""):
    numChannels = len(dataDict[0])
    chLabels    = list(dataDict[0].keys())
    sweepLength = len(dataDict[0][chLabels[0]])

    if 'Time' in chLabels:    
        timeSignal = dataDict[0]['Time']
        chLabels.remove('Time')
    else:
        timeSignal = np.arange(0,sweepLength/2e4,1/2e4)
    
    numPlots = len(chLabels)
    fig,axs     = plt.subplots(numPlots,1,sharex=True)
    
    for sweepData in dataDict.values():
        for i,ch in enumerate(chLabels):
            if ch == 'Cmd':
                axs[i].plot(timeSignal[::5],sweepData[ch][::5],'r')
                axs[i].set_ylabel('Ch#0 Command')
            else:
                axs[i].plot(timeSignal[::5],sweepData[ch][::5],'b')
                axs[i].set_ylabel('Ch# '+str(ch))

    axs[-1].set_xlabel('Time (s)')
    axs[-1].annotate('* Data undersampled for plotting', xy=(1.0, -0.5), xycoords='axes fraction',ha='right',va="center",fontsize=6)
    fig.suptitle(label + ' - ABF Data*')
    plt.show()


def plot_data_from_df(df, data_start_column = 49, signals_to_plot=['Cell','FrameTTL', 'PD', 'Field'], signal_colors=['black','red','cyan','orange'], combine=False, fig=None, ax=None, signal_mapping={}):
    start = data_start_column
    Fs = 2e4
    sweeps = df.shape[0]
    width = int( (df.shape[1] - start - 24) / 4 )
    T = width / Fs
    num_plots = len(signals_to_plot)
    signal_location = {'Cell':slice(start, start+width),
                       'FrameTTL':slice(start+width, start+2*width),
                       'PD':slice(start+2*width, start+3*width),
                       'Field':slice(start+3*width, start+4*width)}
    signalcolors = {}
    for sig in signals_to_plot:
        signalcolors[sig] = signal_colors[signals_to_plot.index(sig)]

    assert len(signals_to_plot) == len(signalcolors)

    # if combine plots is false, draw all the 4 signals separately on 4 subplots
    if combine is False:
        print('Plotting all signals separately')

        # check if fig and ax are supplied:
        if fig is None:
            fig = plt.figure(layout='constrained', figsize=(10, 4), )
        else:
            gridspec = ax.get_subplotspec().get_gridspec()
            ax.remove()
            subfig = fig.add_subfigure(gridspec[:, 0])
        subfigs = fig.subfigures(num_plots,1)
        
        axs = []
        for f in range(num_plots):
            subfig_axs = subfigs[f].subplots(1,1, sharex=True, sharey=True)
            axs.append(subfig_axs)

        time = np.linspace(0, T, num=width, endpoint=False)
        # copy time vector as many times as there are sweeps
        Time = np.tile(time, (sweeps,1) )

        for s, signal in enumerate(signals_to_plot):
            locs = signal_location[signal]
            for i in range(sweeps):
                # start = data_start_column
                trace = df.iloc[i, locs]
                trace = utils.map_range(trace, 0, 5, 0,5)
                axs[s].plot(time, trace, signalcolors[signal], linewidth=1, alpha=0.1)
                axs[s].set_ylabel(signal)
            axs[s].plot(time, df.iloc[:,    locs].mean(axis=0), color=signalcolors[signal], linewidth=1, label=signal)

        axs[-1].set_xlabel('Time (s)')

        return fig, axs
    
    # if combine plots is true, draw all the 4 signals on a single plot
    elif combine is True:
        print('Plotting all 4 signals on a single plot')
        cell_max, cell_min = np.round(np.max(df.iloc[:,49:20049]),2) , np.round(np.min(df.iloc[:,49:20049]),2)
        # print(cell_max, cell_min)
        # cell_max, cell_min = np.round(cell_max, -2), np.sign(cell_min) * (np.remainder(cell_min, 10) - cell_min)
        # print(cell_max, cell_min)
        if not signal_mapping:
            print('remapping to default')
            signal_mapping = {'Cell':[cell_min, cell_max, 2, 4],
                            'FrameTTL': [0, 5, 5, 6],
                            'PD': [0, 1, 4, 5],
                            'Field': [-0.5, 0.5, 0, 2]}

        # check if ax is supplied
        if fig is None and ax is None:
            fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.set_ylim([0,6])
        time = np.linspace(0, T, num=width, endpoint=False)
        # copy time vector as many times as there are sweeps
        Time = np.tile(time, (sweeps,1) )

        for s,signal in enumerate(signals_to_plot):
            locs = signal_location[signal]
            from0, from1, to0, to1 = signal_mapping[signal]
            # print(s, signal, from0, from1, to0, to1)
            for i in range(sweeps):
                # start = data_start_column
                trace = df.iloc[i, locs]                
                trace = utils.map_range(trace, from0, from1, to0, to1)
                ax.plot(time, trace, signal_colors[s], linewidth=1, alpha=0.1)
                ax.set_ylabel(signal)
                trace_average = df.iloc[:,   locs].mean(axis=0)
                trace_average = utils.map_range(trace_average, from0, from1, to0, to1)
                ax.plot(time, trace_average, color=signal_colors[s], linewidth=1, label=signal)
                if signal=='Cell':
                    add_floating_scalebar(ax, scalebar_origin=[0.05, 3.0], xlength=0.1, ylength=2, labelx='', labely=f'{np.round(cell_max, 2)}', unitx='', unity='',
                    fontsize=12, color=signal_colors[s], linewidth=2, pad=0.1, show_labels=True)
                if signal == 'Field':
                    add_floating_scalebar(ax, scalebar_origin=[0.03, 1.2], xlength=0.1, ylength=2, labelx='', labely='1.0', unitx='', unity='mV',
                    fontsize=12, color=signal_colors[s], linewidth=2, pad=0.1, show_labels=False)

                
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Voltage')

        return fig, ax, signal_mapping
    

def plot_grid(spot_locs=[], spot_values=[], grid=[24,24], ax=None, vmin=0, vmax=1, cmap='gray', locs_is_patternID=False, **kwargs):
    '''
    plot a grid of values as a heatmap. input spot_values should have coord column as index
    spot_locs is raw coordinates of the spots, spot_values is the values at those spots
    '''
    if ax is None:
        fig, ax = plt.subplots()


    if len(spot_values) == 0:
        raise ValueError(f'spot_locs and spot_values must be the same length but spot_locs has length {len(spot_locs)} and spot_values has length {len(spot_values)}')
    elif len(spot_values) == 1:
        spot_values = np.repeat(spot_values, len(spot_locs))
    elif len(spot_values) != len(spot_locs):
        raise ValueError(f'spot_locs and spot_values must be the same length, but spot_locs has length {len(spot_locs)} and spot_values has length {len(spot_values)}')


    # make a zero array of the grid size
    grid_array = np.zeros(grid)
    
    # fill the grid array with the spot locations
    for i, spot in enumerate(spot_locs):
        if locs_is_patternID:
            spots = pattern_index.patternID[spot]
            for s in spots:
                locx = (s-1) % grid[0]
                locy = (s-1) // grid[1]
                grid_array[locy, locx] = spot_values[i]
        else:    
            locx = (spot-1) % grid[0]
            locy = (spot-1) // grid[1]
            # print(i, spot, locx, locy)
            grid_array[locy, locx] = spot_values[i]

    ax.imshow(grid_array, cmap=cmap)#, vmin=vmin, vmax=vmax)
    # have the axis scaled
    # ax.axis('scaled')

    ax.set_xlim(0, grid[0])
    ax.set_ylim(0, grid[1])

    # invert the y axis
    ax.invert_yaxis()

    # add the colorbar
    cbar = plt.colorbar(ax.imshow(grid_array, cmap=cmap, ), ax=ax, label='Response',**kwargs)#vmin=vmin, vmax=vmax
    # add colorbar label
    # cbar.set_ylabel('Depolorization (pA)')
    
    ax.set_aspect(1/pattern_index.polygon_frame_properties['aspect_ratio'])

    return locx, locy, ax


def pairwise_draw_and_annotate_line_plot(ax, df, x='', y='', hue='', draw=True, kind='violin', palette='viridis', stat_across='hue', stat=stats.kruskal, skip_first_xvalue=True, annotate_wrt_data=False, offset_btw_star_n_line=0.1, color='grey', coord_system='data', fontsize=12, zorder=10):
    ''' This function takes a dataframe, and makes pairwise comparisons between the groups in the hue column
    for each x value. The function then annotates the line plot with the p-values of the comparisons.'''

    if draw:
        # draw the plots
        if kind == 'violin':
            sns.violinplot(data=df_melt, x=x, y=y, hue=hue, palette=palette, ax=ax, alpha=0.8, split=True, inner='quartile', linewidth=1)
        elif kind == 'strip':
            sns.stripplot(data=df_melt, x=x, y=y, hue=hue, palette=palette, ax=ax, alpha=0.8, dodge=1,)
        elif kind == 'line':
            sns.lineplot(data=df_melt, x=x, y=y, hue=hue, palette=palette, ax=ax, alpha=0.5, errorbar=('sd', 1), err_style='bars', linewidth=3,err_kws={"elinewidth": 3, 'capsize':5})


    hue_values = df[hue].unique() # group labels for each x-axis categorical value
    x_values = df[x].unique() # x-axis categorical value labels

    # get the xticks and xticklabels
    xticks = ax.get_xticks()
    xticklabels = ax.get_xticklabels()

    # get the max value of data across all x and all hue groups
    max_ydata = df[y].max()
    # set ypos to be 0.9*ylim
    ypos = 0.9*ax.get_ylim()[1]


    # for each x-value, get the ygroup values for hue1 and hue2
    for ix, x_val in enumerate(x_values):
        if skip_first_xvalue:
            if ix==0:
                continue
                    
        group_data = df_melt[(df_melt[x]==x_val)].groupby(hue)[y].apply(list)
        # convert all the group data into a list of lists
        group_data = group_data.values.tolist()
        kruskal_statistic, kruskal_pval = stats.kruskal(*group_data)

        # get the location of x_val on the x-axis of ax
        # get x-ticks and x-tick-labels
        xpos = xticks[ix]

        # get the maximum value of y for the given x_val across all the groups, add the offset to get the ypos for annotation
        if annotate_wrt_data:
            ypos = 1.1* np.max(group_data)

        # convert xpos and ypos into axes coordinate system if coord_system=='axes'
        if coord_system=='axes':
            xpos = ax.transAxes.inverted().transform(ax.transData.transform([xpos, ypos]))[0]
            ypos = ax.transAxes.inverted().transform(ax.transData.transform([xpos, ypos]))[1]
        

        
        annotate_stat_stars(ax, kruskal_pval, star_loc=[xpos, ypos], add_line=False, color=color, coord_system=coord_system, fontsize=12, zorder=10)

        # print(ix, x_val, kruskal_statistic, kruskal_pval, xpos, ypos)
    

