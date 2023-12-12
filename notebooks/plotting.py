"""
Created on Fri Mar 24 10:18:22 2023

@author: atownse2
"""
import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import hist as bh

# import seaborn as sns
# sns.set_palette('colorblind')

hep.style.use(hep.style.CMS)

fig_size = (8, 8)
cms_font_size = 20

def plot_stack(hists, labels, title=None,
               setlogy=False):
    fig, ax = plt.subplots()
    fig.set_size_inches(*fig_size)
    if setlogy:
        ax.set_yscale('log')
    hep.histplot(hists, ax=ax, histtype='fill', stack=True, label=labels)
    ax.set_title(title, fontsize=12)
    ax.legend(loc='upper left', frameon=False)
    hep.cms.label(ax=ax, year='2018')

def plot1d(h, sub_title=None,
           cms_label='Work in Progress', is_data=False,
           ylabel=None, setlogy=False, 
           xlim=None,
           fig=None, ax=None, legend=None, label=None,
           density=False,
           save_as=None,):
    if fig is None or ax is None:
        fig, ax = plt.subplots()
    fig.set_size_inches(*fig_size)
    if setlogy:
        ax.set_yscale('log')
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        xlim=get_1D_axis_lims(h)
        ax.set_xlim(xlim)
    if sub_title is not None:
        # place title in top middle under top axis
        ax.text(0.5, 0.9, sub_title, ha='center', va='bottom', transform=ax.transAxes, fontsize=16)
    hep.histplot(h, ax=ax, flow=None, label=label, density=density)
    # ax.errorbar(h.axes[0].centers, h.values(), yerr=np.sqrt(h.values()), fmt='o', markersize=4, linestyle='None', label=label)
    hep.cms.label(ax=ax, year='2018', label=cms_label, data=is_data, fontsize=cms_font_size)

    if legend is not None:
        ax.legend(loc=legend, fontsize=16)
    plt.tight_layout()
    if save_as is not None:
        plt.savefig(save_as)
    return fig, ax

def plot2d(h, title=None, 
            cmap='binary', xticks=None, yticks=None,
            display_bin_vales=False,
            xlim=None, ylim=None,
            setlogy=False, setlogx=False, setlogz=False,
            signal_point=None,
            is_data=False, cms_label='Work in Progress',
            save_as=None):
    
    if isinstance(h, bh.Hist):
        if 'dataset' in h.axes.name:
            h = h[{'dataset':sum}]

    fig, ax = plt.subplots()
    fig.set_size_inches(*fig_size)
    ax.set_title(title)
    if setlogy:
        ax.set_yscale('log')
    if setlogx:
        ax.set_xscale('log')
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)

    hep.hist2dplot(h, ax=ax, cmap=cmap, flow=None, norm=mpl.colors.LogNorm() if setlogz else None)
    if display_bin_vales:
        H, xedges, yedges = h.to_numpy()
        for i in range(len(xedges)-1):
            for j in range(len(yedges)-1):
                ax.text((xedges[i]+xedges[i+1])/2, (yedges[j]+yedges[j+1])/2, np.round(H[i,j], decimals=2), 
                        ha='center', va='center', color='black')
    hep.cms.label(ax=ax, year='2018', label=cms_label, data=is_data, fontsize=cms_font_size)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if signal_point is not None:
        m_Bkk, m_Radion = signal_point
        ax.text(0.05, 0.95, f'$m_{{Bkk}} = {m_Bkk}$ GeV\n$m_{{Radion}} = {m_Radion}$ GeV',
                ha='left', va='top', transform=ax.transAxes, fontsize=12)
        ax.set_xlim(m_Bkk-50, m_Bkk+50)
        ax.set_ylim(m_Radion-2, m_Radion+2)

    fig.tight_layout()
    if save_as is not None:
        plt.savefig(save_as)

def get_2D_axis_lims(hist):
    x_hist = hist[:,sum]
    y_hist = hist[sum,:]

    x_lim = get_1D_axis_lims(x_hist)
    y_lim = get_1D_axis_lims(y_hist)

    return x_lim, y_lim

def get_1D_axis_lims(hists):

    x_min = np.inf
    x_max = -np.inf

    if isinstance(hists, bh.Hist):
        hists = [hists]

    for hist in hists:
        bin_contents = hist.values()
        # Mask to ignore zero bins
        mask = bin_contents != 0

        # Apply the mask and find the minimum non-zero bin
        bin_indicies = np.arange(len(bin_contents))

        min_non_zero_bin_index = bin_indicies[mask][0]
        max_non_zero_bin_index = bin_indicies[mask][-1] + 3

        if max_non_zero_bin_index > len(bin_contents):
            max_non_zero_bin_index = len(bin_contents)

        x_min = min(x_min, hist.axes[0].edges[min_non_zero_bin_index])
        x_max = max(x_max, hist.axes[0].edges[max_non_zero_bin_index])
        

    xlim = (x_min, x_max)
    return xlim

def plot1ds(hists, labels, sub_title=None,
            cms_label='Work in Progress', is_data=False,
            setlogy=False, ylabel=None,
            legend='best',
            xlim=None, xticks=None, 
            density=False,
            save_as=None):
    fig, ax = plt.subplots()

    xlim=get_1D_axis_lims(hists)

    for hist, label in zip(hists,labels):
        plot1d(hist, label=label,
               sub_title=sub_title, cms_label=cms_label, 
               is_data=is_data, 
               setlogy=setlogy, 
               fig=fig, ax=ax,
               legend=legend,
               xlim=xlim,
               density=density)

def plot_residuals(hist, model, labels=None, title=None, setlogy=False, ylabel=None, xlim=None, xticks=None, save_as=None):
    # Plot histograms with residuals underneath
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    ax1.set_title(title)
    if setlogy:
        ax1.set_yscale('log')
    if ylabel is not None:
        ax1.set_ylabel(ylabel)

    x = np.linspace(0,4000,10000)

    hep.histplot(hist, ax=ax1, histtype='fill', stack=True)
    ax1.plot(x, model(x))

    # Residual
    hist_vals, bin_edges = hist.to_numpy()
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2

    model_vals = model(bin_centers)
    residuals = model_vals - hist_vals
    standardized_residuals = residuals / np.sqrt(hist_vals)

    ax2.errorbar(bin_centers, standardized_residuals, fmt='o', color='black', markersize=2, linestyle='None')
    # hep.histplot(hists, ax=ax2, histtype='errorbar', label=labels, stack=True, density=True)

    hep.cms.label(ax=ax1, year='2018')
    ax1.legend(loc='upper left', frameon=False)
    ax1.set_xlim(xlim)
    ax2.set_xlabel(None)
    ax2.set_ylabel('Residual')
    ax2.axhline(1, color='black', linestyle='--')
    if save_as is not None:
        plt.savefig(save_as)
