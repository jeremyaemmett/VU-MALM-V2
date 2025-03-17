""" Diagnostic plot code, to ensure model functionality """

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LogFormatter
from datetime import datetime, timedelta
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import parameters as params
import conversions
import numpy as np
import output

# 3/17/2025


def forcing_plots(fd):

    """ Create 2D plots to represent the time-depth evolution of directly ingested (e.g. temperature, moisture)
     and calculated (e.g. gas equilibrium concentration, DOC, effective diffusivity) forcing fields
     If a given field differs between e.g. chem. species, a separate plot window is used to represent these in a group
    :param fd:
    :return:
    """

    # General forcing fields
    plots_per_row = 4
    fig, axes = plt.subplots()
    fig.patch.set_facecolor('ghostwhite')
    fig.set_size_inches(0.25 * plots_per_row * 45.0, 45.0, forward=True)
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0.25, hspace=0.3)
    keys = list(fd.keys())
    count = 0
    for i in range(0, len(keys)):
        key_dimensions = len(np.shape(fd[keys[i]]))
        if key_dimensions == 2:
            row, col = int(np.floor(count / plots_per_row)), count % plots_per_row
            ax = plt.subplot2grid((plots_per_row, plots_per_row), (row, col), colspan=1, rowspan=1)
            ylbl = 'Depths' if col == 0 else ' '
            title, cmap, col = conversions.key2plotparams(keys[i])
            profile_evolution(fd[keys[i]], 'mkm', fd['tem'], fd['sat'], fd['tim'], fd['dep'], fig, ax, ' ', ylbl,
                              title, 48, cmap, False)
            count += 1
    plt.savefig(params.main_directory + 'output/' + 'composition.png', bbox_inches='tight')
    plt.close(fig)

    # Create plots that show how chemical-dependent forcing parameters vary for each chemical
    plots_per_row = 3
    fn = ['diffusivities.png', 'diffusivities_full.png', 'equilibriums.png',
          'tempfactors.png', 'schmidtnumbers.png', 'transfervelocities.png']
    mk = ['efd', 'efdf', 'eqc', 'ftq', 'schm', 'tvel']
    for m in range(0, len(mk)):

        fig2, axes2 = plt.subplots()
        fig2.patch.set_facecolor('ghostwhite')
        fig2.set_size_inches(0.35 * plots_per_row * 45.0, 45.0, forward=True)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0.25, hspace=0.3)
        sub_keys = list(fd[mk[m]].keys())
        count = 0
        for i in range(0, len(sub_keys)):
            row, col = int(np.floor(count / plots_per_row)), count % plots_per_row
            ax = plt.subplot2grid((plots_per_row, plots_per_row), (row, col), colspan=1, rowspan=1)
            ylbl = 'Depths' if col == 0 else ' '
            title, cmap, col = conversions.key2plotparams(sub_keys[i])
            profile_evolution(fd[mk[m]][sub_keys[i]], mk[m], fd['tem'], fd['sat'], fd['tim'], fd['dep'], fig, ax, ' ',
                              ylbl, title, 48, cmap, False)
            count += 1
        plt.savefig(str(params.main_directory + 'output/' + fn[m]), bbox_inches='tight')
        plt.close(fig2)


def profile_evolution(profiles, mkm, temps, sats, times, depths, fig,
                      axis, xtitle, ytitle, title, fsize, cmap, logflag):

    """ Actual plot commands to visualize the forcing fields, referenced in the function 'forcing_plots' above
    :param profiles:
    :param mkm:
    :param temps:
    :param sats:
    :param times:
    :param depths:
    :param fig:
    :param axis:
    :param xtitle:
    :param ytitle:
    :param title:
    :param fsize:
    :param cmap:
    :param logflag:
    :return:
    """

    levels = np.linspace(np.nanmin(profiles), np.nanmax(profiles), 21)

    if mkm == 'efd' and np.max(profiles) > 0.0:
        levels = np.linspace(np.nanmin(profiles[np.nonzero(profiles)]), np.nanmax(profiles), 21)
    if title == 'sat': levels = np.linspace(0.0, np.nanmax(profiles), 21)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_title(title, fontsize=1.3 * fsize)
    axis.set_ylabel(xtitle, fontsize=fsize)
    axis.set_ylabel(ytitle, fontsize=fsize)
    axis.xaxis.set_tick_params(labelsize=fsize)
    axis.yaxis.set_tick_params(labelsize=fsize)

    axis.xaxis_date()
    date_format = mdates.DateFormatter('%b')
    axis.xaxis.set_major_formatter(date_format)
    # axis.set_xticklabels(['        May', '        Jun', '        Jul', '        Aug', '        Sep', ' ', ' '])
    start, end = datetime(2022, 6, 1), datetime(2022, 11, 1)
    axis.set_xlim([start, end])
    axis.set_ylim(0.63, min(depths))

    datetimes = [conversions.doy2datetime(2022, times[i]) for i in range(0, len(times))]
    X, Y = np.meshgrid(datetimes, depths)
    Z = profiles
    divider = make_axes_locatable(axis)
    cax = divider.append_axes('bottom', size='5%', pad=1.0)
    if levels[0] == 0.0: levels[0] = 1e-20
    if np.nanmax(levels) == np.nanmin(levels): levels = [0.0, np.nanmax(levels)]
    if logflag:
        c3 = axis.contourf(X, Y, Z, cmap=cmap, norm=LogNorm(vmin=levels[0], vmax=levels[-1]))
    else:
        if np.max(Z) != 0.0:
            c3 = axis.contourf(X, Y, Z, levels=levels, cmap=cmap)
            cbar = fig.colorbar(c3, orientation='horizontal', cax=cax)
            cbar.ax.tick_params(axis="both", labelsize=fsize, rotation=45)

    temp_mask = 1.0 + 0.0 * temps
    temp_mask[temps > 0.0] = np.nan
    axis.contourf(X, Y, temp_mask, colors='white', alpha=1.0)
    sats2 = np.copy(sats)
    sats2[temps <= 0] = np.nan
    c = axis.contour(X, Y, sats2, linewidths=25.0, levels=[0.0, 99.0], colors='blue', alpha=0.5)
    c = axis.contour(X, Y, sats2, linewidths=5.0, levels=[0.0, 99.0], colors='cyan', alpha=1.0)

    for i in range(0, len(depths)):
        lstyle, textps = ['-', '-'][i % 2], [start + timedelta(2.5), start + timedelta(8.0)][i % 2]
        axis.plot([start, end], [depths[i], depths[i]], color='white', linewidth=2, linestyle=lstyle)
        if depths[i] < 0.63: axis.text(textps, depths[i] + np.diff(depths)[i] / 1.5, i,
                                       fontsize=max(12, (i / len(depths)) * fsize), color='black')

    return ()