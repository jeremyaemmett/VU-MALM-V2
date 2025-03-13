""" Output plots """

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LogFormatter
from datetime import datetime, timedelta
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from scipy.stats import zscore
import parameters as params
import conversions
import numpy as np
import output

# 3/2/2025


def prepare_plot_grid(plots_per_row):

    fig, axes = plt.subplots()
    fig.patch.set_facecolor('ghostwhite')
    fig.set_size_inches(0.35 * plots_per_row * 45.0, 45.0, forward=True)
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0.25, hspace=0.3)

    return(fig, axes)


def format_subplot(axis, title, fsize, depths, start, end, n_dimensions):

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.set_facecolor('ghostwhite')
    date_format = mdates.DateFormatter('%b')
    axis.set_title(title, fontsize=fsize, fontweight="bold")
    axis.xaxis.set_tick_params(labelsize=fsize)
    axis.yaxis.set_tick_params(labelsize=fsize)
    #axis.grid()
    if n_dimensions == 2:
        axis.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
        axis.set_ylim(0.55, min(depths))
    axis.xaxis_date()
    axis.xaxis.set_major_formatter(date_format)
    axis.set_xlim([start, end])
    plt.minorticks_off()
    divider = make_axes_locatable(axis)
    cax = divider.append_axes('bottom', size='5%', pad=1.5)

    return cax


def plots(dict_name, dict_dict, fd, log_flag, fsize, n_dimensions, flip_y_flag, symmetric_y_flag, weekly_mean_flag):

    temps = fd['tem']

    plots_per_row = 3
    start, end = datetime(2022, 6, 1), datetime(2022, 11, 1)

    fig, axes = prepare_plot_grid(3)

    keys = list(dict_dict[dict_name].keys())
    count = 0
    profiles, depths, datetimes = output.read_output(params.main_directory + 'output/' + dict_name + '.txt', n_dimensions)
    for i in range(0, len(keys)):

        row, col = int(np.floor(count / plots_per_row)), count % plots_per_row
        axis = plt.subplot2grid((plots_per_row, plots_per_row), (row, col), colspan=1, rowspan=1)
        ylbl = 'Depths' if col == 0 else ' '
        title, cmap, clr = conversions.key2plotparams(keys[i])

        profile = profiles[keys[i]]

        cax = format_subplot(axis, title, fsize, depths, start, end, n_dimensions)

        if n_dimensions == 1:

            X = datetimes
            z = profile.flatten()

            z[abs(z) > 3.0 * np.std(z)] = 0.0
            #if weekly_mean_flag: X, z = conversions.resampleDTvals(X, z, '3D')

            maxspan = max(abs(np.min(z)), abs(np.max(z)))
            if flip_y_flag: axis.set_ylim([np.max(z), np.min(z)])
            if symmetric_y_flag:
                axis.set_ylim([-maxspan, maxspan])
                if flip_y_flag: axis.set_ylim([maxspan, -maxspan])
            c = axis.fill_between(X, z, color=clr, alpha=0.35)
            #c2 = axis.scatter(X, z, color=clr, alpha = 0.35)
            axis.plot([X[0], X[-1]], [0.0, 0.0], linestyle = '-', color = 'black')
            l_f = LogFormatter(10, labelOnlyBase=False)
            cbar = fig.colorbar(c, orientation='horizontal', cax=cax, format=l_f)
            cbar.remove()

        if n_dimensions == 2:

            X, Y = np.meshgrid(datetimes, depths)
            z = profile.transpose()
            temps = temps[:, 0:z.shape[1]]
            #z[temps <= 0.0] = np.nan

            if log_flag:
                l_f = LogFormatter(10, labelOnlyBase=False)
                minval, maxval = 3e-5, np.nanmax(profile[np.nonzero(profile)])
                minval, maxval = np.nanmin(profile[np.nonzero(profile)]), np.nanmax(profile[np.nonzero(profile)])
                # Catch the case where everything is 'zero' = 1e-10 (Python plays stupid with contour levels)
                lvls = [1e-20, 1e-10] if maxval < 2e-10 else np.logspace(np.log10(minval), np.log10(maxval), 100).tolist()
                c3 = axis.contourf(X, Y, z, levels=100, cmap=cmap)
                cbar = fig.colorbar(c3, orientation='horizontal', cax=cax, extendfrac=1.5)
                cbar.ax.tick_params(axis="both", labelsize=fsize, rotation=45, length=20, width = 5)

            for i in range(0, len(depths)):
                lstyle, textps = ['-', '-'][i % 2], [start + timedelta(2.5), start + timedelta(8.0)][i % 2]
                axis.plot([start, end], [depths[i], depths[i]], color='white', linewidth=0.5, linestyle = lstyle)
                if depths[i] < 0.63: axis.text(textps, depths[i] + np.diff(depths)[i] / 1.5, i,
                                                   fontsize = max(12, (i / len(depths)) * 36.0), color = 'black')

            # Ice boundary overlay
            temps = fd['tem']
            temp_times = fd['tim']
            temp_datetimes = [conversions.doy2datetime(2022, temp_times[i]) for i in range(0, len(temp_times))]
            X, Y = np.meshgrid(temp_datetimes, depths)
            temp_mask = 1.0 + 0.0 * temps
            temp_mask[temps > 0.0] = np.nan
            axis.contourf(X, Y, temp_mask, colors='white', alpha=1.0)
            c = axis.contour(X, Y, temps, linewidths=10.0, levels=[0.0], colors='gainsboro', alpha=1.0)
            # Water boundary overlay
            sats = fd['sat']
            sat_times = fd['tim']
            sat_datetimes = [conversions.doy2datetime(2022, sat_times[i]) for i in range(0, len(sat_times))]
            X, Y = np.meshgrid(sat_datetimes, depths)
            sats2 = np.copy(sats)
            sats2[temp_mask == 1] = np.nan
            c = axis.contour(X, Y, sats2, linewidths=25.0, levels=[0.0, 99.999],colors='blue',alpha=0.5)
            c = axis.contour(X, Y, sats2, linewidths=5.0, levels=[0.0, 99.999],colors='cyan',alpha=1.0)

            for j in range(0, len(depths)):
                lstyle, textps = ['-', '-'][j % 2], [start + timedelta(2.5), start + timedelta(8.0)][j % 2]
                axis.plot([start, end], [depths[j], depths[j]], color='black', linewidth=2, linestyle=lstyle)
                if depths[i] < 0.63: axis.text(textps, depths[j] + np.diff(depths)[j] / 1.5, i,
                                               fontsize=max(12, (j / len(depths)) * fsize), color='black')
        count += 1

    plt.savefig(params.main_directory + 'output/' + dict_name+ '.png', bbox_inches='tight')
