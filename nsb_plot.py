#!/usr/bin/env python3

import sys
import argparse
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import griddata, Rbf
import matplotlib.tri as tri
from astropy.time import Time
import matplotlib.dates as mdates
import pdb


# setup plotting
plt.rc('text', usetex=False) # change to False for no latex
plt.rc('font', family='serif')
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)


def interpolate_plot(files, min_mag, max_mag, autoscale=False, sqm=False):
    """
    Interpolate over the data and create a contourf plot. We must
    convert to cartesian, create a grid using np.griddata, then convert
    back to polar coordinates before plotting.

    Parameters
    ----------
    files : the list of files

    Return
    ------
    None
    """
    for file in files:

        #data = pd.read_csv(file).dropna(subset = ['nsb']) #Drops rows where the 'nsb' column has nans

        data = ascii.read(file)

        original_values = np.array(data['nsb']) #nsbvalues from original data table that may contain NaNs
        values = original_values[np.isfinite(original_values)]
        azimuth = np.array(data['az'][np.isfinite(original_values)])
        elv = np.array(data['elv'][np.isfinite(original_values)])

        if autoscale:
            min_mag = values.min()
            max_mag = values.max()

        theta = np.radians(azimuth)
        r = elv
        z = values

        x = r * np.cos(theta)
        y = r * np.sin(theta)


        xgrid = np.linspace(x.min(), x.max(), len(x))
        ygrid = np.linspace(y.min(), y.max(), len(y))

        xgrid, ygrid = np.meshgrid(xgrid, ygrid)
        zgrid = griddata((x, y), z, (xgrid, ygrid), method='cubic')


        #zgrid = ma.masked_invalid(zgrid)

        r = np.sqrt(xgrid**2 + ygrid**2)
        theta = np.arctan2(ygrid, xgrid)

        fig = plt.figure(figsize=(8,7))
        ax = plt.subplot(111, projection='polar')
        #ax.set_facecolor('k')
        colors = plt.cm.get_cmap('jet_r')
        levels, steps = np.linspace(min_mag, max_mag, 50, retstep=True)
        ticks = np.linspace(min_mag, max_mag, 9)
        cax = ax.contourf(theta, 90-r, zgrid, levels=levels, cmap=colors)

        cbar = plt.colorbar(cax, fraction=0.046, pad=0.04, ticks=ticks)
        cbar.set_label(r'mag/arcsec$^2$')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_yticks(range(0, 90+10, 10))
        ax.set_xticklabels([r'N', r'NE', r'E', r'SE', r'S', r'SW', r'W', r'NW'])
        yLabel = [r'90$^\circ$', '', '', r'60$^\circ$', '', '', r'30$^\circ$', '', '', '']
        ax.set_yticklabels(yLabel)
        #ax.yaxis.grid(False)
        #ax.xaxis.grid(False)
        ax.grid(alpha=0.3)

        if sqm:
            date = data['sqm_ut'][0].split('T')[0]
            title = 'Telescope SQM Night Sky Brightness'
            filter_used = 'sqm'
        else:
            date = data['ut'][0].split('T')[0]
            filter_used = data['filt'][0]
            if filter_used != 'g' and filter_used != 'z' and filter_used != 'w':
                filter_used = filter_used.upper()
            title = "{}' Filter Night Sky Brightness".format(filter_used)

        plt.title(r"{} {} UT".format(
            title,
            date
        ))

        plt.tight_layout()
        plt.ylim(0,90)
        circle = plt.Circle((0,0), 15, transform=ax.transData._b, color='white')

        ax.add_artist(circle)

        plt.savefig('{}_interpolated_plot.png'.format(filter_used), dpi=600)

def plot(files, min_mag, max_mag, autoscale=False, sqm=False):
    """
    Creates a polar plot using triangle contour plotting. This does not
    require the data to be in a grid format or interpolated. Use this
    plot to show just the raw data. It does not interpolate between
    data points.

    Parameters
    ----------
    files : the list of files

    Return
    ------
    None
    """
    for file in files:

        #data = pd.read_csv(file).dropna(subset = ['nsb']) #Drops rows where the 'nsb' column has nans

        data = ascii.read(file)

        #Grab needed info
        original_values = np.array(data['nsb']) #nsbvalues from original data table that may contain NaNs
        values = original_values[np.isfinite(original_values)]
        azimuth = np.array(data['az'][np.isfinite(original_values)])
        elv = np.array(data['elv'][np.isfinite(original_values)])

        if autoscale:
            min_mag = values.min()
            max_mag = values.max()

        theta = np.radians(azimuth)
        r = elv
        z = values
        fig = plt.figure(figsize=(8,7))
        ax = fig.add_subplot(111, projection='polar')
        #ax.set_facecolor('k')
        colors = plt.cm.get_cmap('jet_r')
        levels, steps = np.linspace(min_mag, max_mag, 35, retstep=True)
        ticks = np.linspace(min_mag, max_mag, 9)
        cax = ax.tricontourf(theta, 90-r, z, levels=levels, cmap=colors)

        cbar = plt.colorbar(cax, fraction=0.046, pad=0.04, ticks=ticks)
        cbar.set_label(r'mag/arcsec$^2$')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_yticks(range(0, 90+10, 10))
        ax.set_xticklabels([r'N', r'NE', r'E', r'SE', r'S', r'SW', r'W', r'NW'])
        yLabel = [r'90$^\circ$', '', '', r'60$^\circ$', '', '', r'30$^\circ$', '', '', '']
        ax.set_yticklabels(yLabel)
        #ax.yaxis.grid(False)
        #ax.xaxis.grid(False)
        ax.grid(alpha=0.3)

        if sqm:
            date = data['sqm_ut'][0].split('T')[0]
            title = 'Telescope SQM Night Sky Brightness'
            filter_used = 'tel_sqm'
        else:
            date = data['ut'][0].split('T')[0]
            filter_used = data['filt'][0]
            if filter_used != 'g' and filter_used != 'z' and filter_used != 'w':
                filter_used = filter_used.upper()

            title = "{}' Filter Night Sky Brightness".format(filter_used)



        plt.title(r"{} {} UT".format(
            title,
            date
        ))

        plt.tight_layout()
        plt.ylim(0,90)
        plt.savefig('{}_plot.png'.format(filter_used), dpi=600)

    return None

def plot_zenith(files, autoscale=True):
    """

    """
    for file in files:
        #data = pd.read_csv(file).dropna(subset = ['nsb']) #Drops rows where the 'nsb' column has nans

        data = ascii.read(file)

        original_values = np.array(data['nsb']) #nsbvalues from original data table that may contain NaNs
        values = original_values[np.isfinite(original_values)]
        date_title = data['sqm_ut'][0].split('T')[0]
        date = Time(data['sqm_ut']).plot_date[np.isfinite(original_values)]
        time_sorted_date = np.sort(date)
        time_sorted_values = values[np.argsort(date)]


        if autoscale:
            min_mag = values.min() - .5
            max_mag = values.max() + .5

        f, ax = plt.subplots(1,1, figsize=(8,8))
        plt.plot_date(date, values, 'k-')

        fmt = mdates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(fmt)
        plt.gcf().autofmt_xdate()
        plt.xlabel(r'Time (UT)')
        plt.ylabel(r'Night Sky Brightness (mag/arcsec$^2$)')
        plt.title(r'Zenith SQM Night Sky Brightness {} UT'.format(date_title))
        plt.ylim(min_mag, max_mag)
        ax.grid(alpha=0.3)

        plt.savefig('zenith_sqm.png', dpi=600)
        plt.clf()

def parse_arguments():
    """
    Create a parser to get the file from command line.

    Parameters
    ----------
    None

    Return
    ------
    args.files : the files in list form
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n',
        '--nsb_data',
        nargs=1,
        help='nsb data file',
    )
    parser.add_argument(
        '-a',
        '--autoscale',
        action='store_true',
        help='autoscale plot by min/max magnitude'
    )
    parser.add_argument(
        '-t',
        '--telescope_sqm',
        nargs=1,
        help='plot the telescope sqm file'
    )
    parser.add_argument(
        '-z',
        '--zenith_sqm',
        nargs=1,
        help='plot the zenith sqm file'
    )

    args = parser.parse_args()

    return args.nsb_data, args.autoscale, args.telescope_sqm, args.zenith_sqm

def main():
    """
    Main function.
    """
    print('\n########## nsb_plot.py')
    print('Creates interpolated plots and polar plots for the nsb data')
    nsb_data, autoscale, telescope_sqm, zenith_sqm = parse_arguments()

    min_mag = 18.
    max_mag = 22.

    if nsb_data:
        print('\n#---------- Plotting Pomenis measurements')
        print('Creating polar plot')
        plot(nsb_data, min_mag, max_mag, autoscale)
        print('\nDone...')

        print('\nCreating interpolated polar plot')
        interpolate_plot(nsb_data, min_mag, max_mag, autoscale)
        print('\nDone...')

    if telescope_sqm:
        print('\n#--------- Plotting telescope SQM measurements')
        print('Creating polar plot')
        plot(telescope_sqm, min_mag, max_mag, autoscale, sqm=True)
        print('\nDone')
        print('\nCreating interpolated polar plot')
        interpolate_plot(telescope_sqm, min_mag, max_mag, autoscale, sqm=True)
        print('\nDone...')

    if zenith_sqm:
        print('\n#--------- Plotting zenith SQM measurements')
        print('Creating zenith sqm line plot')
        plot_zenith(zenith_sqm)
        print('\nDone...')

    print('\nFINISHED\n')
    sys.exit(0)

if __name__ == '__main__':
    main()

##### DUMP
    """
    def tri_data(files):
        for file in files:

            data = ascii.read(file)

            azimuth = np.array(data['az'])
            elv = np.array(data['elv'])
            values = np.array(data['nsb'])

            theta = np.radians(azimuth)
            r = elv
            z = values

            x = r * np.cos(theta)
            y = r * np.sin(theta)

            xgrid = np.linspace(x.min(), x.max(), len(x))
            ygrid = np.linspace(y.min(), y.max(), len(y))

            triang = tri.Triangulation(x, y)
            interpolator = tri.LinearTriInterpolator(triang, z)

            Xi, Yi = np.meshgrid(xgrid, ygrid)
            zi = interpolator(Xi, Yi)


            r = np.sqrt(Xi**2 + Yi**2)
            #r = ma.masked_where(r < 15, r)
            theta = np.arctan2(Yi, Xi)


            ax1 = plt.subplot(111, projection='polar')
            #ax1.contour(theta, 90-r, zi, 14, linewidths=0.5, colors='k')
            cax = ax1.contourf(theta, 90-r, zi, 14, cmap='jet_r')
            cbar = plt.colorbar(cax)
            #cbar.set_clim(19, 22)
            cbar.set_label(r'mag/arcsec$^2$')
            ax1.set_theta_zero_location('N')
            ax1.set_theta_direction(-1)
            plt.ylim(0,90)

            plt.show()



    def rbf_data(files):



        for file in files:
            data = ascii.read(file)


            azimuth = np.array(data['az'])
            elv = np.array(data['elv'])
            values = np.array(data['nsb'])

            theta = np.radians(azimuth)
            r = elv
            z = values

            x = r * np.cos(theta)
            y = r * np.sin(theta)


            rgrid = np.linspace(0, 90, 100)
            tgrid = np.linspace(np.radians(0), np.radians(360), 100)
            xg = rgrid * np.cos(tgrid)
            tg = rgrid * np.sin(tgrid)

            xgrid = np.linspace(xg.min(), xg.max(), len(xg))
            ygrid = np.linspace(tg.min(), tg.max(), len(tg))

            xgrid, ygrid = np.meshgrid(xgrid, ygrid)

            rbf = Rbf(x, y, z, function='multiquadric', smooth=2)
            ZI = rbf(xgrid, ygrid)

            ri = np.sqrt(xgrid*xgrid + ygrid*ygrid)
            ti = np.arctan2(ygrid, xgrid)

            # polar plot
            fig = plt.figure()
            ax = plt.subplot(111, polar=True)
            levels, steps = np.linspace(18.375, 22.375, 35, retstep=True)
            #cax = ax.contour(ti, ri, zi, 10, linewidths=0.5, colors='k')


            cax = ax.contourf(ti, 90-ri, ZI, levels=levels, cmap=plt.cm.cubehelix_r)
            #ax.set_rmax(2)
            plt.ylim(0,90)

            plt.show()
    """
