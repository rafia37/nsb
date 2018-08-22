#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from astropy.io import ascii
from numpy.polynomial.polynomial import polyfit
from astropy.visualization import astropy_mpl_style
from datetime import datetime
import matplotlib.dates as mdates
from astropy.time import Time
from numpy.polynomial.polynomial import polyfit
plt.style.use(astropy_mpl_style)
#plt.style.use('ggplot')
#plt.rc('font', family='serif', serif='Times')
#plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', titlesize=14, labelsize=12)


SQM_TEL_FILE = 'nsb_sqm_tel.txt'
SQM_ZENITH_FILE = 'nsb_sqm_zenith.txt'
#NSB_FILE = 'V_data.csv'


def zenith_plot(data):
    ax = plt.subplot(111)
    x = Time(data['sqm_ut']).plot_date
    plt.plot_date(x, data['sqm_nsb'], 'k-')
    fmt = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(fmt)
    plt.gcf().autofmt_xdate()
    plt.xlabel('Time (UT)')
    plt.ylabel(r'Night Sky Brightness (mag/arcsec$^2$)')
    plt.title('SQM NSB Zenith Measurements for 2018-05-08 UT')
    plt.ylim(20.9,21.5)
    plt.savefig('zenith_sqm_nsb.png', dpi=300)
    plt.clf()

def polar_plot(data, format):
    if format == 'Pomenis':
        color = 'viridis_r'
        filter_used = data['filt'][0]
    elif format == 'SQM':
        color = 'viridis_r'
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    colors = data['sqm_nsb']
    cm = plt.cm.get_cmap(color)

    nsb = ax.scatter(np.radians(data['az']),90-data['elv'], marker='o',edgecolors='k', facecolors='none', c=colors, cmap=cm, s=100)
    #plt.plot(data['az'], 90-data['elv'], 'k.')
    ax.set_yticks(range(0, 90+10, 10))                   # Define the yticks
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    ax.yaxis.grid(True, color='k')
    ax.xaxis.grid(True, color='k')
    #ax.set_xlabel('Azimuth (degrees)')
    #ax.set_ylabel('Elevation (degrees)')
    plt.ylim(0,90)
    plt.title(r'{} NSB Measurements for 2018-05-08 UT'.format(format))
    cbar = plt.colorbar(nsb)
    #cbar.set_ytickslabels(np.arange(15.5,22, .5))
    cbar.set_label(r'mag/arcsec$^2$')
    plt.savefig('{}_nsb.png'.format(format.lower()), dpi=300)
    plt.clf()

def experimental_plot(data, format):

    if format == 'Pomenis':
      color = 'viridis_r'
    elif format == 'SQM':
      color = 'magma_r'

    az = np.radians(data['az'])
    elv = data['elv']
    value = data['nsb']

    x = elv * np.cos(az)
    y = elv * np.sin(az)

    xi = np.linspace(0, 90, 100)
    yi = np.linspace(np.radians(0), np.radians(360.), 20)
    zi = griddata((x, y), value, (xi[None,:], yi[:,None]), method='cubic')


    # try converting xi and yi to polar, i think they in rect

    ax = plt.subplot(111, projection='polar')
    cm = plt.cm.get_cmap(color)
    levels = np.linspace(min(value),max(value), 100)
    cf = plt.contourf(yi,90-xi,zi.T,levels=levels,cmap=cm)
    cbar = plt.colorbar()
    cbar.set_ticks(np.arange(15.5,22, .5))
    cbar.set_label(r'mag/arcsec$^2$')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90+10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    #ax.yaxis.grid(True, alpha=0.5, color='k')
    ax.xaxis.grid(False)
    ax.scatter(az, 90-elv, marker='.', edgecolors='k', facecolor='none')

    plt.ylim(0,90)
    plt.title(
      r'{} INTERPOLATED Night Sky Brightness UT20180508'.format(format)
    )
    plt.savefig('{}_interpolated_nsb.png'.format(format.lower()), dpi=300)
    plt.clf()
    # try converting xi and yi to polar, i think they in rect

"""
    ax = plt.subplot(111, projection='polar')
    cm = plt.cm.get_cmap(color)
    levels = np.linspace(min(value),max(value), 100)
    cf = plt.contourf(theta,r,zi,cmap=cm)
    cbar = plt.colorbar()
    cbar.set_ticks(np.arange(15.5,22, .5))
    cbar.set_label(r'mag/arcsec$^2$')
    #ax.set_theta_zero_location('N')
    #ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90+10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    #ax.yaxis.grid(True, alpha=0.5, color='k')
    ax.xaxis.grid(False)
    ax.scatter(az, 90-elv, marker='.', edgecolors='k', facecolor='none')

    plt.ylim(0,90)
    plt.title(
        r'{} INTERPOLATED Night Sky Brightness UT20180508'.format(format)
    )
    plt.savefig('{}_interpolated_nsb.png'.format(format.lower()), dpi=300)
    plt.clf()
"""
"""
    xi = np.linspace(0, 90, 5)
    yi = np.radians(np.linspace(0, 360., 20))
    zi = griddata((x, y), value, (xi[None,:], yi[:,None]), method='cubic')


    # try converting xi and yi to polar, i think they in rect

    ax = plt.subplot(111, projection='polar')
    cm = plt.cm.get_cmap(color)

    cf = plt.contourf(xi,yi,zi.T,cmap=cm)
    cbar = plt.colorbar()
    cbar.set_ticks(np.arange(15.5,22, .5))
    cbar.set_label(r'mag/arcsec$^2$')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90+10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    #ax.yaxis.grid(True, alpha=0.5, color='k')
    ax.xaxis.grid(False)
    ax.scatter(az, 90-elv, marker='.', edgecolors='k', facecolor='none')

    plt.ylim(0,90)
    plt.title(
        r'{} INTERPOLATED Night Sky Brightness UT20180508'.format(format)
    )
    #plt.savefig('{}_interpolated_nsb.png'.format(format.lower()), dpi=300)
    plt.show()
    plt.clf()
"""
def plot_extinction(data):
    ax = plt.subplot(111)
    x = (1./(np.cos(np.radians(90-data['elv']))))
    plt.plot(x, data['zp'], 'k.')
    b, m = polyfit(x, data['zp'], 1)
    plt.plot(x, b+m*x, 'r-')
    #print(m)
    plt.xlabel('Airmass')
    plt.ylabel('Zeropoint (mag)')
    plt.title('Extinction for 2018-05-08 UT')
    #plt.savefig('zenith_sqm_nsb.pdf', dpi=300)
    #plt.show()
    plt.clf()

if __name__ == '__main__':

    sqm_tel_data = ascii.read(SQM_TEL_FILE)
    sqm_zenith_data = ascii.read(SQM_ZENITH_FILE)
    #nsb_data = ascii.read(NSB_FILE)

    #plot_extinction(nsb_data)
    #polar_plot(nsb_data, 'Pomenis')

    polar_plot(sqm_tel_data, 'SQM')

    #experimental_plot(nsb_data, 'Pomenis')

    #experimental_plot(sqm_tel_data, 'SQM')

    zenith_plot(sqm_zenith_data)




#BACKUP
"""
alt = []
az = []
value = []
with open('nsb_data.txt', 'r') as f:
    for line in f:
        if line[0] == '#':
            pass
        else:
            line = line.strip('\n').split('\t')
            az.append(np.radians(float(line[2])))
            #alt.append(90 - float(line[1]))
            #az.append(float(line[0]))
            alt.append(float(line[3]))
            value.append(float(line[6]))

# r, theta
# r = alt
# theta = az
# x = r x cos(theta)
# y = r x sin(theta)
alt = np.array(alt)
x,y = alt * np.cos(az), alt*np.sin(az)
xi = np.linspace(0, 90, 100)
yi = np.linspace(np.radians(10.), np.radians(350.), 100)
zi = griddata((x, y), value, (xi[None,:], yi[:,None]), method='cubic')


ax = plt.subplot(111, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
cm = plt.cm.get_cmap('jet_r')
levels = np.linspace(min(value),max(value), 50)
cf = plt.contourf(yi,90-xi,zi.T,levels=levels,cmap=cm)

cbar = plt.colorbar()
cbar.set_ticks(np.arange(15.5,22, .5))
cbar.set_label(r'mag/arcsec$^2$')
ax.set_yticks(range(0, 90+10, 10))                   # Define the yticks
yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
ax.set_yticklabels(yLabel)
#ax.grid(True)
ax.yaxis.grid(True, alpha=0.5, color='k')
ax.xaxis.grid(False)
#ax.set_rticks([])
ax.scatter(az,90-alt, marker='.', edgecolors='k', facecolors='none')
#ax.xaxis.grid(False)
#ax.set_rmin(0)
plt.ylim(0,90)
plt.title(r'Night Sky Brightness UT 20180508')
plt.show()
"""
