import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
from mpl_toolkits.mplot3d import Axes3D



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




"""alt = np.array(alt)
x,y = alt * np.cos(az), alt*np.sin(az)
xi = np.linspace(0, 90, 100)
yi = np.linspace(np.radians(10.), np.radians(350.), 100)
zi = griddata((x, y), value, (xi[None,:], yi[:,None]), method='cubic')
ax = plt.subplot(111, projection='polar')
cm = plt.cm.get_cmap('jet_r')
levels = np.linspace(min(value),max(value), 100)
cf = plt.contourf(yi,90-xi,zi.T,levels=levels,cmap=cm)

cbar = plt.colorbar()
cbar.set_ticks(np.arange(15.5,22, .5))
cbar.set_label(r'mag/arcsec$^2$')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.yaxis.grid(False)
#ax.set_rticks([])
ax.scatter(az,90-alt, marker='.', edgecolors='k', facecolors='none')
ax.xaxis.grid(False)
#ax.set_rmin(0)
plt.ylim(0,90)

plt.show()

plt.clf

ax = plt.subplot(111, projection='polar')
ax.pcolormesh(xi,yi,zi)
ax.scatter(az,alt, marker='.', edgecolors='k', facecolors='none')
plt.show()"""
