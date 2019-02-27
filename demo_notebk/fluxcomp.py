#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import leastsq
from astropy.table import Table

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

filename = raw_input("Enter name of csv file [edge_coflux_de20.csv]: "
    ) or 'edge_coflux_de20.csv'

data = Table.read(filename, format='ascii.ecsv')
fig=plt.figure()

xstart, xend = [-25, 500]

# Left panel: Unmasked vs Smoothed
ax1=plt.subplot(1,3,1,aspect='equal')
idc = np.where(np.minimum(data['coSmooth'],data['coNomask']) > xstart)[0]
plt.errorbar(data['coSmooth'], data['coNomask'], xerr=data['coeSmooth'], 
    yerr=data['coeNomask'], capsize=0, ls='none', color='blue', elinewidth=1,
    zorder=2)
plt.xlabel('Flux (smooth mask) [Jy km/s]', labelpad=3, size=8)
plt.ylabel('Flux (unmasked) [Jy km/s]', labelpad=2, size=8)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=8)
ax1.set_xlim(xstart, xend)
ax1.set_ylim(xstart, xend)
#ax1.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
ax1.text(0.04,0.9,'N='+str(len(idc)),fontsize=7,horizontalalignment='left',
    transform=ax1.transAxes)
xmod = np.arange(xstart,xend)
plt.plot(xmod,xmod,color='r',lw=1,zorder=1)
# Label selected galaxies
for xpos, ypos, name in zip(data['coSmooth'], data['coNomask'], data['Name']):
    if xpos > 250:
        plt.annotate(name, (xpos,ypos), xytext=(-3,4), size=5, 
            textcoords='offset points',horizontalalignment='right',zorder=3)

# Middle panel: 2D masked vs Smoothed
ax2=plt.subplot(1,3,2,aspect='equal')
idc = np.where(np.minimum(data['coSmooth'],data['coMask2d']) > xstart)[0]
plt.errorbar(data['coSmooth'], data['coMask2d'], xerr=data['coeSmooth'], 
    yerr=data['coeMask2d'], capsize=0, ls='none', color='blue', elinewidth=1,
    zorder=2)
plt.xlabel('Flux (smooth mask) [Jy km/s]', labelpad=2, size=8)
plt.ylabel('Flux (2D mask) [Jy km/s]', labelpad=3, size=8)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=8)
ax2.set_xlim(xstart, xend)
ax2.set_ylim(xstart, xend)
ax2.text(0.04,0.9,'N='+str(len(idc)),fontsize=7,horizontalalignment='left',
    transform=ax2.transAxes)
#ax2.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
plt.plot(xmod,xmod,color='r',lw=1,zorder=1)

xstart, xend = [1, 10**3]
xmod = np.arange(xstart,xend)

# Right panel: Dilated vs Smoothed, log scale
ax3=plt.subplot(1,3,3,aspect='equal')
idc = np.where(np.minimum(data['coSmooth'],data['coDilated']) > xstart)[0]
plt.errorbar(data['coSmooth'][idc], data['coDilated'][idc], xerr=data['coeSmooth'][idc], 
    yerr=data['coeDilated'][idc], capsize=0, ls='none', color='blue', elinewidth=1,
    zorder=2, marker='o', markersize=2)
plt.xlabel('Flux (smooth mask) [Jy km/s]', labelpad=3, size=8)
plt.ylabel('Flux (dilated mask) [Jy km/s]', labelpad=2, size=8)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=8)
ax3.set_xlim(xstart, xend)
ax3.set_ylim(xstart, xend)
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.text(0.04,0.9,'N='+str(len(idc)),fontsize=7,horizontalalignment='left',
    transform=ax3.transAxes)
#ax2.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
plt.plot(xmod,xmod,color='r',lw=1,zorder=1)
# for xpos, ypos, name in zip(data['coSmooth'], data['coDilated'], data['Name']):
#     if ypos < 7:
#         plt.annotate(name, (xpos,ypos), xytext=(2,-5), size=4, 
#             textcoords='offset points',horizontalalignment='middle',zorder=3)

plt.subplots_adjust(wspace=0.4)

plt.savefig('output/'+filename+'.pdf', bbox_inches='tight')
plt.close()
