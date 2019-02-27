#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import leastsq

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

filename = raw_input("Enter name of csv file [edge_coobs_de10.csv]: "
    ) or 'edge_coobs_de10.csv'

data=np.genfromtxt(filename,delimiter=',',skip_header=1,usecols=[2,4,5,6,7,9,10],
  dtype=float,names=['obstime','bmaj','bmin','rmsjy','rmsk','snrpk','flux'])

boot1=np.genfromtxt('bootfac_e.csv',delimiter=',',skip_header=1,usecols=[3],
  dtype=float,names=['bootfac1'])

boot2=np.genfromtxt('bootfac_d.csv',delimiter=',',skip_header=1,usecols=[3],
  dtype=float,names=['bootfac2'])

boot3=np.genfromtxt('bootfac_pre.csv',delimiter=',',skip_header=1,usecols=[3],
  dtype=float,names=['bootfac3'])

boot=np.concatenate([boot1['bootfac1'],boot2['bootfac2'],boot3['bootfac3']])

fig=plt.figure()

ax1=fig.add_subplot(2,3,1)
n,bins,patches=plt.hist(data['obstime'],20,histtype='bar')
plt.xlabel('obstime [h]', labelpad=3, size=11)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=9)
ax1.text(0.65,0.9,'$\mu$ = %4.2f' % np.mean(data['obstime']),size=10,transform=ax1.transAxes)

ax2=plt.subplot(2,3,2)
beam=np.sqrt(data['bmaj']*data['bmin'])
n,bins,patches=plt.hist(beam,20,histtype='bar')
plt.xlabel('beam [arcsec]', labelpad=3, size=11)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=9)
ax2.text(0.65,0.9,'$\mu$ = %4.2f' % np.mean(beam),size=10,transform=ax2.transAxes)

ax3=plt.subplot(2,3,3)
n,bins,patches=plt.hist(data['rmsjy'],20,histtype='bar')
plt.xlabel('RMS [mJy/bm]', labelpad=3, size=11)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=9)
ax3.text(0.65,0.9,'$\mu$ = %4.1f' % np.mean(data['rmsjy']),size=10,transform=ax3.transAxes)

ax4=plt.subplot(2,3,4)
n,bins,patches=plt.hist(data['rmsk'],20,histtype='bar')
plt.xlabel('RMS [mK]', labelpad=3, size=11)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=9)
ax4.text(0.65,0.9,'$\mu$ = %4.1f' % np.mean(data['rmsk']),size=10,transform=ax4.transAxes)

ax5=plt.subplot(2,3,5)
n,bins,patches=plt.hist(data['snrpk'],bins=np.logspace(0.,2,20),histtype='bar')
plt.xlabel('Peak SNR', labelpad=3, size=11)
plt.gca().set_xscale("log")
ax5.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x)))
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=9)
ax5.text(0.65,0.9,'$\mu$ = %4.1f' % np.mean(data['snrpk']),size=10,transform=ax5.transAxes)

ax6 = plt.subplot(2,3,6)
#fluxes=data['flux']/1000.
#fluxes=fluxes[~np.isnan(fluxes)]
n,bins,patches=plt.hist(boot,20,histtype='bar')
plt.xlabel('Bootstrap factor', labelpad=3, size=11)
plt.tick_params(axis='x', labelsize=8)
plt.tick_params(axis='y', labelsize=9)
# Gaussian fit
fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
errfunc  = lambda p, x, y: (y - fitfunc(p, x))
init  = [1.0, 1.0, 0.5]
xbins = 0.025+bins[:-1]
out   = leastsq( errfunc, init, args=(xbins, n))
c = out[0]
print "A exp[-0.5((x-mu)/sigma)^2]"
print "Fit Coefficients:"
print c[0],c[1],abs(c[2])
plt.plot(xbins, fitfunc(c, xbins))
ax6.text(0.65,0.91,'$\mu$ = %4.2f' % c[1],size=10,transform=ax6.transAxes)
ax6.text(0.65,0.85,'$\sigma$ = %4.2f' % abs(c[2]),size=10,transform=ax6.transAxes)

# ax = fig.add_axes( [0., 0., 1, 1] )
# ax.set_axis_off()
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
# ax.text(0.5, 0.95, filename, horizontalalignment='center', verticalalignment='center')

plt.savefig('output/obshist_six.pdf', bbox_inches='tight')
plt.close()
