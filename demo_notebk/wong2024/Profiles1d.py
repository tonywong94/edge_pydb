#!/usr/bin/env python

import numpy as np 
import datetime
import matplotlib.pyplot as plt 
from astropy.table import Table, join, unique
from matplotlib.ticker import MultipleLocator, LogFormatterExponent
from edge_pydb import EdgeTable
from matplotlib.backends.backend_pdf import PdfPages
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# Display EDGE CO rotation curves at native resolution

bbtable = EdgeTable('bb_natv_fitvd_dilmsk.csv')
bblist  = unique(bbtable.table, keys='Name').columns[0].tolist()
print(len(bblist),' galaxies in bbtable1')

rftable = EdgeTable('rf_CO_natv.csv')
rflist = unique(rftable.table, keys='Name').columns[0].tolist()
print(len(rflist),' galaxies in rftable')
rfpars = EdgeTable('edge_rfpars.csv')
rftable.join(rfpars)

jamtable = EdgeTable('jam_rotcurves.csv')
jamlist = unique(jamtable.table, keys='Name').columns[0].tolist()
print(len(jamlist),' galaxies in jamtable')
jampars = EdgeTable('edge_jampars.csv')
jamtable.join(jampars)

gallist = sorted(jamlist)
print(len(gallist),' galaxies in total')

nx = 6
ny = 9
pages = int(np.ceil(float(len(gallist)) / (nx*ny)))
pp = PdfPages('rotcurves.pdf')

for num in range(0,pages):
    aa = nx*ny*num
    bb = nx*ny+aa
    sublist = gallist[aa:bb]
    nrows = int(np.ceil(float(len(sublist)) / nx))
    fig = plt.figure(0)
    fig.set_size_inches(nx*3, ny*2.5)
    for idx, gal in enumerate(sublist):
        bbrows = bbtable[bbtable['Name']==gal]
        if gal in ['NGC3994','NGC4149']:
            bbrows = bbrows[:-1]
        # Shift to top of bin
        bbrows['radius'] = bbrows['bbRad']+1.5
        rfrows = rftable[rftable['Name']==gal]
        rfrows.sort('radius')
        jamrows = jamtable[jamtable['Name']==gal]
        jamrows.sort('radius')
        row, col = divmod(idx,nx)
        ax = plt.subplot2grid((ny,nx),(row,col))
        # --- for scaling VROT by sin(i)
        scl = [1, 1, 1]
        for i, vec in enumerate([bbrows['bbInc'], rfrows['rfInc'], jamrows['jaIncl']]):
            if len(vec) > 0:
                scl[i] = np.sin(np.radians(vec[0]))
                #vec *= scl[i]
        # --- for choosing the y-limits
        ymax = 0
        for i, vec in enumerate([bbrows['bbVrot'], rfrows['Vrot'], jamrows['Vrot']]):
            if len(vec) > 0:
                ymax = max(ymax, vec.max()*scl[i])
        ymax = min(700,100*np.ceil(0.6+ymax/100))-5
        # --- plot the rotation curves
        ax.errorbar(bbrows['radius'],bbrows['bbVrot']*scl[0],
                    yerr=[abs(bbrows['bbVrot_e1']*scl[0]),bbrows['bbVrot_e2']*scl[0]],
                    c='b',marker='o',label='Bbarolo',zorder=3)
        ax.errorbar(rfrows['radius'],rfrows['Vrot']*scl[1],yerr=rfrows['e_Vrot']*scl[1],
                    c='r',marker='s',label='Levy+18',zorder=2)
        ax.errorbar(jamrows['radius'],jamrows['Vrot']*scl[2],yerr=jamrows['e_Vrot']*scl[2],
                    c='orange',lw=6,elinewidth=2,alpha=0.5,label='Leung+18',zorder=1)
        #ax.errorbar(jamrows['radius'],jamrows['Vrot_bsc'],yerr=jamrows['e_Vrot_bsc'],
        #            c='green',lw=6,alpha=0.5,label='Leung+18bsc')
        ax.set_xlim([0, 45])
        ax.set_ylim([0, ymax])
        ax.yaxis.set_major_locator(MultipleLocator(100))
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(12)
        if row != nrows-1:
            ax.set_xticklabels([])
        if row == 0 and col == 0:
            ax.legend(loc='lower right', fontsize=10)
        ax.text(0.95, 0.9, gal, ha='right', va='center', transform=ax.transAxes,
            fontsize=15)

    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.25)
    fig.text(0.5, 0.08, 'Radius [arcsec]', ha='center', fontsize=22)
    fig.text(0.08, 0.5, '$V_{rot}$ sin(i) [km/s]', va='center', rotation='vertical', 
                fontsize=22)
    pp.savefig(bbox_inches = 'tight', pad_inches=0.1)
    plt.close()

d = pp.infodict()
d['Title'] = 'EDGE Rotation Curves'
d['Author'] = 'Tony Wong'
d['CreationDate'] = datetime.datetime.today()
pp.close()


# Plot the radial CO profiles
rptable = EdgeTable('rprof_smo7_smo.csv')

pages = int(np.ceil(float(len(gallist)) / (nx*ny)))
pp = PdfPages('coprofiles.pdf')

for num in range(0,pages):
    aa = nx*ny*num
    bb = nx*ny+aa
    gals = gallist[aa:bb]
    nrows = int(np.ceil(float(len(gals)) / nx))
    fig = plt.figure(0)
    fig.set_size_inches(nx*3, ny*2.5)
    for i, gal in enumerate(gals):
        idx = np.where(rptable['Name']==gal)[0]
        row,col = divmod(i,nx)
        ax = plt.subplot2grid((ny,nx),(row,col))
        xdata = rptable['r_arcs'][idx]
        ydata = rptable['ico_avg'][idx]
        y_err = rptable['ico_rms'][idx]
        y_lim = rptable['detlim'][idx]
        good = (ydata > y_lim)
        plt.errorbar(xdata[good], ydata[good], yerr=y_err[good], color='b', 
             marker='^', ms=7, ecolor='dimgray', capsize=0, zorder=1, ls='None', lw=1)
        plt.errorbar(xdata[~good], ydata[~good], yerr=y_err[~good], color='darkgray', 
             marker='^', ms=7, ecolor='silver', capsize=0, zorder=1, ls='None', lw=1)
        plt.plot(xdata, y_lim, color='k', lw=2, ls='--', marker=None)
        ax.set_yscale('log')
        ax.set_ylim(10**-2.25, 10**2.75)
        ax.yaxis.set_major_formatter(LogFormatterExponent(labelOnlyBase=True))
        ax.minorticks_off()
        ax.set_xlim(0, 55)
        ax.text(0.95, 0.9, gal, ha='right', va='center', transform=ax.transAxes,
            fontsize=15)
        if col != 0:
            ax.set_yticklabels([])
        if row != nrows-1:
            ax.set_xticklabels([])
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(14)
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.1)
    fig.text(0.5, 0.08, 'Radius [arcsec]', ha='center', fontsize=22)
    fig.text(0.08, 0.5, r'log $I_{CO}$ [K km/s]', va='center', 
        rotation='vertical', fontsize=22)
    pp.savefig(bbox_inches = 'tight', pad_inches=0.1)
    plt.close()

d = pp.infodict()
d['Title'] = 'EDGE Radial Profiles'
d['Author'] = 'Tony Wong'
d['CreationDate'] = datetime.datetime.today()
pp.close()    
