#!/usr/bin/env python

# Plot the radial profiles from the original rprof.txt files.  These should have
# been generated in the subdirectory 'rproftxt' with blanksig=0 and blanksig=1.  
# Horizontal axis is radius in arcsec, kpc, or normalized to $R_{25}$.  Vertical
# axis is average CO intensity or average H$_2$ surface density including helium.  
# Green shading separates the blanksig=0 and blanksig=1 curves; the blanksig=0 curve 
# has the triangles plotted and is the basis for the CO effective radius shown as the 
# vertical red line.  The conservative 3$\sigma$ detection limit is shown as the dashed 
# black curve.  For the surface density vs. kpc plot, the slanted purple line gives the 
# fitted CO scale length and normalization from Alberto's fits in edge_rdist.csv.

import datetime, os, glob
import numpy as np
from astropy.table import Table
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# Function to profile for one galaxy
def pltcoprof(gal, ax, rtype='arcsec', itype='wtmean', imtype='de20_smo', 
              txtdir='rproftxt', ebar=True, mod0=None, modh=None):
    print('{0}...'.format(gal), end='')
    b0file = txtdir + '/' + gal + '.' + imtype + 'b0.rprof.txt'
    b1file = txtdir + '/' + gal + '.' + imtype + 'b1.rprof.txt'
    alphaco = 4.3 * u.solMass * u.s / (u.K * u.km * u.pc**2)
    if rtype == 'kpc':
        rmax = 16.
    elif rtype == 'r25':
        rmax = 1.05
    else:
        rmax = 60.
    if os.path.isfile(b0file):
        b0dat = Table.read(b0file, format='ascii.ecsv')
        b1dat = Table.read(b1file, format='ascii.ecsv')
        if rtype == 'kpc':
            rad = b0dat['r_kpc']
        elif rtype == 'r25':
            rad = b0dat['r_r25']
        else:
            rad = b0dat['radius']
        plt.plot(rad, b0dat[itype], color='b', ls='-', marker='^', ms=7)
        plt.plot(rad, b1dat[itype], color='g', ls='--', marker=None)
        plt.fill_between(rad,b0dat[itype], b1dat[itype], 
            facecolor='green', alpha=0.5)
        if ebar == True:
            plt.errorbar(rad, b0dat[itype], yerr=b0dat['rms'], 
                ecolor='dimgray', capsize=0, 
                zorder=1, marker=None, ls='None', lw=1, label=None)
        if (itype == 'sigmol' and rtype == 'kpc' and np.isfinite(mod0) and modh > 0):
            print('Scale length {0} with normalization {1}'.format(modh,mod0))
            imod = mod0 * np.exp(-rad/modh)
            plt.plot(rad, imod, color='m', ls='-', marker=None, lw=3)
        rco_eff = np.interp(0.5*b0dat['cummass'][-1], b0dat['cummass'], rad)
        if itype == 'wtmean':
            plt.plot(rad, b0dat['detlim'], color='k', lw=2, ls='--', marker=None)
        else:
            plt.plot(rad, b0dat['detlim']*alphaco, color='k', lw=2, ls='--', marker=None)
        ax.axvline(x=rco_eff, lw=4, color='r', alpha=0.5)
    ax.set_yscale('log')
    if itype == 'wtmean':
        ax.set_ylim(10**-1.5, 10**3)
    else:
        ax.set_ylim(10**-0.5, 10**4)
    ax.set_xlim(0,rmax)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(14)
    ax.text(0.95, 0.9, gal, horizontalalignment='right', fontsize=16,
        verticalalignment='center', transform=ax.transAxes)
    return


# Make the plots
dbdir  = '../../../'
rdist = Table.read(dbdir+'dat_glob/external/edge_rdist.csv', format='ascii.ecsv')
if not os.path.isdir('pdf'):
    os.makedirs('pdf')

nx=7
ny=5

dotypes = ['de20_smo', 'smo7_smo']
for imtype in dotypes:

    # From the list of galaxies determine the number of pages required
    flist = glob.glob('rproftxt/*.'+imtype+'b0.rprof.txt')
    glist = [ os.path.basename(x).split('.')[0] for x in flist ]
    pages = int(np.ceil(float(len(glist)) / (nx*ny)))

    for rtype in ['arcsec', 'r25', 'kpc']:
        for itype in ['wtmean', 'sigmol']:
            if rtype == 'r25':
                rlabel = 'Radius / '+r'$R_{25}$'
            else:
                rlabel = 'Radius ['+rtype+']'
            if itype == 'sigmol':
                prefix = 'sig'
                ilabel = r'log $\Sigma_{mol}$ [$M_\odot$ pc$^{-2}$]'
            else:
                prefix = 'i'
                ilabel = r'log $I_{CO}$ [K km/s]'
            pp = PdfPages('pdf/'+prefix+'prof_'+imtype+'_'+rtype+'.pdf')

            for num in range(0,pages):
                print('\n{} vs {}, Page {} of {}'.format(itype,rtype,num+1,pages))
                aa = nx*ny*num
                bb = nx*ny+aa
                gals = glist[aa:bb]
                nrows = int(np.ceil(float(len(gals)) / nx))
                figure = plt.figure(0)
                figure.set_size_inches(nx*4.5, ny*4.1)
                for i, gal in enumerate(gals):
                    idx = np.where(rdist['Name']==gal)[0][0]
                    row,col = divmod(i,nx)
                    ax = plt.subplot2grid((ny,nx),(row,col))
                    pltcoprof(gal, ax, rtype=rtype, itype=itype, imtype=imtype, 
                              ebar=False, mod0=rdist['rdNormMol'][idx], 
                              modh=rdist['rdScaleMol'][idx])
                    if col != 0:
                        ax.set_yticklabels([])
                    if row != nrows-1:
                        ax.set_xticklabels([])
                figure.subplots_adjust(hspace=0.1)
                figure.subplots_adjust(wspace=0.1)
                #if num == 3:
                #    figure.text(0.5, 0.03, rlabel, ha='center', fontsize=24)
                #else:
                figure.text(0.5, 0.06, rlabel, ha='center', fontsize=24)
                figure.text(0.09, 0.5, ilabel, va='center', 
                    rotation='vertical', fontsize=24)
                pp.savefig(bbox_inches = 'tight', pad_inches=0.1)
                plt.close()

            d = pp.infodict()
            d['Title'] = 'EDGE Radial Profiles'
            d['Author'] = 'Tony Wong'
            d['CreationDate'] = datetime.datetime.today()
            pp.close()    

