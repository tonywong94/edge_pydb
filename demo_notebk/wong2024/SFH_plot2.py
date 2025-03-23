#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import datetime
from edge_pydb import EdgeTable

# Plot sSFR profile for all galaxies

nx=1
ny=2

pdfname = 'SSFR'
sspcols   = ['Name','ix','iy','rad_arc','sigstar_sm']
fluxcols  = ['Name','ix','iy','flux_sigsfr_corr_sm','e_flux_sigsfr_corr_sm',
             'flux_sigsfr_adopt_sm']
sfrtab = EdgeTable('edge_carma_allpix.2d_smo7.hdf5',path='flux_elines_sm', cols=fluxcols)
ssptab = EdgeTable('edge_carma_allpix.2d_smo7.hdf5',path='SSP_sm', cols=sspcols)
sfrtab.join(ssptab, keys=['Name', 'ix', 'iy'])
califa = EdgeTable('edge_califa.csv')
califa.add_index('Name')

#gallist = ['NGC3687', 'NGC7738']
gallist = ['NGC3687', 'NGC5614']

print('Plotting SFH for',len(gallist),'galaxies')

# Get default axis limits
pages = int(np.ceil(float(len(gallist)) / (nx*ny)))
if pdfname is not None:
    pp = PdfPages(pdfname+'.pdf')

for num in range(0,pages):
    aa = nx*ny*num
    bb = nx*ny+aa
    thispage = gallist[aa:bb]
    fig = plt.figure(figsize=(3.5,7))
    print('Plotting', thispage[0], 'to', thispage[-1])

    for i in range(0, len(thispage)):
        gname = thispage[i]
        ax = plt.subplot(ny,nx,i+1)

        subtab = Table(sfrtab[sfrtab['Name'] == gname])
        r_re = Column(subtab['rad_arc']/califa.loc[gname]['caRe'], unit='')
        subtab.add_column(r_re, index=1, name='frac_rad')
        ssfr = Column(subtab['flux_sigsfr_adopt_sm']/subtab['sigstar_sm'], unit='1/Gyr')
        subtab.add_column(ssfr, name='ssfr')
        subtab.keep_columns(['Name', 'frac_rad', 'ssfr'])
        print(subtab)

        df = subtab.to_pandas()
        df.dropna(how='any', inplace=True)
        print('Min, max radius for',gname,'is',df['frac_rad'].min(),df['frac_rad'].max())

        # Plot the SSFR for binned ranges in radius.
        radcuts = pd.cut(df['frac_rad'],[0,0.5,1,1.5,2,2.5,3])
        df_radbin = df.groupby(radcuts).mean()
        print(df_radbin)
        colors = plt.cm.rainbow_r(np.linspace(0, 1, len(df_radbin)))
        df_radbin.plot(x=0, logx=False, legend=False, color='r', 
                         title=gname, ax=ax, style='.-')
        if i//nx < (len(thispage)-1)//nx:
            ax.set_xticklabels([])
            ax.xaxis.label.set_visible(False)
        if i % nx != 0:
            ax.set_yticklabels([])
        else:
            ax.set_xlabel("Radius / $R_e$")
            ax.set_ylabel("sSFR (1/Gyr)")

    fig.subplots_adjust(hspace=0.16)
    fig.subplots_adjust(wspace=0.05)
    if pdfname is not None:
        pp.savefig(bbox_inches = 'tight', pad_inches=0.1)
        plt.close()
    else:
        plt.show()

if pdfname is not None:
    d = pp.infodict()
    d['Title'] = 'EDGE Gallery'
    d['Author'] = os.getlogin()
    d['CreationDate'] = datetime.datetime.today()
    pp.close()

