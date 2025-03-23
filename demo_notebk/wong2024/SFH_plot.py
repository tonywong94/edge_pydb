#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import datetime
from edge_pydb import EdgeTable

# Plot mass-based and luminosity-based SFH for all galaxies

nx=1
ny=2

pdfname = 'SFH'
sfh = EdgeTable('edge_carma_allpix.2d_smo7.hdf5',path='SFH_sm')
califa = EdgeTable('edge_califa.csv')
califa.add_index('Name')

#gallist = ['NGC3687', 'NGC7738']
gallist = ['NGC3687', 'NGC5614']

for ftype in ['mass']:
    #gallist = list(np.unique(sfh['Name']))
    print('Plotting SFH for',len(gallist),'galaxies')

    # Get default axis limits
    pages = int(np.ceil(float(len(gallist)) / (nx*ny)))
    if pdfname is not None:
        pp = PdfPages(pdfname+'_'+ftype+'.pdf')

    for num in range(0,pages):
        aa = nx*ny*num
        bb = nx*ny+aa
        thispage = gallist[aa:bb]
        fig = plt.figure(figsize=(4.5,7))
        print('Plotting', thispage[0], 'to', thispage[-1])

        for i in range(0, len(thispage)):
            gname = thispage[i]
            ax = plt.subplot(ny,nx,i+1)

            # Keep only mass/lum fraction columns and galactocentric radius
            if ftype == 'mass':
                subtab = Table(sfh[sfh['Name'] == gname].columns[9+43+43:9+43+43+39])
            else:
                subtab = Table(sfh[sfh['Name'] == gname].columns[9:9+39])
            subtab.add_column(sfh[sfh['Name'] == gname]['rad_arc'], index=0)
            r_re = Column(subtab['rad_arc']/califa.loc[gname]['caRe'], unit='')
            subtab.add_column(r_re, index=1, name='frac_rad')

            df = subtab.to_pandas()
            df.dropna(how='all', subset=df.columns[2:], inplace=True)
            df.columns = df.columns.str.replace(ftype+'frac_age_', '')
            df.columns = df.columns.str.replace('_sm', '')
            print(df.columns)
            print('Min, max radius for',gname,'is',df['frac_rad'].min(),df['frac_rad'].max())
            df_cum = pd.concat([df[['frac_rad']], df.iloc[:,2:].cumsum(axis=1)],axis=1)

            # Plot the cumulative SFH vs time for binned ranges in radius.
            radcuts = pd.cut(df_cum['frac_rad'],[0,0.5,1,1.5,2,2.5,3])
            #radcuts = pd.cut(df_cum['frac_rad'],[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2])
            df_radbin = df_cum.groupby(radcuts).mean()
            radtab = df_radbin.drop(columns=['frac_rad']).T.reset_index(names='Age')
            radtab['Age'] = radtab['Age'].astype(float)
            colors = plt.cm.rainbow_r(np.linspace(0, 1, len(df_radbin)))
            radtab.plot(x=0, logx=False, legend=False, color=colors, 
                             title=gname, ax=ax)
            ax.axhline(0.5, ls='--', alpha=0.5, color='k', zorder=-2)
            
            if i == 0:
                ax.legend(fontsize='small', loc='lower right')
            if i//nx < (len(thispage)-1)//nx:
                ax.set_xticklabels([])
                ax.xaxis.label.set_visible(False)
            if i % nx != 0:
                ax.set_yticklabels([])
            else:
                ax.set_xlabel("Age [Gyr]")
                ax.set_ylabel(ftype+" fraction < Age")

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

