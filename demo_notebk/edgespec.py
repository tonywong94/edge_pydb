#!/usr/bin/env python -i

import numpy as np
import matplotlib.pyplot as plt
#import glob
import os
from astropy.io import ascii
from matplotlib.ticker import MultipleLocator

from fill_between_steps import *

set = raw_input("Enter det, mar, non, or all [all]: ") or 'all'
#type = raw_input("Enter dil, smo, mk2, or str [dil]: ") or 'dil'
dvmin = raw_input("Enter min vel window in km/s [500]: ") or 500.

gallist = np.genfromtxt('class_'+set+'.txt',dtype=None,names=['gname'])
vlims = np.genfromtxt('vlims.csv',delimiter=',',dtype=None,names=True,filling_values=-1)

if set == 'det' or set == 'all':
    nx=7
    ny=5
else:
    nx=4
    ny=3
pages = int(np.ceil(float(len(gallist['gname'])) / (nx*ny)))

for type in ('dil', 'smo', 'mk2', 'str'):
    for num in range(0,pages):
        aa = nx*ny*num
        bb = nx*ny+aa
        sublist = gallist['gname'][aa:bb].tolist()

        figure = plt.figure(0)
        figure.set_size_inches(nx*4.5, ny*4.)

        for gal in sublist:
            index = sublist.index(gal)
            row,col = divmod(index,nx)
            ax = plt.subplot2grid((ny,nx),(row,col))
            ax.axhline(y=0, xmin=0, xmax=1, color = 'black')

            fname = 'data/'+gal+'.co.de20_'+type+'.flux.out'
            if os.path.isfile(fname):
                data = ascii.read(fname)

                vel    = data['col1']
                flux   = data['col2']
                eflux  = data['col3']
                ceflux = data['col4']
                vlimidx = vlims['Galaxy'].tolist().index(gal)
                vsys = vlims['Vsys'][vlimidx]
                vmin = vlims['Vmin'][vlimidx]
                vmax = vlims['Vmax'][vlimidx]
                if vmin == -1 or vmax == -1:
                    vmin = vsys - 0.5*dvmin
                    vmax = vsys + 0.5*dvmin
                elif (vmax-vmin)<dvmin:
                    vmin2 = 0.5*(vmin+vmax-dvmin)
                    vmax2 = 0.5*(vmin+vmax+dvmin)
                    vmin = vmin2
                    vmax = vmax2
                deltv  = vel[1]-vel[0]
                print gal,vmin,vmax

                imax = max(1,np.ceil(0.25+max(2*flux)))
                ax.set_xlim(vel[0], vel[-1])
                ax.set_ylim([-0.15*imax, 0.5*imax])
                plt.hlines(-0.075*imax, xmin=vmin, xmax=vmax, linewidth=3, colors='r')
                plt.vlines(vsys, ymin=-0.1*imax, ymax=-0.05*imax, linewidth=3, colors='r')
                plt.tick_params(axis='x', labelsize='medium')
                plt.tick_params(axis='y', labelsize='medium')

                spec  = ax.step(vel, flux, color = 'b', where='mid')
                start, end = ax.get_xlim()
                if (end-start)>1000:
                    majorLocator = MultipleLocator(400)
                else:
                    majorLocator = MultipleLocator(200)
                ax.xaxis.set_major_locator(majorLocator)
                if type == 'dil' or type == 'smo':
                    fill_between_steps(vel,flux+eflux,y2=flux-eflux,color='white',
                        h_align='mid',ax=ax,facecolor='steelblue',alpha=0.6)
            ax.set_title(gal, fontsize='large')

    #    figure.text(0.5, 0.06, 'VELO-LSR (km/s)', ha='center',fontsize='xx-large')
    #    figure.text(0.1, 0.5, 'Flux(Jy)', va='center', rotation='vertical', fontsize='xx-large')

        plt.savefig(set+'_'+type+'_flux-'+str(num)+'.pdf', bbox_inches='tight')
        plt.close()

    os.system("pdfunite "+set+'_'+type+"_flux-*.pdf spec20_"+set+'_'+type+".pdf")
    os.system("rm -f "+set+'_'+type+"_flux-*.pdf")
