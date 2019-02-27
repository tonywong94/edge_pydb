#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.table import Table, join

def plttwofiles(file1=None, file2=None, outfile=None, dolog=False):
    data1=Table.read(file1,format='ascii.ecsv')
    data2=Table.read(file2,format='ascii.ecsv')
    alldata=join(data1,data2,keys='Name')
    rootname1=os.path.splitext(os.path.basename(file1))[0]
    rootname2=os.path.splitext(os.path.basename(file2))[0]
    fig=plt.figure(figsize=(8,11))
    xstart, xend = [-100, 500]

    # Dilated mask
    ax1=plt.subplot(2,2,1,aspect='equal')
    plt.errorbar(alldata['coDilated_1'], alldata['coDilated_2'], 
        xerr=alldata['coeDilated_1'], yerr=alldata['coeDilated_2'], 
        capsize=0, ls='none', color='blue', elinewidth=1)
    plt.xlabel('Dilated file1', labelpad=3, size=11)
    plt.ylabel('Dilated file2', labelpad=1, size=11)
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    if dolog == True:
        ax1.set_xscale("log")
        ax1.set_yscale("log")
    #ax1.set_xlim(xstart, xend)
    #ax1.set_ylim(xstart, xend)
    #ax1.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
    xmod = np.arange(xstart,xend)
    plt.plot(xmod,xmod,color='r')

    # Smoothed mask
    ax2=plt.subplot(2,2,2,aspect='equal')
    plt.errorbar(alldata['coSmooth_1'], alldata['coSmooth_2'], 
        xerr=alldata['coeSmooth_1'], yerr=alldata['coeSmooth_2'], 
        capsize=0, ls='none', color='blue', elinewidth=1)
    plt.xlabel('Smoothed file1', labelpad=3, size=11)
    plt.ylabel('Smoothed file2', labelpad=1, size=11)
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    if dolog == True:
        ax2.set_xscale("log")
        ax2.set_yscale("log")
    #ax2.set_xlim(xstart, xend)
    #ax2.set_ylim(xstart, xend)
    #ax2.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
    plt.plot(xmod,xmod,color='r')

    # Projected smoothed mask
    ax3=plt.subplot(2,2,3,aspect='equal')
    plt.errorbar(alldata['coMask2d_1'], alldata['coMask2d_2'], 
        xerr=alldata['coeMask2d_1'], yerr=alldata['coeMask2d_2'], 
        capsize=0, ls='none', color='blue', elinewidth=1)
    plt.xlabel('Smoothed 2D mask file1', labelpad=3, size=11)
    plt.ylabel('Smoothed 2D mask file2', labelpad=1, size=11)
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    if dolog == True:
        ax3.set_xscale("log")
        ax3.set_yscale("log")
    #ax2.set_xlim(xstart, xend)
    #ax2.set_ylim(xstart, xend)
    #ax2.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
    plt.plot(xmod,xmod,color='r')

    # No mask
    ax3=plt.subplot(2,2,4,aspect='equal')
    plt.errorbar(alldata['coNomask_1'], alldata['coNomask_2'], 
        xerr=alldata['coeNomask_1'], yerr=alldata['coeNomask_2'], 
        capsize=0, ls='none', color='blue', elinewidth=1)
    plt.xlabel('Unmasked file1', labelpad=3, size=11)
    plt.ylabel('Unmasked file2', labelpad=1, size=11)
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    if dolog == True:
        ax3.set_xscale("log")
        ax3.set_yscale("log")
    #ax2.set_xlim(xstart, xend)
    #ax2.set_ylim(xstart, xend)
    #ax2.xaxis.set_ticks(np.arange(xstart, xend+1, 100))
    plt.plot(xmod,xmod,color='r')

    plt.tight_layout()
    ax = fig.add_axes( [0., 0., 1, 1] )
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.05, 0.95, 'File 1: '+os.path.basename(file1))
    ax.text(0.05, 0.92, 'File 2: '+os.path.basename(file2))
    plt.savefig('output/'+outfile)
    plt.close()
    return

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

plttwofiles(file1='edge_coflux_de20.csv', file2='edge_coflux_e20.csv', 
    outfile='de20_e20_fluxcomp.pdf', dolog=False)
plttwofiles(file1='edge_coflux_de20.csv', file2='edge_coflux_e20.csv', 
    outfile='de20_e20_fluxcomp_log.pdf', dolog=True)
plttwofiles(file1='edge_coflux_de20.csv', file2='edge_coflux_de10.csv', 
    outfile='de20_de10_fluxcomp.pdf', dolog=False)
plttwofiles(file1='edge_coflux_de20.csv', file2='edge_coflux_de20_nosenmsk.csv', 
    outfile='de20_senmsk_fluxcomp.pdf', dolog=False)
plttwofiles(file1='edge_coflux_de20.csv', file2='edge_coflux_de20_smo7.csv', 
    outfile='de20_smo7_fluxcomp.pdf', dolog=False)

