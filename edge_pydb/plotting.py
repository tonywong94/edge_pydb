import os
import datetime
import numpy as np
from scipy import stats
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Prepare a 2D histogram from a scatterplot
def xy2hist(xarr, yarr, log=True, bins=[100,100]):
    if log:
        x = np.log10(xarr)
        y = np.log10(yarr)
    else:
        x = xarr
        y = yarr
    # Histogram the data
    # https://stackoverflow.com/questions/49662964/density-scatter-plot-for-huge-dataset-in-matplotlib
    hh, locx, locy = np.histogram2d(x, y, bins=bins)
    # Get the bin value for each point
    z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    if log:
        z = np.log10(z)
    return x, y, z, hh, locx, locy


# Prepare binned averages from a scatterplot
# If log=True, averaging is done after taking the log
def xy2binned(xarr, yarr, log=True, bins=20, range=None, yval='mean'):
    if log:
        good = (xarr>0) & (yarr>0)
        x = np.log10(xarr[good])
        y = np.log10(yarr[good])
    else:
        x = xarr
        y = yarr
    if yval == 'median':
        stat = lambda y: np.nanmedian(y)
        estat = lambda y: stats.iqr(y, scale='normal', nan_policy='omit')
    else:
        stat = lambda y: np.nanmean(y)
        estat = lambda y: np.nanstd(y)
    if bins > 0:
        ymean, xbinedge, _ = stats.binned_statistic(x, y,
            statistic=stat, bins=bins, range=range)
        ystd, xbinedge, _  = stats.binned_statistic(x, y,
            statistic=estat, bins=bins, range=range)
        xbin = 0.5*(xbinedge[1:]+xbinedge[:-1])
        return xbin, ymean, ystd
    else:
        return [0], [0], [0]


# Prepare a patch collection for a dotplot
# def dotpatch(x, y, colors, size=1, vmin=None, vmax=None, cmap='jet'):
# 
#     patches = []
#     for x1, y1 in zip(x,y):
#         circle = Circle((x1, y1), size)
#         patches.append(circle)
# 
#     p = PatchCollection(patches, cmap=cmap)
#     p.set_array(np.array(colors))
#     if vmax is None:
#         vmax = np.nanmax(colors)
#     if vmin is None:
#         vmin = np.nanmin(colors)
# 
#     p.set_clim([vmin, vmax])
#     xymin = np.min([x.min(),y.min()])
#     xymax = np.max([x.max(),y.max()])
# 
#     return p, xymin, xymax


# Prepare a patch collection for a dotplot
def dotpatch(x, y, imval, blank=None, dotsize=1, axes=None, **kwargs):
    # blank is a boolean array that sets certain pixels to NaN
    if blank is not None:
        imval[blank] = np.nan
    if axes is None:
        axes = plt.gca()
    if 'vmin' in kwargs:
        vmin = kwargs.pop('vmin')
    else:
        vmin = np.nanmin(imval)
    if 'vmax' in kwargs:
        vmax = kwargs.pop('vmax')
    else:
        vmax = np.nanmax(imval)
    patches = []
    for x1, y1 in zip(x,y):
        circle = Circle((x1, y1), dotsize)
        patches.append(circle)
    p = PatchCollection(patches, **kwargs)
    p.set_array(np.array(imval))
    p.set_clim([vmin, vmax])
    xminmax = [x.min(), x.max()]
    yminmax = [y.min(), y.max()]
    img = axes.add_collection(p)
    return img, xminmax, yminmax


# Plot a non-dotpatch image
def imarrayplot(x, y, imval, blank=None, axes=None, **kwargs):
    if blank is not None:
        imval[blank] = np.nan
    if axes is None:
        axes = plt.gca()
    xdim = len(np.unique(x))
    ydim = len(np.unique(y))
    imarray = np.reshape(imval, [ydim,xdim], order='F')
    xminmax = [0, xdim-1]
    yminmax = [0, ydim-1]
    img = axes.imshow(imarray, origin='lower', **kwargs)
    return img, xminmax, yminmax


# Multi-page plots with all galaxies
def gridplot(edgetab=None, gallist=None, column='flux_Halpha_rg', 
            xrange=None, yrange=None, blank=None, plotstyle='image',
            cmap='jet', nx=7, ny=6, dotsize=1, pdfname=None, **kwargs):

    # Plot all galaxies by default
    if gallist is None:
        gallist = list(np.unique(edgetab['Name']))

    pages = int(np.ceil(float(len(gallist)) / (nx*ny)))
    if pdfname is not None:
        pp = PdfPages(pdfname)

    for num in range(0,pages):
        aa = nx*ny*num
        bb = nx*ny+aa
        thispage = gallist[aa:bb]
        fig = plt.figure(figsize=(20,14))
        print('Plotting', thispage[0], 'to', thispage[-1])

        for i in range(0,len(thispage)):
            gname = thispage[i]
            ax = plt.subplot(ny,nx,i+1)
            galtab = edgetab['Name'] == gname

            if plotstyle == 'dot':
                if blank is not None:
                    img, xran, yran = dotpatch(edgetab[galtab]['ix'], 
                                               edgetab[galtab]['iy'],
                                               edgetab[galtab][column], 
                                               blank=blank[galtab], 
                                               dotsize=dotsize, cmap=cmap, 
                                               axes=ax, **kwargs)
                else:
                    img, xran, yran = dotpatch(edgetab[galtab]['ix'], 
                                               edgetab[galtab]['iy'],
                                               edgetab[galtab][column],
                                               dotsize=dotsize, cmap=cmap, 
                                               axes=ax, **kwargs)
            else:
                if blank is not None:
                    img, xran, yran = imarrayplot(edgetab[galtab]['ix'], 
                                                  edgetab[galtab]['iy'],
                                                  edgetab[galtab][column], 
                                                  blank=blank[galtab], 
                                                  cmap=cmap, axes=ax, **kwargs)        
                else:
                    img, xran, yran = imarrayplot(edgetab[galtab]['ix'], 
                                                  edgetab[galtab]['iy'],
                                                  edgetab[galtab][column],
                                                  cmap=cmap, axes=ax, **kwargs)        

            if xrange is None:
                ax.set_xlim(xran)
                if i == 0:
                    print("Default x limits used:",xran)
            else:
                ax.set_xlim(xrange)
            if yrange is None:
                ax.set_ylim(yran)
                if i == 0:
                    print("Default y limits used:",yran)
            else:
                ax.set_ylim(yrange)

            ax.set_aspect('equal')
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.text(0.04,0.9,gname,ha='left',va='center',transform=ax.transAxes,
               bbox=dict(facecolor='white', edgecolor='none', pad=1))
        fig.subplots_adjust(hspace=0.05)
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
    return
