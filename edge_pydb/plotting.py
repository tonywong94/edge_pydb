import os
import datetime
import numpy as np
from scipy import stats
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from edge_pydb.conversion import kewley01, kauffm03, cidfer10


def xy2hist(xarr, yarr, log=True, bins=[100,100]):
    '''
    Prepare a 2D density histogram from a scatterplot.
    Based on a response to a Stack Overflow question:
    stackoverflow.com/questions/49662964/density-scatter-plot-for-huge-dataset-in-matplotlib

    === Parameters ===
    xarr : numpy.array
        The x values in the scatter plot
    yarr : numpy.array
        The y values in the scatter plot
    log : boolean
        If True, take the log of input x and y, and of output density z
    bins : list of int
        The number of bins in x and y for the histogram

    === Returns ===
    x : numpy.array
        x values sorted by increasing density
    y : numpy.array
        y values sorted by increasing density
    z : numpy.array
        histogram count for each point
    hh : numpy.ndarray
        2D histogram values
    locx : numpy.array
        bin edges along first dimension
    locy : numpy.array
        bin edges along first dimension
    '''
    if log:
        x = np.log10(xarr)
        y = np.log10(yarr)
    else:
        x = xarr
        y = yarr
    # Histogram the data
    hh, locx, locy = np.histogram2d(x, y, bins=bins)
    # Get the bin value z for each (x,y) point
    z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
    # Plot the highest density points last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    if log:
        z = np.log10(z)
    return x, y, z, hh, locx, locy


def xy2binned(xarr, yarr, log=True, bins=20, range=None, yval='mean'):
    '''
    Prepare binned averages from a scatterplot.
    If log=True, averaging is done after taking the log of x and y.

    === Parameters ===
    xarr : numpy.array
        The x values in the scatter plot
    yarr : numpy.array
        The y values in the scatter plot
    log : boolean
        If True, take the log of input x and y
    bins : int
        The number of bins in x for the averaging
    range : (float, float)
        The lower and upper range of the bins. If not provided, uses (x.min(), x.max())
    yval : string
        'mean' (default) or 'median'

    === Returns ===
    xbin : numpy.array
        the centers (NOT edges) of the x bins
    ymean : numpy.array
        the binned y-values
    ystd : numpy.array
        standard deviation of y-values
    '''
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


def dotpatch(x, y, imval, blank=None, dotsize=1, axes=None, **kwargs):
    '''
    Generate and plot a patch collection for a dot plot.

    === Parameters ===
    x : numpy.array
        The x values in the scatter plot (typically 'ix')
    y : numpy.array
        The y values in the scatter plot (typically 'iy')
    imval : numpy.array
        The z values which determine the dot colors
    blank : boolean array
        True values in this array are set to NaN
    dotsize : float
        The size of the filled circles
    axes : matplotlib.axes
        Axes for plotting
    **kwargs :
        Additional arguments including vmin, vmax, colormap normalization

    === Returns ===
    img : matplotlib image object
        Plotted image
    xminmax : tuple
        [x.min(), x.max()]
    yminmax : tuple
        [y.min(), y.max()]
    '''
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


def imarrayplot(x, y, imval, blank=None, axes=None, **kwargs):
    '''
    Plot a pixel image from a data column.

    === Parameters ===
    x : numpy.array
        The x values in the image plot (typically 'ix')
    y : numpy.array
        The y values in the image plot (typically 'iy')
    imval : numpy.array
        The z values which determine the dot colors
    blank : boolean array
        True values in this array are set to NaN
    axes : matplotlib.axes
        Axes for plotting
    **kwargs :
        Additional arguments including vmin, vmax, colormap normalization

    === Returns ===
    img : matplotlib image object
        Plotted image
    xminmax : tuple
        [x.min(), x.max()]
    yminmax : tuple
        [y.min(), y.max()]
    '''
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


def gridplot(edgetab=None, gallist=None, column='flux_Halpha_rg', 
            xrange=None, yrange=None, blank=None, plotstyle='image',
            cmap='jet', nx=7, ny=6, dotsize=1, pdfname=None, **kwargs):
    '''
    Plot multiple galaxies on a grid.

    === Parameters ===
    edgetab : EdgeTable
        Table containing the galaxies and data to plot
    gallist : list of strings
        List of galaxy names; default is to plot all galaxies
    column : string
        Name of column in the table that has the image data
    xrange : tuple of float
        x limits applied to each panel (pixels)
    yrange : tuple of float
        y limits applied to each panel (pixels)
    blank : boolean array with same length as 'column'
        True values in this array are set to NaN
    plotstyle : string
        'dot' for dot plot, 'image' for pixel image
    cmap : string
        name of color map
    nx : int
        number of sub-panels in x direction
    ny : int
        number of sub-panels in y direction
    dotsize : float
        size of plot symbol for dot plot
    pdfname : string
        name of output PDF file, otherwise plot to screen
    **kwargs :
        Additional arguments including vmin, vmax, colormap normalization
    '''

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
                    img, xlims, ylims = dotpatch(edgetab[galtab]['ix'], 
                                               edgetab[galtab]['iy'],
                                               edgetab[galtab][column], 
                                               blank=blank[galtab], 
                                               dotsize=dotsize, cmap=cmap, 
                                               axes=ax, **kwargs)
                else:
                    img, xlims, ylims = dotpatch(edgetab[galtab]['ix'], 
                                               edgetab[galtab]['iy'],
                                               edgetab[galtab][column],
                                               dotsize=dotsize, cmap=cmap, 
                                               axes=ax, **kwargs)
            else:
                if blank is not None:
                    img, xlims, ylims = imarrayplot(edgetab[galtab]['ix'], 
                                                  edgetab[galtab]['iy'],
                                                  edgetab[galtab][column], 
                                                  blank=blank[galtab], 
                                                  cmap=cmap, axes=ax, **kwargs)        
                else:
                    img, xlims, ylims = imarrayplot(edgetab[galtab]['ix'], 
                                                  edgetab[galtab]['iy'],
                                                  edgetab[galtab][column],
                                                  cmap=cmap, axes=ax, **kwargs)        

            if xrange is None:
                ax.set_xlim(xlims)
                if i == 0:
                    print("Default x limits used:",xlims)
            else:
                ax.set_xlim(xrange)
            if yrange is None:
                ax.set_ylim(ylims)
                if i == 0:
                    print("Default y limits used:",ylims)
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


def plot_uncertainty_ellipse(xval_n, xval_s, yval_n, yval_s, indices, x_arr, save_to=''):
    '''
    parameters
    xval_n, yval_u : list of nominal values of coordinates (pixels) 
    xval_s, yval_s : list of standard deviation of coordinates (pixels)
    indices : indices of the list of coordinates to plot with
    save_to: file to save the plot to, optional
    '''
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    plt.figure(figsize=(8,8))
    ax = plt.gca()
    for i in indices:
        ax.add_patch(Ellipse(xy=(xval_n[i], yval_n[i]),
                            width=xval_s[i], height=yval_s[i],
                            edgecolor='red', facecolor='none'))
    plt.plot(x_arr, kewley01(x_arr), 'k-.', label="Kewley")
    plt.plot(x_arr, kauffm03(x_arr), 'k--', label="Kauffmann")
    plt.plot(x_arr, cidfer10(x_arr), 'k-', label="Cidfer")
    if save_to:
        plt.savefig(save_to)
    plt.show()
