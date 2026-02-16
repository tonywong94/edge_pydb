import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.patches import Circle, Ellipse
from matplotlib.collections import PatchCollection
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from edge_pydb.conversion import kewley01, kauffm03, cidfer10
from astropy.visualization import PercentileInterval, ImageNormalize
from astropy.visualization import LinearStretch, SqrtStretch, LogStretch

'''
Utility functions for plotting.
'''

def xy2hist(xarr, yarr, log=True, bins=[100,100]):
    '''
    Prepare a 2D density histogram from a scatterplot.
    Based on a response to a Stack Overflow question:
    stackoverflow.com/questions/49662964/density-scatter-plot-for-huge-dataset-in-matplotlib

    Parameters
    ----------
    xarr : numpy.array
        The x values in the scatter plot
    yarr : numpy.array
        The y values in the scatter plot
    log : boolean
        If True, take the log of input x and y, and of output density z
    bins : list of int
        The number of bins in x and y for the histogram

    Returns
    -------
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
        bin edges along second dimension
    '''
    if log:
        valid = (xarr > 0) & (yarr > 0)
        x = np.log10(xarr[valid])
        y = np.log10(yarr[valid])
    else:
        valid = np.isfinite(xarr+yarr)
        x = xarr[valid]
        y = yarr[valid]
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

    Parameters
    ----------
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

    Returns
    -------
    xbin : numpy.array
        the centers (NOT edges) of the x bins
    ymean : numpy.array
        the binned y-values
    ystd : numpy.array
        standard deviation of y-values
    ycnt : numpy.array
        number of y-values in each bin
    '''
    if log:
        good = (xarr>0) & (yarr>0)
        x = np.log10(xarr[good])
        y = np.log10(yarr[good])
    else:
        x = xarr
        y = yarr
    if yval == 'median':
        stat  = lambda y: np.nanmedian(y)
        estat = lambda y: stats.iqr(y, scale='normal', nan_policy='omit')
    else:
        stat  = lambda y: np.nanmean(y)
        estat = lambda y: np.nanstd(y)
    if bins > 0:
        ymean, xbinedge, _ = stats.binned_statistic(x, y,
            statistic=stat, bins=bins, range=range)
        ystd, xbinedge, _  = stats.binned_statistic(x, y,
            statistic=estat, bins=bins, range=range)
        ycnt, xbinedge, _  = stats.binned_statistic(x, y,
            statistic='count', bins=bins, range=range)
        xbin = 0.5*(xbinedge[1:]+xbinedge[:-1])
        return xbin, ymean, ystd, ycnt
    else:
        return [0], [0], [0], [0]


def dotpatch(x, y, imval, blank=None, dotsize=1, clipedge=True, pad=5, axes=None, **kwargs):
    '''
    Generate and plot a patch collection for a dot plot.

    Parameters
    ----------
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
    clipedge : boolean
        True to derive a square bounding box around non-Nan values.
        False to show full image.
    pad : int
        Padding in pixels around edges if clipedge=True
    axes : matplotlib.axes
        Axes for plotting
    **kwargs :
        Additional arguments passed to PatchCollection

    Returns
    -------
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
    if clipedge:
        valid = ~np.isnan(imval)
        xspan = x[valid].max() - x[valid].min()
        yspan = y[valid].max() - y[valid].min()
        if xspan > yspan:
            xminmax = [x[valid].min()-pad, x[valid].max()+pad]
            yminmax = [y[valid].min()-pad-np.floor((xspan-yspan)/2), 
                       y[valid].min()+pad-np.floor((xspan-yspan)/2)+xspan]
        else:
            yminmax = [y[valid].min()-pad, y[valid].max()+pad]
            xminmax = [x[valid].min()-pad-np.floor((yspan-xspan)/2), 
                       x[valid].min()+pad-np.floor((yspan-xspan)/2)+yspan]
    else:
        xminmax = [x.min(), x.max()]
        yminmax = [y.min(), y.max()]
    img = axes.add_collection(p)
    return img, xminmax, yminmax


def imarrayplot(x, y, imval, blank=None, clipedge=True, pad=5, axes=None, **kwargs):
    '''
    Plot a pixel image from a data column.

    Parameters
    ----------
    x : numpy.array
        The x values in the image plot (typically 'ix')
    y : numpy.array
        The y values in the image plot (typically 'iy')
    imval : numpy.array
        The z values which determine the dot colors
    blank : boolean array
        True values in this array are set to NaN
    clipedge : boolean
        True to derive a square bounding box around non-Nan values.
        False to show full image.
    pad : int
        Padding in pixels around edges if clipedge=True
    axes : matplotlib.axes
        Axes for plotting
    **kwargs :
        Additional arguments passed to imshow

    Returns
    -------
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
    sortorder = np.lexsort((y,x))
    imval = imval[sortorder]
    imarray = np.reshape(imval, [ydim,xdim], order='F')
    if clipedge:
        valid = ~np.isnan(imarray)
        xflat = np.flatnonzero(np.any(valid, axis=0))
        yflat = np.flatnonzero(np.any(valid, axis=1))
        xspan = xflat.max() - xflat.min()
        yspan = yflat.max() - yflat.min()
        if xspan > yspan:
            xminmax = [xflat.min()-pad, xflat.max()+pad]
            yminmax = [yflat.min()-pad-np.floor((xspan-yspan)/2), 
                       yflat.min()+pad-np.floor((xspan-yspan)/2)+xspan]
        else:
            yminmax = [yflat.min()-pad, yflat.max()+pad]
            xminmax = [xflat.min()-pad-np.floor((yspan-xspan)/2), 
                       xflat.min()+pad-np.floor((yspan-xspan)/2)+yspan]
    else:
        xminmax = [0, xdim-1]
        yminmax = [0, ydim-1]
    img = axes.imshow(imarray, origin='lower', **kwargs)
    return img, xminmax, yminmax


def gridplot(edgetab=None, gallist=None, columnlist=None, 
            xrange=None, yrange=None, blank=None, plotstyle='image',
            cmap='jet', nx=7, ny=6, dotsize=1, pdfname=None, pct=99,
            allnorm=False, vshow=True, clipedge=False, pad=5, verbose=False, 
            stretch='linear', maxlabel=18, do_cbar=True, **kwargs):
    '''
    Plot one column for multiple galaxies or multiple columns for 
    one galaxy on a grid.

    Parameters
    ----------
    edgetab : EdgeTable
        Table containing the galaxies and data to plot
    gallist : string or list of strings
        List of galaxy names to plot.  If plotting multiple columns,
        only one galaxy should be given.
    columnlist : string or list of strings
        List of columns in the table to plot.  If plotting multiple
        galaxies, only one column should be given.
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
    pct : float
        percentile interval for colormap normalization
    allnorm : boolean
        True to get percentile interval for whole sample rather than each galaxy
    stretch : string
        colormap stretch, can use 'log', 'sqrt', or 'linear'
    pdfname : string
        name of output PDF file, otherwise plot to screen
    vshow : boolean
        True to show vmin and vmax in plotting window
    clipedge : boolean
        True to derive a square bounding box around non-Nan values.
        False to use show full image.  Overridden by xrange, yrange.
    pad : int
        Padding in pixels around edges if clipedge=True
    maxlabel : int
        Units longer than this many characters are not shown when vshow=True
    do_cbar : boolean
        True to plot a colorbar beneath the first subplot
    **kwargs :
        Additional arguments passed to imarrayplot or dotpatch, including 
        vmin, vmax
    '''

    match stretch:
        case 'linear':
             stretch = LinearStretch()
        case 'sqrt':
             stretch = SqrtStretch()
        case 'log':
             stretch = LogStretch()
        case _:
            print('Invalid stretch',stetch,'- assuming linear')
            stretch = LinearStretch()

    if gallist is None and columnlist is None:
        raise TypeError('Either gallist or columnlist must be provided!')
    if isinstance(gallist, str):
        gallist = [gallist]
    if isinstance(columnlist, str):
        columnlist = [columnlist]
    if 'vmin' in kwargs:
        usermin = kwargs.pop('vmin')
    else:
        usermin = None
    if 'vmax' in kwargs:
        usermax = kwargs.pop('vmax')
    else:
        usermax = None

    # Plot mode: multiple galaxies or multiple columns
    if columnlist is not None and len(columnlist) == 1:
        mode = 'onecol'
        # Plot all galaxies by default
        if gallist is None:
            gallist = list(np.unique(edgetab['Name']))
        print('\nPlotting column',columnlist[0],'for',len(gallist),'galaxies')
        pagelist = gallist
        if allnorm:
            if not np.isnan(edgetab[columnlist[0]]).all():
                vmin, vmax = PercentileInterval(pct).get_limits(edgetab[columnlist[0]])
                norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    elif gallist is not None and len(gallist) == 1:
        mode = 'onegal'
        # Plot all non-coordinate columns by default
        if columnlist is None:
            columnlist = edgetab.colnames
            for key in ['Name', 'ix', 'iy', 'ra_abs', 'dec_abs', 'ra_off', 'dec_off']:
                if key in columnlist:
                    columnlist.remove(key)
        print('Plotting',len(columnlist),'columns for galaxy',gallist[0])
        pagelist = columnlist
        if allnorm:  # not sure this works
            vmin, vmax = PercentileInterval(pct).get_limits(
                             edgetab[edgetab['Name']==gallist[0]][columnlist])
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    else:
        raise ValueError('Specify either one galaxy or one column to plot')

    # Get default axis limits
    pages = int(np.ceil(float(len(pagelist)) / (nx*ny)))
    if pdfname is not None:
        pp = PdfPages(pdfname)

    for num in range(0,pages):
        aa = nx*ny*num
        bb = nx*ny+aa
        thispage = pagelist[aa:bb]
        fig = plt.figure(figsize=(18,14))
        if not verbose:
            print('Plotting', thispage[0], 'to', thispage[-1])
        if do_cbar:
            cbar_done = False

        for i in range(0,len(thispage)):
            if mode == 'onecol':
                gname = thispage[i]
                column = columnlist[0]
                label = gname
            else:
                gname = gallist[0]
                column = thispage[i]
                label = column
            if verbose:
                print('Plotting', thispage[i])
            ax = plt.subplot(ny,nx,i+1)
            galtab = (edgetab['Name'] == gname)
            if blank is not None:
                galblank = blank[galtab]
            else:
                galblank = None
            if not np.isnan(edgetab[galtab][column]).all():
                if not allnorm:
                    if usermin is not None and usermax is not None:
                        vmin = usermin
                        vmax = usermax
                    else:
                        vmin, vmax = PercentileInterval(pct).get_limits(edgetab[galtab][column])
                    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
                if plotstyle == 'dot':
                    img, xlims, ylims = dotpatch(edgetab[galtab]['ix'], 
                                                edgetab[galtab]['iy'],
                                                edgetab[galtab][column], 
                                                blank=galblank, clipedge=clipedge,
                                                pad=pad, dotsize=dotsize, cmap=cmap, 
                                                norm=norm, axes=ax, **kwargs)
                else:
                    img, xlims, ylims = imarrayplot(edgetab[galtab]['ix'], 
                                                edgetab[galtab]['iy'],
                                                edgetab[galtab][column], 
                                                blank=galblank, clipedge=clipedge,
                                                pad=pad, cmap=cmap, 
                                                norm=norm, axes=ax, **kwargs)        
                if xrange is None:
                    ax.set_xlim(xlims)
                    if i == 0:
                        print(label, "Default x limits used:",xlims)
                else:
                    ax.set_xlim(xrange)
                if yrange is None:
                    ax.set_ylim(ylims)
                    if i == 0:
                        print(label, "Default y limits used:",ylims)
                else:
                    ax.set_ylim(yrange)
                if vshow:
                    if vmax < 1:
                        labelstr = '[{:.3f} .. {:.3f}]'.format(vmin,vmax)
                    else:
                        labelstr = '[{:.2f} .. {:.2f}]'.format(vmin,vmax)
                    if edgetab[galtab][column].unit is not None:
                        if len(edgetab[galtab][column].unit.to_string()) <= maxlabel:
                            labelstr = f"{labelstr} {edgetab[galtab][column].unit:latex_inline}"
                    plt.text(0.04,0.06,labelstr,ha='left',va='center',size='small',
                        transform=ax.transAxes, bbox=dict(boxstyle='square,pad=0.1',facecolor='white',edgecolor='none'))
                if do_cbar and not cbar_done:
                    cax = ax.inset_axes([0.0, -0.07, 1.0, 0.04]) 
                    cbar = fig.colorbar(img, cax=cax, orientation='horizontal')
                    if not allnorm:
                        cbar.set_ticks([])
                    else:
                        cbar.ax.tick_params(labelsize='x-small')
                        cbar.ax.tick_params(size=0)
                    cbar_done = True
            ax.set_aspect('equal')
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.text(0.04,0.92,label,ha='left',va='center',transform=ax.transAxes,
               bbox=dict(facecolor='white', edgecolor='none', pad=1))
        fig.subplots_adjust(hspace=0.1)
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

def gridplotv2(edgetab=None, gallist=None, columnlist=None,
               global_table=None, pixel_table=None, 
               xrange=None, yrange=None, blank=None, plotstyle='image',
               cmap='jet', nx=7, ny=6, dotsize=1, pdfname=None, pct=99,
               allnorm=False, vshow=True, clipedge=False, pad=5, verbose=False, 
               stretch='linear', maxlabel=18, do_cbar=True,
               show_center=False, show_ellipse=False,
               Re_col=None, PA_col=None, AxIncl_col=None,
               **kwargs):
    
    from matplotlib.patches import Ellipse
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import LogNorm
    from matplotlib.colors import Normalize
    from astropy.visualization import ImageNormalize, PercentileInterval, LinearStretch, SqrtStretch, LogStretch
    from matplotlib.backends.backend_pdf import PdfPages

    Re_col = Re_col or 'Re'
    PA_col = PA_col or 'ledaPA'
    AxIncl_col = AxIncl_col or 'ledaAxIncl'

    def get_galaxy_center_pixel(galname, global_table=None, pixel_table=None):
        if global_table is None or pixel_table is None:
            return None
        galdata = global_table[global_table['Name'] == galname]
        if len(galdata) == 0: return None
        try:
            leda_ra = float(galdata['ledaRA'][0])
            leda_dec = float(galdata['ledaDE'][0])
            leda_coord = SkyCoord(ra=leda_ra*u.deg, dec=leda_dec*u.deg, frame='fk5')
        except Exception:
            return None
        galaxy_rows = pixel_table[pixel_table['Name'] == galname]
        if len(galaxy_rows) == 0: return None
        try:
            pixel_coords = SkyCoord(
                ra=np.array(galaxy_rows['ra_abs'])*u.deg,
                dec=np.array(galaxy_rows['dec_abs'])*u.deg,
                frame='fk5'
            )
            separations = leda_coord.separation(pixel_coords)
            nearest_idx = np.argmin(separations)
            ix = int(galaxy_rows['ix'][nearest_idx])
            iy = int(galaxy_rows['iy'][nearest_idx])
            pixel_ra = float(galaxy_rows['ra_abs'][nearest_idx])
            pixel_dec = float(galaxy_rows['dec_abs'][nearest_idx])
            offset_arcsec = separations[nearest_idx].to(u.arcsec).value
        except Exception:
            return None
        return ix, iy, pixel_ra, pixel_dec, offset_arcsec

    def add_ellipse(ax, galname, global_table, pixel_table):
        galdata = global_table[global_table['Name'] == galname]
        galaxy_rows = pixel_table[pixel_table['Name'] == galname]
        if len(galdata) == 0 or len(galaxy_rows) == 0:
            return
        center_info = get_galaxy_center_pixel(galname, global_table, pixel_table)
        if center_info is None:
            return

        ix_center, iy_center, _, _, _ = center_info

        pix_scale = np.abs(galaxy_rows['dec_abs'][1] - galaxy_rows['dec_abs'][0]) * 3600

        Re_arcsec = float(galdata[Re_col][0])
        ax_incl_deg = float(galdata[AxIncl_col][0])
        pa_deg = float(galdata[PA_col][0])

        a_pix = Re_arcsec / pix_scale
        b_pix = a_pix * np.cos(np.radians(ax_incl_deg))
        theta_deg = pa_deg + 90

        if verbose:
            print(f"\n=== Ellipse debug for {galname} ===")
            print("ix_center_panel:", ix_center)
            print("iy_center_panel:", iy_center)
            print("a_pix:", a_pix)
            print("b_pix:", b_pix)
            print("theta:", theta_deg)

        ell = Ellipse(
            (ix_center, iy_center),
            width=2*a_pix,
            height=2*b_pix,
            angle=theta_deg,
            edgecolor='red',
            facecolor='none',
            lw=2
        )
        ax.add_patch(ell)
        ax.plot(ix_center, iy_center, marker='o', color='red',
                markersize=14, markeredgewidth=3, markeredgecolor='black')

    match stretch:
        case 'linear': stretch = LinearStretch()
        case 'sqrt': stretch = SqrtStretch()
        case 'log': stretch = LogStretch()
        case _: stretch = LinearStretch()

    if gallist is None and columnlist is None:
        raise TypeError('Either gallist or columnlist must be provided!')
    if isinstance(gallist, str): gallist = [gallist]
    if isinstance(columnlist, str): columnlist = [columnlist]
    if 'norm' in kwargs: kwargs.pop('norm')

    if columnlist is not None and len(columnlist) == 1:
        mode = 'onecol'
        if gallist is None:
            gallist = list(np.unique(edgetab['Name']))
        pagelist = gallist
        if allnorm:
            vmin, vmax = PercentileInterval(pct).get_limits(edgetab[columnlist[0]])
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    elif gallist is not None and len(gallist) == 1:
        mode = 'onegal'
        if columnlist is None:
            columnlist = [c for c in edgetab.colnames if c not in 
                          ['Name','ix','iy','ra_abs','dec_abs','ra_off','dec_off']]
        pagelist = columnlist
        if allnorm:
            vmin, vmax = PercentileInterval(pct).get_limits(
                edgetab[edgetab['Name']==gallist[0]][columnlist])
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    else:
        raise ValueError('Specify either one galaxy or one column to plot')

    pages = int(np.ceil(len(pagelist)/(nx*ny)))
    if pdfname: pp = PdfPages(pdfname)

    for num in range(pages):
        aa, bb = nx*ny*num, nx*ny*(num+1)
        thispage = pagelist[aa:bb]
        fig = plt.figure(figsize=(18,14))
        cbar_done = False

        for i, item in enumerate(thispage):
            ax = plt.subplot(ny,nx,i+1)
            if mode == 'onecol':
                gname = item
                column = columnlist[0]
            else:
                gname = gallist[0]
                column = item
            galtab = (edgetab['Name'] == gname)
            galblank = blank[galtab] if blank is not None else None

            if not np.isnan(edgetab[galtab][column]).all():
                if not allnorm:
                    vmin, vmax = PercentileInterval(pct).get_limits(edgetab[galtab][column])
                    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)

                if plotstyle == 'dot':
                    img, xlims, ylims = dotpatch(edgetab[galtab]['ix'],
                                                 edgetab[galtab]['iy'],
                                                 edgetab[galtab][column],
                                                 blank=galblank, clipedge=clipedge,
                                                 pad=pad, dotsize=dotsize, cmap=cmap,
                                                 norm=norm, axes=ax, **kwargs)
                else:
                    img, xlims, ylims = imarrayplot(edgetab[galtab]['ix'],
                                                   edgetab[galtab]['iy'],
                                                   edgetab[galtab][column],
                                                   blank=galblank, clipedge=clipedge,
                                                   pad=pad, cmap=cmap,
                                                   norm=norm, axes=ax, **kwargs)

                ax.set_xlim(xrange if xrange else xlims)
                ax.set_ylim(yrange if yrange else ylims)

                if show_center:
                    center_info = get_galaxy_center_pixel(gname, global_table, pixel_table)
                    if center_info:
                        ix, iy, _, _, _ = center_info
                        ax.plot(ix, iy, marker='x', color='red', markersize=20, markeredgewidth=2.5)

                if show_ellipse:
                    add_ellipse(ax, gname, global_table, pixel_table)

                if vshow:
                    labelstr = '[{:.3f} .. {:.3f}]'.format(vmin,vmax) if vmax < 1 else '[{:.2f} .. {:.2f}]'.format(vmin,vmax)
                    if hasattr(edgetab[galtab][column], 'unit') and edgetab[galtab][column].unit is not None:
                        if len(edgetab[galtab][column].unit.to_string()) <= maxlabel:
                            labelstr += f" {edgetab[galtab][column].unit:latex_inline}"
                    plt.text(0.04,0.06,labelstr,ha='left',va='center',size='small',transform=ax.transAxes,
                             bbox=dict(boxstyle='square,pad=0.1',facecolor='white',edgecolor='none'))

                if do_cbar and not cbar_done:
                    cax = ax.inset_axes([0.0, -0.07, 1.0, 0.04])
                    cbar = fig.colorbar(img, cax=cax, orientation='horizontal')
                    if allnorm:
                        cbar.ax.tick_params(labelsize='x-small')
                        cbar.ax.tick_params(size=0)
                    else:
                        cbar.set_ticks([])
                    cbar_done = True

            ax.set_aspect('equal')
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.text(0.04,0.92, column if mode=='onegal' else gname, ha='left', va='center',
                     transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='none', pad=1))

        fig.subplots_adjust(hspace=0.1, wspace=0.05)
        if pdfname:
            pp.savefig(bbox_inches='tight', pad_inches=0.1)
            plt.close()
        else:
            plt.show()
    if pdfname: pp.close()