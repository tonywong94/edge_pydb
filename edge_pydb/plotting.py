import numpy as np
from scipy import stats
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

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
def dotpatch(x, y, colors, size=1, vmin=None, vmax=None, cmap='jet'):

    patches = []
    for x1, y1 in zip(x,y):
        circle = Circle((x1, y1), size)
        patches.append(circle)

    p = PatchCollection(patches, cmap=cmap)
    p.set_array(np.array(colors))
    if vmax is None:
        vmax = np.nanmax(colors)
    if vmin is None:
        vmin = np.nanmin(colors)

    p.set_clim([vmin, vmax])
    xymin = np.min([x.min(),y.min()])
    xymax = np.max([x.max(),y.max()])

    return p, xymin, xymax
