import numpy as np

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
