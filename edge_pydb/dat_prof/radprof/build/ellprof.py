#!/usr/bin/env python

import numpy as np
from astropy.table import Table
from astropy.io.fits import getdata 

def ellprof(image, pa=0, incline=0, center=None, truncate=False, 
        binsize=1, wts=None, blank=0, fracmask=0, nanmask=False):
    """
    Calculate the azimuthally averaged radial profile in elliptical rings.
    Based on radialprofile.py from Adam Ginsburg: 
    https://github.com/keflavich/image_tools

    INPUTS:
    image - The 2D image
    pa - position angle of major axis in degrees E of N.  Default=0.
    incline - inclination of disk in degrees, must be < 90.  Default=0.
    center - The [x,y] pixel coordinates used as the center. The default is None, 
        which then uses the center of the image (including fractional pixels).
    binsize - size of the averaging bins in pixels.
    wts - can do a weighted average instead of a simple average if this keyword 
        parameter is set.  wts.shape must = image.shape.  
    blank - value to use in place of NaNs, can be 2D image of same shape. Default = 0.
    nanmask - set to True to ignore NaNs in profile. Default = False (use blank value).
    fracmask - minimum fraction of non-NaN pixels required in a ring.  Default=0.

    OUTPUTS: bin centers in pixels, # of pixels per bin, radial profile, cum profile

    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...
    
    """
    # Replace NaNs with blank value
    inds=np.where(np.isnan(image))
    good=np.where(~np.isnan(image))
    if len(inds[0]) > 0:
        print(len(good[0]),'good pixels;',len(inds[0]),'nan pixels')
        if isinstance(blank, np.ndarray):
            image[inds] = blank[inds]
            image[np.where(np.isnan(image))]=0.
        else:
            image[inds] = blank

    # Calculate the indices from the image
    y, x = np.indices(image.shape)
    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    dx = x - center[0]
    dy = y - center[1]
    cospa = np.cos(np.radians(pa))
    sinpa = np.sin(np.radians(pa))
    cosi  = np.cos(np.radians(incline))

    r = np.hypot( dy*cospa-dx*sinpa, (dy*sinpa+dx*cospa)/cosi )
    if len(good[0]) > 0:
        rmax = np.amax(r[good])
        print('The maximum radius for a valid pixel is {:.2f}'.format(rmax))
    else:
        rmax = np.nan

    if wts is None:
        wts = np.ones(image.shape)

    # Initial mask is False for NaN values
    mask = np.ones(image.shape,dtype='bool')
    mask[inds] = False

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    if truncate and np.isfinite(rmax):
        nbins = int(np.round(rmax / binsize)+1)
    else:
        nbins = int(np.round(r.max() / binsize)+1)
    maxbin = nbins * binsize
    bin_edges = np.linspace(0,maxbin,nbins+1)
    # but we're probably more interested in the bin centers than their edges
    bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.0

    # number of pixels per bin
    ngood = np.histogram(r, bin_edges, weights=mask.astype('int'))[0]
    if nanmask==False:
        mask[inds] = True
    nr = np.histogram(r, bin_edges, weights=mask.astype('int'))[0]
    
    # Mask based on fill fraction if requested
    which2d = np.digitize(r,bin_edges)
    for b in range(1,nbins+1):
        if ngood[b-1] < fracmask*nr[b-1]:
            mask[which2d==b] = False

    # rms in each annulus
    # Find out which radial bin each point in the map belongs to
    which1d = np.digitize(r.flat,bin_edges)
    rms = np.array([image.flat[mask.flat*(which1d==b)].std() 
        for b in range(1,nbins+1)])

    # mean in each annulus
    annsum = np.histogram(r, bin_edges, weights=(image*mask))[0]
    cumsum = np.cumsum(annsum)
    wtmean = (np.histogram(r, bin_edges, weights=(image*wts*mask))[0] 
        / np.histogram(r, bin_edges, weights=(mask*wts))[0])

    # create a Table structure
    tab = Table([bin_centers, nr, ngood, wtmean, rms, annsum, cumsum], 
        names=('radius', 'npix', 'ngood', 'wtmean', 'rms', 'annsum', 'cumsum'))
    return tab, r

