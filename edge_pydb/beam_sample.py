from datetime import datetime
import warnings
import numpy as np
from astropy.table import Table, Column
from astropy.convolution import convolve
from astropy import units as u
import radio_beam
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def beam_sample(edgetab, gallist=None, columnlist=None, beam_orig=7, beam_final=20, 
                outfile=None, doplots=False, pdfprefix='beamsample', nx=5, ny=4):

    """
    Process one or more galaxies from an allpix HDF5 file and output central surface
    brightness after convolution to a user-specified Gaussian beam.  The pixel with
    the smallest value of rad_arc is used; this should be recalculated using gc_polr
    beforehand if the default central pixel is not desired.

    === Parameters ===
    edgetab : EdgeTable
        Table, extracted from an allpix HDF5 file
    gallist : string or list of strings
        List of galaxy names to process.  Default is to process all galaxies.
    columnlist : string or list of strings
        List of columns in the table to process.  Default is to process all columns,
        excluding the initial coordinate columns.
        NOTE: Use caution when convolving non-flux columns.
    beam_orig : float
        The FWHM of the original beam in arcsec.  Default 7".
    beam_final : float
        The FWHM of the target beam (to convolve to) in arcsec.  Default 20".
    outfile : string
        The name of the output table.  Default is not to write the table.
    doplots : boolean
        Whether to output plots for debugging purposes
    pdfprefix : string
        Prefix for plots if doplots=True.  The galaxy name is appended to the prefix.
    nx : int
        number of sub-panels in x direction
    ny : int
        number of sub-panels in y direction

    === Returns ===
    tabout : Table
        Astropy table with the results
    """

    # Check the table spacing
    stride = edgetab["iy"][1] - edgetab["iy"][0]
    if stride != 1:
        warnings.warn("### Pixel stride for does not appear to be 1")
    aspp = round(3600. * (edgetab["dec_abs"][1]-edgetab["dec_abs"][0]) / stride, 1)
    print("Assuming pixel size of {} arcseconds\n".format(aspp))

    # Determine the original resolution
    oldbeam = radio_beam.Beam(major = beam_orig * u.arcsec, 
                              minor = beam_orig * u.arcsec, pa = 0 * u.deg)
    newbeam = radio_beam.Beam(major = beam_final * u.arcsec, 
                              minor = beam_final * u.arcsec, pa = 0 * u.deg)
    conv_beam = newbeam.deconvolve(oldbeam, failure_returns_pointlike=True)
    if conv_beam.major > 0:
        conv_kernel = conv_beam.as_kernel(aspp*u.arcsec)
    else:
        print("### Beam deconvolution failed - sampling original images")
    
    # Process all galaxies and columns by default
    if gallist is None:
        gallist = list(np.unique(edgetab['Name']))
    elif isinstance(gallist, str):
        gallist = [gallist]
    if columnlist is None:
        columnlist = edgetab.colnames
        for key in ['Name', 'ix', 'iy', 'ra_abs', 'dec_abs', 
                    'ra_off', 'dec_off', 'rad_arc', 'azi_ang']:
            if key in columnlist:
                columnlist.remove(key)
    elif isinstance(columnlist, str):
        columnlist = [columnlist]

    # Create the output table
    namecol = Column(name='Name', data=gallist)
    tabout = Table()
    tabout.add_column(namecol)
    for colname in columnlist:
        emptyarray = np.full(len(gallist), np.nan)
        tabout.add_column(emptyarray, name=colname)
    
    # Loop over galaxies
    for i_gal, gname in enumerate(gallist):
        galtab = edgetab[edgetab['Name'] == gname]
        print('\nWorking on galaxy {}'.format(gname))
        if doplots:
            pdfname = pdfprefix+'_'+gname+'.pdf'
            pp = PdfPages(pdfname)
            fig = plt.figure(figsize=(18,14))
        # Loop over columns
        for i_col, colname in enumerate(columnlist):
            print('Working on {} with units {}'.format(colname,galtab[colname].unit))
            if not np.isnan(galtab[colname]).all():
                xdim = len(np.unique(galtab['ix']))
                ydim = len(np.unique(galtab['iy']))
                imarray = np.reshape(galtab[colname], [ydim,xdim], order='F')
                if conv_beam.major > 0:
                    conv_image = convolve(imarray, conv_kernel, normalize_kernel=True,
                                      preserve_nan=True, boundary='wrap')
                else:
                    conv_image = imarray
                if doplots:
                    ax = plt.subplot(ny,nx,i_col+1)
                    ax.imshow(conv_image, origin='lower')
                    ax.set_aspect('equal')
                    ax.xaxis.set_ticks([])
                    ax.yaxis.set_ticks([])
                    ax.text(0.04,0.92,colname,transform=ax.transAxes,ha='left',
                       va='center',bbox=dict(facecolor='white', edgecolor='none', pad=1))
                galtab.add_index('rad_arc')
                ix = galtab.iloc['rad_arc',0]['ix']
                iy = galtab.iloc['rad_arc',0]['iy']
                print('Extracting pixel at [x,y]=[{},{}]'.format(ix,iy))
                print('The original value was',imarray[iy, ix])
                pickval = conv_image[iy, ix]
                print('The convolved value is',pickval)
                if hasattr(pickval, 'unit'):
                    pickval = pickval.value
                tabout[i_gal][colname] = pickval
                tabout[colname].unit = galtab[colname].unit
                tabout[colname].description = galtab[colname].description
        if doplots:
            fig.subplots_adjust(hspace=0.05)
            fig.subplots_adjust(wspace=0.05)
            pp.savefig(bbox_inches = 'tight', pad_inches=0.1)
            pp.close()
            print('\nPDF output to file',pdfname)

    tabout.pprint()
    tabout.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    tabout.meta['comments'] = ('Central brightness after convolution with '
                'original beam {} and final beam {}'.format(beam_orig, beam_final))
    if outfile is not None:
        tabout.write(outfile, format='ascii.ecsv', delimiter=',', overwrite=True)
        print('\nCSV output to file',outfile)

    return tabout
