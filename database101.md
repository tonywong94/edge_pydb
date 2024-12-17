# Database 101
This document aims to provide a detailed "under-the-hood" guide to the `edge_pydb` infrastructure for constructing matched resolution CO-IFU databases.

## Input global parameters
For the overall galaxy sample, a CSV file (referred to as `ortpar` in the code) is required which contains columns of basic information including:
* Center coordinates `coln_ra` and `coln_dc`
* Disk inclination `coln_inc` (in degrees or as an axial ratio b/a)
* Disk position angle `coln_pa`
Another CSV file (referred to as `distpar` in the code, which could be the same file as `ortpar`) should contain the galaxy distance that was assumed by Pipe3D as `coln_dmpc`.

Note that additional global parameters (such as total stellar mass and SFR) are useful for scienctific analysis but not needed for database generation.

## Input CO data
The CO FITS files (all 2-D images) are converted into tables using `do_comom.py` (import command: `from edge_pydb.img_comom import do_comom`).  It is assumed that the CO cubes have already been processed into moment maps using the [maskmoment](https://github.com/tonywong94/maskmoment) package.  Up to three masking types are supported per galaxy:
* Straight (`str`): no masking, only moment-0 and its uncertainty are saved.
* Dilated (`dil`): dilated mask without smoothing, all moments and a peak S/N ratio image are saved.
* Smoothed (`smo`): dilated mask with smoothing, all moments are saved.

The FITS files should all exist within a single directory with file names obeying the following convention: `GNAME.LINELBL.SEQ_MASK.FTYPE.fits.gz`, where GNAME is the galaxy name, LINELBL is a string identifying the spectral line, SEQ is a string identifying the spatial resolution, MASK is one of the three mask identifiers in the bulleted list above, and FTYPE describes the type of data product, such as `mom0` for moment-0 or `emom1` for the uncertainty in moment-1.  For example, `UGC10710.co.smo7_smo.emom2.fits.gz`.

To produce the `sigmol` column (units of solMass/pc2), a nominal value of `alphaco` is required.  The default value is 4.3, corresponding to the standard Galactic value for CO(1-0).  Any variable &alpha;<sub>CO</sub> is a scaling factor relative to this default value.  Variable conversion factors for the 2-1 line are not yet implemented.

## Input IFU data
It is assumed that the IFU data have been processed by some version of the Pipe3D package.  Typically Pipe3D produces the following five outputs for a single galaxy, all on the same astrometric grid.  Each output is a pseudo-cube in which each plane is a map of a particular quantity.
1. `ELINES`: parametric fits to 9 bright emission lines, plus the centroid velocity and velocity dispersion for H&alpha;.  For CALIFA, there are 20 planes including uncertainties for the line fluxes.  For MaNGA, the uncertainties are omitted and there are 11 planes.
2. `flux_elines`: non-parametric fits to many more emission lines.  For each line, 8 planes (corresponding to a value and uncertainty for flux, velocity, dispersion, and equivalent width) are generated.  The total number of planes can range from 240 (AMUSING) to 456 (MaNGA), but `edge_pydb` only keeps data for 10 emission lines (the nine used for ELINES plus [OI]6300).
3. `indices`: Nine spectral indices and their uncertainties, providing a total of 18 planes.
4. `SFH`: Luminosity fractions in different age and metallicity bins, resulting from the single stellar population (SSP) analysis.  For CALIFA there are 39 age bins and 4 metallicity bins; for MaNGA there are 39 age bins and 7 metallicity bins.  Mass fractions are also calculated using the mass-to-light ratios for each SSP.  Both luminosity and mass fractions, aggregated by either age or metallicity, are saved by `edge_pydb`. An "SSP-based" star formation rate is estimated by taking the mass fraction of stars less than 33 Myr old, scaling by the stellar mass, and dividing by 33 Myr.
5. `SSP`: Average properties of the fitted stellar populations, including stellar mass, age, metallicity, velocity, and velocity dispersion.  For CALIFA the `GSD156` stellar library is used, whereas for MaNGA the `MaStar_sLOG` stellar library is used (as of DR17).  A Salpeter (1955) IMF is assumed throughout.  20 planes are provided for CALIFA, with an additional plane (uncertainty in stellar mass) provided for MaNGA.

In older Pipe3D runs (`packed=False`), the five pseudo-cubes are separate FITS files with specific naming conventions.  In more recent runs using pyPipe3D (`packed=True`), the five pseudo-cubes are provided as multiple extensions in a single FITS file, with extension 0 having no data array but providing the astrometry in the header.

The Pipe3D outputs are converted into tables by `do_pipe3d.py` (import command: `from edge_pydb.img_califa import do_pipe3d`).  The following important issues should be noted that may differ between data sets.

* __Masking/Blanking:__ In the CALIFA DR3 (v2.2) processing outputs, the bad pixel masks are a separate set of FITS images with values of 0 at the location of foreground stars.  These are currently applied before running `do_pipe3d.py`.  Aside from applying the mask files, one should also mask any zero values as these normally occur outside the field of view of the IFU.  For the MaNGA and CALIFA eDR (v2.3) processing, the masks are incorporated into the packed FITS files as two additional extensions.  The `GAIA_MASK` is set to 1 at the locations of foreground stars, whereas `SELECT_REG` is set to 0 at locations not observed or where the S/N ratio is low.  Note that this is the S/N ratio of the stellar continuum, so applying this selection to the emission lines may not be appropriate.

* __Extinction:__ The gas extinction is estimated from the Balmer decrement ratio between H&alpha; and H&beta;.  An intrinsic ratio of 2.86 and the Cardelli et al. (1989) extinction curve is assumed.  When the ratio is less than 2.86 (formally negative extinction), the extinction is assumed to be 0.  Also, very high extinctions (>6 mag by default) are blanked.  The noise in the Balmer decrement can be reduced by pre-smoothing the H&alpha; and H&beta; images.  This technique is used to generate a "smooth" A(Ha) value (`AHa_smooth`) and an "adopted" SFR surface density (`sigsfr_adopt`) which should be compared with the (extinction) "corrected" SFR surface density without additional smoothing (`sigsfr_corr`) to evaluate the significance of any improvement.

* __Metallicity:__ Three metallicity calculations are included, namely the O3N2 calibration from Marino et al. (2013), O3N2 calibration from Pettini & Pagel (2004), and the N2 calibration from Marino et al. (2013).  No S/N cuts on the individual line fluxes are used, though there is a requirement that the BPT type of the spaxel is star-forming (below the Kauffmann line with H&alpha; EW > 6), requiring non-negligible flux in all 4 emission lines.

* __Dezonification:__ The stellar mass map in the `SSP` extension is spatially binned, so should be multiplied by the dezonification map (the ratio of the original continuum to the spatially binned continuum) to yield a map closer to the native resolution of the data.  The `sigstar` maps in the output database have been dezonified.

## Sampling
The current code supports outputting every pixel into the database (for `allpix=True`) or outputting every 3rd pixel in RA and Dec (for `allpix=False`).  Hexagonal grid sampling is also supported (`hexgrid=True`), though not well tested at this point.

## Regridding
Images are regridded using the `reproject_interp` function of the `reproject` package, using bilinear interpolation by default (`interp_order=1`).  To regrid the CO maps to match the IFU grid, provide the `p3d_dir` and `p3dtempl` parameters to `do_comom.py`, which is run after `do_pipe3d.py`.  This is advisable when the CO and IFU resolutions are not matched to begin with and the IFU data have substantially better resolution: one is better off regridding the lower resolution image because regridding the higher resolution image will degrade its resolution more noticeably.  More typically, the FITS files are already resolution matched, and one will regrid the IFU maps to match the CO grid, since the CO maps will cover a larger area.  This is accomplished by providing the `comomdir` and `cotempl` parameters to `do_pipe3d.py`, which is run after `do_comom.py`.

## Resolution matching
Originally the resolution matching was performed directly on the FITS images, before any processing using `edge_pydb`, and often before even running Pipe3D.  Now it is possible for some basic resolution matching to be performed in `do_pipe3d.py`.  This should be done cautiously as convolving the output of Pipe3D is not the same as convolving the original data and re-running Pipe3D.  The IFU data are always matched to the resolution of the CO data.  This works best when both the IFU and CO PSFs are Gaussian (i.e. MaNGA), since the convolution kernel is determined by `create_matching_kernel` in `photutils` which does not check that the output beam is actually larger than the input beam.  To activate this capability, provide `matchres=True` in `do_pipe3d.py` and ensure that the `comomdir` and `cotempl` parameters are also set.  The IFU data will be matched to the CO resolution but _not_ placed on the CO grid.  Then `do_comom.py` should be run as the second step, with parameters `p3d_dir` and `p3dtempl` set to ensure regridding to the IFU grid.

