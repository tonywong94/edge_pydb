The 8 txt files are astropy tables in ascii.ecsv format, corresponding to four different runs

* gal_nrad8 contains the result of Bbarolo fitting with 8 * 3 (arcsec) radii
* if the name of the file contains 'smolist', the corresponding run uses SMOOTH masks; otherwise, dilated mask
* if the name of the file contains 'vdisp8', the corresponding run fixes velocity dispersion at 8 km/s' otherwise, set as free variable

Densprof contains the density profile of the corresponding run. they are required to generate the plots of parameters.

GenPlot provides a function that plots the fitting results of the parameters, including both stages. See example in genPlot.
