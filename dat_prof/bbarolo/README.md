Results of Bbarolo fitting are in self-documenting ascii.ecsv format, corresponding to 3 x 2 x 2 = 12 different runs

* 'natv', 'smo5', 'smo7' indicates cubes are at native resolution or common resolution of 5" or 7".
* 'fitvd' or 'fixvd' indicates velocity dispersion is fit as a function of radius or fixed at 8 km/s.
* 'bbmsk' or 'dilmsk' indicates whether Bbarolo's SMOOTH masking was used or the IDL dilated mask.

For quick inspection, bb_natv_fitvd_dilmsk.csv is recommended.

An initial run is always performed with PA and VSYS fit per ring.  The results of this run are in the *_freepa.csv files.
