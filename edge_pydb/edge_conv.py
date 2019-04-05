import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy import units as u

# Convert Halpha intensity to A_V-corrected SFR surface density
def sfr_ha(flux_ha, flux_hb):
    # Extinction curve from Cardelli+(1989).
    K_Ha = 2.53
    K_Hb = 3.61
    # Eq(1) from Catalan-Torrecilla+(2015). 
    A_Ha = K_Ha/(-0.4*(K_Ha-K_Hb)) * np.log10((flux_ha/flux_hb)/2.86)
    A_Ha[A_Ha < 0] = 0.
    flux_ha_cor = flux_ha * 10**(0.4*A_Ha)
    # line flux is actually flux per arcsec2
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    sb_ha  = flux_ha_cor * sterad   # flux per steradian
    lsd_ha = 4*np.pi * sb_ha
    # Eq(4) from Catalan-Torrecilla+(2015).
    lumcon = 5.5e-42 * (u.solMass/u.yr) / (u.erg/u.s)
    sigsfr = (lumcon * lsd_ha).to(u.solMass/(u.pc**2*u.Gyr))
    return sigsfr
    
# Convert CO intensity to H2(+He) surface density
def msd_co(sb_co, alphaco=4.3):
    convfac = alphaco * (u.solMass/u.pc**2) / (u.K*u.km/u.s)
    sig_mol = (convfac*sb_co).to(u.solMass/u.pc**2)
    return sig_mol
    
# Convert units for stellar surface density
def stmass_pc2(stmass_as2, dist=10*u.Mpc):
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    pxarea = (dist**2/sterad).to(u.pc**2)
    stmass_pc2 = 10**stmass_as2 * u.solMass / pxarea
    return stmass_pc2

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
