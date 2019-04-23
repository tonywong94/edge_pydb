import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy import units as u
from uncertainties import ufloat
import uncertainties.unumpy as unp 

# To-do: Return Column object only if inputs are Columns

# Convert Halpha intensity to A_V-corrected SFR surface density
def sfr_ha(flux_ha, flux_hb, name=None):
    # Extinction curve from Cardelli+(1989).
    K_Ha = 2.53
    K_Hb = 3.61
    # Eq(1) from Catalan-Torrecilla+(2015). 
    A_Ha = K_Ha/(-0.4*(K_Ha-K_Hb)) * np.log10((flux_ha/flux_hb)/2.86)
    A_Ha[A_Ha < 0] = 0.
    flux_ha_cor = flux_ha * 10**(0.4*A_Ha)
    # input line flux is actually flux per arcsec^2
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    sb_ha  = flux_ha_cor * sterad   # flux per steradian
    lsd_ha = 4*np.pi * sb_ha
    # Eq(4) from Catalan-Torrecilla+(2015).
    lumcon = 5.5e-42 * (u.solMass/u.yr) / (u.erg/u.s)
    sigsfr = (lumcon * lsd_ha).to(u.solMass/(u.pc**2*u.Gyr))
    return Column(sigsfr, name=name)
    
# Convert CO intensity to H2(+He) surface density
def msd_co(sb_co, alphaco=4.3, name=None):
    convfac = alphaco * (u.solMass/u.pc**2) / (u.K*u.km/u.s)
    sig_mol = (convfac*sb_co).to(u.solMass/u.pc**2)
    return Column(sig_mol, name=name)
    
# Convert units for stellar surface density
def stmass_pc2(stmass_as2, dist=10*u.Mpc, name=None):
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    pxarea = (dist**2/sterad).to(u.pc**2)
    stmass_pc2 = 10**stmass_as2 * u.solMass / pxarea
    return Column(stmass_pc2, name=name)

# return log10(a/b), taking care of -InF values from the logarithm. 
def ulogratio(a, b, ae = 0.0, be = 0.0): 
    try: 
        ua = ufloat(a, ae)
        ub = ufloat(b, be)
        logr = unp.log10(ua) - unp.log10(ub)
        return logr.n, logr.s
    except:
        return np.nan, np.nan

# BPT classification, see Husemann et al. (2013) Figure 7.
def bpt_type(flux_nii, flux_oiii, flux_ha, flux_hb, ew_ha):
    n2ha, _ = ulogratio(flux_nii,  flux_ha)
    o3hb, _ = ulogratio(flux_oiii, flux_hb)    

    kewley01 = lambda nii: 1.19 + 0.61/(nii - 0.47)
    kauffm03 = lambda nii: 1.30 + 0.61/(nii - 0.05)
    cidfer10 = lambda nii: 0.48 + 1.01*nii

    try:
        if n2ha > -1.5 and n2ha < -0.1 and o3hb < kewley01(n2ha) and abs(ew_ha) > 6.0:
            BPT = -1  # Starforming
        elif n2ha < 0.3 and o3hb < kauffm03(n2ha):
            BPT = 0   # Intermediate
        elif o3hb > -1  and o3hb < cidfer10(n2ha):
            BPT = 1   # LINER
        elif o3hb > -1:
            BPT = 2   # Seyfert
        else:
            BPT = np.nan
        return BPT
    except:
        return np.nan    

# Metallicity derived from O3N2 line ratio. 
def ZOH_o3n2(flux_nii, flux_oiii, flux_ha, flux_hb, ew_ha, err_nii, err_oiii, 
             err_ha, err_hb, err=False):
    try:
        uNIIF = ufloat(flux_nii, err_nii)
        uOIIIF = ufloat(flux_oiii, err_oiii)
        uHaF = ufloat(flux_ha, err_ha)
        uHbF = ufloat(flux_hb, err_hb)
    
        uO3N2 = unp.log10(uOIIIF) - unp.log10(uHbF) - (unp.log10(uNIIF) 
                    - unp.log10(uHaF))  
    
        BPT = bpt_type(flux_nii, flux_oiii, flux_ha, flux_hb, ew_ha)
                
        if BPT == -1:
            uOH_M13  = 8.533 - 0.214*uO3N2 # Eq(2) from Marino+2013
            if err: 
                return uOH_M13.s
            else:            
                return uOH_M13.n
        else:
            return np.nan
    except:
        return np.nan
 
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
