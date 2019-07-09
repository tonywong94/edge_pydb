import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy import units as u
from astropy.coordinates import SkyCoord


# Calculate galactocentric polar coordinates 
# (radius in arcsec, azangle in degrees from receding majaxis)
# Inputs should all be in degrees
def gc_polr(ra, dec, ra_gc, dec_gc, pa, inc):
    ctr = SkyCoord(ra_gc, dec_gc, unit="deg")
    pos = SkyCoord(ra, dec, unit="deg")
    # Polar vector in sky plane
    t_sky = ctr.position_angle(pos).degree
    r_sky = ctr.separation(pos).arcsecond
    #print(t_sky,r_sky)
    # Convert to galaxy plane
    x_sky = r_sky * np.cos(np.radians(t_sky - pa))
    y_sky = r_sky * np.sin(np.radians(t_sky - pa))
    x_gal = x_sky
    y_gal = y_sky/np.cos(np.radians(inc))
    # Reject very high inclinations (cannot be deprojected)
    y_gal[inc>89] = np.nan
    radius = np.sqrt(x_gal**2 + y_gal**2)
    azang  = np.degrees(np.arctan2(y_gal, x_gal))
    return radius, azang


# Convert Halpha intensity to A_V-corrected SFR surface density
def sfr_ha(flux_ha, flux_hb, name='sig_sfr'):
    # TODO: Return A_Ha column as well.
    # Extinction curve from Cardelli+(1989).
    K_Ha = 2.53
    K_Hb = 3.61
    # Eq(1) from Catalan-Torrecilla+(2015). 
    A_Ha = K_Ha/(-0.4*(K_Ha-K_Hb)) * np.log10((flux_ha/flux_hb)/2.86)
    # Do not apply negative extinction.
    A_Ha[A_Ha < 0] = 0.
    flux_ha_cor = flux_ha * 10**(0.4*A_Ha)
    # input line flux is actually flux per arcsec^2
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    sb_ha  = flux_ha_cor * sterad   # flux per steradian
    lsd_ha = 4*np.pi * sb_ha
    # Eq(4) from Catalan-Torrecilla+(2015).
    lumcon = 5.5e-42 * (u.solMass/u.yr) / (u.erg/u.s)
    sig_sfr = (lumcon * lsd_ha).to(u.solMass/(u.pc**2*u.Gyr))
    if isinstance(flux_ha, Column) and isinstance(flux_hb, Column):
        return Column(sig_sfr, name=name)
    else:
        return sig_sfr
    

# Convert CO intensity to H2(+He) surface density
def msd_co(sb_co, alphaco=4.3, name='sig_mol'):
    convfac = alphaco * (u.solMass/u.pc**2) / (u.K*u.km/u.s)
    sig_mol = (convfac*sb_co).to(u.solMass/u.pc**2)
    if isinstance(sb_co, Column):
        return Column(sig_mol, name=name)
    else:
        return sig_mol
    

# Convert units for stellar surface density
def stmass_pc2(stmass_as2, dist=10*u.Mpc, name='sig_star'):
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    pxarea = (dist**2/sterad).to(u.pc**2)
    if isinstance(stmass_as2, Column):
        stmass_pc2 = 10**np.array(stmass_as2) * u.solMass / pxarea
        stmass_pc2[~np.isfinite(stmass_pc2)] = np.nan
        return Column(stmass_pc2, name=name)
    else:
        stmass_pc2 = 10**stmass_as2 * u.solMass / pxarea
        stmass_pc2[~np.isfinite(stmass_pc2)] = np.nan
        return stmass_pc2


# return log10(a/b), taking care of -InF values from the logarithm. 
# def ulogratio(a, b, ae = 0.0, be = 0.0): 
#     try: 
#         ua = ufloat(a, ae)
#         ub = ufloat(b, be)
#         logr = unp.log10(ua) - unp.log10(ub)
#         return logr.n, logr.s
#     except:
#         return np.nan, np.nan


# BPT classification, see Husemann et al. (2013A&A...549A..87H) Figure 7.
def bpt_type(flux_nii, flux_oiii, flux_ha, flux_hb, ew_ha):

    good = (flux_nii>0) & (flux_oiii>0) & (flux_ha>0) & (flux_hb>0) & (~np.isnan(ew_ha))
    n2ha = np.full(len(flux_nii), np.nan)
    n2ha[good] = np.log10(flux_nii[good])  - np.log10(flux_ha[good])
    o3hb = np.full(len(flux_oiii), np.nan)
    o3hb[good] = np.log10(flux_oiii[good]) - np.log10(flux_hb[good])    

    kewley01 = lambda nii: 1.19 + 0.61/(nii - 0.47) # Eq. 5 of 2001ApJ...556..121K
    kauffm03 = lambda nii: 1.30 + 0.61/(nii - 0.05) # Eq. 1 of 2003MNRAS.346.1055K
    cidfer10 = lambda nii: 0.48 + 1.01*nii          # Eq. 3 of 2010MNRAS.403.1036C

    BPT = np.full(len(n2ha), np.nan)
    sf = (n2ha > -1.5) & (n2ha < -0.1) & (o3hb < kauffm03(n2ha)) & (abs(ew_ha) > 6.0)
    BPT[sf] = -1
    inter = (~sf) & (n2ha < 0.3) & (o3hb < kewley01(n2ha))
    BPT[inter] = 0
    liner = (~sf) & (~inter) & (o3hb > -1) & (o3hb < cidfer10(n2ha))
    BPT[liner] = 1
    seyfert = (~sf) & (~inter) & (~liner) & (o3hb > -1)
    BPT[seyfert] = 2

#     BPT = np.zeros(len(n2ha))
#     for i in range(len(n2ha)):
#         if (n2ha[i] > -1.5 and n2ha[i] < -0.1 and 
#                 o3hb[i] < kauffm03(n2ha[i]) and abs(ew_ha[i]) > 6.0):
#             BPT[i] = -1  # Starforming, below Kauffmann line
#         elif n2ha[i] < 0.3 and o3hb[i] < kewley01(n2ha[i]):
#             BPT[i] = 0   # Intermediate, between Kewley & Kauffmann
#         elif o3hb[i] > -1  and o3hb[i] < cidfer10(n2ha[i]):
#             BPT[i] = 1   # LINER
#         elif o3hb[i] > -1:
#             BPT[i] = 2   # Seyfert
#         else:
#             BPT[i] = np.nan
    return BPT


# Metallicity derived from Marino+13 calibration.
# use method='o3n2' or method='n2'
# Require star-forming in BPT diagram.
# Input is a table containing the appropriate columns.
def ZOH_M13(fluxtab, method='o3n2', name='ZOH', err=False):

    N2F = fluxtab['flux_[NII]6583']
    O3F = fluxtab['flux_[OIII]5007']
    HaF = fluxtab['flux_Halpha']
    HbF = fluxtab['flux_Hbeta']
    if method == 'o3n2':
        good = (N2F>0) & (O3F>0) & (HaF>0) & (HbF>0)
    elif method == 'n2':
        good = (N2F>0) & (HaF>0)
    else:
        raise Exception('Method {} is not recognized'.format(method))
    nelt = len(N2F)

    BPT = bpt_type(N2F, O3F, HaF, HbF, fluxtab['EW_Halpha'])

    if err == False:
        O3N2 = np.full(nelt, np.nan)
        O3N2[good] = (np.log10(O3F[good]) - np.log10(HbF[good]) 
                    - (np.log10(N2F[good]) - np.log10(HaF[good])))                    
        N2 = np.full(nelt, np.nan)
        N2[good] = np.log10(N2F[good]) - np.log10(HaF[good])                   
    else:
        try:
            from uncertainties import ufloat
            import uncertainties.unumpy as unp 
        except:
            print('uncertainties package required for err=True')

        uN2F = unp.uarray(N2F, fluxtab['e_flux_[NII]6583'])
        uO3F = unp.uarray(O3F, fluxtab['e_flux_[OIII]5007'])
        uHaF = unp.uarray(HaF, fluxtab['e_flux_Halpha'])
        uHbF = unp.uarray(HbF, fluxtab['e_flux_Hbeta'])
        O3N2 = unp.uarray(np.full(nelt, np.nan),np.full(nelt, np.nan))
        O3N2[good] = (unp.log10(uO3F[good]) - unp.log10(uHbF[good]) 
                    - (unp.log10(uN2F[good]) - unp.log10(uHaF[good])))
        N2 = unp.uarray(np.full(nelt, np.nan),np.full(nelt, np.nan))
        N2[good] = unp.log10(uN2F[good]) - unp.log10(uHaF[good])

    if method == 'o3n2':
        ZOH_M13 = 8.533 - 0.214 * O3N2 # Eq(2) from Marino+2013
    else:
        ZOH_M13 = 8.743 + 0.462 * N2   # Eq(4) from Marino+2013

    if err == False: 
        ZOH_M13[BPT != -1] = np.nan
        return Column(ZOH_M13, name=name)
    else:            
        ZOH_M13[BPT != -1] = ufloat(np.nan, np.nan)
        return (Column(unp.nominal_values(ZOH_M13), name=name), 
                Column(unp.std_devs(ZOH_M13), name=name+'_err'))


# Metallicity derived from N2 line ratio, Marino+13 calibration.
# Input is a table containing the appropriate columns.
# def ZOH_n2(fluxtab, err=False):
#     nelt = len(fluxtab['flux_[NII]6583'])
#     uN2F = unp.uarray(fluxtab['flux_[NII]6583'], fluxtab['e_flux_[NII]6583'])
#     uHaF = unp.uarray(fluxtab['flux_Halpha'], fluxtab['e_flux_Halpha'])
# 
#     uN2 = unp.uarray(np.full(nelt, np.nan),np.full(nelt, np.nan))
#     good = ((unp.nominal_values(uN2F)>0) & (unp.nominal_values(uHaF)>0))
#     uN2[good] = (unp.log10(uN2F[good]) - unp.log10(uHaF[good]))
# 
#     BPT = bpt_type(fluxtab['flux_[NII]6583'], fluxtab['flux_[OIII]5007'], 
#             fluxtab['flux_Halpha'], fluxtab['flux_Halpha'], fluxtab['EW_Halpha'])
#             
#     uOH_N2 = 8.743 + 0.462*uN2  
#     uOH_N2[BPT != -1] = ufloat(np.nan, np.nan)
# 
#     if err: 
#         return unp.std_devs(uOH_N2)
#     else:            
#         return unp.nominal_values(uOH_N2)
