import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy import units as u
from astropy.coordinates import SkyCoord
#from astropy.modeling.functional_models import Gaussian2D
from scipy.optimize import fsolve
from scipy import ndimage
from scipy.stats import multivariate_normal as ndNormal


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
    # NEW: Returns A_Ha column as well.
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
        return Column(sig_sfr, name=name, dtype='f4',
            description='BD corrected SFR surface density'), Column(A_Ha, 
            name='AHa_'+name, dtype='f4', unit='mag', 
            description='Ha extinction from BD')
    else:
        return sig_sfr, A_Ha
    

# Convert CO intensity to H2(+He) surface density
def msd_co(sb_co, alphaco=4.3, name='sig_mol'):
    convfac = alphaco * (u.solMass/u.pc**2) / (u.K*u.km/u.s)
    sig_mol = (convfac*sb_co).to(u.solMass/u.pc**2)
    if isinstance(sb_co, Column):
        return Column(sig_mol, name=name, dtype='f4',
                      description='mol gas surface density')
    else:
        return sig_mol
    

# Convert units for stellar surface density
def stmass_pc2(stmass_as2, dist=10*u.Mpc, name='sig_star'):
    # Assume Mpc units if not given
    try:
        unit = dist.unit
    except:
        dist = dist * u.Mpc
    sterad = (u.sr/u.arcsec**2).decompose()   # 206265^2
    pxarea = (dist**2/sterad).to(u.pc**2)
    if isinstance(stmass_as2, Column):
        stmass_ary = np.array(stmass_as2)
        stmass_ary[~np.isfinite(stmass_ary)] = np.nan
        stmass_pc2 = 10**stmass_ary * u.solMass / pxarea
        stmass_pc2[~np.isfinite(stmass_pc2)] = np.nan
        return Column(stmass_pc2, name=name, dtype='f4',
                      description='stellar mass surface density')
    else:
        stmass_as2[~np.isfinite(stmass_as2)] = np.nan
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


# TODO not a formal position to import 
#from astropy.modeling.functional_models import Gaussian2D


def findIntersec(fun1, fun2, x0):
    return fsolve(lambda x:fun1(x)-fun2(x), x0)


def bpt_prob(n2ha_u, o3hb_u, bpt_type, grid_size=None):
    '''
    @Parameters
    criterions: a list of functions in BPT diagram
    x_u, y_u: data point with uncertainty
    grid_size: the size of the square grid where normal dist constructed
    
    @output
    probability of data points in regions
    '''
    try:
        from uncertainties import unumpy, umath
    except:
        print("uncertainties package is required for this function")
    x = unumpy.nominal_values(n2ha_u)
    y = unumpy.nominal_values(o3hb_u)
    x_std = unumpy.std_devs(n2ha_u)
    y_std = unumpy.std_devs(o3hb_u)
    if not grid_size:
        grid_size = 5
    x_arr, y_arr = np.meshgrid(np.linspace(x - x_std, x + x_std, grid_size),
                 np.linspace(y - y_std, y + y_std, grid_size))
    pos = np.dstack((x_arr, y_arr))
    kewley01 = lambda nii: 1.19 + 0.61/(nii - 0.47) # Eq. 5 of 2001ApJ...556..121K
    kauffm03 = lambda nii: 1.30 + 0.61/(nii - 0.05) # Eq. 1 of 2003MNRAS.346.1055K
    cidfer10 = lambda nii: 0.48 + 1.01*nii          # Eq. 3 of 2010MNRAS.403.1036C
    kew = kewley01(x_arr[0])
    kau = kauffm03(x_arr[0])
    cid = cidfer10(x_arr[0])
    grid = np.zeros((grid_size, grid_size))   
    # HII star forming
    sf = np.logical_and(y_arr[:, 0] < kau, x_arr[0] < -0.1)
    grid[:, sf] = -1 
    # Composite
    inter = np.logical_and.reduce([y_arr[:, 0] > kau, y_arr[:, 0] < kew, x_arr[0] < 0.3, ~sf])
    grid[:, inter] = 0 
    # LINER
    liner = np.logical_and.reduce([~sf, ~inter, y_arr[:, 0] > -1, y_arr[:, 0] < cid ])
    grid[:, liner] = 1 
    # Seyfert
    seyfert = np.logical_and.reduce([~sf, ~inter, ~liner, y_arr[:, 0] > -1])
    grid[:, seyfert] = 2
    
#     cid_region = np.where(x_arr[0, :]>  -0.19935)[0]
#     print(cid[cid_region])
#     print(cid_region)
#     print(y_arr[cid_region])
#     grid[x_arr[0][cid_region], y_arr[cid_region, 0] > cid[cid_region]] = 4
#     grid[x_arr[0][cid_region], y_arr[cid_region, 0] < cid[cid_region] and y_arr[cid_region, 0] > kew[cid_region]] = 5
    gauss2d = ndNormal(mean=(x, y), cov=(x_std, y_std))
    pdf = gauss2d.pdf(pos)
    normal_prob = np.zeros((grid_size, grid_size))
    delta_x = x_arr[0][1] - x_arr[0][0]
    delta_y = y_arr[1][0] - y_arr[0][0]
    total = 0
    for i in range(grid_size):
        for j in range(grid_size):
            normal_prob[i, j] = pdf[i, j] * delta_x * delta_y
            total += normal_prob[i, j]
            # Normalization constant
    return ndimage.sum(normal_prob * 1/total, grid, index=bpt_type)
    
        
# BPT classification, see Husemann et al. (2013A&A...549A..87H) Figure 7.
# Input is a flux_elines table.
def bpt_type(fluxtab, ext='', name='BPT', prob=False, grid_size=5):

    flux_nii  = fluxtab['flux_[NII]6583'+ext]
    flux_oiii = fluxtab['flux_[OIII]5007'+ext]
    flux_ha   = fluxtab['flux_Halpha'+ext]
    flux_hb   = fluxtab['flux_Hbeta'+ext]
    ew_ha     = fluxtab['EW_Halpha'+ext]

    good = (flux_nii>0) & (flux_oiii>0) & (flux_ha>0) & (flux_hb>0) & (~np.isnan(ew_ha))
    n2ha = np.full(len(flux_nii), np.nan)
    n2ha[good] = np.log10(flux_nii[good])  - np.log10(flux_ha[good])
    o3hb = np.full(len(flux_oiii), np.nan)
    o3hb[good] = np.log10(flux_oiii[good]) - np.log10(flux_hb[good])   

    kewley01 = lambda nii: 1.19 + 0.61/(nii - 0.47) # Eq. 5 of 2001ApJ...556..121K
    kauffm03 = lambda nii: 1.30 + 0.61/(nii - 0.05) # Eq. 1 of 2003MNRAS.346.1055K
    cidfer10 = lambda nii: 0.48 + 1.01*nii          # Eq. 3 of 2010MNRAS.403.1036C

    kewley_start = findIntersec(kewley01, kauffm03, -1)[0]  # -1.2805
    cidfer_start = findIntersec(kewley01, cidfer10, 0)[0]   # -0.1993
    kewley_end = findIntersec(kewley01, lambda x: -4, 0)[0] #  0.352466
    kauffm_end = findIntersec(kauffm03, lambda x: -4, 0)[0] # -0.06509

    BPT = np.full(len(n2ha), np.nan)
    # Star forming: below Kauffmann line and EW > 6
    sf = (n2ha < kauffm_end) & (o3hb < kauffm03(n2ha)) & (abs(ew_ha) > 6.0)
    BPT[sf] = -1
    # Intermediate: below Kewley line and not star-forming
    inter = (~sf) & (n2ha < kewley_end) & (o3hb < kewley01(n2ha))
    BPT[inter] = 0
    # LINER: above Kewley line and below Cid Fernandes line
    liner = (~sf) & (~inter) & (o3hb < cidfer10(n2ha))
    BPT[liner] = 1
    # Seyfert: above Kewley line and above Cid Fernandes line
    seyfert = (~sf) & (~inter) & (~liner) & good
    BPT[seyfert] = 2
    bpt_col = Column(BPT, name=name, dtype='f4', format='.1f',
                description='BPT type (-1=SF 0=inter 1=LINER 2=Sy)')

    if prob:
        try:
            from uncertainties import unumpy, umath
        except:
            print("uncertainties package is required for prob calculation")
       
        BPT_prob = np.full(len(n2ha), np.nan)
        eN2 = fluxtab['e_flux_[NII]6583'+ext]
        eO3 = fluxtab['e_flux_[OIII]5007'+ext] 
        eHa = fluxtab['e_flux_Halpha'+ext]
        eHb = fluxtab['e_flux_Hbeta'+ext]

        Ha_u = unumpy.uarray(np.array(flux_ha), np.array(eHa))
        Hb_u = unumpy.uarray(np.array(flux_hb), np.array(eHb))
        N2_u = unumpy.uarray(np.array(flux_nii), np.array(eN2))
        O3_u = unumpy.uarray(np.array(flux_oiii), np.array(eO3))
        t1 = [umath.log10(N2_u[good][i]) - umath.log10(Ha_u[good][i]) for i in range(len(Ha_u[good]))]
        t2 = [umath.log10(O3_u[good][i]) - umath.log10(Hb_u[good][i]) for i in range(len(Hb_u[good]))]
        n2ha_u = unumpy.uarray(np.full(len(Ha_u), np.nan), np.full(len(Ha_u), np.nan))
        o3hb_u = unumpy.uarray(np.full(len(Hb_u), np.nan), np.full(len(Hb_u), np.nan))
        n2ha_u[good] = unumpy.uarray(unumpy.nominal_values(t1), unumpy.std_devs(t1))
        o3hb_u[good] = unumpy.uarray(unumpy.nominal_values(t2), unumpy.std_devs(t2)) 
        for i in np.where(sf)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], -1, grid_size)
        for i in np.where(inter)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], 0, grid_size)
        for i in np.where(liner)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], 1, grid_size)
        for i in np.where(seyfert)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], 2, grid_size)
        prob_col = Column(BPT_prob, name='prob', dtype='f4', description='BPT probability')
        return bpt_col, prob_col
    else:
        return bpt_col


# Metallicity derived from Marino+13 calibration.
# use method='o3n2' or method='n2'
# Require star-forming in BPT diagram.
# Input is a flux_elines table.
def ZOH_M13(fluxtab, ext='', method='o3n2', name='ZOH', err=False):

    N2F = fluxtab['flux_[NII]6583'+ext]
    O3F = fluxtab['flux_[OIII]5007'+ext]
    HaF = fluxtab['flux_Halpha'+ext]
    HbF = fluxtab['flux_Hbeta'+ext]
    if method == 'o3n2':
        good = (N2F>0) & (O3F>0) & (HaF>0) & (HbF>0)
    elif method == 'n2':
        good = (N2F>0) & (HaF>0)
    else:
        raise Exception('Method {} is not recognized'.format(method))
    nelt = len(N2F)

    BPT = bpt_type(fluxtab, ext=ext)

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

        uN2F = unp.uarray(N2F, fluxtab['e_flux_[NII]6583'+ext])
        uO3F = unp.uarray(O3F, fluxtab['e_flux_[OIII]5007'+ext])
        uHaF = unp.uarray(HaF, fluxtab['e_flux_Halpha'+ext])
        uHbF = unp.uarray(HbF, fluxtab['e_flux_Hbeta'+ext])
        O3N2 = unp.uarray(np.full(nelt, np.nan),np.full(nelt, np.nan))
        O3N2[good] = (unp.log10(uO3F[good]) - unp.log10(uHbF[good]) 
                    - (unp.log10(uN2F[good]) - unp.log10(uHaF[good])))
        N2 = unp.uarray(np.full(nelt, np.nan),np.full(nelt, np.nan))
        N2[good] = unp.log10(uN2F[good]) - unp.log10(uHaF[good])

    if method == 'o3n2':
        ZOH_M13 = 8.533 - 0.214 * O3N2 # Eq(2) from Marino+2013
    else:
        ZOH_M13 = 8.743 + 0.462 * N2   # Eq(4) from Marino+2013

    desc = '12+log(O/H) using {} method in Marino+13'.format(method)
    if err == False: 
        ZOH_M13[BPT != -1] = np.nan
        return Column(ZOH_M13, name=name, unit='dex', dtype='f4', description=desc)
    else:            
        ZOH_M13[BPT != -1] = ufloat(np.nan, np.nan)
        return (Column(unp.nominal_values(ZOH_M13), name=name, dtype='f4',
                       unit='dex', description=desc), 
                Column(unp.std_devs(ZOH_M13), name='e_'+name, dtype='f4',
                       unit='dex', description='error in '+desc))


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
