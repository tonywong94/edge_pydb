import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, join
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import fsolve
from scipy import ndimage
from scipy.stats import multivariate_normal as ndNormal

try:
    from uncertainties import unumpy, umath
except (ImportError, ModuleNotFoundError) as error:
    # ? should we require the uncertainties package in the requirement.txt? 
    print(error.__class__.__name__ + ": " + str(error))
    print("Could not import uncertainties package, \
        please try to install it or do not input the error as argument.")
except Exception as exception: 
    print(exception.__class__.__name__ + ": " + str(exception))


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


def uarray_to_list(target):
    '''
    break array of uncertainties values into 2 lists, 
    retval[0] -> nominal values
    retval[1] -> standard deviation
    if no value for std, then return None
    '''
    result = map(lambda x: (x.n, x.s) if not (isinstance(x, float) 
                    or isinstance(x, int)) else (x, None), target)
    return list(map(list, zip(*list(result))))


# Convert Halpha intensity to A_V-corrected SFR surface density
def sfr_ha(flux_ha, flux_hb=None, e_flux_ha=None, e_flux_hb=None, 
            name='sig_sfr', column=True, filter_bad=True,
            imf='kroupa'):
    '''
    Note that both e_flux_ha and e_flux_hb have to be not None
    in order to propagate the error
    '''
    # input line flux is actually flux per arcsec^2
    sterad = (u.sr/u.arcsec**2).decompose().scale # 206265^2
    
    # Calzetti+(2010ApJ...714.1256C), Kroupa IMF, 1 Gyr old pop
    lumcon = 5.45e-42 * (u.solMass/u.yr) / (u.erg/u.s)
    if imf == 'salpeter':
        lumcon *= 1.51

    # If only flux_ha exists, just do a scaling
    if flux_hb is None:
        sb_ha  = flux_ha * sterad  # flux per steradian
        lsd_ha = 4 * np.pi * sb_ha
        sig_sfr = (lumcon * lsd_ha).to(u.solMass/(u.pc**2*u.Gyr))
        if e_flux_ha is not None:
            e_sig_sfr = e_flux_ha / flux_ha * sig_sfr
            return sig_sfr, e_sig_sfr
        else:
            return sig_sfr

    def apply_extinction(flux_ha, flux_hb, log10):
        good = True
        if filter_bad:
            good = (flux_ha > 0) & (flux_hb > 0)
        # Extinction curve from Cardelli+(1989).
        K_Ha = 2.53
        K_Hb = 3.61
        # Get A_Ha using Eq(1) from Catalan-Torrecilla+(2015). 
        # By default A_Ha is NaN when flux_ha is NaN.
        A_Ha = flux_ha * 0.
        A_Ha[good] = K_Ha/(-0.4*(K_Ha-K_Hb)) * log10((flux_ha[good]/flux_hb[good])/2.86)
        # Do not apply negative extinction.
        A_Ha[A_Ha < 0] = 0.
        flux_ha_cor = flux_ha * 10**(0.4*A_Ha)
        # TODO make it as a separate function and do the convert
        sb_ha  = flux_ha_cor * sterad   # flux per steradian
        lsd_ha = 4 * np.pi * sb_ha
        sig_sfr = (lumcon * lsd_ha).to(u.solMass/(u.pc**2*u.Gyr))
        return sig_sfr, A_Ha

    if e_flux_ha is not None:
        if e_flux_hb is None:
            e_flux_hb = np.zeros_like(flux_hb)
        u_flux_ha = unumpy.uarray(flux_ha, e_flux_ha)
        u_flux_hb = unumpy.uarray(flux_hb, e_flux_hb)
        sig_sfr, A_Ha = apply_extinction(u_flux_ha, u_flux_hb, unumpy.log10)
        # * sig_sfr is a astropy Quantity, and A_Ha is still a Column
        sig_sfr_out = uarray_to_list(sig_sfr.value)
        A_Ha_out = uarray_to_list(A_Ha.data)
        if column:
            return Column(sig_sfr_out[0], name=name, dtype='f4', unit=sig_sfr.unit,
                description='BD corrected SFR surface density'), \
                Column(A_Ha_out[0], name='AHa_'+name, dtype='f4', unit='mag', 
                description='Ha extinction from BD'), \
                Column(sig_sfr_out[1], name='e_'+name, dtype='f4', unit=sig_sfr.unit,
                description='error of BD corrected SFR surface density'), \
                Column(A_Ha_out[1], name='e_AHa_'+name, dtype='f4', unit='mag', 
                description='error of Ha extinction from BD')
        else:
            return sig_sfr_out[0], A_Ha_out[0], sig_sfr_out[1], A_Ha_out[1]
    else:
        sig_sfr, A_Ha = apply_extinction(flux_ha, flux_hb, np.log10)
        if column:
            return Column(sig_sfr, name=name, dtype='f4',
                description='BD corrected SFR surface density'), \
                Column(A_Ha, name='AHa_'+name, dtype='f4', unit='mag', 
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
# dz = dezonification file from Pipe3D
def stmass_pc2(stmass_as2, dz=None, dist=10*u.Mpc, name='sig_star'):
    # Assume Mpc units if not given
    try:
        unit = dist.unit
    except:
        dist = dist * u.Mpc
    sterad = (u.sr/u.arcsec**2).decompose().scale   # 206265^2
    pxarea = (dist**2/sterad).to(u.pc**2)

    def convert_stmass(stmass_in, dz):
        stmass_in[~np.isfinite(stmass_in)] = np.nan
        stmass_pc2 = 10**stmass_in * u.solMass / pxarea
        if dz is not None:
            stmass_pc2 *= np.array(dz)
        stmass_pc2[~np.isfinite(stmass_pc2)] = np.nan
        return stmass_pc2

    if isinstance(stmass_as2, Column):
        stmass_ary = np.array(stmass_as2)
        stmass_pc2 = convert_stmass(stmass_ary, dz)
        return Column(stmass_pc2, name=name, dtype='f4',
                      description='stellar mass surface density')
    else:
        stmass_pc2 = convert_stmass(stmass_as2, dz)
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


def findIntersec(fun1, fun2, x0):
    return fsolve(lambda x:fun1(x)-fun2(x), x0)


def kewley01(nii):
    # Eq. 5 of 2001ApJ...556..121K
    return 1.19 + 0.61/(nii - 0.47) 


def kauffm03(nii):
    # Eq. 1 of 2003MNRAS.346.1055K
    return 1.30 + 0.61/(nii - 0.05)


def cidfer10(nii):
    # Eq. 3 of 2010MNRAS.403.1036C
    return 0.48 + 1.01*nii    


def bpt_region(n2ha, o3hb, good=True):
    # just to make sure ew > 6 to take all the data if ew_ha is not available
    # ? might be a good idea to store these values, as they are constants
    kewley_start = findIntersec(kewley01, kauffm03, -1)[0]  # -1.2805
    cidfer_start = findIntersec(kewley01, cidfer10, 0)[0]   # -0.1993
    kewley_end = findIntersec(kewley01, lambda x: -4, 0)[0] #  0.352466
    kauffm_end = findIntersec(kauffm03, lambda x: -4, 0)[0] # -0.06509
    # Star forming: below Kauffmann line
    sf = (n2ha < kauffm_end) & (o3hb < kauffm03(n2ha))
    # Intermediate: below Kewley line and not star-forming
    inter = (~sf) & (n2ha < kewley_end) & (o3hb < kewley01(n2ha))
    # LINER: above Kewley line and below Cid Fernandes line
    liner = (~sf) & (~inter) & (o3hb < cidfer10(n2ha))
    # Seyfert: above Kewley line and above Cid Fernandes line
    seyfert = (~sf) & (~inter) & (~liner) & good
    return sf, inter, liner, seyfert     


def bpt_prob(n2ha_u, o3hb_u, bpt_type, grid_size=None):
    '''
    Calculate the probability of a single point in the **bpt_type** region of the BPT plot.
    This is done by first constructing a grid with **grid_size** resolution. Then, the grid is marked by 
    three dividing lines (Kewley, Kauffmann, Cid Fernandes)

    Parameters
    ----------
    n2ha_u, o3hb_u : one data point with uncertainty
    bpt_type : the bpt_type we are considering, -1 for star forming, 0 for intermediate, 1 for liner, 2 for seyfert
    grid_size : the size of the square grid where normal dist constructed
    
    Returns
    -------
    probability of data points in regions
    '''
    x = unumpy.nominal_values(n2ha_u)
    y = unumpy.nominal_values(o3hb_u)
    x_std = unumpy.std_devs(n2ha_u)
    y_std = unumpy.std_devs(o3hb_u)
    if not grid_size:
        grid_size = 5
    x_arr, y_arr = np.meshgrid(np.linspace(x - x_std, x + x_std, grid_size),
                 np.linspace(y - y_std, y + y_std, grid_size))
    pos = np.dstack((x_arr, y_arr))
    grid = np.zeros((grid_size, grid_size))   
    sf, inter, liner, seyfert = bpt_region(x_arr[0], y_arr[:, 0])
    grid[:, sf] = -1 
    grid[:, inter] = 0 
    grid[:, liner] = 1 
    grid[:, seyfert] = 2
    
    gauss2d = ndNormal(mean=(x, y), cov=(x_std, y_std))
    pdf = gauss2d.pdf(pos)
    normal_prob = np.zeros((grid_size, grid_size))
    delta_x = x_arr[0][1] - x_arr[0][0]
    delta_y = y_arr[1][0] - y_arr[0][0]
    total = 0
    # the integrand
    normal_prob = pdf * delta_x * delta_y
    # Normalization constant, count how many cells we have
    total = np.sum(normal_prob)
    # do the sum reduction
    return ndimage.sum(normal_prob * 1/total, grid, index=bpt_type)



# BPT classification, see Husemann et al. (2013A&A...549A..87H) Figure 7.
# Input is a flux_elines table.
def bpt_type(fluxtab, ext='', name='BPT', prob=False, grid_size=5):

    flux_nii  = fluxtab['flux_[NII]6583'+ext]
    flux_oiii = fluxtab['flux_[OIII]5007'+ext]
    flux_ha   = fluxtab['flux_Halpha'+ext]
    flux_hb   = fluxtab['flux_Hbeta'+ext]
    ew_ha     = fluxtab['EW_Halpha'+ext]

    good = (flux_nii>0) & (flux_oiii>0) & (flux_ha>0) & (flux_hb>0) #& (~np.isnan(ew_ha))
    n2ha = np.full(len(flux_nii), np.nan)
    n2ha[good] = np.log10(flux_nii[good])  - np.log10(flux_ha[good])
    o3hb = np.full(len(flux_oiii), np.nan)
    o3hb[good] = np.log10(flux_oiii[good]) - np.log10(flux_hb[good])   

    sf, inter, liner, seyfert = bpt_region(n2ha, o3hb, good)

    BPT = np.full(len(n2ha), np.nan)
    BPT[sf] = -1
    BPT[inter] = 0
    BPT[liner] = 1
    BPT[seyfert] = 2
    bpt_col = Column(BPT, name=name, dtype='f4', format='.1f',
                description='BPT type (-1=SF 0=inter 1=LINER 2=Sy)')
    bpt_sf = (BPT == -1) & (abs(ew_ha)> 6.0)
    bpt_sfcol = Column(bpt_sf, name='SF_' + name, dtype='?',
                description='True if star forming (BPT=-1 and EW_Ha>6)')
    if prob:
        BPT_prob = np.full(len(n2ha), np.nan)
        eN2 = fluxtab['e_flux_[NII]6583'+ext]
        eO3 = fluxtab['e_flux_[OIII]5007'+ext] 
        eHa = fluxtab['e_flux_Halpha'+ext]
        eHb = fluxtab['e_flux_Hbeta'+ext]
        print("Calculating the probability")
        Ha_u = unumpy.uarray(np.array(flux_ha), np.array(eHa))
        Hb_u = unumpy.uarray(np.array(flux_hb), np.array(eHb))
        N2_u = unumpy.uarray(np.array(flux_nii), np.array(eN2))
        O3_u = unumpy.uarray(np.array(flux_oiii), np.array(eO3))
        t1 = [unumpy.log10(N2_u[good][i]) - unumpy.log10(Ha_u[good][i]) for i in range(len(Ha_u[good]))]
        t2 = [unumpy.log10(O3_u[good][i]) - unumpy.log10(Hb_u[good][i]) for i in range(len(Hb_u[good]))]
        n2ha_u = unumpy.uarray(np.full(len(Ha_u), np.nan), np.full(len(Ha_u), np.nan))
        o3hb_u = unumpy.uarray(np.full(len(Hb_u), np.nan), np.full(len(Hb_u), np.nan))
        n2ha_u[good] = unumpy.uarray(unumpy.nominal_values(t1), unumpy.std_devs(t1))
        o3hb_u[good] = unumpy.uarray(unumpy.nominal_values(t2), unumpy.std_devs(t2)) 
        print("Working on SFR")
        for i in np.where(sf)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], BPT[sf][0], grid_size)
        print("Working on intermediate region, composite")
        for i in np.where(inter)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], BPT[inter][0], grid_size)
        print("Working on liner")
        for i in np.where(liner)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], BPT[liner][0], grid_size)
        print("Working on Seyfert")
        for i in np.where(seyfert)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], BPT[seyfert][0], grid_size)
        prob_col = Column(BPT_prob, name='p_'+name, dtype='f4', description='BPT probability')
        return bpt_col, bpt_sfcol, prob_col
    else:
        return bpt_col, bpt_sfcol


def plot_uncertainty_ellipse(xval_u, yval_u, indices, x_arr, save_to=''):
    '''
    parameters
    xval_u, yval_u : list of coordinates with uncertainty 
    indices : indices of the list of coordinates to plot with
    save_to: file to save the plot to, optional
    '''
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    plt.figure(figsize=(8,8))
    ax = plt.gca()
    for i in indices:
        ax.add_patch(Ellipse(xy=(xval_u[i].n, yval_u[i].n),
                            width=xval_u[i].s, height=yval_u[i].s,
                            edgecolor='red', fc='white'))
    plt.plot(x_arr, kewley01(x_arr), 'k-.', label="Kewley")
    plt.plot(x_arr, kauffm03(x_arr), 'k--', label="Kauffmann")
    plt.plot(x_arr, cidfer10(x_arr), 'k-', label="Cidfer")
    if save_to:
        plt.savefig(save_to)
    plt.show()


# Metallicity derived from Marino+13 calibration.
# use method='o3n2' or method='n2'
# Require star-forming in BPT diagram.
# Input is a flux_elines table.
def ZOH_M13(fluxtab, ext='', method='o3n2', name='ZOH', err=True):

    N2F = fluxtab['flux_[NII]6583'+ext]
    O3F = fluxtab['flux_[OIII]5007'+ext]
    HaF = fluxtab['flux_Halpha'+ext]
    HbF = fluxtab['flux_Hbeta'+ext]
    BPT_sf = fluxtab['SF_BPT'+ext]

    if method == 'o3n2':
        good = (N2F>0) & (O3F>0) & (HaF>0) & (HbF>0) & BPT_sf
    elif method == 'n2':
        good = (N2F>0) & (HaF>0) & BPT_sf
    else:
        raise Exception('Method {} is not recognized'.format(method))
    nelt = len(N2F)

    #BPT = bpt_type(fluxtab, ext=ext)

    if err == False:
        O3N2 = np.full(nelt, np.nan)
        O3N2[good] = (np.log10(O3F[good]) - np.log10(HbF[good]) 
                    - (np.log10(N2F[good]) - np.log10(HaF[good])))                    
        N2 = np.full(nelt, np.nan)
        N2[good] = np.log10(N2F[good]) - np.log10(HaF[good])                   
    else:
        from uncertainties import ufloat
        import uncertainties.unumpy as unp 

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
        #ZOH_M13[BPT != -1] = np.nan
        return Column(ZOH_M13, name=name, unit='dex', dtype='f4', description=desc)
    else:            
        #ZOH_M13[BPT != -1] = ufloat(np.nan, np.nan)
        return (Column(unp.nominal_values(ZOH_M13), name=name, dtype='f4',
                       unit='dex', description=desc), 
                Column(unp.std_devs(ZOH_M13), name='e_'+name, dtype='f4',
                       unit='dex', description='error in '+desc))

