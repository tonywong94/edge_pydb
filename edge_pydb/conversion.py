import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import fsolve
from scipy import ndimage
from scipy.stats import multivariate_normal as ndNormal

try:
    import uncertainties.unumpy as unp
except (ImportError, ModuleNotFoundError) as error:
    print(error.__class__.__name__ + ": " + str(error))
    print("Could not import uncertainties package, \
        please try to install it or do not input the error as argument.")
except Exception as exception: 
    print(exception.__class__.__name__ + ": " + str(exception))


# BPT codes are defined here
STAR_FORMING = -1
INTERMEDIATE =  0
LINER        =  1
SEYFERT      =  2


def gc_polr(ra, dec, ra_gc, dec_gc, pa, inc, reject=89):
    '''
    Returns galactocentric polar coordinates (radius in arcsec, 
    azimuthal angle in degrees from receding major axis).

    === Parameters ===
    ra : float or numpy.array
        The RA coordinate(s) to transform, in degrees
    dec : float or numpy.array
        The DEC coordinate(s) to transform, in degrees
    ra_gc : float
        The RA of the galaxy center, in degrees
    dec_gc : float
        The DEC of the galaxy center, in degrees
    pa : float
        The position angle of the galaxy's receding major axis, in degrees
    inc : float
        The inclination of the galaxy disk, in degrees (0 = face-on)
    reject : float
        Return NaNs for inclinations larger than this (default=89)

    === Returns ===
    pair of floats or numpy.array [radius, azang]
    '''
    ctr = SkyCoord(ra_gc, dec_gc, unit="deg")
    pos = SkyCoord(ra, dec, unit="deg")
    # Polar vector in sky plane 
    t_sky = ctr.position_angle(pos).degree
    r_sky = ctr.separation(pos).arcsecond
    # Convert to galaxy plane
    x_sky = r_sky * np.cos(np.radians(t_sky - pa))
    y_sky = r_sky * np.sin(np.radians(t_sky - pa))
    x_gal = x_sky
    y_gal = y_sky/np.cos(np.radians(inc))
    # Reject very high inclinations (cannot be deprojected)
    y_gal[inc>reject] = np.nan
    radius = np.sqrt(x_gal**2 + y_gal**2)
    azang  = np.degrees(np.arctan2(y_gal, x_gal))
    return radius, azang


def uarray_to_list(target):
    '''
    Break array of uncertainties values into 2 separate arrays.

    === Parameters ===
    target : unumpy.uarray
        The array of values with uncertainties

    === Returns ===
    retval[0] -> nominal values
    retval[1] -> standard deviation
    '''
    return [unp.nominal_values(target), unp.std_devs(target)]


def get_AHa(flux_ha, flux_hb, log10):
    '''
    Returns Halpha extinction in magnitudes given flux_ha and flux_hb.

    === Parameters ===
    flux_ha : astropy.table.Column
        The H-alpha flux values
    flux_hb : astropy.table.Column
        The H-beta flux values
    log10 : function
        The log10 function, normally np.log10 but use unumpy.log10 for error propagation.

    === Returns ===
    numpy.ndarray representing AHa in magnitudes
    '''
    A_Ha = flux_ha.data * np.nan
    if flux_hb is not None:
        # Need positive flux to take the log.
        good = (flux_ha > 0) & (flux_hb > 0)
        # Extinction curve from Cardelli+ (1989).
        K_Ha = 2.53
        K_Hb = 3.61
        # Get A_Ha using Eq(1) from Catalan-Torrecilla+ (2015). 
        A_Ha[good] = K_Ha/(-0.4*(K_Ha-K_Hb)) * log10((flux_ha[good]/flux_hb[good])/2.86)
    return A_Ha


def sfr_ha(flux_ha, flux_hb=None, e_flux_ha=None, e_flux_hb=None, 
            name='sigsfr', column=True, imf='kroupa'):
    '''
    Convert Halpha intensity to SFR surface density, optionally
    with extinction estimates and corrections (if flux_hb is provided).
    Note that both e_flux_ha and e_flux_hb must not be None
    in order to propagate the error.  Otherwise the SFR is computed
    without dust correction (if flux_hb=None) or without error estimation. 

    === Parameters ===
    flux_ha : astropy.table.Column
        The H-alpha flux values
    flux_hb : astropy.table.Column
        The H-beta flux values.  If not provided, no dust correction is made.
    e_flux_ha : astropy.table.Column
        Errors in H-alpha flux values.  If not provided, error propagation skipped.
    e_flux_hb : astropy.table.Column
        Errors in H-beta flux values.
    name : string
        Name of the output Column; error column will be prepended with 'e_'
    column : boolean
        True to return astropy Column objects, otherwise return arrays
    imf : string
        'salpeter' to scale result by 1.51 (default uses Kroupa IMF)

    === Returns ===
    several columns depending on input (see code)
    '''
    # input line flux is actually flux per arcsec^2
    sterad = (u.sr/u.arcsec**2).decompose().scale # 206265^2
    
    # Calzetti+(2010ApJ...714.1256C), Kroupa IMF, 1 Gyr old pop
    lumcon = 5.45e-42 * (u.solMass/u.yr) / (u.erg/u.s)
    if imf == 'salpeter':
        lumcon *= 1.51

    # Case with error estimation
    if e_flux_ha is not None:
        u_flux_ha = unp.uarray(flux_ha, e_flux_ha)
        sig_sfr = (lumcon * 4*np.pi * u_flux_ha * sterad).to(u.solMass/(u.pc**2*u.Gyr))
        if flux_hb is None:   # no Balmer decrement
            sig_sfr_out = uarray_to_list(sig_sfr.value)
            if column:
                return Column(sig_sfr_out[0], name=name, dtype='f4', unit=sig_sfr.unit,
                    description='SFR surface density no extinction'), \
                    Column(sig_sfr_out[1], name='e_'+name, dtype='f4', unit=sig_sfr.unit,
                    description='error of uncorrected SFR surface density')
            else:
                return sig_sfr_out[0], sig_sfr_out[1]
        else:   # with Balmer decrement
            if e_flux_hb is None:
                e_flux_hb = np.zeros_like(flux_hb)
            u_flux_hb = unp.uarray(flux_hb, e_flux_hb)
            A_Ha = get_AHa(u_flux_ha, u_flux_hb, unp.log10)
            sig_sfr = sig_sfr * 10**(0.4*A_Ha)
            sig_sfr_out = uarray_to_list(sig_sfr.value)
            A_Ha_out = uarray_to_list(A_Ha.data)
            if column:
                return Column(sig_sfr_out[0], name=name, dtype='f4', unit=sig_sfr.unit,
                    description='BD corrected SFR surface density'), \
                    Column(A_Ha_out[0], name=name.replace('sigsfr','AHa'), 
                    dtype='f4', unit='mag', description='Ha extinction from BD'), \
                    Column(sig_sfr_out[1], name='e_'+name, dtype='f4', unit=sig_sfr.unit,
                    description='error of BD corrected SFR surface density'), \
                    Column(A_Ha_out[1], name='e_'+name.replace('sigsfr','AHa'), 
                    dtype='f4', unit='mag', description='error of Ha extinction from BD')
            else:
                return sig_sfr_out[0], A_Ha_out[0], sig_sfr_out[1], A_Ha_out[1]

    # Case with no error estimation
    else:
        sig_sfr = (lumcon * 4*np.pi * flux_ha * sterad).to(u.solMass/(u.pc**2*u.Gyr))
        if flux_hb is None:   # no Balmer decrement
            if column:
                return Column(sig_sfr, name=name, dtype='f4',
                    description='SFR surface density no extinction')
            else:
                return sig_sfr
        else:   # with Balmer decrement
            A_Ha = get_AHa(flux_ha, flux_hb, np.log10)
            sig_sfr = sig_sfr * 10**(0.4*A_Ha)
            if column:
                return Column(sig_sfr, name=name, dtype='f4',
                    description='BD corrected SFR surface density'), \
                    Column(A_Ha, name=name.replace('sigsfr','AHa'), dtype='f4', 
                    unit='mag', description='Ha extinction from BD')
            else:
                return sig_sfr, A_Ha


def msd_co(sb_co, alphaco=4.3, name='sigmol'):
    '''
    Convert CO intensity to H_2 (+ He) mass surface density.

    === Parameters ===
    sb_co : numpy.array or astropy.table.Column
        CO integrated intensity values, in K km/s
    alphaco : float
        CO to H2 conversion factor, from K km/s to Msol/pc2
    name : string
        The name of the output Column, if input is a Column

    === Returns ===
    numpy.array or astropy.table.Column, matching input
    '''
    convfac = alphaco * (u.solMass/u.pc**2) / (u.K*u.km/u.s)
    sig_mol = (convfac*sb_co).to(u.solMass/u.pc**2)
    if isinstance(sb_co, Column):
        return Column(sig_mol, name=name, dtype='f4',
                      description='mol gas surface density')
    else:
        return sig_mol
    

def stmass_pc2(stmass_as2, dz=None, dist=10*u.Mpc, name='sigstar'):
    '''
    Convert units for stellar surface density to Msol/pc2, optionally 
    applying dezonification image.
    
    === Parameters ===
    stmass_as2 : numpy.array or astropy.table.Column
        stellar surface density in Pipe3D units
    dz : numpy.array or astropy.table.Column
        dezonification image from Pipe3D (optional)
    dist : float or astropy.Quantity
        The distance of the galaxy, assumed to be in Mpc if no units
    name : string
        The name of the output Column, if input is a Column

    === Returns ===
    numpy.array or astropy.table.Column, depending on input
    '''
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
    '''
    Determine BPT classification based on [NII]/Ha and [OIII]/Hb.  Does 
    not impose a requirement on EW(Ha).

    === Parameters ===
    n2ha : numpy.ndarray
        log10([NII]/Halpha)
    o3hb : numpy.ndarray
        log10([OIII]/Hbeta)
    good : boolean
        mask to apply to output arrays    

    === Returns ===
    boolean index arrays for SF, intermediate, LINER, and Seyfert.
    '''
    # These are constant values determined by the curve equations
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
    retval = np.full_like(n2ha, fill_value=np.nan)
    retval[sf]      = STAR_FORMING
    retval[inter]   = INTERMEDIATE
    retval[liner]   = LINER
    retval[seyfert] = SEYFERT
    return retval


def bpt_prob(n2ha_u, o3hb_u, bpt_type, grid_size=5):
    '''
    Calculate the probability that a measurement with uncertainty is in a 
    particular region of the BPT diagram, as specified by the parameter bpt_type.  
    This is done by constructing a grid and measuring the normalized Gaussian 
    area within the chosen BPT region.

    === Parameters ===
    n2ha_u : ufloat
        log10([NII]/Halpha) and uncertainty
    o3hb_u : ufloat
        log10([OIII]/Hbeta) and uncertainty
    bpt_type : float
        BPT type to measure prob for, -1=SF, 0=intermediate, 1=LINER, 2=Seyfert
    grid_size : float
        The number of grid points placed between +/- one sigma of the mean
    
    === Returns ===
    probability that measurement is in the specified region
    '''
    x = unp.nominal_values(n2ha_u)
    y = unp.nominal_values(o3hb_u)
    x_std = unp.std_devs(n2ha_u)
    y_std = unp.std_devs(o3hb_u)
    x_arr, y_arr = np.meshgrid( np.linspace(x - x_std, x + x_std, grid_size),
                                np.linspace(y - y_std, y + y_std, grid_size) )
    pos = np.dstack((x_arr, y_arr))
    grid = np.zeros((grid_size, grid_size))  
    grid = bpt_region(x_arr, y_arr)
    gauss2d = ndNormal(mean=(x, y), cov=(x_std**2, y_std**2))
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


def bpt_type(fluxtab, ext='', name='BPT', sf=True, prob=False, grid_size=5):
    '''
    Adds BPT classification to the flux_elines table. BPT==-1 means in the
    star-forming region of the diagram.  An additional column SF_BPT is
    set to True when BPT==-1 and Halpha EW > 6. 
    This function should be run before ZOH_M13.
    For more information see Husemann et al. (2013A&A...549A..87H) Figure 7.

    === Parameters ===
    fluxtab : astropy.Table
        flux_elines table extracted from Pipe3D output
    ext : string
        suffix for selected column names, e.g. '_rg' or '_sm'
    name : string
        name of the output column
    sf : boolean
        True to also calculate the SF_BPT column using Halpha EW.
    prob : boolean
        True to also calculate the BPT probabilities as an additional column
    grid_size : float
        The size of the square grid where BPT probabilities constructed
    
    === Returns ===
    if sf==False and prob==False:
        one astropy.table.Column [BPT]
    if sf==False and prob==True:
        two astropy.table.Column [BPT, p_BPT]
    if sf==True and prob==False:
        two astropy.table.Column [BPT, SF_BPT]
    if sf==True and prob==True:
        three astropy.table.Column [BPT, SF_BPT, p_BPT]
    '''
    flux_nii  = fluxtab['flux_[NII]6583'+ext]
    flux_oiii = fluxtab['flux_[OIII]5007'+ext]
    flux_ha   = fluxtab['flux_Halpha'+ext]
    flux_hb   = fluxtab['flux_Hbeta'+ext]

    good = (flux_nii>1e-5) & (flux_oiii>1e-5) & (flux_ha>1e-5) & (flux_hb>1e-5)
    n2ha = np.full(len(flux_nii), np.nan)
    n2ha[good] = np.log10(flux_nii[good])  - np.log10(flux_ha[good])
    o3hb = np.full(len(flux_oiii), np.nan)
    o3hb[good] = np.log10(flux_oiii[good]) - np.log10(flux_hb[good])   

    BPT = bpt_region(n2ha, o3hb, good)
    bpt_col = Column(BPT, name=name, dtype='f4', format='.1f',
                description='BPT type (-1=SF 0=inter 1=LINER 2=Sy)')
    if sf:
        ew_ha = fluxtab['EW_Halpha'+ext]
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
        Ha_u = unp.uarray(np.array(flux_ha), np.array(eHa))
        Hb_u = unp.uarray(np.array(flux_hb), np.array(eHb))
        N2_u = unp.uarray(np.array(flux_nii), np.array(eN2))
        O3_u = unp.uarray(np.array(flux_oiii), np.array(eO3))
        t1 = [unp.log10(N2_u[good][i]) - unp.log10(Ha_u[good][i]) for i in range(len(Ha_u[good]))]
        t2 = [unp.log10(O3_u[good][i]) - unp.log10(Hb_u[good][i]) for i in range(len(Hb_u[good]))]
        n2ha_u = unp.uarray(np.full(len(Ha_u), np.nan), np.full(len(Ha_u), np.nan))
        o3hb_u = unp.uarray(np.full(len(Hb_u), np.nan), np.full(len(Hb_u), np.nan))
        n2ha_u[good] = unp.uarray(unp.nominal_values(t1), unp.std_devs(t1))
        o3hb_u[good] = unp.uarray(unp.nominal_values(t2), unp.std_devs(t2)) 
        print("Working on SFR")
        for i in np.where(BPT == STAR_FORMING)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], STAR_FORMING, grid_size)
        print("Working on intermediate region, composite")
        for i in np.where(BPT == INTERMEDIATE)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], INTERMEDIATE, grid_size)
        print("Working on liner")
        for i in np.where(BPT == LINER)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], LINER, grid_size)
        print("Working on Seyfert")
        for i in np.where(BPT == SEYFERT)[0]:
            BPT_prob[i] = bpt_prob(n2ha_u[i], o3hb_u[i], SEYFERT, grid_size)
        prob_col = Column(BPT_prob, name='p_'+name, dtype='f4', description='BPT probability')

    if sf and prob:
        return bpt_col, bpt_sfcol, prob_col
    elif sf:
        return bpt_col, bpt_sfcol
    elif prob:
        return bpt_col, prob_col
    else:
        return bpt_col


def ZOH_M13(fluxtab, ext='', method='o3n2', name='ZOH', err=True):
    '''
    Adds gas-phase metallicity from Marino+13 calibration to flux_elines table.
    This function should be run after bpt_type since it requires spaxel to be
    star-forming in the BPT diagram.
    Both O3N2 and N2 methods are supported.

    === Parameters ===
    fluxtab : astropy.Table
        flux_elines table extracted from Pipe3D output
    ext : string
        suffix for selected column names, e.g. '_rg' or '_sm'
    method : string
        choose 'o3n2' (default) or 'n2'
    err : boolean
        True to calculate uncertainty as an additional column
    
    === Returns ===
    if err==False:
        astropy.table.Column ZOH
    if err==True:
        two astropy.table.Column [ZOH, e_ZOH]
    '''
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
    
    # Require SF in BPT if available.
    if 'SF_BPT'+ext in fluxtab.colnames:
        BPT_sf = fluxtab['SF_BPT'+ext]
        good = good & BPT_sf
    else:
        print('ZOH_M13 warning: The SF_BPT{} column is missing'.format(ext))

    nelt = len(N2F)

    if err == False:
        O3N2 = np.full(nelt, np.nan)
        O3N2[good] = (np.log10(O3F[good]) - np.log10(HbF[good]) 
                    - (np.log10(N2F[good]) - np.log10(HaF[good])))                    
        N2 = np.full(nelt, np.nan)
        N2[good] = np.log10(N2F[good]) - np.log10(HaF[good])                   
    else:
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
        return Column(ZOH_M13, name=name, unit='dex', dtype='f4', description=desc)
    else:            
        return (Column(unp.nominal_values(ZOH_M13), name=name, dtype='f4',
                       unit='dex', description=desc), 
                Column(unp.std_devs(ZOH_M13), name='e_'+name, dtype='f4',
                       unit='dex', description='error in '+desc))

