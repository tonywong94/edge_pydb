#!/usr/bin/env python

from astropy.coordinates import SkyCoord
from astropy.table import Table, join
import numpy as np
#import numpy.ma as ma

# Which tables to rewrite with this execution of the script
#rewrite=['wise','rdist','rfpars']
rewrite=['ned']

# Alberto's master table, last updated April 14, 2017
t = Table.read('inputs/DETableFinal.csv', format='ascii.csv')

# CALIFA Quality Tables
tq5 = Table.read('inputs/QCflags_std_V500_DR3.csv', format='ascii.csv', 
            names=('CALIFAID','Name','FLAG_OBS_SKYMAG','FLAG_OBS_EXT',
                'FLAG_OBS_AM','FLAG_RED_STRAYLIGHT','FLAG_RED_DISP',
                'FLAG_RED_CDISP','FLAG_RED_SKYLINES','FLAG_RED_LIMSB',
                'FLAG_RED_ERRSPEC','FLAG_CAL_SPECPHOTO','FLAG_CAL_WL',
                'FLAG_CAL_IMGQUAL','FLAG_CAL_SPECQUAL','FLAG_CAL_FLATSDSS',
                'FLAG_CAL_REGISTRATION','FLAG_RELEASE','NOTES'),
            data_start=0, comment='#')
tq12 = Table.read('inputs/QCflags_std_V1200_DR3.csv', format='ascii.csv', 
            names=('CALIFAID','Name','FLAG_OBS_SKYMAG','FLAG_OBS_EXT',
                'FLAG_OBS_AM','FLAG_RED_STRAYLIGHT','FLAG_RED_DISP',
                'FLAG_RED_CDISP','FLAG_RED_SKYLINES','FLAG_RED_LIMSB',
                'FLAG_RED_ERRSPEC','FLAG_CAL_WL',
                'FLAG_CAL_IMGQUAL','FLAG_CAL_SPECQUAL','FLAG_CAL_FLATSDSS',
                'FLAG_CAL_REGISTRATION','FLAG_RELEASE','NOTES'),
            data_start=0, comment='#')
ca_qc = join(tq5, tq12, keys='Name', join_type='left', table_names=['V500','V1200'])

# NED positions
nedt = Table.read('inputs/nedpos.txt', format='ascii.csv', delimiter='\t')

# ---------------------------------------------------------------------------------

# Write the LEDA table:
t['ledaRA'].unit = 'hourangle'
t['ledaRA'].description = 'RA J2000 from LEDA <celpos>'
t['ledaDE'].unit = 'deg'
t['ledaDE'].description = 'DEC J2000 from LEDA <celpos>' 
t['ledaA_Bgal'].unit = 'mag'
t['ledaA_Bgal'].description = 'Galactic A_B from LEDA <ag>'
t['ledaType'].description = 'Morphological type from LEDA <t>'
t['ledaD25'].unit = 'arcmin'
t['ledaD25'].description = 'Apparent B diameter from LEDA <logd25> linearized' 
t['ledaAxRatio'].description = 'Maj/min axis ratio from LEDA <logr25> linearized' 
t['ledaPA'].unit = 'deg'
t['ledaPA'].description = 'PA from LEDA <pa>, N to E' 
t['ledaIncl'].unit = 'deg'
t['ledaIncl'].description = 'Morph inclination from LEDA <incl>' 
t['ledaVrad'].unit = 'km / s'
t['ledaVrad'].description = 'cz from radio data from LEDA <vrad>' 
t['ledaVmaxg'].unit = 'km / s'
t['ledaVmaxg'].description = 'HI max v_rot not corr for incl from LEDA <vmaxg>' 
t['ledaVrot'].unit = 'km / s'
t['ledaVrot'].description = 'HI max v_rot corr for incl from LEDA <vrot>' 
t['ledaMorph'].description = 'Hubble type from LEDA <type>' 
t['ledaBar'].description = 'B = bar present from LEDA <bar>' 
t['ledaRing'].description = 'R = ring present from LEDA <ring>' 
t['ledaMultiple'].description = 'M = multiple system from LEDA <multiple>' 
t['ledaBt'].unit = 'mag'
t['ledaBt'].description = 'Apparent B total magnitude from LEDA <bt>' 
t['ledaIt'].unit = 'mag'
t['ledaIt'].description = 'Apparent I total magnitude from LEDA <it>' 
t['ledaMfir'].unit = 'mag'
t['ledaMfir'].description = 'FIR flux as magnitude from LEDA <mfir>' 
t['ledaM21'].unit = 'mag'
t['ledaM21'].description = 'HI line flux as magnitude from LEDA <m21>' 
t['ledaVvir'].unit = 'km / s'
t['ledaVvir'].description = 'Virgo infall corr cz from LEDA <vvir>' 
t['ledaModz'].unit = 'mag'
t['ledaModz'].description = 'Dist modulus from LEDA <modz> based on <vvir>; NGC2880 NGC4211 and UGC05498 substituted with NED scaled to same cosmology' 
t['ledaDistMpc'].unit = 'Mpc'
t['ledaDistMpc'].description = 'Distance in Mpc corresponding to ledaModz' 
if 'leda' in rewrite:
    outcols=['Name']
    for cname in t.colnames:
        if "leda" in cname:
            outcols.append(cname)
    newt = t[outcols]
    newt.write('edge_leda.csv', format='ascii.ecsv', delimiter=',', overwrite=True)


# Write the NED table:
if 'ned' in rewrite:
    sc = SkyCoord(nedt['RA'],nedt['Dec']) #convert to degrees
    nedRA=[]
    nedDE=[]
    for obj in sc:
        nedRA.append(float(obj.to_string(precision=5).split(' ')[0]))
        nedDE.append(float(obj.to_string(precision=5).split(' ')[1]))
    newt=Table()
    newt['Name'] = nedt['Name']
    newt['nedRA'] = nedRA
    newt['nedRA'].unit = 'deg'
    newt['nedRA'].description = 'J2000 RA from NED'
    newt['nedDE'] = nedDE
    newt['nedDE'].unit = 'deg'
    newt['nedDE'].description = 'J2000 DEC from NED'
    newt['nedVopt']=nedt['cz']
    newt['nedVopt'].unit='km / s'
    newt['nedVopt'].description='cz from NED'
    newt.write('edge_ned.csv', format='ascii.ecsv', delimiter=',', overwrite=True)


# Write the WISE table:
t['W1'].unit = 'mag'
t['W1'].description = '3.4 um Vega magnitude from image photometry from table by Bitsakis. NGC0598 and NGC4676A are from w1mpro in allWISE catalog'
t['eW1'].unit = 'mag'
t['eW1'].description = 'Error in W1 in Vega magnitudes from table by Bitsakis. NGC0598 and NGC4676A are from w1mpro+2.5*log(w1snr) in allWISE catalog'
t['W2'].unit = 'mag'
t['W2'].description = '4.6 um Vega magnitude from Bitsakis'
t['eW2'].unit = 'mag'
t['eW2'].description = 'Error in W2 in Vega magnitudes from Bitsakis'
t['W3'].unit = 'mag'
t['W3'].description = '12 um Vega magnitude from Bitsakis'
t['eW3'].unit = 'mag'
t['eW3'].description = 'Error in W3 in Vega magnitudes from Bitsakis'
t['W4'].unit = 'mag'
t['W4'].description = '22 um Vega magnitude from Bitsakis'
t['eW4'].unit = 'mag'
t['eW4'].description = 'Error in W4 in Vega magnitudes from Bitsakis'
t['W1lum'].unit = 'dex(erg / s)'
t['W1lum'].description = 'Luminosity in W1 from caDistMpc and W1 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W2lum'].unit = 'dex(erg / s)'
t['W2lum'].description = 'Luminosity in W2 from caDistMpc and W2 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W3lum'].unit = 'dex(erg / s)'
t['W3lum'].description = 'Luminosity in W3 from caDistMpc and W3 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W4lum'].unit = 'dex(erg / s)'
t['W4lum'].description = 'Luminosity in W4 from caDistMpc and W4 magnitude using the zero point, frequency, and bandwidth from Jarrett+11'
t['W4SFR'].unit = 'solMass / yr'
t['W4SFR'].description = 'SFR from W4lum using Catalan-Torrecilla+14 calibration'
if 'wise' in rewrite:
    outcols=['Name']
    clist=['W1', 'eW1', 'W2', 'eW2', 'W3', 'eW3', 'W4', 'eW4', 'W1lum', 'W2lum', 'W3lum', 'W4lum', 'W4SFR']
    for cname in clist:
        outcols.append(cname)
    for cname in ['W1lum', 'W2lum', 'W3lum', 'W4lum']:
        logged = np.around(np.log10(t[cname]),decimals=4)
        t[cname][:] = logged
    newt = t[outcols]
    newt.write('edge_wise.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the NSA table:
t['nsaZdist'].description = 'Redshift from NASA-Sloan Atlas using pecular velocity model of Willick+ 97'
t['nsaAu'].unit = 'mag'
t['nsaAu'].description = 'Galactic extinction in u from Schlegel+ 97'
t['nsaAg'].unit = 'mag'
t['nsaAg'].description = 'Galactic extinction in g from Schlegel+ 97'
t['nsaAr'].unit = 'mag'
t['nsaAr'].description = 'Galactic extinction in r from Schlegel+ 97'
t['nsaAi'].unit = 'mag'
t['nsaAi'].description = 'Galactic extinction in i from Schlegel+ 97'
t['nsaAz'].unit = 'mag'
t['nsaAz'].description = 'Galactic extinction in z from Schlegel+ 97'
if 'nsa' in rewrite:
    outcols=['Name']
    for cname in t.colnames:
        if cname.startswith("nsa"):
            outcols.append(cname)
    newt = t[outcols]
    newt.write('edge_nsa.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the CALIFA table (edge_califa.csv):
t['caMass'].unit = 'dex(solMass)'
t['caMass'].description = 'Stellar mass from col 149 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caeMass'].unit = 'dex(solMass)'
t['caeMass'].description = 'Error in stellar mass from col 150 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caSFR'].unit = 'dex(solMass / yr)'
t['caSFR'].description = 'SFR from col 151 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caeSFR'].unit = 'dex(solMass / yr)'
t['caeSFR'].description = 'Error in SFR from col 152 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caOH'].unit = 'dex'
t['caOH'].description = 'Oxygen abundance as 12+log(O/H) from col 153 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caeOH'].unit = 'dex'
t['caeOH'].description = 'Error in oxygen abundance from col 154 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caAvgas'].unit = 'mag'
t['caAvgas'].description = 'Nebular extinction as Av from col 155 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caeAvgas'].unit = 'mag'
t['caeAvgas'].description = 'Error in nebular extinction from col 156 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caAvstars'].unit = 'mag'
t['caAvstars'].description = 'Stellar extinction as Av from col 157 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['caeAvstars'].unit = 'mag'
t['caeAvstars'].description = 'Error in stellar extinction from col 158 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
t['Su'].unit = 'mag'
t['Su'].description = 'SDSS u magnitude from col 4 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['Sg'].unit = 'mag'
t['Sg'].description = 'SDSS g magnitude from col 8 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['Sr'].unit = 'mag'
t['Sr'].description = 'SDSS r magnitude from col 12 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['Si'].unit = 'mag'
t['Si'].description = 'SDSS i magnitude from col 16 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caB'].unit = 'mag'
t['caB'].description = 'B magnitude from col 20 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caV'].unit = 'mag'
t['caV'].description = 'V magnitude from col 24 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caR'].unit = 'mag'
t['caR'].description = 'R magnitude from col 28 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caRe'].unit = 'arcsec'
t['caRe'].description = 'Equivalent radius Re from col 34 in get_mag_cubes_v2.2.csv'
t['caeRe'].unit = 'arcsec'
t['caeRe'].description = 'Error in equivalent radius from col 35 in get_mag_cubes_v2.2.csv'
t['caEllipticity'].description = 'Ellipticity sqrt(1-b^2/a^2) from col 38 in get_mag_cubes_v2.2.csv'
t['caPA'].unit = 'deg'
t['caPA'].description = 'PA from col 39 in get_mag_cubes_v2.2.csv'
t['caR50'].unit = 'arcsec'
t['caR50'].description = 'R50 from col 40 in get_mag_cubes_v2.2.csv'
t['caeR50'].unit = 'arcsec'
t['caeR50'].description = 'Error in R50 from col 41 in get_mag_cubes_v2.2.csv'
t['caR90'].unit = 'arcsec'
t['caR90'].description = 'R90 from col 42 in get_mag_cubes_v2.2.csv'
t['caeR90'].unit = 'arcsec'
t['caeR90'].description = 'Error in R90 from col 43 in get_mag_cubes_v2.2.csv'
t['caOH_O3N2'].unit = 'dex'
t['caOH_O3N2'].description = 'O3N2-based metallicity from <OH_O3N2> col 4 in get_proc_elines_CALIFA.csv'
t['caZgas'].description = 'Redshift for gas lines, from <z_gas> col 14 in get_proc_elines_CALIFA.csv'
t['caZstars'].description = 'Redshift for stars, from <z_stars> col 15 in get_proc_elines_CALIFA.csv'
t['caAge'].unit = 'dex(Gyr)'
t['caAge'].description = 'Mean stellar age from <log_age_mean_LW> col 37 in get_proc_elines_CALIFA.csv'
t['caeAge'].unit = 'dex(Gyr)'
t['caeAge'].description = 'Error in mean stellar age from <s_log_age_mean_LW> col 38 in get_proc_elines_CALIFA.csv'
t['caFHa'].unit = 'dex(1e-16 erg / (cm2 s))'
t['caFHa'].description = 'Log of Halpha flux, from <log_F_Ha> column 145 in get_proc_elines_CALIFA.csv'
t['caFHacorr'].unit = 'dex(1e-16 erg / (cm2 s))'
t['caFHacorr'].description = 'Log of Halpha flux, extinction corrected, from <log_F_Ha_cor> column 146 in get_proc_elines_CALIFA.csv'
t['caLHacorr'].unit = 'dex(erg / s)'
t['caLHacorr'].description = 'Log of Halpha luminosity, extinction corrected, from <log_L_Ha_cor> column 147 in get_proc_elines_CALIFA.csv'
t['caMstars'].unit = 'dex(solMass)'
t['caMstars'].description = 'Log of stellar mass, from column 73 log_Mass in get_proc_elines_CALIFA.csv'
t['caDistMpc'].unit = 'Mpc'
t['caDistMpc'].description = 'Luminosity distance in Mpc computed from caZgas assuming Ho=70, Om=0.27, Ol=0.73'
if 'califa' in rewrite:
    outcols=['Name']
    for cname in t.colnames:
        if cname.startswith("ca") or cname.startswith("S"):
            outcols.append(cname)
    newt = t[outcols]
    joint = join(t,ca_qc,keys='Name',join_type='left')
    newt.add_columns([joint['FLAG_CAL_WL_V500'], 
                      joint['FLAG_CAL_REGISTRATION_V500'], 
                      joint['FLAG_CAL_IMGQUAL_V500']], 
                      names=['caFlgWav5', 'caFlgReg5', 'caFlgImg5'])
    newt.add_columns([joint['FLAG_CAL_WL_V1200'], 
                      joint['FLAG_CAL_REGISTRATION_V1200'], 
                      joint['FLAG_CAL_IMGQUAL_V1200']], 
                      names=['caFlgWav12', 'caFlgReg12', 'caFlgImg12'])
    newt['caFlgWav5'].description  = 'Flag (-1/0/1/2=NA/good/minor/bad) for wavelength calibration V500'
    newt['caFlgWav12'].description = 'Flag (-1/0/1/2=NA/good/minor/bad) for wavelength calibration V1200'
    newt['caFlgReg5'].description  = 'Flag (-1/0/1/2=NA/good/minor/bad) for 2D registration rel to SDSS V500'
    newt['caFlgReg12'].description = 'Flag (-1/0/1/2=NA/good/minor/bad) for 2D registration rel to SDSS V1200'
    newt['caFlgImg5'].description  = 'Flag (-1/0/1/2=NA/good/minor/bad) for reconstructed image quality V500'
    newt['caFlgImg12'].description = 'Flag (-1/0/1/2=NA/good/minor/bad) for reconstructed image quality V1200'
    newt.write('edge_califa.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the radial distributions table (edge_rdist.csv):
t['coScaleMol'].unit = 'kpc'
t['coScaleMol'].description = 'Exponential scale length for CO disk derived by filling in undetected values in annulii with 1-sigma'
t['coeScaleMol'].unit = 'kpc'
t['coeScaleMol'].description = 'Statistical error from fit to exponential scale length'
t['coScaleMolHi'].unit = 'kpc'
t['coScaleMolHi'].description = 'Upper limit to exponential scale length by filling in annuli with 2-sigma values'
t['coScaleMolLo'].unit = 'kpc'
t['coScaleMolLo'].description = 'Lower limit to exponential scale length by filling in annulii with zeros'
t['coNormMol'].unit = 'solMass / pc2'
t['coNormMol'].description = 'Normalization of CO exponential disk profile, i.e. density at R=0'
t['coeNormMol'].unit = 'solMass / pc2'
t['coeNormMol'].description = 'Error in normalization of CO exponential disk profile'
t['coScaleSt'].unit = 'kpc'
t['coScaleSt'].description = 'Exponential scale length for the mass of the stellar disk'
t['coeScaleSt'].unit = 'kpc'
t['coeScaleSt'].description = 'Formal error in stellar scale length fit' 
t['coNormSt'].unit = 'solMass / pc2'
t['coNormSt'].description = 'Normalization of stellar exponential disk profile, i.e. density at R=0'
t['coeNormSt'].unit = 'solMass / pc2'
t['coeNormSt'].description = 'Formal error in stellar disk normalization'
t['coR50Mol'].unit = 'kpc'
t['coR50Mol'].description = 'Radius enclosing 50% of the molecular mass'
t['coeR50Mol'].unit = 'kpc'
t['coeR50Mol'].description = 'Error in radius enclosing 50% of the molecular mass, including beam size'
t['coR50St'].unit = 'kpc'
t['coR50St'].description = 'Radius enclosing 50% of the stellar mass'
t['coeR50St'].unit = 'kpc'
t['coeR50St'].description = 'Error in radius enclosing 50% of the stellar mass, including beam size'
t['coScaleSFR'].unit = 'kpc'
t['coScaleSFR'].description = 'Exponential scale length for SFR from extinction corrected Ha'
t['coeScaleSFR'].unit = 'kpc'
t['coeScaleSFR'].description = 'Error in exponential scale length for SFR from extinction corrected Ha'
if 'rdist' in rewrite:
    outcols=['Name']
    clist = ['coScaleMol', 'coeScaleMol', 'coScaleMolHi', 'coScaleMolLo', 'coNormMol', 'coeNormMol', 'coScaleSt', 'coeScaleSt', 'coNormSt', 'coeNormSt', 'coR50Mol', 'coeR50Mol', 'coR50St', 'coeR50St', 'coScaleSFR', 'coeScaleSFR']
    newt = t[outcols+clist]
    for cname in clist:
        newt[cname].name = t[cname].name.replace('co','rd')
    newt.write('edge_rdist.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# Write the RINGFIT kinematics table (edge_rfpars.csv) - now moved to derived
# t['coVsys'].unit = 'km / s'
# t['coVsys'].description = 'Systemic galaxy velocity (relativistic convention), derived from CO rotation curve or taken from Ha or LEDA if no curve'
# t['coVsysFlag'].description = 'Origin of coVsys, 0=CO kinematics, 1=Ha kinematics, 3=LEDA'
# t['coVrotmax'].unit = 'km / s'
# t['coVrotmax'].description = 'Maximum rotational velocity of CO, determined using the Lelli+16 algorithm'
# t['coeVrotmax'].unit = 'km / s'
# t['coeVrotmax'].description = 'Error in Vrotmax'
# t['coRflat'].unit = 'arcsec'
# t['coRflat'].description = 'Radius at which the rotation curve flattens, according to Lelli+16 algorithm'
# t['coPA'].unit = 'deg'
# t['coPA'].description = 'Galaxy PA determined from CO/Ha rotation curve by minimizing the radial component, or if not available from outer isophotes or LEDA'
# t['coInc'].unit = 'deg'
# t['coInc'].description = 'Galaxy Inc determined from best fit CO/Ha rotation curve, or if not available from outer isophotes or LEDA'
# t['coPAFlag'].description = 'Origin of coPA, 0=CO kinematics, 1=Ha kinematics, 2=outer optical isophotes (GvdV), 3=LEDA'
# t['coIncFlag'].description = 'Origin of coInc, 0=CO kinematics, 1=Ha kinematics, 2=outer optical isophotes (GvdV), 3=LEDA'
# t['coKinXoff'].unit = 'arcsec'
# t['coKinXoff'].description = 'Offset of CO kinematic center from LEDA RA, RA_kc=RA-xoff'
# t['coKinYoff'].unit = 'arcsec'
# t['coKinYoff'].description = 'Offset of CO kinematic center from LEDA DE, DE_kc=DE-yoff'
# if 'rfpars' in rewrite:
#     outcols=['Name']
#     clist=['coVsys','coVsysFlag','coVrotmax','coeVrotmax','coRflat','coPA','coInc','coPAFlag','coIncFlag','coKinXoff','coKinYoff']
#     newt = t[outcols+clist]
#     for cname in clist:
#         newt[cname].name = t[cname].name.replace('co','rf')
#     newt.write('edge_rfpars.csv', format='ascii.ecsv', delimiter=',', overwrite=True)
