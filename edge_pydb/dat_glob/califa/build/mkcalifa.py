#!/usr/bin/env python

from astropy.table import Table, join
import numpy as np
from datetime import datetime
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
#import numpy.ma as ma

# Write the CALIFA table (edge_califa.csv):
# This combines information from three separate CSV tables.
# Info from bitsakis and califa_inter subdirectories not yet folded in.

fixheader = False

if fixheader:
    # 1064 galaxies
    with open('get_mag_cubes_v2.2.csv', 'r') as infile:
        with open('get_mag_cubes_v2.2.fixheader.csv', 'w') as outfile:
            header_names = []
            for line in infile:
                if line[0] == '#':
                    # remove the number and parentheses from the label
                    header_names.append(line.split(') ')[-1].strip())
                    preamble = True
                else:
                    if preamble:
                        outfile.write(",".join(header_names)+'\n')
                        preamble = False
                    outfile.write(line)

    # 1096 galaxies; this already has the header row
    with open('get_proc_elines_CALIFA.csv', 'r') as infile:
        with open('get_proc_elines_CALIFA.fixheader.csv', 'w') as outfile:
            for line in infile:
                if line.startswith('# HEADER '):
                    hdrs = line.replace("# HEADER |","").split()
                    print('There are {} columns in the old header.'.format(len(hdrs)))
                    # --- Problem 1: Some column names are repeated
                    print('Duplicated column names:')
                    print(set([x for x in hdrs if hdrs.count(x) > 1]))
                    uniquehdr = []
                    for i, v in enumerate(hdrs):
                        totalcount = hdrs.count(v)
                        count = hdrs[:i].count(v)
                        if totalcount > 1:
                            print('Renaming column {} to {}'.format(v,v+'_'+str(count+1)))
                            uniquehdr.append(v+'_'+str(count+1))
                        else:
                            uniquehdr.append(v)
                    print('There are {} columns in the new header.'.format(len(uniquehdr)))
                    # --- Problem 2: The OH_S_ALL_error column is corrupted because the
                    # comma separating it from the next column was omitted.  Remove it from
                    # the header.  Note this results in 'NH_Re_fit_R' column being corrupted.
                    uniquehdr.remove('OH_S_ALL_error')
                    outfile.write(",".join(uniquehdr)+'\n')
                elif line[0] != '#':
                    # --- Problem 3: Some entries are 'BAD' instead of 'nan', cant read as float
                    outfile.write(line.replace('BAD', 'nan'))

    # 657 galaxies
    with open('Pipe3D_NSA_CALIFA-DR3_candidates.csv', 'r') as infile:
        with open('Pipe3D_NSA_CALIFA-DR3_candidates.fixheader.csv', 'w') as outfile:
            header_names = []
            for line in infile:
                if line[0] == '#':
                    header_names.append((line.split(' ',2)[2]).strip())
                    preamble = True
                else:
                    if preamble:
                        outfile.write(",".join(header_names)+'\n')
                        preamble = False
                    # Some rows have exactly 24 missing columns after column 113 'O2FLUX'
                    if line.count(",") == 133:
                        parts= line.split(',', 114)
                        idx = len(line)-len(parts[-1])-1
                        line=line[:idx]+',,,,,,,,,,,,,,,,,,,,,,,,'+line[idx:]
                    outfile.write(line)


# Open csv tables                
pipe3d_candidates = Table.read('Pipe3D_NSA_CALIFA-DR3_candidates.fixheader.csv',
                               format='ascii.csv')
get_mag_cubes     = Table.read('get_mag_cubes_v2.2.fixheader.csv',
                               format='ascii.csv')
get_proc_elines   = Table.read('get_proc_elines_CALIFA.fixheader.csv',
                               format='ascii.csv')

# Quality control flags
tq5 = Table.read('QCflags_std_V500_DR3.csv', format='ascii.csv', 
            names=('CALIFAID','Name','FLAG_OBS_SKYMAG','FLAG_OBS_EXT',
                'FLAG_OBS_AM','FLAG_RED_STRAYLIGHT','FLAG_RED_DISP',
                'FLAG_RED_CDISP','FLAG_RED_SKYLINES','FLAG_RED_LIMSB',
                'FLAG_RED_ERRSPEC','FLAG_CAL_SPECPHOTO','FLAG_CAL_WL',
                'FLAG_CAL_IMGQUAL','FLAG_CAL_SPECQUAL','FLAG_CAL_FLATSDSS',
                'FLAG_CAL_REGISTRATION','FLAG_RELEASE','NOTES'),
            data_start=0, comment='#')
tq12 = Table.read('QCflags_std_V1200_DR3.csv', format='ascii.csv', 
            names=('CALIFAID','Name','FLAG_OBS_SKYMAG','FLAG_OBS_EXT',
                'FLAG_OBS_AM','FLAG_RED_STRAYLIGHT','FLAG_RED_DISP',
                'FLAG_RED_CDISP','FLAG_RED_SKYLINES','FLAG_RED_LIMSB',
                'FLAG_RED_ERRSPEC','FLAG_CAL_WL',
                'FLAG_CAL_IMGQUAL','FLAG_CAL_SPECQUAL','FLAG_CAL_FLATSDSS',
                'FLAG_CAL_REGISTRATION','FLAG_RELEASE','NOTES'),
            data_start=0, comment='#')
ca_qc = join(tq5, tq12, keys='Name', join_type='left', table_names=['V500','V1200'])

# Create the V500 sample
t = Table()
t['ID']   = tq5['CALIFAID'] #[tq5['FLAG_RELEASE']==1]
t['Name'] = tq5['Name'] #[tq5['FLAG_RELEASE']==1]
t['ID'].description = 'CALIFA ID'
t['Name'].description = 'CALIFA Name'
print(t)

# qc = Table.read("QCflags_std_V500_DR3.csv", comment='#', delimiter=',',
#                 format='ascii.no_header')
# gallist = sorted(tq5['col2'])
# #gallist = ['UGC12494NOTES01']
# print(gallist)

# Set Pipe3D_NSA_CALIFA-DR3_candidates.fixheader.csv as the table t 
# t = Table()
# t['ID'] = pipe3d_candidates['CALIFAID']
# t['Name'] = pipe3d_candidates['CALIFANAME']
# t['caMass'] = pipe3d_candidates['log(Mass) Pipe3D']
# t['caeMass'] = pipe3d_candidates['error log(Mass) Pipe3D']
# t['caSFR'] = pipe3d_candidates['log(SFR) Pipe3D']
# t['caeSFR'] = pipe3d_candidates['error log(SFR) Pipe3D']
# t['caOH'] = pipe3d_candidates['12+log(O/H) Pipe3D']
# t['caeOH'] = pipe3d_candidates['error 12+log(O/H) Pipe3D']
# t['caAvgas'] = pipe3d_candidates['Av_gas Pipe3D']
# t['caeAvgas'] = pipe3d_candidates['error Av_gas Pipe3D']
# t['caAvstars'] = pipe3d_candidates['Av_stars Pipe3D']
# t['caeAvstars'] = pipe3d_candidates['error Av_stars Pipe3D']
# t = t.filled(np.nan)
# for col in ['caMass','caeMass','caSFR','caeSFR','caOH','caeOH','caAvgas','caeAvgas','caAvstars','caeAvstars']:
#     #t['ledaPA'].mask.nonzero()
#     #t[col][np.where(t[col]==' ')] = np.nan
#     t[col][np.where(t[col]==0)] = np.nan
# t['ID'].description = 'CALIFA ID'
# t['Name'].description = 'CALIFA Name'
# t['caMass'].unit = 'dex(solMass)'
# t['caMass'].description = 'Stellar mass from col 149 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caeMass'].unit = 'dex(solMass)'
# t['caeMass'].description = 'Error in stellar mass from col 150 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caSFR'].unit = 'dex(solMass / yr)'
# t['caSFR'].description = 'SFR from col 151 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caeSFR'].unit = 'dex(solMass / yr)'
# t['caeSFR'].description = 'Error in SFR from col 152 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caOH'].unit = 'dex'
# t['caOH'].description = 'Oxygen abundance as 12+log(O/H) from col 153 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caeOH'].unit = 'dex'
# t['caeOH'].description = 'Error in oxygen abundance from col 154 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caAvgas'].unit = 'mag'
# t['caAvgas'].description = 'Nebular extinction as Av from col 155 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caeAvgas'].unit = 'mag'
# t['caeAvgas'].description = 'Error in nebular extinction from col 156 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caAvstars'].unit = 'mag'
# t['caAvstars'].description = 'Stellar extinction as Av from col 157 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'
# t['caeAvstars'].unit = 'mag'
# t['caeAvstars'].description = 'Error in stellar extinction from col 158 in Pipe3D_NSA_CALIFA-DR3_candidates.csv'

# Columns from get_mag_cubes_v2.2.fixheader.csv
for str_name in ['caSu','caSg','caSr','caSi','caB','caV','caR','caRe','caeRe',
                 'caEllipticity','caPA','caR50','caeR50','caR90','caeR90']:
    t.add_column(np.nan,name=str_name)
for i in range(len(t)):
    if len(np.where(get_mag_cubes['name-obj'] == t['Name'][i])[0]) != 0:
        ind = np.where(get_mag_cubes['name-obj'] == t['Name'][i])[0][0]
        t['caSu'][i] = get_mag_cubes['u band mag'][ind]
        t['caSg'][i] = get_mag_cubes['g band mag'][ind]
        t['caSr'][i] = get_mag_cubes['r band mag'][ind]
        t['caSi'][i] = get_mag_cubes['i band mag'][ind]
        t['caB'][i] = get_mag_cubes['B band mag'][ind]
        t['caV'][i] = get_mag_cubes['V band mag'][ind]
        t['caR'][i] = get_mag_cubes['R band mag'][ind]
        t['caRe'][i] = get_mag_cubes['Re (arcsec)'][ind]
        t['caeRe'][i] = get_mag_cubes['error Re (arcsec)'][ind]
        t['caEllipticity'][i] = get_mag_cubes['ellipticy'][ind] #column name is ellipticy, not ellipticity
        t['caPA'][i] = get_mag_cubes['Pa (deg)'][ind]
        t['caR50'][i] = get_mag_cubes['R50 (arcsec)'][ind]
        t['caeR50'][i] = get_mag_cubes['error R50 (arcsec)'][ind]
        t['caR90'][i] = get_mag_cubes['R90 (arcsec)'][ind]
        t['caeR90'][i] = get_mag_cubes['error R90 (arcsec)'][ind]
t['caSu'].unit = 'mag'
t['caSu'].description = 'SDSS u magnitude from col 4 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caSg'].unit = 'mag'
t['caSg'].description = 'SDSS g magnitude from col 8 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caSr'].unit = 'mag'
t['caSr'].description = 'SDSS r magnitude from col 12 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
t['caSi'].unit = 'mag'
t['caSi'].description = 'SDSS i magnitude from col 16 in get_mag_cubes_v2.2.csv. From CALIFA synthetic photometry corrected for foreground extinction'
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

# Columns from get_proc_elines_CALIFA.fixheader.csv
for str_name in ['cazgas','cazstars','caAge','caeAge','caFHa','caFHacorr','caLHacorr',
                 'caMstars','caeMstars','caSFR','caeSFR','caOH','caeOH','caAvgas',
                 'caeAvgas','caAvstars','caeAvstars','caDistP3d','caDistMpc']:
    t.add_column(np.nan,name=str_name)
cosmo = FlatLambdaCDM(H0=70, Om0=0.27)
for i in range(len(t)):
    if len(np.where(get_proc_elines['name'] == t['Name'][i])[0]) != 0:
        ind = np.where(get_proc_elines['name'] == t['Name'][i])[0][0]
        t['cazgas'][i]     = get_proc_elines['z_gas'][ind]
        t['cazstars'][i]   = get_proc_elines['z_stars'][ind]
        t['caAge'][i]      = get_proc_elines['log_age_mean_LW'][ind]
        t['caeAge'][i]     = get_proc_elines['s_log_age_mean_LW'][ind]
        t['caFHa'][i]      = get_proc_elines['log_F_Ha'][ind]
        t['caFHacorr'][i]  = get_proc_elines['log_F_Ha_cor'][ind]
        t['caLHacorr'][i]  = get_proc_elines['log_L_Ha_cor'][ind]
        t['caMstars'][i]   = get_proc_elines['log_Mass'][ind]
        t['caeMstars'][i]  = get_proc_elines['error_Mass'][ind]
        t['caSFR'][i]      = get_proc_elines['lSFR'][ind]
        t['caeSFR'][i]     = get_proc_elines['e_lSFR'][ind]
        t['caOH'][i]       = get_proc_elines['OH_O3N2'][ind]
        t['caeOH'][i]      = get_proc_elines['e_OH_O3N2'][ind]
        t['caAvgas'][i]    = get_proc_elines['Av_gas_LW_Re'][ind]
        t['caeAvgas'][i]   = get_proc_elines['e_Av_gas_LW_Re'][ind]
        t['caAvstars'][i]  = get_proc_elines['Av_ssp_stats_mean'][ind]
        t['caeAvstars'][i] = get_proc_elines['Av_ssp_stats_stddev'][ind]
        t['caDistP3d'][i]  = get_proc_elines['DL'][ind]
    if t['cazgas'][i] != np.nan:
        t['caDistMpc'][i] = cosmo.luminosity_distance(t['cazgas'][i]).value
t['cazgas'].description = 'Redshift for gas lines from <z_gas> in get_proc_elines_CALIFA.csv'
t['cazstars'].description = 'Redshift for stars from <z_stars> in get_proc_elines_CALIFA.csv'
t['caAge'].unit = 'dex(Gyr)'
t['caAge'].description = 'Mean stellar age from <log_age_mean_LW> in get_proc_elines_CALIFA.csv'
t['caeAge'].unit = 'dex(Gyr)'
t['caeAge'].description = 'Error in mean stellar age from <s_log_age_mean_LW> in get_proc_elines_CALIFA.csv'
t['caFHa'].unit = 'dex(1e-16 erg / (cm2 s))'
t['caFHa'].description = 'Log of Halpha flux from <log_F_Ha> in get_proc_elines_CALIFA.csv'
t['caFHacorr'].unit = 'dex(1e-16 erg / (cm2 s))'
t['caFHacorr'].description = 'Log of extinction corrected Halpha flux from <log_F_Ha_cor> in get_proc_elines_CALIFA.csv'
t['caLHacorr'].unit = 'dex(erg / s)'
t['caLHacorr'].description = 'Log of extinction corrected Halpha luminosity from <log_L_Ha_cor> in get_proc_elines_CALIFA.csv'
t['caMstars'].unit = 'dex(solMass)'
t['caMstars'].description = 'Log of stellar mass from <log_Mass> in get_proc_elines_CALIFA.csv'
t['caeMstars'].unit = 'dex(solMass)'
t['caeMstars'].description = 'Error in log of stellar mass from <error_Mass> in get_proc_elines_CALIFA.csv'
t['caSFR'].unit = 'dex(solMass / yr)'
t['caSFR'].description = 'SFR from <lSFR> in get_proc_elines_CALIFA.csv'
t['caeSFR'].unit = 'dex(solMass / yr)'
t['caeSFR'].description = 'Error in SFR from <e_lSFR> in get_proc_elines_CALIFA.csv'
t['caOH'].unit = 'dex'
t['caOH'].description = 'Oxygen abundance as 12+log(O/H) from <OH_O3N2> in get_proc_elines_CALIFA.csv'
t['caeOH'].unit = 'dex'
t['caeOH'].description = 'Error in oxygen abundance from <e_OH_O3N2> in get_proc_elines_CALIFA.csv'
t['caAvgas'].unit = 'mag'
t['caAvgas'].description = 'Nebular extinction as Av from <Av_gas_LW_Re> in get_proc_elines_CALIFA.csv'
t['caeAvgas'].unit = 'mag'
t['caeAvgas'].description = 'Error in nebular extinction from <e_Av_gas_LW_Re> in get_proc_elines_CALIFA.csv'
t['caAvstars'].unit = 'mag'
t['caAvstars'].description = 'Stellar extinction as Av from <Av_ssp_stats_mean> in get_proc_elines_CALIFA.csv'
t['caeAvstars'].unit = 'mag'
t['caeAvstars'].description = 'Error in stellar extinction <Av_ssp_stats_stddev> in get_proc_elines_CALIFA.csv'
t['caDistP3d'].unit = 'cm'
t['caDistP3d'].convert_unit_to('Mpc')
t['caDistP3d'].description = 'Luminosity distance in Mpc from <DL> in get_proc_elines_CALIFA.csv'
t['caDistMpc'].unit = 'Mpc'
t['caDistMpc'].description = 'Luminosity distance in Mpc computed from cazgas assuming Ho=70, Om=0.27, Ol=0.73'

# CALIFA Quality Tables
for str_name in ['caFlgWav5','caFlgWav12','caFlgReg5','caFlgReg12','caFlgImg5','caFlgImg12']:
    t.add_column(np.nan,name=str_name)
for i in range(len(t)):
    if len(np.where(ca_qc['Name'] == t['Name'][i])[0]) != 0:
        ind = np.where(ca_qc['Name'] == t['Name'][i])[0][0]
        t['caFlgWav5'][i] = ca_qc['FLAG_CAL_WL_V500'][ind]
        t['caFlgReg5'][i] = ca_qc['FLAG_CAL_REGISTRATION_V500'][ind]
        t['caFlgImg5'][i] = ca_qc['FLAG_CAL_IMGQUAL_V500'][ind]
        t['caFlgWav12'][i] = ca_qc['FLAG_CAL_WL_V1200'][ind]
        t['caFlgReg12'][i] = ca_qc['FLAG_CAL_REGISTRATION_V1200'][ind]
        t['caFlgImg12'][i] = ca_qc['FLAG_CAL_IMGQUAL_V1200'][ind]
t['caFlgWav5'].description  = 'Flag (-1/0/1/2=NA/good/minor/bad) for wavelength calibration V500'
t['caFlgWav12'].description = 'Flag (-1/0/1/2=NA/good/minor/bad) for wavelength calibration V1200'
t['caFlgReg5'].description  = 'Flag (-1/0/1/2=NA/good/minor/bad) for 2D registration rel to SDSS V500'
t['caFlgReg12'].description = 'Flag (-1/0/1/2=NA/good/minor/bad) for 2D registration rel to SDSS V1200'
t['caFlgImg5'].description  = 'Flag (-1/0/1/2=NA/good/minor/bad) for reconstructed image quality V500'
t['caFlgImg12'].description = 'Flag (-1/0/1/2=NA/good/minor/bad) for reconstructed image quality V1200'


t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
t.meta['comments'] = ('Galaxy properties determined from CALIFA')
print(t.meta)
t.sort('ID')
t.write('edge_califa.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

