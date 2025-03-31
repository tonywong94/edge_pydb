#!/usr/bin/env python
# coding: utf-8

from astropy.io import fits
from astropy.table import Table, join
from astropy import units as u
from astropy import constants as const
import numpy as np
from datetime import datetime

# DR17 Pipe3D analysis from Sanchez+22.
# Download from https://data.sdss.org/sas/dr17/manga/spectro/pipe3d/v3_1_1/3.1.1/
# For data model see https://data.sdss.org/datamodel/files/MANGA_PIPE3D/MANGADRP_VER/PIPE3D_VER/SDSS17Pipe3D.html
p3dtab = Table.read('SDSS17Pipe3D_v3_1_1.fits', hdu=1)

# Final summary FITS binary table for the MaNGA Data Reduction Pipeline
# Download from https://www.sdss4.org/dr17/manga/manga-data/catalogs/
# For data model see https://data.sdss.org/datamodel/files/MANGA_SPECTRO_REDUX/DRPVER/drpall.html
drptab = Table.read('drpall-v3_1_1.fits', hdu=1)
# From documentation for nsa_zdist: Distance estimate using peculiar velocity model of Willick et al. (1997); mulitply by c/Ho for Mpc
drptab['nsa_z_dMpc'] = drptab['nsa_zdist'] * (const.c/(70*u.km/(u.s*u.Mpc))).to(u.Mpc)

# Join the tables
gtab = join(p3dtab, drptab, join_type='left', 
            keys=['mangaid','plateifu','objra','objdec'])

keepcols = ['plateifu', 'mangaid', 'objra', 'objdec', 
            'drp3qual', 'seemed', 'gfwhm', 'rfwhm', 'ifwhm', 'zfwhm', 'QCFLAG',
            'nsa_redshift', 'nsa_iauname', 'nsa_sersic_ba', 'nsa_sersic_phi',
            'nsa_inclination', 'nsa_z_dMpc', 'r_band_abs_mag', 'g-r', 'DL', 'DA', 
            'Re_arc', 'Re_kpc', 'PA', 'ellip', 
            'log_Mass', 'e_log_Mass', 'log_SFR_Ha', 'e_log_SFR_Ha',
            'log_SFR_ssp', 'log_NII_Ha_cen', 'e_log_NII_Ha_cen',
            'log_OIII_Hb_cen', 'e_log_OIII_Hb_cen', 'log_SII_Ha_cen',
            'e_log_SII_Ha_cen', 'log_OII_Hb_cen', 'e_log_OII_Hb_cen',
            'EW_Ha_cen', 'e_EW_Ha_cen', 'Age_LW_Re_fit', 'e_Age_LW_Re_fit',
            'Age_MW_Re_fit', 'e_Age_MW_Re_fit',
            'vel_sigma_Re', 'log_SFR_SF', 'log_SFR_D_C', 
            'OH_O3N2_cen', 'e_OH_O3N2_cen', 'OH_N2_cen', 'e_OH_N2_cen',
            'Ha_Hb_cen', 'e_Ha_Hb_cen', 'log_NII_Ha_Re', 'e_log_NII_Ha_Re',
            'log_OIII_Hb_Re', 'e_log_OIII_Hb_Re',
            'log_SII_Ha_Re', 'e_log_SII_Ha_Re', 'log_OII_Hb_Re', 'e_log_OII_Hb_Re',
            'log_OI_Ha_Re', 'e_log_OI_Ha_Re', 'EW_Ha_Re', 'e_EW_Ha_Re',
            'Ha_Hb_Re', 'e_Ha_Hb_Re', 'log_NII_Ha_ALL', 'e_log_NII_Ha_ALL',
            'log_OIII_Hb_ALL', 'e_log_OIII_Hb_ALL',
            'log_SII_Ha_ALL', 'e_log_SII_Ha_ALL', 'log_OII_Hb_ALL', 'e_log_OII_Hb_ALL',
            'log_OI_Ha_ALL', 'e_log_OI_Ha_ALL', 'EW_Ha_ALL', 'e_EW_Ha_ALL',
            'Ha_Hb_ALL', 'Sigma_Mass_cen', 'e_Sigma_Mass_cen',
            'Sigma_Mass_Re', 'e_Sigma_Mass_Re', 'Sigma_Mass_ALL', 'e_Sigma_Mass_ALL',
            'vel_disp_Ha_cen', 'vel_disp_ssp_cen',
            'vel_disp_Ha_1Re', 'vel_disp_ssp_1Re', 
            'Av_gas_Re', 'e_Av_gas_Re',
            'Av_ssp_Re', 'e_Av_ssp_Re',
            'flux_Halpha6562.85_Re_fit', 'e_flux_Halpha6562.85_Re_fit',
            'flux_Halpha6562.85_alpha_fit','e_flux_Halpha6562.85_alpha_fit',
            'OH_Mar13_N2_Re_fit', 'e_OH_Mar13_N2_Re_fit',
            'OH_Mar13_N2_alpha_fit', 'e_OH_Mar13_N2_alpha_fit',
            'OH_Mar13_O3N2_Re_fit', 'e_OH_Mar13_O3N2_Re_fit',
            'OH_Mar13_O3N2_alpha_fit', 'e_OH_Mar13_O3N2_alpha_fit',
            'OH_Pet04_O3N2_Re_fit', 'e_OH_Pet04_O3N2_Re_fit',
            'OH_Pet04_O3N2_alpha_fit', 'e_OH_Pet04_O3N2_alpha_fit',
            'OH_Pil16_S_Re_fit', 'e_OH_Pil16_S_Re_fit',
            'OH_Pil16_S_alpha_fit', 'e_OH_Pil16_S_alpha_fit']

gtab = gtab[keepcols]
gtab.meta.clear()

gtab.rename_column('plateifu','Name')
gtab['Name'].description = 'MaNGA Name'
gtab.rename_column('QCFLAG','QC_flag')
gtab['QC_flag'].description= 'QC flag 0=good 2=bad >2 warning'

# Flag negative values for Re
gtab['Re_arc'][np.where(gtab['Re_arc'] < 0)] = np.nan
gtab['Re_kpc'][np.where(gtab['Re_kpc'] < 0)] = np.nan

# Set the units
gtab['objra'].unit = 'deg'
gtab['objdec'].unit = 'deg'
gtab['seemed'].unit = 'arcsec'
gtab['gfwhm'].unit = 'arcsec'
gtab['rfwhm'].unit = 'arcsec'
gtab['ifwhm'].unit = 'arcsec'
gtab['zfwhm'].unit = 'arcsec'
gtab['nsa_inclination'].unit = 'deg'
gtab['r_band_abs_mag'].unit = 'mag'
gtab['g-r'].unit = 'mag'
gtab['Re_arc'].unit = 'arcsec'
gtab['DL'].unit = 'cm'
gtab['DL'].convert_unit_to(u.Mpc)
gtab['DA'] *= 206.265
gtab['DA'].unit = 'Mpc'
gtab['PA'].unit = 'deg'
gtab['Re_kpc'].unit = 'kpc'
gtab['log_Mass'].unit = 'dex(solMass)'
gtab['e_log_Mass'].unit = 'dex(solMass)'
gtab['log_SFR_Ha'].unit = 'dex(solMass/yr)'
gtab['e_log_SFR_Ha'].unit = 'dex(solMass/yr)'
gtab['log_SFR_ssp'].unit = 'dex(solMass/yr)'
gtab['EW_Ha_cen'].unit = 'Angstrom'
gtab['e_EW_Ha_cen'].unit = 'Angstrom'
gtab['Age_LW_Re_fit'].unit = 'dex(yr)'
gtab['e_Age_LW_Re_fit'].unit = 'dex(yr)'
gtab['Age_MW_Re_fit'].unit = 'dex(yr)'
gtab['e_Age_MW_Re_fit'].unit = 'dex(yr)'
gtab['log_SFR_SF'].unit = 'dex(solMass/yr)'
gtab['log_SFR_D_C'].unit = 'dex(solMass/yr)'

for name in (gtab.colnames):
    if name.startswith(('OH','e_OH')):
        gtab[name].unit = 'dex'
    elif name.startswith(('flux','e_flux')):
        gtab[name].unit = '1e-16 erg/(cm2 s)'
    elif name.startswith(('Av','e_Av')):
        gtab[name].unit = 'mag'
    elif name.startswith(('EW','e_EW')):
        gtab[name].unit = 'Angstrom'
    elif name.startswith('vel_disp'):
        gtab[name].unit = 'km/s'
    elif name.startswith(('Sigma_Mass','e_Sigma_Mass')):
        gtab[name].unit = 'dex(solMass/pc2)'

# Fill masked values with NaN
for coln in gtab.colnames[4:]:
    if hasattr(gtab[coln], 'mask'):
        gtab[coln] = gtab[coln].filled(fill_value=np.nan)
    if gtab[coln].dtype.name.startswith('float'):
        gtab.replace_column(coln, gtab[coln].astype(np.float32))

# Write the table
gtab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
gtab.meta['comments'] = ('Galaxy properties determined from Pipe3D for MaNGA')
print(gtab.meta)
gtab.write('manga_global.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

