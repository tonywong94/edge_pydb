#!/usr/bin/env python
# coding: utf-8

from astropy.io import fits
from astropy.table import Table, join
from astropy import units as u
import numpy as np
from datetime import datetime

# Properties determined outside pyPipe3D in galaxies_properties.fits

gtab = Table.read('galaxies_properties.fits.gz')
gtab.keep_columns(['ID','cubename','RA','DEC','z','type','r-abs','g-r','Re',
                   'QC_flag','QC_reason','multiple'])
gtab.rename_column('cubename','Name')
gtab.rename_column('r-abs','r_band_abs_mag')

gtab['ID'].description = 'CALIFA ID'
gtab['Name'].description = 'CALIFA Name'
gtab['RA'].unit = 'deg'
gtab['DEC'].unit = 'deg'
gtab['r_band_abs_mag'].unit = 'mag'
gtab['g-r'].unit = 'mag'
gtab['Re'].unit = 'arcsec'
gtab['QC_flag'].description = '0:OK;1:BAD;2:WARNING'
gtab['QC_reason'].description = '0:OK,1:partially covered/too large,2:bad quality,3:not good quality,4:too small,5:repeated,6:Crowded with field stars,7:Bright single field star,8:Scale factor is clearly wrong,9:Wrong convolution,10:Integrated and resolved masses do not match,11:pyPipe3D fitting issues'
gtab['multiple'].description = '0:single,1:multiple'
gtab.meta.clear()

# Properties determined from pyPipe3D in eCALIFA.pyPipe3D.fits

ptab = Table.read('eCALIFA.pyPipe3D.fits.gz')
ptab.sort('ID')
ptab.meta.clear()

names1 = ['cubename', 'DL', 'DA', 'PA', 'ellip', 'Re_kpc', 'log_Mass', 'e_log_Mass', 
          'log_SFR_Ha', 'e_log_SFR_Ha', 'log_SFR_ssp', 'e_log_SFR_ssp',
          'log_NII_Ha_cen', 'e_log_NII_Ha_cen', 'log_OIII_Hb_cen',
          'e_log_OIII_Hb_cen', 'log_SII_Ha_cen', 'e_log_SII_Ha_cen',
          'log_OII_Hb_cen', 'e_log_OII_Hb_cen', 'EW_Ha_cen',
          'e_EW_Ha_cen','Age_LW_Re_fit', 'e_Age_LW_Re_fit','Age_MW_Re_fit',
          'e_Age_MW_Re_fit', 'vel_sigma_Re', 'log_SFR_SF', 'log_SFR_D_C']

names2 = ['OH_O3N2_cen', 'e_OH_O3N2_cen', 'OH_N2_cen', 'e_OH_N2_cen', 'Ha_Hb_cen',
          'e_Ha_Hb_cen', 'log_NII_Ha_Re', 'e_log_NII_Ha_Re',
          'log_OIII_Hb_Re', 'e_log_OIII_Hb_Re', 'log_SII_Ha_Re',
          'e_log_SII_Ha_Re', 'log_OII_Hb_Re', 'e_log_OII_Hb_Re',
          'log_OI_Ha_Re', 'e_log_OI_Ha_Re', 'EW_Ha_Re', 'e_EW_Ha_Re',
          'Ha_Hb_Re', 'e_Ha_Hb_Re', 'log_NII_Ha_ALL', 'e_log_NII_Ha_ALL',
          'log_OIII_Hb_ALL', 'e_log_OIII_Hb_ALL', 'log_SII_Ha_ALL',
          'e_log_SII_Ha_ALL', 'log_OII_Hb_ALL', 'e_log_OII_Hb_ALL',
          'log_OI_Ha_ALL', 'e_log_OI_Ha_ALL', 'EW_Ha_ALL', 'e_EW_Ha_ALL',
          'Ha_Hb_ALL', 'Sigma_Mass_cen', 'e_Sigma_Mass_cen',
          'Sigma_Mass_Re', 'e_Sigma_Mass_Re', 'Sigma_Mass_ALL',
          'e_Sigma_Mass_ALL', 'vel_disp_Ha_cen', 'vel_disp_ssp_cen',
          'vel_disp_Ha_1Re', 'vel_disp_ssp_1Re', 'Av_gas_Re',
          'e_Av_gas_Re', 'Av_ssp_Re', 'e_Av_ssp_Re']

names3 = ['flux_Halpha6562.85_Re_fit', 'e_flux_Halpha6562.85_Re_fit',
          'flux_Halpha6562.85_alpha_fit',
          'e_flux_Halpha6562.85_alpha_fit','OH_Mar13_N2_Re_fit',
          'e_OH_Mar13_N2_Re_fit', 'OH_Mar13_N2_alpha_fit',
          'e_OH_Mar13_N2_alpha_fit', 'OH_Mar13_O3N2_Re_fit',
          'e_OH_Mar13_O3N2_Re_fit', 'OH_Mar13_O3N2_alpha_fit',
          'e_OH_Mar13_O3N2_alpha_fit', 'OH_Pet04_O3N2_Re_fit',
          'e_OH_Pet04_O3N2_Re_fit', 'OH_Pet04_O3N2_alpha_fit',
          'e_OH_Pet04_O3N2_alpha_fit', 'OH_Pil16_S_Re_fit',
          'e_OH_Pil16_S_Re_fit', 'OH_Pil16_S_alpha_fit',
          'e_OH_Pil16_S_alpha_fit']

allnames = names1+names2+names3
ptab.keep_columns(allnames)
ptab.rename_column('cubename','Name')

# Fix distance columns (not actually Mpc)
ptab['DL'].unit = 'cm'
ptab['DL'].convert_unit_to(u.Mpc)
ptab['DA'] *= 206.265
ptab['DA'].unit = 'Mpc'

ptab['PA'].unit = 'deg'
ptab['Re_kpc'].unit = 'kpc'
ptab['log_Mass'].unit = 'dex(solMass)'
ptab['e_log_Mass'].unit = 'dex(solMass)'
ptab['log_SFR_Ha'].unit = 'dex(solMass/yr)'
ptab['e_log_SFR_Ha'].unit = 'dex(solMass/yr)'
ptab['log_SFR_ssp'].unit = 'dex(solMass/yr)'
ptab['e_log_SFR_ssp'].unit = 'dex(solMass/yr)'
ptab['EW_Ha_cen'].unit = 'Angstrom'
ptab['e_EW_Ha_cen'].unit = 'Angstrom'
ptab['Age_LW_Re_fit'].unit = 'dex(yr)'
ptab['e_Age_LW_Re_fit'].unit = 'dex(yr)'
ptab['Age_MW_Re_fit'].unit = 'dex(yr)'
ptab['e_Age_MW_Re_fit'].unit = 'dex(yr)'
ptab['log_SFR_SF'].unit = 'dex(solMass/yr)'
ptab['log_SFR_D_C'].unit = 'dex(solMass/yr)'

for name in (names2+names3):
    if name.startswith(('OH','e_OH')):
        ptab[name].unit = 'dex'
    elif name.startswith(('flux','e_flux')):
        ptab[name].unit = '1e-16 erg/(cm2 s)'
    elif name.startswith(('Av','e_Av')):
        ptab[name].unit = 'mag'
    elif name.startswith(('EW','e_EW')):
        ptab[name].unit = 'Angstrom'
    elif name.startswith('vel_disp'):
        ptab[name].unit = 'km/s'
    elif name.startswith(('Sigma_Mass','e_Sigma_Mass')):
        ptab[name].unit = 'dex(solMass/pc2)'

# Fix incorrect column dtypes and fill masked values with NaN
for coln in ptab.colnames[1:]:
    if ptab[coln].dtype.name.startswith('bytes'):
        if hasattr(ptab[coln], 'mask'):
            ptab[coln] = ptab[coln].filled(fill_value=np.nan).astype(float)
        else:
            ptab[coln] = ptab[coln].astype(float)
    else:
        if hasattr(ptab[coln], 'mask'):
            ptab[coln] = ptab[coln].filled(fill_value=np.nan)
    ptab[coln].format = '{:.4f}'

# Join the tables and output
ctab = join(gtab, ptab, keys=['Name'])
ctab.sort('ID')
ctab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
ctab.meta['comments'] = ('Galaxy properties determined from eCALIFA')
print(ctab.meta)
ctab.write('ecalifa_global.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

