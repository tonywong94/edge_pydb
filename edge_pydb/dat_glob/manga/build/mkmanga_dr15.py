#!/usr/bin/env python
# coding: utf-8

from astropy.io import fits
from astropy.table import Table, Column, join
from astropy import units as u
from astropy import constants as const
import numpy as np
from datetime import datetime

# DR15 Pipe3D analysis from Sanchez+16.
# Download from https://www.sdss4.org/dr15/manga/manga-data/manga-pipe3d-value-added-catalog/
# For data model see https://data.sdss.org/datamodel/files/MANGA_PIPE3D/MANGADRP_VER/PIPE3D_VER/manga.Pipe3D.html
p3dtab = Table.read('/Users/tonywong/Scratch2/almaquest/pipe3d_dr15/manga.Pipe3D-v2_4_3.fits')
p3dtab.remove_columns(['objra','objdec','mangaid'])

# Final summary FITS binary table for the MaNGA Data Reduction Pipeline
# Download from https://www.sdss4.org/dr15/manga/manga-data/catalogs/
drptab = Table.read('/Users/tonywong/Scratch2/almaquest/manga/drpall-v2_4_3.fits', hdu=1)
# From documentation for nsa_zdist: Distance estimate using peculiar velocity model of Willick et al. (1997); mulitply by c/Ho for Mpc
drptab['nsa_z_dMpc'] = drptab['nsa_zdist'] * (const.c/(70*u.km/(u.s*u.Mpc))).to(u.Mpc)

# Join the tables
gtab = join(p3dtab, drptab, join_type='left', keys=['plateifu'])
print(gtab.colnames)

keepcols = ['plateifu', 'mangaid', 'objra', 'objdec', 'seemed', 'gfwhm',
            'rfwhm', 'ifwhm', 'zfwhm', 'nsa_z', 'nsa_sersic_ba',
            'nsa_sersic_phi', 'nsa_sersic_n', 'nsa_elpetro_ba',
            'nsa_elpetro_phi', 'nsa_elpetro_th50_r', 'nsa_z_dMpc', 'dl',
            're_arc', 're_kpc', 'pa', 'ellip', 'log_mass', 'e_log_mass',
            'log_sfr_ha', 'e_log_sfr_ha', 'log_sfr_ssp', 'e_log_sfr_ssp',
            'log_nii_ha_cen', 'e_log_nii_ha_cen', 'log_oiii_hb_cen',
            'e_log_oiii_hb_cen', 'log_sii_ha_cen', 'e_log_sii_ha_cen',
            'log_oii_hb_cen', 'e_log_oii_hb_cen', 'ew_ha_cen', 'e_ew_ha_cen',
            'age_lw_re_fit', 'e_age_lw_re_fit', 'age_mw_re_fit',
            'e_age_mw_re_fit', 'vel_sigma_re', 'e_vel_sigma_re', 'sigma_cen',
            'e_sigma_cen', 'sigma_cen_ha', 'e_sigma_cen_ha', 'av_gas_re',
            'e_av_gas_re', 'av_ssp_re', 'e_av_ssp_re', 'oh_re_fit_n2',
            'e_oh_re_fit_n2', 'alpha_oh_re_fit_n2', 'e_alpha_oh_re_fit_n2',
            'oh_re_fit_o3n2', 'e_oh_re_fit_o3n2', 'alpha_oh_re_fit_o3n2',
            'e_alpha_oh_re_fit_o3n2' ]

gtab = gtab[keepcols]
gtab.meta.clear()

gtab.rename_column('plateifu','Name')
gtab['Name'].description = 'MaNGA Name'

# Flag negative values for Re
gtab['re_arc'][np.where(gtab['re_arc'] < 0)] = np.nan
gtab['re_kpc'][np.where(gtab['re_kpc'] < 0)] = np.nan

# Set the units
gtab['objra'].unit  = 'deg'
gtab['objdec'].unit = 'deg'
gtab['seemed'].unit = 'arcsec'
gtab['gfwhm'].unit  = 'arcsec'
gtab['rfwhm'].unit  = 'arcsec'
gtab['ifwhm'].unit  = 'arcsec'
gtab['zfwhm'].unit  = 'arcsec'
gtab['nsa_sersic_phi'].unit = 'deg'
gtab['nsa_elpetro_phi'].unit = 'deg'
gtab['nsa_elpetro_th50_r'].unit = 'arcsec'
gtab['dl'].unit = 'Mpc'
gtab['re_arc'].unit = 'arcsec'
gtab['re_kpc'].unit = 'kpc'
gtab['pa'].unit = 'deg'
gtab['log_mass'].unit = 'dex(solMass)'
gtab['e_log_mass'].unit = 'dex(solMass)'
gtab['log_sfr_ha'].unit = 'dex(solMass/yr)'
gtab['e_log_sfr_ha'].unit = 'dex(solMass/yr)'
gtab['log_sfr_ssp'].unit = 'dex(solMass/yr)'
gtab['e_log_sfr_ssp'].unit = 'dex(solMass/yr)'
gtab['ew_ha_cen'].unit = 'Angstrom'
gtab['e_ew_ha_cen'].unit = 'Angstrom'
gtab['age_lw_re_fit'].unit = 'dex(yr)'
gtab['e_age_lw_re_fit'].unit = 'dex(yr)'
gtab['age_mw_re_fit'].unit = 'dex(yr)'
gtab['e_age_mw_re_fit'].unit = 'dex(yr)'

for name in (gtab.colnames):
    if name.startswith(('oh','e_oh')):
        gtab[name].unit = 'dex'
    elif name.startswith(('av','e_av')):
        gtab[name].unit = 'mag'
    elif name.startswith(('ew','e_ew')):
        gtab[name].unit = 'Angstrom'
    elif name.startswith(('sigma','e_sigma')):
        gtab[name].unit = 'km/s'

# Fill masked values with NaN
for coln in gtab.colnames[4:]:
    if hasattr(gtab[coln], 'mask'):
        gtab[coln] = gtab[coln].filled(fill_value=np.nan)
    if gtab[coln].dtype.name.startswith('float'):
        gtab.replace_column(coln, gtab[coln].astype(np.float32))

# Add the calculated inclination
min_ratio = 0.13
ratio_ba = np.clip(gtab['nsa_elpetro_ba'], a_min=min_ratio, a_max=None)
cos_theta_inc = np.sqrt((ratio_ba**2-min_ratio**2)/(1.0-min_ratio**2))
incdeg = Column(np.degrees(np.arccos(cos_theta_inc)), name='nsa_inclination', 
    description='Inclination from nsa_elpetro_ba with q_min of '+str(min_ratio),
    unit='degree')
gtab.add_column(incdeg, index=16)

# Write the table
gtab.meta['date'] = datetime.today().strftime('%Y-%m-%d')
gtab.meta['comments'] = ('Galaxy properties determined from DRP and Pipe3D for MaNGA DR15')
print(gtab.meta)
gtab.write('manga_global_dr15.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

