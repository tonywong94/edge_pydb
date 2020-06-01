#!/usr/bin/env python

from astropy.table import Table, join
import numpy as np
from datetime import datetime

# Process coobs tables into ECSV format

# D-array, 20 km/s channels
t = Table.read('edge_coobs_d20.csv', format='ascii.csv')
t['Name'].description = 'Galaxy Name'
t['coVsys'].description = 'LSR Radio Velocity'
t['coVsys'].unit = 'km / s'
t['coObstim'].name = 'coDobstim'
t['coDobstim'].description = 'D array observing time'
t['coDobstim'].unit = 'hour'
t['coNpt'].description = 'No. of CARMA pointings'
t['Bmaj'].name = 'coDbmaj'
t['coDbmaj'].description = 'Native D-array beam major axis'
t['coDbmaj'].unit = 'arcsec'
t['Bmin'].name = 'coDbmin'
t['coDbmin'].description = 'Native D-array beam minor axis'
t['coDbmin'].unit = 'arcsec'
t['coRMS_mJybm'].description = 'Noise in 20 km/s channel'
t['coRMS_mJybm'].unit = 'mJy / beam'
t['coRMS_mK'].description = 'Noise in 20 km/s channel'
t['coRMS_mK'].unit = 'mK'
t['coTpk_mK'].description = 'Peak brightness in 20 km/s channel'
t['coTpk_mK'].unit = 'mK'
t['coSNRpeak'].description = 'Peak brightness in SNR units'
t.remove_column('co_mmom0_flux_Kkmsas^2')
t['Last_updated'].name = 'ImagingDate'
t['ImagingDate'].description = 'Date of imaging'
t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
print(t.meta)
t.write('edge_coobs_D.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# E-array, 20 km/s channels
t = Table.read('edge_coobs_e20.csv', format='ascii.csv')
t['Name'].description = 'Galaxy Name'
t['Dname'].name = 'coDname'
t['coDname'].description = 'Name repeated if also observed in D'
t['coVsys'].description = 'LSR Radio Velocity'
t['coVsys'].unit = 'km / s'
t['coObstim'].name = 'coEobstim'
t['coEobstim'].description = 'E array observing time'
t['coEobstim'].unit = 'hour'
t['coNpt'].description = 'No. of CARMA pointings'
t['Bmaj'].name = 'coEbmaj'
t['coEbmaj'].description = 'Native E-array beam major axis'
t['coEbmaj'].unit = 'arcsec'
t['Bmin'].name = 'coEbmin'
t['coEbmin'].description = 'Native E-array beam minor axis'
t['coEbmin'].unit = 'arcsec'
t['coRMS_mJybm'].description = 'Noise in 20 km/s channel'
t['coRMS_mJybm'].unit = 'mJy / beam'
t['coRMS_mK'].description = 'Noise in 20 km/s channel'
t['coRMS_mK'].unit = 'mK'
t['coTpk_mK'].description = 'Peak brightness in 20 km/s channel'
t['coTpk_mK'].unit = 'mK'
t['coSNRpeak'].description = 'Peak brightness in SNR units'
t.remove_column('co_mmom0_flux_Kkmsas^2')
t['Last_updated'].name = 'ImagingDate'
t['ImagingDate'].description = 'Date of imaging'
t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
print(t.meta)
t.write('edge_coobs_E.csv', format='ascii.ecsv', delimiter=',', overwrite=True)

# D + E array, 10 and 20 km/s channels
t1 = Table.read('edge_coobs_de10.csv', format='ascii.csv')
t1.remove_column('coRMS_mJybm')
t1.remove_column('co_mmom0_flux_Kkmsas^2')
t2 = Table.read('edge_coobs_de20.csv', format='ascii.csv')
for col in ['coVsys', 'coObstim', 'coNpt', 'Bmaj', 'Bmin',
            'coRMS_mJybm', 'co_mmom0_flux_Kkmsas^2']:
    t2.remove_column(col)
t = join(t1,t2,keys='Name')
t['Name'].description = 'Galaxy Name'
t['coVsys'].description = 'LSR Radio Velocity'
t['coVsys'].unit = 'km / s'
t['coObstim'].name = 'coDEobstim'
t['coDEobstim'].description = 'D+E array observing time'
t['coDEobstim'].unit = 'hour'
t['coNpt'].description = 'No. of CARMA pointings'
t['Bmaj'].name = 'coDEbmaj'
t['coDEbmaj'].description = 'Native D+E-array beam major axis'
t['coDEbmaj'].unit = 'arcsec'
t['Bmin'].name = 'coDEbmin'
t['coDEbmin'].description = 'Native D+E-array beam minor axis'
t['coDEbmin'].unit = 'arcsec'
for i in ['1','2']:
    t['coRMS_mK_'+i].name = 'coRMS_'+i+'0'
    t['coRMS_'+i+'0'].description = 'Noise in '+i+'0 km/s channel'
    t['coRMS_'+i+'0'].unit = 'mK'
    t['coTpk_mK_'+i].name = 'coTpk_'+i+'0'
    t['coTpk_'+i+'0'].description = 'Peak brightness in '+i+'0 km/s channel'
    t['coTpk_'+i+'0'].unit = 'mK'
    t['coSNRpeak_'+i].name = 'coSNRpeak_'+i+'0'
    t['coSNRpeak_'+i+'0'].description = 'Peak brightness at '+i+'0 km/s in SNR units'
    t['Last_updated_'+i].name = 'ImagingDate_'+i+'0'
    t['ImagingDate_'+i+'0'].description = 'Date of '+i+'0 km/s imaging'
t.meta['date'] = datetime.today().strftime('%Y-%m-%d')
print(t.meta)
t.write('edge_coobs_DE.csv', format='ascii.ecsv', delimiter=',', overwrite=True)
