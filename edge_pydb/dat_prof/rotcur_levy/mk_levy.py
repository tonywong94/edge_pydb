#!/usr/bin/env python

import glob
import os
from astropy.table import Table,Column,vstack,join
from datetime import datetime

sets = ['CO_natv', 'CO_smo6', 'HA_natv', 'HA_smo6']

for set in sets:

    files=glob.glob(set+'/*RC.txt')

    tablelist=[]
    for file in files:
        txttab = Table.read(file, format='ascii', comment='%', 
                  names=('radius','Vrot','e_Vrot','Vrad','e_Vrad','dVsys','e_dVsys'))
        gal = os.path.basename(file).split('.')[0]
        # --- Add the Galaxy Name
        nrows = len(txttab)
        if nrows > 0:
            gname = Column([gal]*nrows, name='Name', description='Galaxy Name')
            txttab.add_column(gname, index=0)
            tablelist.append(txttab)

    t_rfcurve = vstack(tablelist)
    t_rfcurve.pprint()
    t_rfcurve['radius'].unit = 'arcsec'
    t_rfcurve['Vrot'].description = 'Rotation curve'
    t_rfcurve['Vrot'].unit = 'km / s'
    t_rfcurve['e_Vrot'].description = 'Error in rotation curve'
    t_rfcurve['e_Vrot'].unit = 'km / s'
    t_rfcurve['Vrad'].description = 'Expansion velocity'
    t_rfcurve['e_Vrad'].unit = 'km / s'
    t_rfcurve['e_Vrad'].description = 'Error in expansion velocity'
    t_rfcurve['e_Vrad'].unit = 'km / s'
    t_rfcurve['dVsys'].description = 'Systemic velocity offset'
    t_rfcurve['dVsys'].unit = 'km / s'
    t_rfcurve['e_dVsys'].description = 'Error in systemic velocity offset'
    t_rfcurve['e_dVsys'].unit = 'km / s'
    t_rfcurve.meta['comments'] = (set+' Rotation Curves from 2018ApJ...860...92L')
    t_rfcurve.meta['date'] = datetime.today().strftime('%Y-%m-%d')
    t_rfcurve.write('rf_'+set+'.csv', overwrite=True, delimiter=',', format='ascii.ecsv')
