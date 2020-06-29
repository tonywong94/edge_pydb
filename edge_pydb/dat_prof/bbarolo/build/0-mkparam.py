#!/usr/bin/env python

import numpy as np
from astropy.table import Table, join
from astropy.io.fits import getheader
from astropy.wcs import WCS
import os.path
import glob
from edge_pydb import EdgeTable

def write_par(gallist, run='nrad8', template='edge_bb.par', edgedir='../../carmaedge', db=None):
    # Prepare output directory
    if not os.path.exists(run):
        os.makedirs(run)
    else:
        for f in glob.glob(run+'/*.par'):
            os.remove(f)

    # List of galaxies with override positions
    if os.path.exists('gallist_fitctr.txt'):
        with open('gallist_fitctr.txt') as f:
            namelist_fitctr = f.read().splitlines()
        gallist_fitctr = [gal for gal in namelist_fitctr if not gal.startswith("#")]
    else:
        gallist_fitctr = []

    # Input template
    with open(template) as p:
        paramlist = p.read()

    for gal in gallist:
        if run.startswith('smo5'):
            fitsin = edgedir+'/comb_de_10/smo5/'+gal+'.co.smo5norm.fits.gz'
            maskin = 'FILE('+edgedir+'/comb_de_10/smo5/mom/'+gal+'.co.de10_smo5_dil.mask.fits.gz)'
        elif run.startswith('smo7'):
            fitsin = edgedir+'/comb_de_10/smo7/'+gal+'.co.smo7norm.fits.gz'
            maskin = 'FILE('+edgedir+'/comb_de_10/smo7/mom/'+gal+'.co.de10_smo7_dil.mask.fits.gz)'
        else:
            fitsin = edgedir+'/comb_de_10/native/'+gal+'.co.cmnorm.fits'
            maskin = 'FILE('+edgedir+'/comb_de_10/native/mom/'+gal+'.co.de10_dil.mask.fits.gz)'
        if not os.path.exists(fitsin):
            print('Galaxy {} not found'.format(gal))
            continue
        if 'bbmsk' in run: 
            mask = 'SMOOTH'
        else:
            mask = maskin
        if 'fixvd' in run:
            free = 'VROT VSYS PA'
        else:
            free = 'VROT VDISP VSYS PA'
        hdr = getheader(fitsin)
        w = WCS(hdr)
        # --- Number of rings
        if gal in ['NGC5908','NGC6060','NGC6361']:
            nrad = 12
        elif gal in ['UGC04132']:
            nrad = 11
        elif gal in ['IC0944','IC2487','NGC2410','NGC5980','NGC6478','UGC05111']:
            nrad = 10
        elif gal in ['NGC3994','NGC4149','NGC5953','UGC09067']:
            nrad = 7
        elif gal in ['NGC5657','NGC7819']:
            nrad = 6
        elif gal in ['NGC0447','NGC4676A','UGC05108']:
            nrad = 5
        else:
            nrad = 8
        # --- Set default parameters for VROT and VDISP
        if gal in ['NGC5784']:
            vrot = 600.
        elif gal in ['NGC2639']:
            vrot = 400.
        elif gal in ['NGC0496','NGC0551','NGC4210','NGC5480']:
            vrot = 150
        elif gal in ['NGC4961','NGC5016','NGC5520','NGC6155','UGC04461','UGC09542']:
            vrot = 100
        else:
            vrot = 200.
        vdisp = 8.
        if len(np.where(db['Name'] == gal)[0])==0:
            print('Skipping {} because it is not in db'.format(gal))
            continue
        i = np.where(db['Name'] == gal)[0][0]
        vsys = db['coVsys'][i]
        # --- Get center position from NED, convert to pixel units
        if gal in ['ARP220','NGC0496','NGC0523','NGC6155','UGC03973','UGC10043','UGC10123']:
            ractr = db['nedRA'][i]
            dcctr = db['nedDE'][i]
        else:
            ractr = db['ledaRA'][i]
            dcctr = db['ledaDE'][i]
        xposdg = hdr['crval1']
        yposdg = hdr['crval2']
        xpospx = hdr['crpix1']
        ypospx = hdr['crpix2']
        print('Working on galaxy {}'.format(gal))
        pixcrd = np.array([[xpospx, ypospx, 1, 1]])
        refpos = w.wcs_pix2world(pixcrd, 1)
        nedpx  = w.wcs_world2pix([[ractr,dcctr,refpos[0][2],refpos[0][3]]],0)
        print('  LEDA pixel (0-based) is {:.2f} {:.2f}'.format(nedpx[0][0],nedpx[0][1]))
        tstpx  = w.wcs_world2pix([[refpos[0][0],refpos[0][1],refpos[0][2],refpos[0][3]]],0)
        print('  Test: Ref  pixel (0-based) is {:.2f} {:.2f} should be {:.2f} {:.2f}'.
              format(tstpx[0][0],tstpx[0][1],xpospx-1,ypospx-1))
        # --- Use RA and DEC from NED
        xpos   = nedpx[0][0]
        ypos   = nedpx[0][1]
        # --- Use RA and DEC from 2D Gaussian fitting for certain galaxies
        if gal in gallist_fitctr:
            gal_table = Table.read('fitctr/'+gal+'_fitctr.txt',format='ascii.csv')
            ractr_fit = gal_table['ractr_fit'][0]
            dcctr_fit = gal_table['dcctr_fit'][0]
            fitpx = w.wcs_world2pix([[ractr_fit,dcctr_fit,refpos[0][2],refpos[0][3]]],0)
            xpos = fitpx[0][0]
            ypos = fitpx[0][1]
            print('  Using fitted center of {:.2f} {:.2f}'.format(fitpx[0][0],fitpx[0][1]))
        # --- Get PA and INC from LEDA
        inc = db['ledaAxIncl'][i]
        if (inc > 88):
            inc=88.0
        pa = (db['rfPA'][i] + 180) % 360
        # --- Use default values for distance and thickness = 100 pc
        dmpc = db['caDistMpc'][i]
        z0 = 206265*100/(dmpc*1e6)  # 100 pc thickness, fixed
        print('  Assumed INC, PA, Z0: {:.2f} {:.2f} {:.2f}'.format(inc,pa,z0))
        gal_param = paramlist % (fitsin, nrad, vsys, xpos, ypos, vrot, inc, pa, z0, free, mask)
        file = open(run+'/param_'+gal+'.par','w')
        file.write(gal_param)
        file.close()
    print (run+' Done')
    return

# CALIFA table: source for DISTANCE
db = EdgeTable('edge_califa.csv', cols=['Name', 'caDistMpc'])
# NED table: source for CENTER RA & DEC
ned = EdgeTable('edge_ned.csv', cols=['Name', 'nedRA', 'nedDE'])
db.join(ned)
# LEDA table: source for INC, CENTER RA & DEC
leda = EdgeTable('edge_leda.csv', cols=['Name', 'ledaRA', 'ledaDE', 'ledaPA', 'ledaAxIncl'])
leda['ledaRA'].convert_unit_to('deg')
leda['ledaRA'].format = '.5f'
db.join(leda)
# CO observations table: source for VSYS
coobs = EdgeTable('edge_coobs_DE.csv', cols=['Name', 'coVsys', 'coTpk_10'])
db.join(coobs)
# Becca's fits: source for PA
rfpars = EdgeTable('edge_rfpars.csv', cols=['Name', 'rfKinRA', 'rfKinDecl', 'rfPA', 'rfInc'])
db.join(rfpars)
print(db.keys())

# Get the list of galaxies to work on
# listfile = 'detected.txt'
listfile = 'resolved.txt'
with open(listfile) as f:
    namelist = f.read().splitlines()
gallist = [gal for gal in namelist if not gal.startswith("#")]
print (gallist)

#masks = ['dilmsk', 'bbmsk']
masks = ['dilmsk']
fits  = ['fitvd', 'fixvd']
#sets  = ['natv', 'smo5', 'smo7']
sets  = ['natv', 'smo7']
runs  = []
for set in sets:
    for fit in fits:
        for mask in masks:
            runs.append(set+'_'+fit+'_'+mask)
print(runs)

edgedir = os.path.normpath(os.getcwd()+'/..')
for run in runs:
    write_par(gallist, run=run, template='edge_bb.par', edgedir=edgedir, db=db)    

