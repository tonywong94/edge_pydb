import numpy as np
from astropy.table import Table, join
from astropy.io.fits import getheader
from astropy.wcs import WCS
import os.path
import glob

def write_par(gallist, run='nrad8', template='edge_bb.par', edgedir='../../carmaedge', db=None):
    # Prepare output directory
    if not os.path.exists('gal_'+run):
        os.makedirs('gal_'+run)
    else:
        for f in glob.glob('gal_'+run+'/*.par'):
            os.remove(f)

    # List of galaxies with override positions
    with open('gallist_fitctr.txt') as f:
        namelist_fitctr = f.read().splitlines()
    gallist_fitctr = [gal for gal in namelist_fitctr if not gal.startswith("#")]

    # Input template
    with open('templates/'+template) as p:
        paramlist = p.read()

    for gal in gallist:
        if '20' in run:
            fitsin = edgedir+'comb_de/cmnorm_sub/'+gal+'.co.cmnormsub.fits.gz'
            maskdir = edgedir+'comb_de/reprojmask/'
        else:
            fitsin = edgedir+'comb_de_10/cmnorm_sub/'+gal+'.co.cmnormsub.fits.gz'
            maskdir = edgedir+'comb_de_10/reprojmask/'
        if 'smolist' in run: 
            mask = 'SMOOTH'
        else:
            mask = 'FILE('+maskdir+gal+'.co.de10_dil.masksub.fits.gz)'
        if 'vdisp' in run:
            free = 'VROT VSYS PA'
        else:
            free = 'VROT VDISP VSYS PA'
        hdr = getheader(fitsin)
        hdr['specsys'] = 'LSRK'
        hdr['ctype3']  = 'VRAD'
        hdr['velref']  = 257
        w = WCS(hdr)
        # --- Set default parameters for VROT and VDISP
        vrot = 200.
        vdisp = 8.
        if len(np.where(db['Name'] == gal)[0])==0:
            print('Skipping {} because it is not in db'.format(gal))
            continue
        i = np.where(db['Name'] == gal)[0][0]
        vsys = db['coVsys'][i]
        # --- Get center position from NED, convert to pixel units
        ractr = db['nedRA'][i]
        dcctr = db['nedDE'][i]
        xposdg = hdr['crval1']
        yposdg = hdr['crval2']
        xpospx = hdr['crpix1']
        ypospx = hdr['crpix2']
        print('Working on galaxy {}'.format(gal))
        pixcrd = np.array([[xpospx, ypospx, 1, 1]])
        refpos = w.wcs_pix2world(pixcrd, 1)
        nedpx  = w.wcs_world2pix([[ractr,dcctr,refpos[0][2],refpos[0][3]]],0)
        print('  NED  pixel (0-based) is {:.2f} {:.2f}'.format(nedpx[0][0],nedpx[0][1]))
        tstpx  = w.wcs_world2pix([[refpos[0][0],refpos[0][1],refpos[0][2],refpos[0][3]]],0)
        print('  Test: Ref  pixel (0-based) is {:.2f} {:.2f} should be {:.2f} {:.2f}'.
              format(tstpx[0][0],tstpx[0][1],xpospx-1,ypospx-1))
        # --- Use RA and DEC from NED
        xpos   = nedpx[0][0]
        ypos   = nedpx[0][1]
        # --- Use RA and DEC from 2D Gaussian fitting for certain galaxies
        if gal in namelist_fitctr:
            gal_table = Table.read('fitctr/'+gal+'_fitctr.txt',format='ascii.csv')
            ractr_fit = gal_table['ractr_fit'][0]
            dcctr_fit = gal_table['dcctr_fit'][0]
            fitpx = w.wcs_world2pix([[ractr_fit,dcctr_fit,refpos[0][2],refpos[0][3]]],0)
            xpos = fitpx[0][0]
            ypos = fitpx[0][1]
            print('  Using fitted center of {:.2f} {:.2f}'.format(fitpx[0][0],fitpx[0][1]))
        # --- Get PA and INC from Becca's fits
        inc = db['rfInc'][i]
        if (inc > 88):
            inc=88.0
        pa = (db['rfPA'][i] + 180) % 360
        # --- Use default values for distance and thickness = 100 pc
        dmpc = db['caDistMpc'][i]
        z0 = 206265*100/(dmpc*1e6)  # 100 pc thickness, fixed
        print('  Assumed INC, PA, Z0: {:.2f} {:.2f} {:.2f}'.format(inc,pa,z0))
        gal_param = paramlist % (fitsin, vsys, xpos, ypos, inc, pa, z0, free, mask)
        file = open('./gal_'+run+'/param_'+gal+'.par','w')
        file.write(gal_param)
        file.close()
    print (run+' Done')
    return

#edgedir = '/home/jovyan/carmaedge/'
edgedir = '/Volumes/Scratch2/tonywong/EDGE/'
basedir = edgedir+'edge-sql-base/'
# CALIFA table: source for distance
db_ca = Table.read(basedir+'global_values/external/edge_califa.csv', format='ascii.ecsv')
# CO observations table: source for VSYS
db_co = Table.read(basedir+'global_values/obs/edge_coobs_de20.csv', format='ascii.csv')
# Becca's fits: source for PA, INC
db_rf = Table.read(basedir+'global_values/derived/edge_rfpars.csv', format='ascii.ecsv')
# NED table: source for XPOS, YPOS
db_nd = Table.read(basedir+'global_values/external/edge_ned.csv', format='ascii.ecsv')

db12 = join(db_ca, db_co, keys='Name')
db123 = join(db12, db_rf, keys='Name')
db = join(db123, db_nd, keys='Name')
print (db.keys())

# Get the list of galaxies to work on
listfile = 'detected.txt'
with open(listfile) as f:
    namelist = f.read().splitlines()
gallist = [gal for gal in namelist if not gal.startswith("#")]
print (gallist)

runs = ['nrad8', 'nrad8_vdisp8']
for run in runs:
    write_par(gallist, run=run, template='edge_bb.par', edgedir=edgedir, db=db)    

