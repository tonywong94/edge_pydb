#!/usr/bin/env python

# This script takes the galaxies which failed twice and changes the MASK
# to SMOOTH and runs Bbarolo again with both varying and fixed VDISP.
# Galaxies which fail are listed in nrad8_fail3.txt and nrad8_vdisp8_fail3.txt

import os
import glob
import subprocess
import shutil

os.system('/Users/tonywong/Work/bin/bbarolo/BBarolo -v')

tout = 180
runs = ['nrad8', 'nrad8_vdisp8']

for run in runs:
    with open(run+'_fail2.txt') as f:
        namelist = f.read().splitlines()
    noPlot_txt = open(run+'_fail3.txt','w')
    noPlot=[]
    if os.path.isdir('./gal_'+run+'_smolist'):
        shutil.rmtree('./gal_'+run+'_smolist')
    os.makedirs('./gal_'+run+'_smolist')
    os.chdir('./gal_'+run+'_smolist')
    for gal in namelist:
        print(gal)
        oldpfile = '../gal_'+run+'/param_'+gal+'.par'
        pfile = 'param_'+gal+'.par'
        # Change the MASK to SMOOTH
        with open(oldpfile) as p:
            with open(pfile,'w') as p_smo:
                for line in p:
                    rest = line.split('MASK        ')
                    if len(rest) == 2:
                        p_smo.write(line.replace(rest[1],'SMOOTH'))
                    else:
                        p_smo.write(line)
        os.makedirs('output/'+gal[:8])
        with open('output/'+gal[:8]+'/'+gal+'.bblog.txt', "wb") as outfile:
            try:
                subprocess.run(['/Users/tonywong/Work/bin/bbarolo/BBarolo','-p',pfile], 
                            timeout=tout, stdout=outfile, stderr=outfile)
            except subprocess.TimeoutExpired:
                print('Timeout after {} seconds'.format(tout))
        plot_paths=['/plot_pv_local.pdf','/plot_parameters.pdf','/plot_maps_local.pdf','/plot_chanmaps_local.pdf']
        for plot_path in plot_paths:
            if not os.path.isfile('output/'+gal[:8]+plot_path):
                noPlot.append(gal)
                print(gal+' in '+ run + ' is missing plots')
                break
    os.chdir('..')
    for obj in noPlot:
        noPlot_txt.write(obj+'\n')
    noPlot_txt.close()
    print (run+' Done')
 
