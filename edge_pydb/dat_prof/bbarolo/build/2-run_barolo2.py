#!/usr/bin/env python

# This script tries repeating the runs which failed the first time.
# Galaxies which still fail are listed in nrad8_fail2.txt and nrad8_vdisp8_fail2.txt

import os
import glob
import subprocess
import shutil

os.system('/Users/tonywong/Work/bin/bbarolo/BBarolo -v')

tout = 180
masks = ['dilmsk', 'bbmsk']
fits = ['fitvd', 'fixvd']
sets = ['natv', 'smo5', 'smo7']
runs = []
for set in sets:
    for fit in fits:
        for mask in masks:
            runs.append(set+'_'+fit+'_'+mask)
print(runs)
#runs = ['nrad8', 'nrad8_vdisp8', 'nrad8_smolist', 'nrad8_vdisp8_smolist']

for run in runs:
    with open(run+'_fail1.txt') as f:
        namelist = f.read().splitlines()
    noPlot_txt = open(run+'_fail2.txt','w')
    noPlot=[]
    os.chdir(run)
    for gal in namelist:
        print(gal)
        pfile = 'param_'+gal+'.par'
        if os.path.isdir('output/'+gal[:8]):
            shutil.rmtree('output/'+gal[:8])
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
 
