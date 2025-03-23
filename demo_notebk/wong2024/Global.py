#!/usr/bin/env python

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.table import Table, join
from astropy import units as u
from edge_pydb import EdgeTable
from edge_pydb.beam_sample import beam_sample

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# Calculate molecular mass for 1 Jy km/s at 1 Mpc.
# These serve as fiducial distance and flux which we will scale later.
sco  = 1 * u.Jy * u.km / u.s
dmpc = (1 * u.Mpc).to(u.pc)
alphaco = 4.3 * u.solMass * u.s / (u.K * u.km * u.pc**2) # Bolatto+13
freq = 115.271204 * u.GHz
kperjy = (1*u.Jy/u.sr).to(u.K, equivalencies=u.brightness_temperature(freq))
molmass = alphaco * kperjy * (sco/u.Jy) * dmpc**2
print('At 1 Mpc, 1 Jy km/s of CO flux is {} of H_2'.format(molmass))

# Measured CO fluxes
db = EdgeTable('edge_coflux_smo7.csv')
# CALIFA params incl. distance
califa = EdgeTable('edge_califa.csv')
db.join(califa)
print(db.colnames)
print('Stellar mass units:', db['caMstars'].unit)
print('SFR units:', db['caSFR'].unit)

fig, (ax1,ax2) = plt.subplots(1,2, figsize=(7,4))

# Define the subsets
valid_ssfr = (~np.isnan(db['caMstars'])) & (~np.isnan(db['caSFR']))
valid_mgas = (~np.isnan(db['caMstars'])) & (~np.isnan(db['caSFR']))
det  = valid_ssfr & (~np.isnan(db['coDilated_smo7']))
print('Number of CO detections:',np.count_nonzero(det))
ndet = valid_ssfr & (np.isnan(db['coDilated_smo7']))
print('Number of CO non-detections:',np.count_nonzero(ndet))

# Plot the data
ax1.scatter(db['caMstars'][det],db['caSFR'][det],
            c='tab:blue', alpha=0.5, label='CO det')
ax1.scatter(db['caMstars'][ndet],db['caSFR'][ndet], marker='s',
            c='tab:red', alpha=0.5, s=20, label='CO ndet')

# Main sequence from Cano-Diaz et al. (2016ApJ...821L..26C)
x_ms = np.linspace(9,12,num=50)
y_ms = 0.81*x_ms-8.34
ax1.plot(x_ms,y_ms,'m--', zorder=-2)
ax1.legend(loc='lower right', handletextpad=0.01, fontsize=9)

# Axis labels
ax1.set_xlabel('log($M_*$ [$M_\odot$])', fontsize=14)
ax1.set_ylabel('log(SFR [$M_\odot$ yr$^{-1}$])', fontsize=14)
ax1.set_aspect('equal')
ax1.set_xlim(8.75,12.25)
ax1.set_ylim(-3,2)

# Annotations
ax1.text(0.04,0.94, '(a)', size=12, ha='left',
         va='center', transform=ax1.transAxes)
ax1.text(0.92,0.92, 'MS', size=12, ha='center', color='m',
         va='center', transform=ax1.transAxes)

# Scale noise estimate from unmasked moment map
print('Median unmasked velocity width is',np.nanmedian(db['coNomaskDv_smo7']))
print('Median masked velocity width is',np.nanmedian(db['coDilatedDv_smo7']))
nsefactor = np.sqrt(np.nanmedian(db['coDilatedDv_smo7'])/np.nanmedian(db['coNomaskDv_smo7']))
print('Noise estimates will be scaled down by',nsefactor)

# Plot the data
delsfr  = db['caSFR'] - (0.81*db['caMstars']-8.34)
mgas = np.log10(molmass.value * db['coDilated_smo7'] * db['caDistMpc']**2)
mgas[ndet] = np.log10(molmass.value * 2*nsefactor*db['coeNomask_smo7'][ndet]
                      * db['caDistMpc'][ndet]**2)
fgas = mgas - db['caMstars']
ax2.scatter(fgas[det],delsfr[det],c='tab:blue', alpha=0.5)

# Plot the upper limits
uplims = np.zeros(fgas[ndet].shape)
uplims[:] = True
ax2.errorbar(fgas[ndet], delsfr[ndet], xuplims=uplims, 
             xerr=0.2, ls='none', color='tab:red', zorder=-2, alpha=0.7)

# Axis labels
ax2.set_ylabel(r'log(SFR) - log(SFR$_{\rm MS}$)', fontsize=14)
ax2.set_xlabel('log($M(H_2)/M_*$)', fontsize=14)
ax2.set_aspect('equal')
ax2.set_xlim(-3,0)
ax2.set_ylim(-3,2)
ax2.text(0.05,0.94, '(b)', size=12, ha='left',
         va='center', transform=ax2.transAxes)

plt.subplots_adjust(wspace=0.1)
plt.savefig('sfms.pdf', bbox_inches='tight')
