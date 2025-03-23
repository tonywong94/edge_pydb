import numpy as np
from edge_pydb import EdgeTable
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from edge_pydb.plotting import dotpatch, imarrayplot

ig, axs = plt.subplots(1, 3, figsize=(12,4), sharex=True, sharey=True)

# Left panel
alltab = EdgeTable('edge_carma_allpix.2d_smo7.hdf5', path='flux_elines_sm')
sel = (alltab['Name'] == 'NGC3687')
npix = np.count_nonzero(~np.isnan(alltab[sel]['flux_Halpha_sm']))
xmid = alltab[sel]['ix'][np.argmin(alltab[sel]['rad_arc'])]
ymid = alltab[sel]['iy'][np.argmin(alltab[sel]['rad_arc'])]
alltab['ix'] = xmid - alltab['ix']
alltab['iy'] = alltab['iy'] - ymid
img, xymin, xymax = dotpatch(alltab[sel]['ix'], alltab[sel]['iy'], alltab[sel]['flux_Halpha_sm'], 
                             axes=axs[0], vmin=0, vmax=2.5, cmap='jet', clipedge=False)
print('Center at',xmid,ymid)
axs[0].set_title('All pixels ('+str(npix)+')')
axs[0].set_xlim([40,-40])
axs[0].set_ylim([-35,45])
axs[0].set_aspect('equal')
axs[0].set_xlabel(r'$\Delta\alpha$ (")', fontsize='large')
axs[0].set_ylabel(r'$\Delta\delta$ (")', fontsize='large')
axs[0].yaxis.set_major_locator(ticker.MultipleLocator(20))

# Middle panel
samtab = EdgeTable('edge_carma.2d_smo7.hdf5', path='flux_elines_sm')
sel = (samtab['Name'] == 'NGC3687')
npix = np.count_nonzero(~np.isnan(samtab[sel]['flux_Halpha_sm']))
xmid = samtab[sel]['ix'][np.argmin(samtab[sel]['rad_arc'])]
ymid = samtab[sel]['iy'][np.argmin(samtab[sel]['rad_arc'])]
samtab['ix'] = xmid - samtab['ix']
samtab['iy'] = samtab['iy'] - ymid
img, xymin, xymax = dotpatch(samtab[sel]['ix'], samtab[sel]['iy'], samtab[sel]['flux_Halpha_sm'], 
                             axes=axs[1], vmin=0, vmax=2.5, cmap='jet', clipedge=False)
print('Center at',xmid,ymid)
axs[1].set_title('Square grid ('+str(npix)+')')
axs[1].set_aspect('equal')
axs[1].set_xlabel(r'$\Delta\alpha$ (")', fontsize='large')

# Right panel
hextab = EdgeTable('edge_carma_hex.2d_smo7.hdf5', path='flux_elines_sm')
sel = (hextab['Name'] == 'NGC3687')
npix = np.count_nonzero(~np.isnan(hextab[sel]['flux_Halpha_sm']))
xmid = hextab[sel]['ix'][np.argmin(hextab[sel]['rad_arc'])]
ymid = hextab[sel]['iy'][np.argmin(hextab[sel]['rad_arc'])]
hextab['ix'] = xmid - hextab['ix']
hextab['iy'] = hextab['iy'] - ymid
img, xymin, xymax = dotpatch(hextab[sel]['ix'], hextab[sel]['iy'], hextab[sel]['flux_Halpha_sm'], 
                             axes=axs[2], vmin=0, vmax=2.5, cmap='jet', clipedge=False)
print('Center at',xmid,ymid)
axs[2].set_title('Hexagonal grid ('+str(npix)+')')
axs[2].set_aspect('equal')
axs[2].set_xlabel(r'$\Delta\alpha$ (")', fontsize='large')
#plt.colorbar(img, ax=axs[2])

plt.savefig('compare_grids.pdf', bbox_inches='tight')
